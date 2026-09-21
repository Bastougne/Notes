function krig = krigeage()
% KRIGEAGE  Le relevé et les modèles d'observation construits dessus.
%
%   krig = krigeage();
%   survey = krig.releve(map, par);
%   model  = krig.modele(kind, map, survey, par);
%
% Un modèle répond à un seul appel :
%
%   [z_hat, R_hat, n_used, n_modes, cout] = model.query(X, ctx)
%
% X est un état par colonne, dont seul le bloc position sert ; ctx est ce que le filtre
% sait du nuage, que les modèles adaptatifs lisent et que les statiques ignorent. R_hat
% est la variance déclarée, et c'est là tout l'enjeu : une interpolation rend un nombre
% sans incertitude, un processus gaussien rend les deux.
%
% CINQ MODÈLES :
%
%   'carte'       le champ vrai, la borne. R est la variance du capteur seule.
%   'bilineaire'  le relevé bruité lu bilinéairement, variance constante. La ligne de base
%                 sans krigeage.
%   'ok'          krigeage ordinaire sur tout le relevé, hyperparamètres ajustés une fois.
%                 L'estimateur est alors un couple fixe de fonctions de la position, donc
%                 tabulé une fois et relu bilinéairement.
%   'ak'          les mêmes équations sur les échantillons proches de là où le filtre se
%                 croit, choisis à chaque pas, hyperparamètres réajustés sur la fenêtre.
%   'cak'         le même, une fois par mode du nuage. Une postérieure bimodale a sa
%                 moyenne entre les deux modes, là où l'avion n'est pas.
%
% Krigeage ordinaire et non simple partout : sa covariance porte le terme qui paie
% l'estimation de la moyenne, celui qui élargit l'intervalle là où les échantillons
% manquent.
    krig = struct('releve', @build_survey, 'modele', @build_model, 'lire', @grid_read);
end

%% ------------------------------------------------------------------- relevé

function survey = build_survey(map, par)
% Le relevé : une grille régulière sur la carte, lue sur le champ vrai et dégradée par le
% bruit de relevé, plus les hyperparamètres ajustés dessus.
    lo    = map.step;                                   % un pas depuis le bord
    hi    = ([size(map.h, 2), size(map.h, 1)] - 2) .* map.step;

    % Le relevé se décrit de deux façons. Par son nombre de points, comme les sauvegardes.
    % Ou par sa MAILLE, et c'est la bonne pour la figure du chapitre : c'est la maille qui
    % est en abscisse, la fixer à des valeurs rondes et régulièrement espacées se lit sans
    % explication, et le nombre de points s'en déduit. La grille ne remplit alors plus
    % exactement la carte ; le reliquat est laissé en marge, ce qui est sans conséquence
    % puisque le corridor du vol est loin des bords.
    if isfield(par.krig, 'pas') && ~isempty(par.krig.pas)
        pas = [1, 1] * par.krig.pas;
        n_s = min(floor((hi - lo) ./ pas));
    else
        n_s = par.krig.n_side;
        pas = (hi - lo) ./ n_s;
    end
    reste = (hi - lo) - n_s * pas;

    % La graine du relevé est distincte de celle des essais : c'est elle que le pilote
    % fait varier pour que l'erreur de carte soit moyennée et non fixée.
    rng(graine_releve(par));

    % Échantillonnage au centre des cellules plutôt qu'aux bords : les n_s points laissent
    % alors une cellule entière de jeu, dans laquelle la grille peut coulisser sans sortir
    % de la carte. u est cette position, en fractions de maille.
    %
    % La tirer au sort décide de la PHASE du relevé par rapport au plan de vol, et c'est
    % nécessaire : à phase fixe, une maille commensurable avec l'espacement des branches
    % met les trois branches au même endroit de la cellule, le biais d'interpolation s'y
    % répète au lieu de se moyenner, et le balayage de maille lit une résonance de géométrie
    % là où il croit lire une finesse de relevé. La graine du relevé étant ce que le pilote
    % fait varier, une phase par tirage suffit à moyenner l'effet.
    if isfield(par.krig, 'phase_alea') && par.krig.phase_alea
        u = rand(1, 2);
    else
        u = [0.5, 0.5];                                 % centré, donc reproductible
    end
    ax    = lo(1) + reste(1) / 2 + ((0:n_s-1)' + u(1)) * pas(1);
    ay    = lo(2) + reste(2) / 2 + ((0:n_s-1)' + u(2)) * pas(2);
    [SX, SY] = meshgrid(ax, ay);
    survey.X     = [SX(:), SY(:)]';                     % 2 x n, un échantillon par colonne

    % Un relevé réel n'est pas un réseau exact : les lignes de vol qui l'ont produit ont
    % leur propre erreur de navigation. Secouer chaque point d'une fraction de maille est
    % donc le modèle honnête — et c'est surtout ce qui tue la résonance.
    %
    % Sur un réseau exact, l'erreur d'interpolation est une fonction PÉRIODIQUE de la
    % position dans la cellule : un biais déterministe, que deux branches distantes d'un
    % multiple entier de la maille subissent à l'identique. Aucune translation de la grille
    % n'y change rien, puisqu'une translation conserve les phases relatives. Secoué, le
    % champ d'erreur redevient aléatoire et se décorrèle en une maille : il se moyenne le
    % long de chaque branche au lieu de se répéter d'une branche à l'autre.
    if isfield(par.krig, 'jitter') && par.krig.jitter > 0
        survey.X = survey.X + par.krig.jitter * pas(:) .* (2 * rand(2, n_s^2) - 1);
        survey.X = min(max(survey.X, lo(:)), hi(:));    % le secouage ne sort pas de la carte
    end
    survey.z     = grid_read(map.h, map.step, survey.X) ...
                 + par.krig.sigma_map * randn(n_s^2, 1);
    survey.n      = n_s^2;
    survey.n_side = n_s;                                % déduit quand la maille est donnée
    survey.pitch  = pas;
    survey.phase  = u;
    survey.jitter = secouage(par);

    % Sur un sous-échantillon, le bruit tenu à sa valeur déclarée plutôt qu'ajusté :
    % laissé libre, fitrgp explique toute la rugosité du champ comme du bruit et rend un
    % ell très en dessous de la maille. La vraisemblance marginale a de plus un second
    % optimum en ell -> 0, où les optimiseurs tombent depuis leurs défauts — d'où le point
    % de départ imposé.
    i_ml = randperm(survey.n, min(par.krig.n_ml, survey.n));
    survey.hyp = fit_hyperparameters(survey.X(:, i_ml)', survey.z(i_ml), par);

    fprintf(['Relevé %d x %d = %d points, maille %.2f x %.2f km, phase %.2f/%.2f,' ...
             ' secouage %.0f %%, bruit %.0f %s\n'], ...
            n_s, n_s, survey.n, survey.pitch / 1e3, u, 100 * survey.jitter, ...
            par.krig.sigma_map, par.unit);
    fprintf('   sigma_f = %.0f %s, ell = %.2f km\n', ...
            survey.hyp.sigma_f, par.unit, survey.hyp.ell / 1e3);
end

function n = noyau(par)
% La famille de noyau, exponentielle quadratique par defaut. Les noms sont ceux de
% fitrgp, pour que l'ajustement et la resolution parlent de la meme chose.
    if isfield(par.krig, 'noyau')
        n = par.krig.noyau;
    else
        n = 'squaredexponential';
    end
end

function j = secouage(par)
% L'amplitude du secouage du relevé, en fractions de maille. Nulle par défaut : le réseau
% exact reste accessible, c'est le témoin qui fait apparaître la résonance.
    if isfield(par.krig, 'jitter'), j = par.krig.jitter;
    else, j = 0; end
end

function s = graine_releve(par)
% La graine du relevé, celle des essais par défaut.
    if isfield(par.krig, "survey_seed")
        s = par.krig.survey_seed;
    else
        s = par.run.seed;
    end
end

function hyp = fit_hyperparameters(X, z, par)
    if ~isempty(which('fitrgp'))
        gp = fitrgp(X, z, 'KernelFunction', noyau(par), ...
                    'BasisFunction', 'constant', 'Sigma', par.krig.sigma_map, ...
                    'ConstantSigma', true, 'Standardize', false, ...
                    'KernelParameters', [par.krig.ell_init; std(z)]);
        hyp.ell     = gp.KernelInformation.KernelParameters(1);
        hyp.sigma_f = gp.KernelInformation.KernelParameters(2);
        hyp.nom     = noyau(par);
    else
        hyp.ell     = par.krig.ell_init;      % pas de toolbox : la valeur initiale, et on le dit
        hyp.sigma_f = std(z);
        hyp.nom     = noyau(par);
        warning('krigeage:fitrgp', 'Pas de fitrgp, les hyperparamètres restent la valeur initiale.');
    end
end

%% ------------------------------------------------------------------ modèles

function model = build_model(kind, map, survey, par)
    R_obs = par.sigma_obs^2;

    switch kind
        case 'carte'
            model = struct('name', 'carte connue', 'kind', kind, 'n_train', 0, ...
                           'query', @(X, ctx) query_carte(map, R_obs, X(1:2, :)));

        case 'ok'
            % Tabulé sur le corridor du vol et non sur toute la carte.
            tab = tabulate_ok(map, survey, par);
            model = struct('name', 'OK carte totale', 'kind', kind, ...
                           'n_train', survey.n, ...
                           'query', @(X, ctx) query_table(tab, X(1:2, :), survey.n));

        case 'bilineaire'
            % La ligne de base sans krigeage. Elle déclare une seule variance, partout la
            % même, faute de savoir son erreur d'interpolation : celle-ci croît en
            % maille^2, donc le modèle devient sur-confiant quand la maille grossit.
            %
            % Deux lectures pour un même modèle. Sur un relevé au réseau exact, la lecture
            % bilinéaire de la grille. Sur un relevé secoué, les échantillons ne sont plus
            % sur le réseau nominal, et les y forcer ferait payer au bilinéaire une erreur
            % que le krigeage ne paie pas, puisque lui connaît les vraies positions : la
            % comparaison serait faussée en sa faveur. L'interpolation linéaire sur la
            % triangulation des points est alors la ligne de base honnête — la même
            % information, vue par un estimateur plus simple.
            R_bil = R_obs + par.krig.sigma_map^2;
            if isfield(survey, 'jitter') && survey.jitter > 0
                F    = scatteredInterpolant(survey.X(1, :)', survey.X(2, :)', ...
                                            survey.z, 'linear', 'nearest');
                lire = @(P) F(P(1, :)', P(2, :)');
            else
                n_s    = survey.n_side;
                grille = struct('h', reshape(survey.z, n_s, n_s), ...
                                'step', [survey.pitch; survey.X(:, 1)']);
                lire   = @(P) grid_read(grille.h, grille.step, P);
            end
            model = struct('name', 'bilinéaire', 'kind', kind, 'n_train', survey.n, ...
                           'query', @(X, ctx) query_bilineaire(lire, R_bil, ...
                                     X(1:2, :), survey.n));

        case "ok_exact"
            % Le même estimateur que "ok", évalué exactement a chaque pas au lieu
            % d'être tabulé. Sert a verifier que la tabulation est invisible ; hors de
            % cette verification il est inutilisable, une résolution contre les n~ points
            % du relevé par particule et par pas.
            model = struct("name", "OK exact", "kind", kind, "n_train", survey.n, ...
                           "query", @(X, ctx) query_ok_exact(survey, X(1:2, :), par));

        case "ak"
            model = struct("name", "OK adaptatif", 'kind', kind, 'n_train', survey.n, ...
                           'query', @(X, ctx) query_ak(survey, X(1:2, :), ctx, par));

        case 'ak_alea'
            % Contrôle négatif de la SÉLECTION : le relevé éclairci au hasard, figé pour
            % la mission. Il sépare ce que l'adaptativité doit au nombre de points de ce
            % qu'elle doit à leur place, et il montre qu'un sous-échantillonnage ne résout
            % que le coût du krigeage statique.
            %
            % ÉCLAIRCISSEMENT DE BERNOULLI, et non un effectif fixe. À probabilité p sur
            % une maille m, la densité moyenne devient celle d'une maille m/sqrt(p) : c'est
            % une ÉQUIVALENCE DE DENSITÉ avec un relevé plus grossier que la campagne A a
            % déjà mesuré, et donc une comparaison. Un effectif fixe n'a pas d'équivalent
            % de ce genre et ne compare rien d'une maille à l'autre.
            %
            % À densité moyenne égale, le réseau garantit une distance maximale de
            % m/sqrt(2) au plus proche échantillon. L'éclairci, lui, laisse des trous : un
            % bloc 3x3 vide arrive avec probabilité p^9, soit une douzaine de fois sur
            % 6241 points à p = 1/2, et son rayon atteint alors la longueur de corrélation.
            % C'est là que l'éclairci perd, et c'est ce qui sépare « moins de points » de
            % « les points au bon endroit ».
            %
            % TABULÉ comme le statique, et pour la même raison : le jeu d'entraînement est
            % fixe, donc l'estimateur est un couple fixe de fonctions de la position. Sans
            % tabulation, l'éclairci d'une maille fine demande 9000 Gflop par essai et le
            % lot ne tient plus en heures. Les deux côtés de chaque paire sont ainsi
            % traités de la même façon, ce que la comparaison exige.
            % LES HYPERPARAMÈTRES VIENNENT DU RELEVÉ COMPLET. C'est ce qui fait de ce
            % modèle un raccourci de coût pour le krigeage statique, et non un relevé plus
            % maigre : on possède tout, on ne peut pas se permettre de kriger dessus. Les
            % ajuster sur l'éclairci mélangerait deux effets — de moins bons
            % hyperparamètres, et une moins bonne prédiction — alors que l'argument ne
            % porte que sur le second.
            %
            % DEUX MODES, selon que le sous-ensemble est figé ou retiré à chaque pas.
            %
            %   figé      l'erreur de carte est une fonction fixe de la position, donc
            %             corrélée le long de la trajectoire : l'avion retrouve le même
            %             trou à chaque passage et le filtre converge vers la position que
            %             ce champ d'erreur lui dicte. Tabulable, donc gratuit par pas.
            %   par pas   l'erreur change d'un pas à l'autre et se moyenne sur la mission.
            %             Il ne reste que le rétrécissement vers la moyenne dans les trous,
            %             bien plus petit. Le retirage DÉCORRÈLE, il ne débiaise pas — la
            %             nuance compte, l'espérance sur les tirages restant rétrécie.
            %             Interdit la tabulation, donc N n^2 par pas.
            %
            % Comparer les deux attribue le dégât : si le retiré converge et le figé non, à
            % même p, le coupable est la corrélation et non la rareté des points.
            if isfield(par.krig, 'alea_par_pas') && par.krig.alea_par_pas
                fprintf('Sous-échantillonnage retiré à chaque pas, p = %.2f\n', par.krig.p_alea);
                model = struct('name', 'OK éclairci par pas', 'kind', kind, ...
                               'n_train', survey.n, ...
                               'query', @(X, ctx) query_alea_pas(survey, X(1:2, :), ctx, par));
            else
                rs        = RandStream('threefry', 'Seed', par.krig.survey_seed + 7);
                idx_a     = find(rand(rs, 1, survey.n) < par.krig.p_alea);
                srv_a     = survey;
                srv_a.X   = survey.X(:, idx_a);
                srv_a.z   = survey.z(idx_a);
                srv_a.n   = numel(idx_a);
                srv_a.p_alea = par.krig.p_alea;   % distingue sa tabulation de celle du plein
                fprintf(['Relevé éclairci à p = %.2f : %d points sur %d, densité équivalente ' ...
                         'à une maille de %.2f km\n'], par.krig.p_alea, srv_a.n, survey.n, ...
                        survey.pitch(1) / sqrt(par.krig.p_alea) / 1e3);
                tab_a = tabulate_ok(map, srv_a, par);
                model = struct('name', 'OK éclairci', 'kind', kind, 'n_train', srv_a.n, ...
                               'query', @(X, ctx) query_table(tab_a, X(1:2, :), srv_a.n));
            end

        case 'ak_boules'
            % Contrôle négatif du FENÊTRAGE ET DE LA CLUSTERISATION : la réunion des
            % sélections faites autour de chaque particule, sans plancher et sans
            % partition. C'est la sélection « évidente », celle qu'un lecteur proposerait,
            % et elle échoue par les deux bouts.
            %
            % L'échec qui porte l'argument est le premier : un SEUL ajustement
            % d'hyperparamètres et une SEULE moyenne estimés sur une réunion qui enjambe
            % les modes. Le krigeage ordinaire estime une moyenne ; si la réunion couvre
            % deux modes distants de dizaines de kilomètres, cette moyenne unique est
            % fausse pour les deux. C'est exactement l'incohérence que le CA-OK supprime en
            % ajustant par mode, et elle ne se voit qu'en régime multimodal.
            %
            % Le second échec motive n_min : quand le nuage se resserre sous une maille,
            % toutes les particules partagent les mêmes voisins et la réunion tombe à un ou
            % deux points.
            %
            % L'arbre est construit une fois ; la recherche de voisinage reste un coût par
            % pas que l'AK ne paie pas, et elle est comptée.
            arbre = KDTreeSearcher(survey.X');
            model = struct('name', 'OK par réunion', 'kind', kind, 'n_train', survey.n, ...
                           'query', @(X, ctx) query_boules(survey, arbre, X(1:2, :), ...
                                     ctx, par));

        case 'cak'
            model = struct('name', 'CA-OK', 'kind', kind, 'n_train', survey.n, ...
                           'query', @(X, ctx) query_cak(survey, X(1:2, :), ctx, par));

        otherwise
            error('krigeage:modele', 'Modèle inconnu %s.', kind);
    end
end

function [z, R, n, n_modes, cout] = query_carte(map, R_obs, P)
% Le champ vrai, variance du capteur seule. Un mode par construction.
    z = grid_read(map.h, map.step, P);
    R = R_obs * ones(size(P, 2), 1);
    n = 0;  n_modes = 1;  cout = [0; 0; 0; 0];
end

function [z, R, n, n_modes, cout] = query_table(tab, P, n_train)
% L'estimateur tabulé, relu bilinéairement. La variance est plancherée : l'arrondi peut
% rendre négative l'interpolation d'une quantité qui ne l'est pas, et une vraisemblance
% divisée par une variance n'est pas l'endroit pour le découvrir.
    z = grid_read(tab.z, tab.step, P);
    R = max(grid_read(tab.R, tab.step, P), eps);
    n = n_train;  n_modes = 1;  cout = [0; 0; 0; 0];   % tout le coût est dans la tabulation
end

function [z, R, n, n_modes, cout] = query_bilineaire(lire, R_fixe, P, n_train)
% Le relevé bruité lu par interpolation linéaire, variance constante. Aucun calcul par pas
% au-delà de la lecture, et aucune information sur sa propre erreur.
    z = lire(P);
    R = R_fixe * ones(size(P, 2), 1);
    n = n_train;  n_modes = 1;  cout = [0; 0; 0; 0];
end

function [z, R, n, n_modes, cout] = query_ok_exact(survey, P, par)
% Tout le relevé, hyperparamètres globaux, aucune tabulation.
    [z, R] = ok_solve(survey.X, survey.z, P, survey.hyp, par);
    n = survey.n;  n_modes = 1;  cout = [n; n^3; 0; 0];
end

function [z, R, n, n_modes, cout] = query_ak(survey, P, ctx, par)
% Une fenêtre : sélectionner, réajuster, kriger. C'est le selection_fenetrage_adaptatif
% de 2A, écrit sur le contexte que le filtre calcule déjà.
    idx = window_select(survey.X, ctx, par);
    hyp = window_hyperparameters(survey, idx, par);
    [z, R] = ok_solve(survey.X(:, idx), survey.z(idx), P, hyp, par);
    n = numel(idx);  n_modes = 1;  cout = [n; n^3; 1; 0];
end
function [z, R, n, n_modes, cout] = query_alea_pas(survey, P, ctx, par)
% Un sous-ensemble de Bernoulli retiré à chaque pas, krigé avec les hyperparamètres du
% relevé complet.
%
% LE FLUX DÉPEND DU PAS ET DU TIRAGE DE RELEVÉ, pas de l'essai : les sous-ensembles sont
% une propriété de l'algorithme et non du hasard de la mission, et les cent tirages de
% relevé donnent déjà cent séquences différentes. Reproductible, et indépendant des flux
% du filtre.
    rs  = RandStream('threefry', 'Seed', par.krig.survey_seed + 7 + ctx.k);
    idx = find(rand(rs, 1, survey.n) < par.krig.p_alea);
    n   = numel(idx);
    n_modes = 1;
    if n < 2
        z    = mean(survey.z) * ones(size(P, 2), 1);
        R    = (survey.hyp.sigma_f^2 + par.krig.sigma_map^2 + par.sigma_obs^2) ...
               * ones(size(P, 2), 1);
        cout = [n; 0; 0; 0];
        return
    end
    [z, R] = ok_solve(survey.X(:, idx), survey.z(idx), P, survey.hyp, par);
    cout   = [n; n^3; 0; 0];      % zéro ajustement : les hyperparamètres sont ceux du plein
end

function [z, R, n, n_modes, cout] = query_boules(survey, arbre, P, ctx, par)
% La réunion des sélections faites autour de chaque particule, sans plancher ni partition.
% Le cas dégénéré — moins de deux points retenus — rend la loi a priori du relevé plutôt
% que de lever : l'échec doit être mesurable, pas fatal, sans quoi la campagne s'arrête au
% lieu de documenter l'effondrement.
    % DEUX SÉLECTIONS SOUS UN SEUL PARAMÈTRE. Un rayon positif donne la réunion des
    % boules ; un rayon nul donne les k_voisins plus proches échantillons de chaque
    % particule. Le second est préféré : un rayon fixe attrape trois points par particule
    % sur une maille de 3,5 km et moins d'un sur une maille de 7,5, si bien que
    % l'effondrement mesuré serait celui du rayon et non celui de la méthode.
    %
    % QUATRE VOISINS, ET NON UN. Sur un réseau, les quatre plus proches échantillons d'un
    % point sont les coins de sa cellule : c'est exactement le pochoir du bilinéaire, donc
    % la sélection qu'un lecteur proposerait, généralisée au krigeage. À un seul voisin la
    % méthode devient caricaturale — la réunion tombe à un ou deux points dès que le nuage
    % se resserre — et on lui ferait un procès trop facile.
    %
    % Elle reste néanmoins maigre en poursuite, toutes les particules tenant alors dans
    % quelques cellules : c'est l'effondrement que n_min corrige, et il reste mesurable
    % sans être fabriqué.
    if par.krig.rayon_boule > 0
        voisins = rangesearch(arbre, ctx.X(1:2, :)', par.krig.rayon_boule);
        idx     = unique([voisins{:}]);
    else
        idx     = unique(knnsearch(arbre, ctx.X(1:2, :)', 'K', par.krig.k_voisins))';
    end
    n       = numel(idx);
    n_modes = 1;
    n_dist  = size(ctx.X, 2) * log2(max(survey.n, 2));   % la recherche, en distances

    if n < 2
        z    = mean(survey.z) * ones(size(P, 2), 1);
        R    = (survey.hyp.sigma_f^2 + par.krig.sigma_map^2 + par.sigma_obs^2) ...
               * ones(size(P, 2), 1);
        cout = [n; 0; 0; n_dist];
        return
    end

    hyp    = window_hyperparameters(survey, idx, par);
    [z, R] = ok_solve(survey.X(:, idx), survey.z(idx), P, hyp, par);
    cout   = [n; n^3; 1; n_dist];
end

function [z, R, n, n_modes, cout] = query_cak(survey, P, ctx, par)
% Une fenêtre par mode du nuage. Les points de requête sont affectés au centre le plus
% proche plutôt que d'hériter des étiquettes des particules : l'APF interroge des points
% qui ne sont pas ceux dont le contexte est issu, et un centre est un lieu auquel les deux
% ensembles peuvent être rattachés.
    [centres, n_dist] = mean_shift(ctx.X(1:2, :), ctx.w, par.krig.bandwidth);
    n_modes = size(centres, 2);
    if n_modes == 1
        % Le mean-shift a tourne pour rien, et il faut quand meme le facturer : c'est le
        % prix de la clusterisation dans la phase de poursuite, ou le nuage n'a qu'un mode.
        [z, R, n, ~, cout] = query_ak(survey, P, ctx, par);
        cout(4) = n_dist;
        return
    end

    label_p = nearest_centre(ctx.X(1:2, :), centres);
    label_q = nearest_centre(P, centres);
    z = zeros(size(P, 2), 1);
    R = zeros(size(P, 2), 1);
    used = false(1, size(survey.X, 2));
    cout = [0; 0; 0; n_dist];

    for c = 1:size(centres, 2)
        rows = label_q == c;
        if ~any(rows)
            continue
        end
        ctx_c = cluster_context(ctx, label_p == c, centres(:, c));
        idx   = window_select(survey.X, ctx_c, par);
        used(idx) = true;
        cout  = cout + [numel(idx); numel(idx)^3; 1; 0];
        [z(rows), R(rows)] = ok_solve(survey.X(:, idx), survey.z(idx), P(:, rows), ...
                                      window_hyperparameters(survey, idx, par), par);
    end
    n = nnz(used);            % the union, which is what the step actually touched
end

%% ------------------------------- fenêtre, clusters, équations du krigeage

function idx = window_select(X, ctx, par)
% Un disque centré sur la position prédite, de rayon alpha_cov fois le demi-grand axe de
% l'ellipse à 3 sigma, contenant entre n_min et n_max échantillons.
%
% Le plancher est le minimum_echantillons de 2A : il empêche un filtre devenu
% sur-confiant de kriger sur deux points. Le plafond n'y est pas, et il le faut : le
% rayon est proportionnel à l'incertitude déclarée, donc un nuage encore large de 15 km
% réclame un disque de 90 km, c'est-à-dire tout le relevé — l'AK devient alors un
% krigeage ordinaire avec une factorisation de 6400 points par pas et par cluster.
% Plafonner en nombre et non en rayon garde la sélection valide que les échantillons
% soient denses ou clairsemés.
    centre = ctx.x_pred(1:2);
    P_pos  = ctx.P_pred(1:2, 1:2);
    radius = par.krig.window_factor * 3 * sqrt(max(eig((P_pos + P_pos') / 2)));

    dist   = vecnorm(X - centre, 2, 1);
    n      = numel(dist);
    sorted = sort(dist);
    radius = max(radius, sorted(min(par.krig.n_min, n)));
    radius = min(radius, sorted(min(par.krig.n_max, n)));
    idx    = find(dist <= radius);
end

function hyp = window_hyperparameters(survey, idx, par)
% 'fenetre' les réestime sur la fenêtre, ce qui rend la méthode adaptative et non
% seulement moins chère : une fenêtre est assez petite pour être presque stationnaire, là
% où un noyau global est un compromis entre les cuvettes lisses et les crêtes rugueuses
% d'un champ réel. 'global' garde ceux du relevé et isole ce que la sélection seule vaut.
%
% n_ml_win plafonne l'échantillon de l'ajustement, pas celui du krigeage : ell et sigma_f
% se lisent aussi bien sur une partie de la fenêtre, et le coût est en n^3 par itération
% interne, donc 2000 points coûtent 64 fois 500.
    switch par.krig.refit
        case 'global'
            hyp = survey.hyp;
        case 'fenetre'
            if isfield(par.krig, 'n_ml_win') && numel(idx) > par.krig.n_ml_win
                sub = idx(randperm(numel(idx), par.krig.n_ml_win));
            else
                sub = idx;
            end
            hyp = fit_hyperparameters(survey.X(:, sub)', survey.z(sub), par);
        otherwise
            error('krigeage:refit', 'par.krig.refit inconnu : %s.', par.krig.refit);
    end
end

function [z_hat, R_hat] = ok_solve(X, z, P, hyp, par)
% Krigeage ordinaire en un appel, sans former l'inverse. R_hat porte aussi le bruit du
% capteur : le filtre le lit comme la variance de l'observation sachant le relevé, et ne
% doit donc pas l'ajouter une seconde fois.
    n = size(X, 2);
    C = ones(n, 1);

    G      = se_kernel(X, X, hyp) + par.krig.sigma_map^2 * eye(n);
    [L, p] = chol(G, 'lower');
    if p > 0
        L = chol(G + 1e-8 * trace(G) / n * eye(n), 'lower');
    end

    sol  = L' \ (L \ [z, C]);
    Minv = C' * sol(:, 2);
    m_ok = Minv \ (C' * sol(:, 1));                  % the estimated mean
    a    = sol(:, 1) - sol(:, 2) * m_ok;

    % Les requêtes par blocs : K est N x n et se_kernel en tient cinq copies, soit plus
    % d'un demi-gigaoctet à n = 2000 et N = 5000. Ses lignes sont indépendantes, donc le
    % découpage ne change pas un bit du résultat.
    N     = size(P, 2);
    z_hat = zeros(N, 1);
    R_hat = zeros(N, 1);
    bloc  = max(1, floor(5e5 / max(n, 1)));   % 5e5 et non 2e6 depuis le 7 septembre : trois
    % plantages de suite ont eu lieu dans se_kernel, faute de memoire, sur une machine ou
    % VS Code tenait cinq gigaoctets sur seize. A 2e6 le noyau demandait un bloc contigu de
    % 16 Mo et ses cinq copies 80 ; a 5e5 c'est 4 et 20. Le resultat est inchange — les
    % lignes de K sont independantes — et la boucle tourne quatre fois plus, ce qui ne se
    % voit pas : le cout est dans l'algebre, pas dans l'iteration.
    for i = 1:bloc:N
        j = min(i + bloc - 1, N);
        K = se_kernel(P(:, i:j), X, hyp);
        z_hat(i:j) = m_ok + K * a;

        V          = L \ K';
        R_sk       = hyp.sigma_f^2 + par.sigma_obs^2 - sum(V.^2, 1)';
        U          = 1 - K * sol(:, 2);              % the cost of estimating the mean
        R_hat(i:j) = max(U.^2 / Minv + R_sk, eps);
    end
end

function K = se_kernel(A, B, hyp)
% Le noyau, un point par colonne. hyp.nom choisit la famille, et les formes sont celles
% de fitrgp — mêmes conventions de ell, donc un ell ajusté par fitrgp s'emploie ici tel
% quel. 1A et 2A omettent le facteur deux de l'exponentielle quadratique tout en lisant
% ell dans fitrgp : leur noyau est plus étroit que celui qu'ils ont ajusté, d'un sqrt(2).
%
% La régularité du noyau à l'origine décide de ce que la reconstruction peut représenter.
% L'exponentielle quadratique est infiniment dérivable, donc incapable de porter la
% structure courte échelle d'un champ d'anomalie, dont le spectre suit une loi de
% puissance. Les Matérn de nu demi-entier sont dérivables nu - 1/2 fois seulement, et
% l'exponentielle pas du tout.
    D2 = max(sum(A.^2, 1)' - 2 * (A' * B) + sum(B.^2, 1), 0);
    if isfield(hyp, 'nom'), nom = hyp.nom; else, nom = 'squaredexponential'; end

    switch nom
        case 'squaredexponential'
            K = exp(-D2 / (2 * hyp.ell^2));
        case 'exponential'                                    % Matérn nu = 1/2
            K = exp(-sqrt(D2) / hyp.ell);
        case 'matern32'
            r = sqrt(3 * D2) / hyp.ell;
            K = (1 + r) .* exp(-r);
        case 'matern52'
            r = sqrt(5 * D2) / hyp.ell;
            K = (1 + r + r.^2 / 3) .* exp(-r);
        otherwise
            error('krigeage:noyau', 'Noyau inconnu : %s.', nom);
    end
    K = hyp.sigma_f^2 * K;
end

function ctx_c = cluster_context(ctx, members, centre)
% Moyenne et covariance prédites d'un mode, sur les particules qu'il a gardées. Un mode
% qui n'en a gardé aucune peut rester le plus proche d'un point de requête, et son centre
% est alors tout ce dont il dispose.
    ctx_c = ctx;
    if ~any(members)
        ctx_c.x_pred(1:2) = centre;
        return
    end
    w = ctx.w(members);  w = w / sum(w);
    X = ctx.X(:, members);
    ctx_c.X      = X;
    ctx_c.w      = w;
    ctx_c.x_pred = X * w';
    Xc           = X - ctx_c.x_pred;
    ctx_c.P_pred = (Xc .* w) * Xc';
end

function [centres, n_dist] = mean_shift(X, w, bandwidth)
% Mean shift à noyau plat sur le nuage pondéré : chaque germe marche vers le mode des
% points dans sa bande, et deux modes distants de moins d'une demi-bande n en font qu'un.
% Germes tirés d'un sous-échantillon : les modes de dix mille particules sont ceux de deux
% cents d'entre elles, et c'est la marche qui coûte.
    % La marche porte sur un nuage éclairci : les modes de dix mille particules sont ceux
    % de cinq cents d'elles, et chaque germe balaie tous les points à chaque itération —
    % c'est là tout le coût de la méthode.
    if size(X, 2) > 500
        keep = round(linspace(1, size(X, 2), 500));
        X = X(:, keep);  w = w(keep);
    end
    n_dist = 0;                  % evaluations de distance, pour le compte de flops
    n_seed = min(size(X, 2), 60);
    seeds  = X(:, round(linspace(1, size(X, 2), n_seed)));
    modes  = zeros(size(seeds));

    for s = 1:n_seed
        m = seeds(:, s);
        for it = 1:100
            in = sum((X - m).^2, 1) < bandwidth^2;
            n_dist = n_dist + size(X, 2);
            if ~any(in)
                break
            end
            wi = w(in);
            m_new = (X(:, in) * wi') / sum(wi);
            if norm(m_new - m) < 1e-3 * bandwidth
                m = m_new;
                break
            end
            m = m_new;
        end
        modes(:, s) = m;
    end

    centres = modes(:, 1);
    for s = 2:n_seed
        if all(vecnorm(centres - modes(:, s), 2, 1) > bandwidth / 2)
            centres(:, end + 1) = modes(:, s);  %#ok<AGROW>
        end
    end
end

function label = nearest_centre(X, centres)
    d2 = zeros(size(centres, 2), size(X, 2));
    for c = 1:size(centres, 2)
        d2(c, :) = sum((X - centres(:, c)).^2, 1);
    end
    [~, label] = min(d2, [], 1);
end

%% ------------------------------------------------------------------ tabulation

function tab = tabulate_ok(map, survey, par)
% L'estimateur statique, évalué une fois sur le corridor du vol, puis mis en cache.
%
% C'est la tabulation qui rend le modèle statique abordable : le relevé et les
% hyperparamètres étant fixes, le krigeage ordinaire est un couple fixe de fonctions de la
% position, et l'évaluer en une particule revient à les lire. Elle coûte une interpolation
% de l'estimateur, dont l'erreur est mesurée plus bas plutôt que supposée.
    box  = par.krig.corridor;                     % [x_lo x_hi y_lo y_hi], metres
    step = par.krig.tab_pitch;
    ax   = (box(1):step:box(2))';
    ay   = (box(3):step:box(4))';
    [TX, TY] = meshgrid(ax, ay);
    Q    = [TX(:), TY(:)]';

    key  = struct('map', par.map.name, 'n_side', survey.n_side, ...
                  'sigma_map', par.krig.sigma_map, 'sigma_obs', par.sigma_obs, ...
                  'ell_init', par.krig.ell_init, 'n_ml', par.krig.n_ml, ...
                  "seed", graine_releve(par), "box", box, "step", step, ...
                  "noyau", noyau(par), "phase", survey.phase, ...
                  "jitter", survey.jitter, "pas", survey.pitch);
    % La maille est dans la clé ET dans le nom : deux relevés de même nombre de points
    % peuvent avoir des mailles différentes dès que le pas est le paramètre premier, et
    % sans elle le second relit la tabulation du premier.
    suff = "";
    if ~strcmp(noyau(par), "squaredexponential"), suff = "_" + string(noyau(par)); end
    if isfield(par.krig, 'phase_alea') && par.krig.phase_alea, suff = suff + "_ph"; end
    if secouage(par) > 0, suff = suff + sprintf("_j%d", round(100 * secouage(par))); end
    file = fullfile(par.path.cache, sprintf("ok_%dpts_%dm_p%d_g%d%s.mat", survey.n, ...
                    round(step), round(survey.pitch(1)), graine_releve(par), suff));

    % UN RELEVÉ ÉCLAIRCI NE PARTAGE PAS LA TABULATION DU PLEIN. Le champ n'existe que
    % dans ce cas, donc les six cents tabulations déjà en cache gardent leur clé et
    % restent valides — ajouter un champ inconditionnellement les aurait toutes invalidées.
    if isfield(survey, 'p_alea')
        key.p_alea = survey.p_alea;
        suff = suff + sprintf("_b%d", round(100 * survey.p_alea));
        file = fullfile(par.path.cache, sprintf("ok_%dpts_%dm_p%d_g%d%s.mat", survey.n, ...
                        round(step), round(survey.pitch(1)), graine_releve(par), suff));
    end

    if exist(file, 'file')
        c = load(file);
        if isequal(c.key, key)
            fprintf('Tabulation relue depuis %s\n', file);
            tab = c.tab;
            return
        end
    end

    fprintf('Tabulation du krigeage sur %d points, %.1f Gflop...\n', ...
            numel(TX), 2 * numel(TX) * survey.n^2 / 1e9);
    t = tic;
    [z, R] = ok_solve_chunked(survey, Q, survey.hyp, par);
    tab = struct('z', reshape(z, numel(ay), numel(ax)), ...
                 'R', reshape(R, numel(ay), numel(ax)), ...
                 'step', [step, step; box(1), box(3)]);
    fprintf('   %.0f s, %.1f Gflop/s\n', toc(t), ...
            2 * numel(TX) * survey.n^2 / toc(t) / 1e9);

    % Ce que la tabulation coûte, mesuré hors de sa propre grille.
    chk = [box(1) + (box(2) - box(1)) * rand(1, 500)
           box(3) + (box(4) - box(3)) * rand(1, 500)];
    [ze, Re] = ok_solve_chunked(survey, chk, survey.hyp, par);
    zi = grid_read(tab.z, tab.step, chk);
    fprintf('   erreur de tabulation : %.2f %s au pire, soit %.1f %% d''un sigma\n', ...
            max(abs(ze - zi)), par.unit, 100 * max(abs(ze - zi)) / mean(sqrt(Re)));

    % Ce que la reconstruction vaut vraiment, contre ce qu'elle annonce. Le rapport est
    % un ANEES du modèle d'observation seul, avant qu'aucun filtre ne le lise : au-dessus
    % de un le krigeage est prudent, en dessous il sous-déclare, et c'est la première
    % chose à regarder quand un modèle statique se met à diverger.
    [GX, GY] = meshgrid(ax, ay);
    vrai = grid_read(map.h, map.step, [GX(:), GY(:)]');
    err  = sqrt(mean((tab.z(:) - vrai).^2));
    decl = sqrt(max(mean(tab.R(:)) - par.sigma_obs^2, 0));
    fprintf('   erreur de carte %.1f %s, incertitude déclarée %.1f %s, rapport %.2f\n', ...
            err, par.unit, decl, par.unit, decl / max(err, eps));
    tab.err_carte = err;  tab.sigma_declare = decl;

    if ~exist(par.path.cache, 'dir')
        mkdir(par.path.cache);
    end
    save(file, 'key', 'tab');
end

function [z, R] = ok_solve_chunked(survey, Q, hyp, par)
% Les équations de ok_solve, les requêtes prises par blocs pour ne jamais tenir k(x, x~)
% en entier, et l'inverse formé une fois puisqu'il sert à tous les blocs — ce qui change
% une résolution triangulaire par bloc en un produit matriciel, que le BLAS exécute à une
% bien meilleure fraction du pic.
    X = survey.X;  n = size(X, 2);
    C = ones(n, 1);
    G = se_kernel(X, X, hyp) + par.krig.sigma_map^2 * eye(n);
    L = chol(G, 'lower');
    Ginv = L' \ (L \ eye(n));  Ginv = (Ginv + Ginv') / 2;
    clear G L

    Ginv_C = Ginv * C;
    Minv   = C' * Ginv_C;
    m_ok   = Minv \ (C' * (Ginv * survey.z));
    Ginv_r = Ginv * (survey.z - C * m_ok);

    z = zeros(size(Q, 2), 1);
    R = zeros(size(Q, 2), 1);
    for b = 1:par.krig.n_chunk:size(Q, 2)
        j  = b:min(b + par.krig.n_chunk - 1, size(Q, 2));
        K  = se_kernel(Q(:, j), X, hyp);
        z(j) = m_ok + K * Ginv_r;
        KG   = K * Ginv;
        U    = 1 - K * Ginv_C;
        R(j) = max(U.^2 / Minv + hyp.sigma_f^2 + par.sigma_obs^2 - sum(KG .* K, 2), eps);
    end
end

%% ----------------------------------------------------------- lecture de grille

function v = grid_read(g, step, P)
% GRID_READ  Lecture bilinéaire d'un champ tabulé sur une grille régulière.
%
% Le nœud (ligne, colonne) porte la valeur en (origine_x + (col-1) pas_x, origine_y +
% (lig-1) pas_y) : les lignes indexent y et les colonnes x, ce que produit meshgrid. P
% donne un point de requête par colonne.
%
% STEP est soit une ligne [pas_x, pas_y], l'origine valant alors un pas — convention des
% sauvegardes — soit les deux lignes [pas_x, pas_y; origine_x, origine_y] pour une grille
% qui commence ailleurs, ce dont une tabulation sur corridor a besoin.
%
% LES DEUX AXES SONT TRAITÉS PAREIL. Le bilineaire.m de 1A ne l'est pas : il indexe les
% colonnes depuis zéro et les lignes depuis un. L'asymétrie est invisible tant qu'une
% seule grille est lue de la même façon partout, et décale un estimateur tabulé d'une
% maille contre le champ dont il est issu dès qu'il y en a deux — ici 96 nT d erreur, onze
% fois l'écart-type du krigeage lui-même.
%
% Les requêtes hors grille sont ramenées au bord plutôt que rejetées : une particule
% partie hors carte lit la valeur tabulée la plus proche, son poids tranche, et aucun
% appelant n'a de cas particulier à porter.
    if size(step, 1) > 1
        origin = step(2, :);
        step   = step(1, :);
    else
        origin = step;
    end

    lam = (P(1, :) - origin(1)) / step(1);
    ix  = floor(lam);
    lam = lam - ix;

    mu = (P(2, :) - origin(2)) / step(2);
    iy = floor(mu);
    mu = mu - iy;

    ny = size(g, 1);
    ix = min(max(ix, 0), size(g, 2) - 2);
    iy = min(max(iy, 0), ny - 2);

    i = ny * ix + iy + 1;
    v = ((1 - mu) .* (1 - lam) .* g(i)          + ...
         (1 - mu) .*      lam  .* g(i + ny)     + ...
              mu  .* (1 - lam) .* g(i + 1)      + ...
              mu  .*      lam  .* g(i + ny + 1))';
end
