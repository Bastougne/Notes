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
    n_s   = par.krig.n_side;
    lo    = map.step;                                   % un pas depuis le bord
    hi    = ([size(map.h, 2), size(map.h, 1)] - 2) .* map.step;
    ax    = linspace(lo(1), hi(1), n_s)';
    ay    = linspace(lo(2), hi(2), n_s)';
    [SX, SY] = meshgrid(ax, ay);

    rng(par.run.seed);
    survey.X     = [SX(:), SY(:)]';                     % 2 x n, un échantillon par colonne
    survey.z     = grid_read(map.h, map.step, survey.X) ...
                 + par.krig.sigma_map * randn(n_s^2, 1);
    survey.n     = n_s^2;
    survey.pitch = [ax(2) - ax(1), ay(2) - ay(1)];

    % Sur un sous-échantillon, le bruit tenu à sa valeur déclarée plutôt qu'ajusté :
    % laissé libre, fitrgp explique toute la rugosité du champ comme du bruit et rend un
    % ell très en dessous de la maille. La vraisemblance marginale a de plus un second
    % optimum en ell -> 0, où les optimiseurs tombent depuis leurs défauts — d'où le point
    % de départ imposé.
    i_ml = randperm(survey.n, min(par.krig.n_ml, survey.n));
    survey.hyp = fit_hyperparameters(survey.X(:, i_ml)', survey.z(i_ml), par);

    fprintf('Relevé %d x %d = %d points, maille %.2f x %.2f km, bruit %.0f %s\n', ...
            n_s, n_s, survey.n, survey.pitch / 1e3, par.krig.sigma_map, par.unit);
    fprintf('   sigma_f = %.0f %s, ell = %.2f km\n', ...
            survey.hyp.sigma_f, par.unit, survey.hyp.ell / 1e3);
end

function hyp = fit_hyperparameters(X, z, par)
    if ~isempty(which('fitrgp'))
        gp = fitrgp(X, z, 'KernelFunction', 'squaredexponential', ...
                    'BasisFunction', 'constant', 'Sigma', par.krig.sigma_map, ...
                    'ConstantSigma', true, 'Standardize', false, ...
                    'KernelParameters', [par.krig.ell_init; std(z)]);
        hyp.ell     = gp.KernelInformation.KernelParameters(1);
        hyp.sigma_f = gp.KernelInformation.KernelParameters(2);
    else
        hyp.ell     = par.krig.ell_init;      % pas de toolbox : la valeur initiale, et on le dit
        hyp.sigma_f = std(z);
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
            n_s  = par.krig.n_side;
            grille = struct('h', reshape(survey.z, n_s, n_s), ...
                            'step', [survey.pitch; survey.X(:, 1)']);
            model = struct('name', 'bilinéaire', 'kind', kind, 'n_train', survey.n, ...
                           'query', @(X, ctx) query_bilineaire(grille, ...
                                     R_obs + par.krig.sigma_map^2, X(1:2, :), survey.n));

        case 'ak'
            model = struct('name', 'OK adaptatif', 'kind', kind, 'n_train', survey.n, ...
                           'query', @(X, ctx) query_ak(survey, X(1:2, :), ctx, par));

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
    n = 0;  n_modes = 1;  cout = [0; 0; 0];
end

function [z, R, n, n_modes, cout] = query_table(tab, P, n_train)
% L'estimateur tabulé, relu bilinéairement. La variance est plancherée : l'arrondi peut
% rendre négative l'interpolation d'une quantité qui ne l'est pas, et une vraisemblance
% divisée par une variance n'est pas l'endroit pour le découvrir.
    z = grid_read(tab.z, tab.step, P);
    R = max(grid_read(tab.R, tab.step, P), eps);
    n = n_train;  n_modes = 1;  cout = [0; 0; 0];   % tout le coût est dans la tabulation
end

function [z, R, n, n_modes, cout] = query_bilineaire(grille, R_fixe, P, n_train)
% Le relevé bruité lu bilinéairement, variance constante. Aucun calcul par pas au-delà
% de la lecture, et aucune information sur sa propre erreur.
    z = grid_read(grille.h, grille.step, P);
    R = R_fixe * ones(size(P, 2), 1);
    n = n_train;  n_modes = 1;  cout = [0; 0; 0];
end

function [z, R, n, n_modes, cout] = query_ak(survey, P, ctx, par)
% Une fenêtre : sélectionner, réajuster, kriger. C'est le selection_fenetrage_adaptatif
% de 2A, écrit sur le contexte que le filtre calcule déjà.
    idx = window_select(survey.X, ctx, par);
    hyp = window_hyperparameters(survey, idx, par);
    [z, R] = ok_solve(survey.X(:, idx), survey.z(idx), P, hyp, par);
    n = numel(idx);  n_modes = 1;  cout = [n; n^3; 1];
end

function [z, R, n, n_modes, cout] = query_cak(survey, P, ctx, par)
% Une fenêtre par mode du nuage. Les points de requête sont affectés au centre le plus
% proche plutôt que d'hériter des étiquettes des particules : l'APF interroge des points
% qui ne sont pas ceux dont le contexte est issu, et un centre est un lieu auquel les deux
% ensembles peuvent être rattachés.
    centres = mean_shift(ctx.X(1:2, :), ctx.w, par.krig.bandwidth);
    n_modes = size(centres, 2);
    if n_modes == 1
        [z, R, n, ~, cout] = query_ak(survey, P, ctx, par);
        return
    end

    label_p = nearest_centre(ctx.X(1:2, :), centres);
    label_q = nearest_centre(P, centres);
    z = zeros(size(P, 2), 1);
    R = zeros(size(P, 2), 1);
    used = false(1, size(survey.X, 2));
    cout = [0; 0; 0];

    for c = 1:size(centres, 2)
        rows = label_q == c;
        if ~any(rows)
            continue
        end
        ctx_c = cluster_context(ctx, label_p == c, centres(:, c));
        idx   = window_select(survey.X, ctx_c, par);
        used(idx) = true;
        cout  = cout + [numel(idx); numel(idx)^3; 1];
        [z(rows), R(rows)] = ok_solve(survey.X(:, idx), survey.z(idx), P(:, rows), ...
                                      window_hyperparameters(survey, idx, par), par);
    end
    n = nnz(used);            % the union, which is what the step actually touched
end

%% ------------------------------- fenêtre, clusters, équations du krigeage

function idx = window_select(X, ctx, par)
% Un disque centré sur la position prédite, de rayon alpha_W fois le demi-grand axe de
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
    bloc  = max(1, floor(2e6 / max(n, 1)));
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
% Exponentiel quadratique, un point par colonne. Le facteur deux de l'exposant est celui
% de fitrgp ; 1A et 2A l'omettent tout en lisant ell dans fitrgp, donc le noyau avec
% lequel ils krigent est plus étroit que celui qu'ils ont ajusté, d'un facteur sqrt(2).
    D2 = sum(A.^2, 1)' - 2 * (A' * B) + sum(B.^2, 1);
    K  = hyp.sigma_f^2 * exp(-max(D2, 0) / (2 * hyp.ell^2));
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

function centres = mean_shift(X, w, bandwidth)
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
    n_seed = min(size(X, 2), 60);
    seeds  = X(:, round(linspace(1, size(X, 2), n_seed)));
    modes  = zeros(size(seeds));

    for s = 1:n_seed
        m = seeds(:, s);
        for it = 1:100
            in = sum((X - m).^2, 1) < bandwidth^2;
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

    key  = struct('map', par.map.name, 'n_side', par.krig.n_side, ...
                  'sigma_map', par.krig.sigma_map, 'sigma_obs', par.sigma_obs, ...
                  'ell_init', par.krig.ell_init, 'n_ml', par.krig.n_ml, ...
                  'seed', par.run.seed, 'box', box, 'step', step);
    file = fullfile(par.path.cache, sprintf('ok_%dpts_%dm.mat', survey.n, round(step)));

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
