function sl = slam()
% SLAM  Les deux algorithmes du chapitre 5, qui estiment la carte en même temps que l'état.
%
%   sl   = slam();
%   base = sl.base(kind, par, domaine);
%   out  = sl.run(kind, base, z_obs, scen, opt);
%
% KIND vaut 'rbpf' ou 'akkf'. OUT a les mêmes champs que celui de navigation.m — erreur,
% NEES, estimée, coût — plus `h_eval`, la carte reconstruite aux points opt.P_eval, qui
% est le second objet estimé et la moitié de ce qu'une campagne du chapitre 5 mesure.
%
% CE QUI SÉPARE CE FICHIER DE navigation.m. Là, la carte est une donnée : le filtre
% interroge un modèle d'observation construit une fois pour toutes sur un relevé, et rien
% de ce qu'il mesure ne revient dans la carte. Ici il n'y a pas de relevé, la carte fait
% partie de l'inconnue, et les deux algorithmes ne se distinguent que par la façon dont
% ils la portent : une par particule pour le RBPF, une seule pour l'AKKF.
%
% LA MÊME REPRÉSENTATION POUR LES DEUX. Le champ est écrit h(x) = Psi(x)' alpha sur une
% base finie, et seules la base et l'estimation de alpha changent :
%
%   'laplace'  les fonctions propres du laplacien de Dirichlet sur un rectangle, pondérées
%              par la densité spectrale du noyau. Grille de fréquences déterministe, base
%              nulle au bord du domaine.
%   'rff'      des cosinus de fréquences tirées selon cette même densité spectrale. Rien
%              ne s'annule nulle part et aucun domaine n'est déclaré, au prix d'une erreur
%              de Monte Carlo là où la grille a une erreur de troncature.
%
% LE CHAMP EST SCALAIRE (d_z = 1), ce qui retire du RBPF tout ce que le cas magnétique
% vectoriel impose : pas de quaternion, pas de rotation du repère capteur vers le repère
% carte, pas de gradient de potentiel. La base est scalaire, l'innovation est un nombre,
% et la covariance d'innovation aussi — ce qui permet de vectoriser la correction de
% Kalman des N cartes en divisions terme à terme.
    sl = struct('base', @build_base, 'run', @run_slam, 'densite', @densite_spectrale);
end

%% ---------------------------------------------------------------------- bases

function base = build_base(kind, par, domaine)
% BUILD_BASE  La base finie sur laquelle la carte est écrite.
%
%   base.eval(P)  r x M, les fonctions de base aux points P (2 x M)
%   base.Lambda   r x 1, la variance a priori de chaque coefficient
%   base.r        le rang
%
% DOMAINE vaut [x_min, x_max, y_min, y_max] et n'est lu que par la base de Laplace.
    switch kind
        case 'laplace'
            base = base_laplace(par, domaine);
        case 'rff'
            base = base_rff(par);
        otherwise
            error('slam:base', 'Base inconnue : %s.', kind);
    end
end

function base = base_laplace(par, domaine)
% Les fonctions propres du laplacien de Dirichlet sur le rectangle `domaine`, indexées par
% un multi-indice j et pondérées par la densité spectrale du noyau évaluée aux fréquences
% qu'elles portent.
%
% ON EN CONSTRUIT PLUS QU'ON N'EN GARDE. Le rang demandé est r, mais les r plus grands
% poids ne sont pas les r premiers multi-indices dans l'ordre lexicographique : une
% fréquence élevée sur un axe court pèse moins qu'une fréquence moyenne sur les deux. La
% grille est donc prise assez large, triée par poids décroissant, et tronquée.
    L  = [domaine(2) - domaine(1), domaine(4) - domaine(3)];
    lo = [domaine(1), domaine(3)];
    r  = par.slam.r_kl;

    % Assez de modes par axe pour que la troncature choisisse réellement : le carré de
    % côté ceil(sqrt(r)) donnerait exactement r modes et le tri ne servirait à rien.
    n_axe = ceil(sqrt(r) * 1.6);
    [J1, J2] = ndgrid(1:n_axe, 1:n_axe);
    J = [J1(:), J2(:)]';                               % 2 x n_axe^2

    W = pi * J ./ L(:);                                % les fréquences portées
    S = densite_spectrale(W, par);
    [S, ord] = sort(S, 'descend');
    J = J(:, ord(1:r));
    S = S(1:r);

    % NORMALISATION. La densité spectrale n'est définie qu'à une constante près selon la
    % convention de transformée retenue, et cette constante se lit directement sur la
    % variance a priori du champ : sum_j lambda_j psi_j(x)^2 doit valoir sigma_f^2. On la
    % fixe au centre du domaine, là où la base de Dirichlet n'est pas encore contrainte
    % par le bord. Prendre la formule analytique de la constante marcherait aussi, mais
    % elle dépend du noyau et de la dimension, et elle ne dirait pas si la troncature a
    % mangé une part notable de la variance — ici le facteur obtenu le dit.
    phi_c  = fonctions(J, L, lo, mean(reshape(domaine, 2, 2), 1)');
    facteur = par.slam.sigma_f^2 / sum(S .* phi_c.^2);

    base = struct('name', 'laplace', 'r', r, 'Lambda', facteur * S, ...
                  'eval', @(P) fonctions(J, L, lo, P), 'J', J, 'domaine', domaine, ...
                  'variance_tronquee', facteur);
end

function Phi = fonctions(J, L, lo, P)
% Les fonctions propres elles-mêmes, r x M. Le produit sur les axes se fait par diffusion
% plutôt qu'en boucle : r et M valent quelques centaines chacun.
%
% Les requêtes hors du rectangle ne sont pas ramenées au bord, contrairement à la lecture
% de grille de krigeage.m : le prolongement du sinus est ce qu'il est, et une particule
% sortie du domaine doit voir une carte qui ne veut rien dire plutôt qu'une valeur de bord
% plausible. C'est son poids qui tranche.
    Phi = sqrt(4 / (L(1) * L(2))) ...
          * sin(pi * J(1, :)' .* (P(1, :) - lo(1)) / L(1)) ...
          .* sin(pi * J(2, :)' .* (P(2, :) - lo(2)) / L(2));
end

function base = base_rff(par)
% Les cosinus de Rahimi-Recht. Les fréquences sont tirées selon la densité spectrale
% normalisée du noyau, les phases uniformément, et le facteur sqrt(2/r) fait de la somme
% des r produits l'estimateur de Monte Carlo de la covariance.
    r = par.slam.r_rff;
    W = tirage_spectral(r, par);
    b = 2 * pi * rand(1, r);

    base = struct('name', 'rff', 'r', r, ...
                  'Lambda', par.slam.sigma_f^2 * ones(r, 1), ...
                  'eval', @(P) sqrt(2 / r) * cos(W' * P + b'), 'W', W, ...
                  'variance_tronquee', 1);
end

function S = densite_spectrale(W, par)
% La densité spectrale du noyau, aux fréquences W (2 x m), à une constante multiplicative
% près — celle-ci est reprise par la normalisation de la base.
    q = sum(W.^2, 1)';
    switch par.slam.noyau
        case 'se'
            S = exp(-par.slam.ell^2 * q / 2);
        case 'matern32'
            S = (3 / par.slam.ell^2 + q).^(-(1.5 + 1));       % nu + d/2, d = 2
        case 'matern52'
            S = (5 / par.slam.ell^2 + q).^(-(2.5 + 1));
        otherwise
            error('slam:noyau', 'Noyau inconnu : %s.', par.slam.noyau);
    end
end

function W = tirage_spectral(r, par)
% Un tirage de r fréquences selon la densité spectrale normalisée.
%
% Pour l'exponentielle quadratique c'est une gaussienne d'écart-type 1/ell. Pour un
% Matérn de régularité nu c'est une Student à 2 nu degrés de liberté et même échelle :
% sa densité est proportionnelle à (1 + ell^2 ||w||^2 / (2 nu))^{-(nu + d/2)}, qui est la
% densité spectrale à la constante près. La queue lourde est ce qui distingue les deux, et
% c'est elle qui porte la structure courte échelle d'un champ d'anomalie.
    switch par.slam.noyau
        case 'se'
            W = randn(2, r) / par.slam.ell;
        case 'matern32'
            W = randn(2, r) / par.slam.ell ./ sqrt(chi2rnd(3, 1, r) / 3);
        case 'matern52'
            W = randn(2, r) / par.slam.ell ./ sqrt(chi2rnd(5, 1, r) / 5);
        otherwise
            error('slam:noyau', 'Noyau inconnu : %s.', par.slam.noyau);
    end
end

%% ------------------------------------------------------------------ dispatch

function out = run_slam(kind, base, z_obs, scen, opt)
    switch kind
        case 'rbpf'
            out = run_rbpf(base, z_obs, scen, opt);
        case 'akkf'
            out = run_akkf(base, z_obs, scen, opt);
        case 'inertiel'
            out = run_inertiel(z_obs, scen, opt);
        otherwise
            error('slam:run', 'Algorithme inconnu : %s.', kind);
    end
end

%% ------------------------------------------------------------------ inertiel

function out = run_inertiel(z_obs, scen, opt)
% LA BORNE BASSE : le même nuage, propagé sans jamais être corrigé. C'est ce que donne la
% navigation à l'estime seule, et c'est à quoi un SLAM doit être comparé.
%
% Sans elle un chiffre de SLAM ne se lit pas. Le chapitre 4 se compare à la carte connue,
% qui dit ce qu'on perd à reconstruire ; le chapitre 5 a besoin des deux bouts, car ce
% qu'il mesure n'est pas une précision mais une CROISSANCE ÉVITÉE — et une erreur qui
% double quand celle de l'estime est multipliée par dix est un résultat, pas un échec.
%
% Aucune carte n'est apprise, donc h_eval est nulle et l'erreur de carte rendue est celle
% du champ lui-même, ce qui est exactement ce que vaut de ne rien estimer.
    [d_x, n_iter] = deal(size(scen.x_true, 1), numel(z_obs));
    n  = opt.n_part;
    Lq = chol(scen.Q, 'lower');
    Lp = chol(scen.P0, 'lower');

    X = scen.x_true(:, 1) + Lp * randn(d_x, 1) + Lp * randn(d_x, n);
    w = ones(1, n) / n;

    out = sortie_vide(n_iter);
    for k = 1:n_iter
        X = scen.F * X + scen.Bu(:, k) + Lq * randn(d_x, n);

        x_hat = X * w';
        Xc    = X - x_hat;
        P_hat = (Xc .* w) * Xc';

        e = scen.x_true(:, k + 1) - x_hat;
        out = consigner(out, k, e, P_hat, x_hat, false, [n; 0; 0; 0]);
    end

    out.alpha  = 0;
    out.h_eval = zeros(size(opt.P_eval, 2), 1);
    out.n_eff  = n;
end

%% ---------------------------------------------------------------------- RBPF

function out = run_rbpf(base, z_obs, scen, opt)
% Le filtre Rao-Blackwellisé : chaque particule porte sa trajectoire ET sa carte, et la
% seconde est corrigée analytiquement conditionnellement à la première.
%
% TROIS ÉTAPES SEULEMENT. L'identification du chapitre annule la pré-correction, la carte
% n'agissant pas sur le véhicule, et la prédiction de Kalman, le champ ne changeant pas
% avec le temps. Restent la propagation des particules, la mise à jour des poids sous la
% vraisemblance marginale, et la correction des cartes.
%
% L'ORDRE COMPTE. La correction vient après le rééchantillonnage, pour qu'une particule
% éliminée ne paie pas la correction d'une carte qu'on s'apprête à jeter — et surtout la
% régularisation déplace les particules, donc la base doit être réévaluée avant de
% corriger, sans quoi la carte serait mise à jour à une position que la particule n'occupe
% plus.
    [d_x, n_iter] = deal(size(scen.x_true, 1), numel(z_obs));
    n  = opt.n_part;
    r  = base.r;
    Lq = chol(scen.Q, 'lower');
    Lp = chol(scen.P0, 'lower');
    R  = opt.R_obs;

    verifier_memoire(r, n, opt);

    X = scen.x_true(:, 1) + Lp * randn(d_x, 1) + Lp * randn(d_x, n);
    w = ones(1, n) / n;

    % Une carte par particule. Le prior est celui de la base, identique pour toutes.
    %
    % opt.alpha_0 est le TÉMOIN DE CARTE CONNUE : les coefficients de la projection du
    % champ vrai sur la base, avec une covariance écrasée. Le filtre navigue alors sur une
    % carte qu'il n'a pas apprise, et tout ce qui le sépare du filtre du chapitre 4 est
    % l'erreur de représentation de la base. Sans ce témoin, un mauvais résultat ne dit
    % pas si c'est l'apprentissage ou le filtre qui échoue.
    [alpha, P] = prior_carte(base, n, opt);
    Pbuf = zeros(r, r, n);              % le second tableau du rééchantillonnage
    bloc = max(1, floor(opt.bloc_octets / (8 * r^2)));

    out = sortie_vide(n_iter);

    for k = 1:n_iter
        X = scen.F * X + scen.Bu(:, k) + Lq * randn(d_x, n);

        [S, innov] = innovation(base, X, alpha, P, z_obs(k), R, opt.sigma_f2);
        w = normalise(log(w') - 0.5 * (innov.^2 ./ S + log(S)));

        x_hat = X * w';
        Xc    = X - x_hat;
        P_hat = (Xc .* w) * Xc';

        resampled = 1 / sum(w.^2) < opt.n_th;
        if resampled
            A   = chol_psd(P_hat);
            idx = resample(w, n);
            X   = X(:, idx);
            alpha = alpha(:, idx);              % la carte suit sa particule

            % LES COVARIANCES CHANGENT DE PLACE PAR BLOCS, DANS UN SECOND TABLEAU ALLOUÉ
            % UNE FOIS. P(:,:,idx) d'un coup construirait un tableau de la taille de P
            % avant de remplacer l'ancien, donc le double au même instant, et un tableau
            % de cette taille demande de la mémoire contiguë qu'une session de travail
            % n'a pas toujours. Les deux tableaux s'échangent ensuite sans rien copier.
            for i0 = 1:bloc:n
                ii = i0:min(i0 + bloc - 1, n);
                Pbuf(:, :, ii) = P(:, :, idx(ii));
            end
            [P, Pbuf] = deal(Pbuf, P);
            w = ones(1, n) / n;

            % La régularisation ne touche que l'état. Perturber aussi les coefficients
            % demanderait un noyau sur la carte, et celle-ci n'a pas de raison d'être
            % dispersée : deux particules issues du même parent portent la même carte
            % parce qu'elles ont vu la même chose, et c'est la suite de leurs trajectoires
            % qui les séparera.
            h_reg = (8 * (d_x + 4) * (2 * sqrt(pi))^d_x / ...
                     (pi^(d_x/2) / gamma(d_x/2 + 1)))^(1/(d_x + 4)) * n^(-1/(d_x + 4));
            X = X + opt.alpha_reg * h_reg * A * epanechnikov(d_x, n);
        end

        % Correction de Kalman des N cartes, vectorisée sur la troisième dimension. Le
        % champ étant scalaire, S est un nombre par particule et le gain est un simple
        % quotient : aucune inversion n'est faite ici, ce qui est tout l'intérêt d'un
        % champ scalaire pour un filtre qui porte N covariances de rang r.
        %
        % La base est réévaluée après le rééchantillonnage, qui a déplacé les particules :
        % corriger une carte à une position que sa particule n'occupe plus la décalerait
        % du bruit de régularisation à chaque pas.
        [S, innov, Psi] = innovation(base, X, alpha, P, z_obs(k), R, opt.sigma_f2);
        PPsi = pagemtimes(P, Psi);                              % r x 1 x n
        K    = PPsi ./ reshape(S, 1, 1, n);
        alpha = alpha + reshape(K, r, n) .* innov';

        % UNE PARTICULE À LA FOIS, et non pagemtimes sur le nuage entier. Celui-ci alloue
        % un second tableau r x r x N avant d'en soustraire quoi que ce soit, donc trois
        % fois la taille de P au même instant, et le faire par blocs ne supprime que le
        % facteur, pas le pic. La mise à jour étant de rang un — le champ est scalaire —
        % chaque covariance se corrige sur place par un produit externe, sans qu'aucun
        % tableau intermédiaire de la taille du nuage n'existe jamais.
        %
        % Le coût est celui d'une boucle MATLAB de N tours par pas, soit quelques dizaines
        % de secondes sur une campagne entière, contre un plantage aléatoire selon ce que
        % la machine a de libre au moment où la campagne passe.
        Kf    = reshape(K, r, n);
        PPsif = reshape(PPsi, r, n);
        for i = 1:n
            P(:, :, i) = P(:, :, i) - Kf(:, i) * PPsif(:, i)';
        end

        e = scen.x_true(:, k + 1) - x_hat;
        % Le coût d'un pas : N corrections de rang r, chacune quadratique en r. Linéaire
        % en le nombre de particules, quadratique en le rang de la carte — l'inverse de
        % l'AKKF, dont le rang n'entre pas dans le coût du filtre.
        out = consigner(out, k, e, P_hat, x_hat, resampled, [n; 3 * n * r^2; r; 0]);
    end

    % La carte rendue est la moyenne des cartes des particules sous les poids finaux,
    % c'est-à-dire l'espérance a posteriori des coefficients. Prendre celle de la particule
    % la plus lourde donnerait une carte cohérente avec une trajectoire, mais une seule.
    out.alpha  = alpha * w';
    out.h_eval = (base.eval(opt.P_eval)' * out.alpha);
    out.n_eff  = 1 / sum(w.^2);
end

function [S, innov, Psi3] = innovation(base, X, alpha, P, z, R, sigma_f2)
% La covariance d'innovation et l'innovation de chaque particule, sous sa propre carte.
% C'est la vraisemblance marginale du chapitre : la carte est intégrée, pas fixée.
%
% LE TROISIÈME TERME EST LA VARIANCE QUE LA BASE NE REPRÉSENTE PAS. Une base de rang r ne
% porte qu'une part de la variance a priori du champ, sigma_f^2 - Psi' Lambda Psi, et le
% reste n'est nulle part : ni dans la carte, ni dans le bruit du capteur. L'omettre revient
% à déclarer que tout ce que la carte n'explique pas est du bruit de mesure — cinq nT ici,
% quand le résidu de troncature en vaut six le long du vol et trente près du bord du
% domaine. Le filtre devient alors sur-confiant là où sa base est la plus pauvre, c'est-à-
% dire précisément près du bord.
%
% Au premier pas, P vaut Lambda et la somme des deux premiers termes vaut sigma_f^2 : la
% covariance d'innovation est la variance a priori du champ plus celle du capteur, ce
% qu'elle doit être. Quand la carte est apprise, P tend vers zéro et il reste le bruit du
% capteur PLUS le résidu de troncature, qui lui ne s'apprend pas.
    n     = size(X, 2);
    Psi   = base.eval(X(1:2, :));                              % r x n
    Psi3  = reshape(Psi, size(Psi, 1), 1, n);
    porte = sum(base.Lambda .* Psi.^2, 1)';                    % ce que la base porte
    S     = reshape(pagemtimes(Psi3, 'transpose', pagemtimes(P, Psi3), 'none'), n, 1) ...
            + R + max(sigma_f2 - porte, 0);
    innov = z - sum(Psi .* alpha, 1)';
end

function [alpha, P] = prior_carte(base, n, opt)
    r = base.r;
    if isfield(opt, 'alpha_0') && ~isempty(opt.alpha_0)
        alpha = repmat(opt.alpha_0(:), 1, n);
        P     = repmat(diag(base.Lambda * opt.alpha_0_facteur), 1, 1, n);
    else
        alpha = zeros(r, n);
        P     = repmat(diag(base.Lambda), 1, 1, n);
    end
end

function verifier_memoire(r, n, opt)
% N covariances r x r en double, plus la copie que le rééchantillonnage en fait — celle-ci
% est inévitable, les indices pouvant se répéter. La mise à jour, elle, passe par blocs et
% ne coûte plus rien. Un dépassement se solde par un swap de plusieurs minutes par pas
% plutôt que par une erreur franche, d'où le garde-fou explicite.
    octets = 2 * 8 * r^2 * n;
    if octets > opt.mem_max
        error('slam:memoire', ...
              ['%d particules portant un rang %d demandent %.1f Go, la limite est %.1f Go. ' ...
               'Baisser slam.r_kl ou mc.n_part.'], n, r, octets / 2^30, opt.mem_max / 2^30);
    end
end

%% ---------------------------------------------------------------------- AKKF

function out = run_akkf(base, z_obs, scen, opt)
% Le filtre de Kalman à noyau adaptatif : une seule carte, et un nuage dont ce sont les
% POIDS qui sont corrigés par une récursion de Kalman écrite dans le RKHS.
%
% DEUX DIFFÉRENCES VISIBLES AVEC UN FILTRE PARTICULAIRE. Les poids ne sont pas des
% probabilités — la correction résout un problème de moindres carrés régularisés dans
% l'espace de Hilbert, et rien ne les contraint à être positifs — donc ni taille effective
% ni rééchantillonnage conditionnel. Le nuage est retiré à chaque pas d'une gaussienne
% portant l'estimée courante, ce qui est le prix à payer : une postérieure multimodale est
% ramenée à sa moyenne à chaque pas.
%
% LA CARTE EST APPRISE PAR MOINDRES CARRÉS RÉCURSIFS sur les couples (estimée, observation
% réelle), et c'est elle qui fournit les observations simulées du pas suivant. Les deux
% défauts que le chapitre analyse sont ici, et sont voulus : l'entrée de la régression est
% bruitée, et elle a été produite par une correction qui a déjà lu cette observation.
    [d_x, n_iter] = deal(size(scen.x_true, 1), numel(z_obs));
    n  = opt.n_part;
    r  = base.r;
    Lq = chol(scen.Q, 'lower');
    Lp = chol(scen.P0, 'lower');

    ech_x = opt.akkf.echelle_x(:);           % le noyau d'état travaille sur X ./ ech_x
    lam_x = opt.akkf.lambda_x;
    lam_z = opt.akkf.lambda_z;

    X  = scen.x_true(:, 1) + Lp * randn(d_x, 1) + Lp * randn(d_x, n);
    Xp = X;                                                  % le nuage de proposition
    Gp = gram(Xp ./ ech_x, Xp ./ ech_x, opt.akkf.gamma_x);
    w  = ones(n, 1) / n;
    C  = eye(n) / n;

    alpha = zeros(r, 1);
    Pa    = opt.akkf.echelle_p * eye(r);     % covariance des coefficients, cf. RLS

    out = sortie_vide(n_iter);

    for k = 1:n_iter
        X = scen.F * Xp + scen.Bu(:, k) + Lq * randn(d_x, n);

        % Prédiction dans le RKHS : ce que la propagation a perdu en représentativité du
        % nuage précédent, lu sur l'écart entre la projection régularisée et l'identité.
        D = (Gp + lam_x * eye(n)) \ Gp - eye(n);
        C = C + (D * D') / n;

        % Les observations simulées viennent de la carte courante, la vraie étant inconnue.
        z_sim = (base.eval(X(1:2, :))' * alpha)' + sqrt(opt.R_obs) * randn(1, n);

        ech_z = echelle_obs(z_sim, z_obs(k), opt);
        Gz    = gram(z_sim, z_sim, ech_z);
        gz    = gram(z_sim, z_obs(k), ech_z);

        K = C / (Gz * C + lam_z * eye(n));
        w = w + K * (gz - Gz * w);
        C = C - K * Gz * C;

        x_hat = X * w;

        % LE NUAGE EST CENTRÉ AVANT DE FORMER LA COVARIANCE, ce que le code de référence
        % ne fait pas — il écrit X C X'. Les deux coïncident tant que les colonnes de C
        % somment à zéro, ce qu'impose la contrainte sur les poids, mais C est initialisée
        % à I/N qui ne la vérifie pas : X C X' vaut alors le moment d'ordre deux NON
        % centré, soit la position quadratique moyenne. Sur un champ où les coordonnées
        % valent quelques unités cela passe inaperçu ; sur une carte où elles valent cent
        % kilomètres, la covariance déclarée au premier pas vaut (100 km)^2 et le nuage
        % retiré à partir d'elle couvre la carte entière. Centrer rend à P_hat la
        % dispersion du nuage, qui est aussi ce que le NEES du chapitre 4 mesure.
        Xc    = X - x_hat;
        P_hat = symetriser(Xc * C * Xc');

        % La carte, par moindres carrés récursifs sur (estimée, observation réelle).
        phi = base.eval(x_hat(1:2));
        g   = Pa * phi / (opt.akkf.oubli + phi' * Pa * phi);
        alpha = alpha + g * (z_obs(k) - phi' * alpha);
        Pa    = (Pa - g * phi' * Pa) / opt.akkf.oubli;

        % Rééchantillonnage inconditionnel, puis reprojection des poids et de leur
        % covariance sur le nouveau nuage. P_hat peut n'être que semi-définie, les poids
        % étant signés, d'où le facteur de Cholesky tronqué plutôt que mvnrnd.
        Xp = x_hat + chol_psd(P_hat) * randn(d_x, n);
        Gp = gram(Xp ./ ech_x, Xp ./ ech_x, opt.akkf.gamma_x);
        T  = (Gp + lam_x * eye(n)) \ gram(Xp ./ ech_x, X ./ ech_x, opt.akkf.gamma_x);
        w  = T * w;
        C  = symetriser(T * C * T');

        e = scen.x_true(:, k + 1) - x_hat;
        % Trois systèmes n x n par pas, plus les produits qui les accompagnent : c'est le
        % coût de l'AKKF, cubique en le nombre de particules là où celui du RBPF est
        % linéaire. Le rang de la carte n'y entre pas.
        out = consigner(out, k, e, P_hat, x_hat, true, [n; 4 * n^3; r; ech_z]);
    end

    out.alpha  = alpha;
    out.h_eval = base.eval(opt.P_eval)' * alpha;
    out.n_eff  = NaN;                        % sans objet : les poids ne sont pas des masses
end

function G = gram(A, B, ech)
% La matrice de Gram du noyau du RKHS. Gaussien, seul cas utilisé ici : les noyaux
% polynomiaux du code de référence ne portent pas la même notion de proximité et n'ont pas
% d'échelle qui se lise en unités du problème.
%
% ATTENTION À CE QUE ECH MESURE. Pour l'état, les colonnes ont déjà été divisées par leurs
% échelles respectives — sans quoi une position en mètres et une vitesse en mètres par
% seconde entreraient dans la même distance euclidienne — et ech est alors sans dimension.
% Pour l'observation, ech est en nT.
    G = exp(-pdist2(A', B').^2 / (2 * ech^2));
end

function ech = echelle_obs(z_sim, z_obs, opt)
% L'échelle du noyau d'observation, prise sur les écarts ENTRE LES OBSERVATIONS SIMULÉES
% ET L'OBSERVATION RÉELLE, et non sur la seule dispersion des premières.
%
% C'est un écart au code de référence, et il est nécessaire. Aux premiers pas la carte est
% vide, donc toutes les observations simulées valent zéro à quelques nT près quand la
% vraie en vaut cent cinquante : une échelle lue sur leur dispersion vaut alors le bruit
% du capteur, le noyau évalué entre les deux vaut exp(-450), et la correction s'applique à
% un vecteur numériquement nul. Elle ne rend pas des poids uniformes, comme la
% renormalisation d'un filtre particulaire le ferait, mais des poids arbitraires — la
% correction du RKHS est linéaire, pas multiplicative, et rien ne la ramène à l'identité
% quand plus rien ne ressemble à rien. Le filtre part alors dès le premier pas.
%
% Prendre la médiane des écarts au réel donne le bruit du capteur quand la carte est
% bonne, ce qui est l'échelle voulue, et une échelle large tant qu'elle ne l'est pas, ce
% qui laisse le noyau classer les particules au lieu de les annuler toutes.
    if strcmp(opt.akkf.echelle_z, 'mediane')
        ech = max(median(abs(z_sim - z_obs)), sqrt(opt.R_obs));
    else
        ech = opt.akkf.echelle_z;
    end
end

%% ------------------------------------------------------------------ communs

function out = sortie_vide(n_iter)
    out = struct('err', zeros(n_iter, 1), 'd2', zeros(n_iter, 1), ...
                 'n_used', zeros(n_iter, 1), 'n_modes', ones(n_iter, 1), ...
                 'cout', zeros(4, n_iter), 'track', zeros(2, n_iter), 'n_res', 0);
end

function out = consigner(out, k, e, P_hat, x_hat, resampled, cout)
    out.err(k)      = norm(e(1:2));
    out.d2(k)       = nees(e, P_hat);
    out.n_used(k)   = cout(1);
    out.cout(:, k)  = cout;
    out.track(:, k) = x_hat(1:2);
    out.n_res       = out.n_res + resampled;
end

function q = nees(e, P)
    [L, p] = chol(symetriser(P), 'lower');
    if p > 0
        q = Inf;
    else
        q = sum((L \ e).^2);
    end
end

function M = symetriser(M)
    M = (M + M') / 2;
end

function w = normalise(log_w)
    log_w(~isfinite(log_w)) = -Inf;
    w = exp(log_w - max(log_w))';
    s = sum(w);
    if isfinite(s) && s > 0
        w = w / s;
    else
        w = ones(1, numel(log_w)) / numel(log_w);
    end
end

function idx = resample(w, n)
% Rééchantillonnage systématique, identique à celui de navigation.m.
    r        = rand;
    cdf      = cumsum(w(:))';
    cdf(end) = 1;
    hi       = ceil(n * cdf - r);
    idx      = repelem(1:numel(w), max(hi - [0, hi(1:end-1)], 0))';
    if numel(idx) ~= n
        idx = [idx; repmat(idx(end), n - numel(idx), 1)];
        idx = idx(1:n);
    end
end

function E = epanechnikov(d, n)
    E    = zeros(d, n);
    todo = 1:n;
    while ~isempty(todo)
        m    = numel(todo);
        U    = randn(d, m);
        U    = U ./ sqrt(sum(U.^2, 1));
        rr   = rand(1, m).^(1 / d);
        keep = rand(1, m) <= 1 - rr.^2;
        E(:, todo(1, keep)) = U(:, keep) .* rr(1, keep);
        todo = todo(1, ~keep);
    end
end

function A = chol_psd(P)
    [V, D] = eig(symetriser(P));
    A      = V * diag(sqrt(max(diag(D), 0)));
end
