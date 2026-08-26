function nav = navigation()
% NAVIGATION  Un essai d'un filtre sur un modèle d'observation.
%
%   nav = navigation();
%   out = nav.run(kind, model, z_obs, scen, opt);
%
% KIND vaut 'PF', 'RPF' ou 'APF'. MODEL est un modèle de krigeage.m. OUT donne, par pas,
% l'erreur de position, le NEES, l'estimée, le nombre de points utilisés et le coût.
%
% Le modèle n'est atteint que par model.query : rien ici ne sait si la carte est vraie,
% interpolée ou krigée, et c'est ce qui permet à la campagne de croiser les deux listes.
    nav = struct('run', @run_filter, 'nees', @nees);
end

function out = run_filter(kind, model, z_obs, scen, opt)
    d_x    = size(scen.x_true, 1);
    n      = opt.n_part;
    n_iter = numel(z_obs);
    Lq     = chol(scen.Q, 'lower');
    Lp     = chol(scen.P0, 'lower');

    % Un point tiré de N(x_0, P_0), puis le nuage autour de lui avec le même P_0, comme
    % les sauvegardes : l'erreur initiale a donc pour covariance 2 P_0 quand le filtre
    % déclare P_0.
    X = scen.x_true(:, 1) + Lp * randn(d_x, 1) + Lp * randn(d_x, n);
    w = ones(1, n) / n;

    out.err     = zeros(n_iter, 1);
    out.d2      = zeros(n_iter, 1);
    out.n_used  = zeros(n_iter, 1);
    out.n_modes = zeros(n_iter, 1);
    out.cout    = zeros(3, n_iter);   % somme des fenêtres, somme de leurs cubes, ajustements
    out.track   = zeros(2, n_iter);
    out.n_res   = 0;

    for k = 1:n_iter
        [X, w, info] = step(kind, X, w, z_obs(k, :), scen.F, scen.Bu(:, k), Lq, model, opt);

        x_hat = X * w';
        Xc    = X - x_hat;
        P_hat = (Xc .* w) * Xc';

        e = scen.x_true(:, k + 1) - x_hat;
        out.err(k)      = norm(e(1:2));
        out.d2(k)       = nees(e, P_hat);
        out.n_used(k)   = info.n_used;
        out.n_modes(k)  = info.n_modes;
        out.cout(:, k)  = info.cout;
        out.track(:, k) = x_hat(1:2);
        out.n_res       = out.n_res + info.resampled;
    end
end

function [X, w, info] = step(kind, X, w, z, F, Bu, Lq, model, opt)
% PF  propager, pondérer, rééchantillonner quand le nuage s'appauvrit.
% RPF le même, plus un noyau qui disperse le nuage rééchantillonné.
% APF classer les parents sur la transition, rééchantillonner là-dessus, puis propager et
%     pondérer contre ce classement. Il rééchantillonne à chaque pas et lit la carte deux
%     fois.
    [d, n] = size(X);

    if strcmp(kind, 'APF')
        Mu  = F * X + Bu;                              % moyenne de la transition
        ctx = context(Mu, w);
        [z_mu, R_mu, n_mu] = model.query(Mu, ctx);
        log_lam = log_likelihood(z, z_mu, R_mu);
        lambda  = normalise(log(w') + log_lam);

        idx = resample(lambda, n);
        X   = F * X(:, idx) + Bu + Lq * randn(d, n);

        % Le contexte de la seconde correction. Par défaut celui du nuage propagé, comme
        % le main.m de 2A : le rééchantillonnage auxiliaire a concentré le nuage, donc les
        % centres d'avant ne décrivent plus où sont les particules — et pour un modèle
        % clusterisé c'est toute la partition qui devient fausse. 'reuse' réemploie celui
        % de la passe auxiliaire, une requête de moins par pas.
        if isfield(opt, 'apf_ctx') && strcmp(opt.apf_ctx, 'reuse')
            ctx2 = ctx;
        else
            ctx2 = context(X, ones(1, n) / n);   % rééchantillonnées, donc équipondérées
        end
        [z_x, R_x, n_x, n_modes, cout] = model.query(X, ctx2);
        w = normalise(log_likelihood(z, z_x, R_x) - log_lam(idx));

        info = struct('n_used', max(n_mu, n_x), 'n_modes', n_modes, 'resampled', true, ...
                      'cout', 2 * cout);   % deux lectures de carte par pas
        return
    end

    X   = F * X + Bu + Lq * randn(d, n);
    ctx = context(X, w);
    [z_x, R_x, n_used, n_modes, cout] = model.query(X, ctx);
    w   = normalise(log(w') + log_likelihood(z, z_x, R_x));

    resampled = 1 / sum(w.^2) < opt.n_th;
    if resampled
        Xc = X - X * w';
        A  = chol_psd((Xc .* w) * Xc');                % P = A A'
        X  = X(:, resample(w, n));
        w  = ones(1, n) / n;

        if strcmp(kind, 'RPF')
            h = (8 * (d + 4) * (2 * sqrt(pi))^d / (pi^(d/2) / gamma(d/2 + 1)))^(1/(d+4)) ...
                * n^(-1/(d+4));
            X = X + opt.alpha_reg * h * A * epanechnikov(d, n);
        end
    end
    info = struct('n_used', n_used, 'n_modes', n_modes, 'resampled', resampled, 'cout', cout);
end

function ctx = context(X, w)
% Ce que le modèle d'observation a le droit de savoir du nuage. Les modèles statiques
% l'ignorent ; c'est toute l'entrée des modèles adaptatifs.
    ctx.X      = X;
    ctx.w      = w;
    ctx.x_pred = X * w';
    Xc         = X - ctx.x_pred;
    ctx.P_pred = (Xc .* w) * Xc';
end

function l = log_likelihood(z, z_hat, R_hat)
% Le terme en log R n'est pas une constante que la normalisation efface : R varie d'une
% particule à l'autre dès que la carte est reconstruite, et l'omettre ferait préférer au
% filtre les régions où la carte est la moins connue. La variance du capteur est déjà
% dans R_hat pour tous les modèles.
    l = -0.5 * ((z - z_hat).^2 ./ R_hat + log(R_hat));
end

function w = normalise(log_w)
% Poids depuis les log-poids, décalés de leur maximum. Un nuage que la carte contredit
% partout tombe à zéro ; l'essai est perdu de toute façon, et le repli uniforme évite de
% répandre des NaN sur la campagne.
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
% Rééchantillonnage systématique : un tirage, puis un peigne régulier sur les poids
% cumulés. Les descendances se lisent sur une différence de plafonds, sans parcours.
    r        = rand;
    cdf      = cumsum(w(:))';
    cdf(end) = 1;
    hi       = ceil(n * cdf - r);
    idx      = repelem(1:numel(w), max(hi - [0, hi(1:end-1)], 0))';

    if numel(idx) ~= n                    % seul un vecteur de poids dégénéré arrive ici
        idx = [idx; repmat(idx(end), n - numel(idx), 1)];
        idx = idx(1:n);
    end
end

function q = nees(e, P)
% L'erreur au carré dans les unités de la covariance déclarée. Infinie quand celle-ci est
% singulière, ce que donne un nuage effondré : le filtre est certain, et faux.
    [L, p] = chol((P + P') / 2, 'lower');
    if p > 0
        q = Inf;
    else
        q = sum((L \ e).^2);
    end
end

function E = epanechnikov(d, n)
% n tirages du noyau d'Epanechnikov : un point uniforme dans la boule unité de R^d, gardé
% avec probabilité 1 - ||x||^2. Le taux d'acceptation vaut 2 / (d + 2), donc la boucle
% recomplète les colonnes rejetées au lieu de tout retirer.
    E    = zeros(d, n);
    todo = 1:n;
    while ~isempty(todo)
        m    = numel(todo);
        U    = randn(d, m);
        U    = U ./ sqrt(sum(U.^2, 1));
        r    = rand(1, m).^(1 / d);
        keep = rand(1, m) <= 1 - r.^2;
        E(:, todo(1, keep)) = U(:, keep) .* r(1, keep);
        todo = todo(1, ~keep);
    end
end

function A = chol_psd(P)
% Un facteur A avec P = A A', pour une covariance que l'appauvrissement peut avoir rendue
% singulière.
    [V, D] = eig((P + P') / 2);
    A      = V * diag(sqrt(max(diag(D), 0)));
end
