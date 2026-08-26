% CAMPAGNE_CHAPITRE_4
%
% Navigation magnétique par filtrage particulaire quand la carte n'est pas connue mais
% reconstruite depuis un relevé. Une campagne croise une liste de filtres et une liste de
% modèles d'observation, et c'est la seconde qui compte.
%
% TROIS FICHIERS : navigation.m tient les filtres, krigeage.m le relevé et les modèles,
% celui-ci les paramètres, la vérité, la campagne, les tableaux et les figures.
%
% DEUX SÉLECTEURS, chacun un switch plus bas.
%
%   par.scenario  le point de fonctionnement. 'reference' est celui du chapitre ; les
%                 quatre autres sont gardés parce que des résultats y sont attachés.
%                 Attention : 'eusipco' est le tableau I de l'article et '2A' le code qui
%                 l'accompagne, et les deux ne portent PAS les mêmes réglages.
%   par.campagne  ce qu'on mesure : les listes à croiser, le nombre d'essais, et le
%                 paramètre balayé.
%
% Avec une carte reconstruite, R_krig est le terme dominant et spatialement variable de
% sigma_eff, et c'est là qu'est la substance de la comparaison : une interpolation déclare
% une variance partout la même, un processus gaussien déclare celle à laquelle il a droit,
% et les variantes adaptatives la tirent des échantillons réellement proches.

clearvars -except over; close all; clc;   % `over` survit : voir surcharge() en fin de fichier
here = fileparts(mfilename('fullpath'));
addpath(here);

%% Paramètres
par.compare.filters = {'RPF'};
par.compare.models  = {'carte', 'ok', 'ak', 'cak'};

par.scenario   = 'reference';       % 'reference', 'acquisition', 'eusipco', '2A', 'balayage'
par.campagne   = 'libre';           % 'A', 'B', 'C', 'D' ou 'libre' — voir le second switch
par.run.seed   = 123456789;
par.mc.n_runs  = 10;

% Ce qui ne dépend pas du scénario. Les cas ci-dessous redéfinissent ce qui leur est
% propre, et la surcharge appliquée après eux a le dernier mot sur les deux.
par.unit    = 'nT';
par.filt.alpha_reg = 0.3;
par.krig.n_side       = 80;          % relevé n_side x n_side, celui de 2A
par.krig.ell_init     = 10e3;        % m, départ imposé à fitrgp, voir krigeage
par.krig.n_ml         = 500;         % points de l'ajustement global
par.krig.n_chunk      = 2000;        % points krigés d'un coup
par.krig.refit        = 'fenetre';   % 'fenetre' ou 'global', voir krigeage
par.krig.window_factor = 2;          % alpha de l'article, coefficient de dilatation
par.krig.n_min        = 25;          % plancher sur les points retenus (2A)
par.krig.bandwidth    = 1200;        % m, celle de 2A ; 'reference' la ramène à 1000
par.map.name   = 'carte_magnetometrie_anomalie_mexique.mat';
par.map.file   = fullfile(here, '..', '..', 'MATLAB 2A - sauvegarde', 'Cartes', par.map.name);
par.path.cache = fullfile(here, 'cache');
par.path.out   = fullfile(here, 'resultats');

par = surcharge(par);               % pour que `over.scenario` soit lu avant le switch

switch par.scenario
    % Le scénario du chapitre. Un seul vol, deux régimes : le nuage est multimodal les
    % trente premiers pas, unimodal ensuite — la comparaison des krigeages et celle des
    % filtres se lisent sur la première phase, le balayage de maille sur la seconde.
    %
    % Les demi-tours débordent des branches d'un rayon, donc le point le plus bas du vol
    % est à `marge` du bord et non à `marge + R` : c'est lui qui limite sigma_0.
    case 'reference'
        par.mc.n_steps = 150;            % la forme en offre 160
        par.mc.n_part  = 5000;
        par.dt         = 5;                      % s, donc des pas de 2,50 km ronds
        par.sigma_obs  = 5;                      % nT, magnétomètre compensé
        par.x0         = [0; 0; 0; 500];         % la tondeuse place elle-même le départ,
                                                 % seule la vitesse de 500 m/s est lue
        par.traj.forme      = 'tondeuse';
        par.traj.n_branches = 3;
        par.traj.bank       = 45;        % deg, 1,41 g, d'où un rayon de 25,5 km
        par.traj.marge      = 35e3;      % m
        par.sigma_0    = [12000, 12000, 2, 2];   % 2,9 sigma de marge au bord : au-delà, le
                                                 % nuage initial déborde et lit une valeur
                                                 % ramenée au bord
        par.sigma_q    = [100, 100, 1, 1];
        par.krig.sigma_map = 10;         % nT
        par.krig.tab_pitch = 250;
        par.krig.n_max     = 2000;
        par.krig.bandwidth = 1000;       % m

    % Le tableau I de l'article. Trois écarts avec le code de 2A — relevé à 80 nT et non
    % 10, Q à 40 m et non 10, P_0 à 3 km et non 2 — qui vont tous dans le sens d'un nuage
    % plus dispersé.
    case 'eusipco'
        par.mc.n_steps = 600;
        par.mc.n_part  = 5000;
        par.dt         = 0.5;                    % s
        par.sigma_obs  = 80;                     % nT, R_k
        par.x0         = [80e3; 15e3; 50; 500];  % 502 m/s
        par.traj.forme = 'ligne';
        par.sigma_0    = [3000, 3000, 3, 3];     % P_0
        par.sigma_q    = [40, 40, 0.4, 0.4];     % Q_k
        par.krig.sigma_map = 80;                 % nT, R tilde
        par.krig.tab_pitch = 500;
        par.krig.n_max     = 6400;       % pas de plafond : l'article mesure justement la
                                         % queue de distribution où l'AK ramasse les points
                                         % situés entre les modes, et un plafond l'effacerait

    % Le prédécesseur de reference : même physique, mais en ligne droite sur 60 pas.
    case 'acquisition'
        par.mc.n_steps = 60;
        par.mc.n_part  = 5000;
        par.dt         = 5;                      % s
        par.sigma_obs  = 5;                      % nT
        par.x0         = [80e3; 25e3; 50; 500];  % assez rentré pour un nuage de 15 km
        par.traj.forme = 'ligne';
        par.sigma_0    = [15000, 15000, 2, 2];
        par.sigma_q    = [126.5, 126.5, 1.265, 1.265];
        par.krig.sigma_map = 10;
        par.krig.tab_pitch = 250;
        par.krig.n_max     = 2000;       % assez haut pour que la fenêtre de l'AK enfle
                                         % comme elle doit — c'est le mécanisme à mesurer —
                                         % mais borné : sans borne elle demande les 6400
                                         % points et une factorisation 6400^3 par pas

    case '2A'
        par.mc.n_steps = 600;
        par.mc.n_part  = 5000;
        par.dt         = 0.5;                    % s
        par.sigma_obs  = 100;                    % nT, le magnétomètre non compensé
        par.x0         = [80e3; 15e3; 50; 500];  % 502 m/s
        par.traj.forme = 'ligne';
        par.sigma_q    = [10, 10, 0.4, 0.4];
        par.sigma_0    = [2000, 2000, 2, 2];
        par.krig.sigma_map = 10;                 % nT, ce que le code de 2A déclare
        par.krig.n_max     = 200;
        par.krig.tab_pitch = 500;        % le vol en ligne fait 151 km, donc un corridor
                                         % quatre fois plus long qu'en tondeuse. À 250 m
                                         % la tabulation coûte 10600 Gflop pour gagner
                                         % 5.8 nT, quand le capteur en apporte 100

    case 'balayage'
        par.mc.n_steps = 400;            % plafond : la forme en donne ce qu'elle peut
        par.mc.n_part  = 10000;
        par.dt         = 5;                      % s
        par.sigma_obs  = 5;                      % nT, magnétomètre compensé
        par.x0         = [80e3; 5e3; 50; 500];   % le cap de 2A, 502 m/s
        par.traj.forme      = 'tondeuse';
        par.traj.n_branches = 4;         % espacées d'un diamètre de virage, 51 km
        par.traj.bank       = 45;        % deg, ce qui fixe le rayon à 25.7 km
        par.traj.marge      = 10e3;      % m, distance minimale au bord de carte
        par.sigma_q    = [126.5, 126.5, 1.265, 1.265];
        par.sigma_0    = [2000, 2000, 2, 2];
        par.krig.sigma_map = 10;
        par.krig.n_max     = 200;        % plafond : sans lui, un nuage large demande un
                                         % disque de 90 km, donc tout le relevé, et l'AK
                                         % redevient un OK avec une factorisation par pas
        par.krig.tab_pitch = 250;        % à 500 m l'interpolation de l'estimateur laissait
                                         % jusqu'à 5.8 nT, contre les 8 nT d'incertitude
                                         % qu'il déclare

    otherwise
        error('campagne_chapitre_4:scenario', 'Scénario inconnu : %s.', par.scenario);
end

% LES CAMPAGNES DU CHAPITRE. Chacune fixe les deux listes à croiser, le nombre d'essais,
% et le paramètre qu'elle balaie — le pilote refait alors le relevé, les modèles et le
% croisement pour chaque valeur, et empile tout dans un seul res.
%
% L'ordre compte : A donne la maille et le filtre que B et D emploieront, et C dit si
% l'ajustement peut être plafonné, ce dont dépend la faisabilité de B.
switch par.campagne
    case 'A'    % maille x filtre, sur les trois modèles qui ne coûtent rien par pas.
                % Sort la maille de référence — celle où le bilinéaire lâche mais où le
                % krigeage tient — et le filtre.
        par.compare.models  = {'carte', 'bilineaire', 'ok'};
        par.compare.filters = {'PF', 'RPF', 'APF'};
        par.mc.n_runs       = 100;
        par.balayage = struct('champ', 'krig.n_side', 'valeurs', [20 30 40 60 80]);

    case 'B'    % les cinq modèles, à maille et filtre fixés. C'est LA comparaison.
        par.compare.models  = {'carte', 'bilineaire', 'ok', 'ak', 'cak'};
        par.compare.filters = {'RPF'};
        par.mc.n_runs       = 100;
        par.balayage = struct('champ', '', 'valeurs', 0);

    case 'C'    % le plafond d'ajustement des hyperparamètres. 6400 vaut « aucun » :
                % la fenêtre ne dépasse jamais n_max, donc le plafond ne mord pas.
        par.compare.models  = {'ak', 'cak'};
        par.compare.filters = {'RPF'};
        par.mc.n_runs       = 30;
        par.balayage = struct('champ', 'krig.n_ml_win', 'valeurs', [250 500 1000 1500 6400]);

    case 'D'    % à coût égal : les particules balayées, pour comparer AK et CA-OK à
                % temps de calcul égal plutôt qu'à paramètres égaux.
        par.compare.models  = {'ak', 'cak'};
        par.compare.filters = {'RPF'};
        par.mc.n_runs       = 100;
        par.balayage = struct('champ', 'mc.n_part', 'valeurs', [1250 2500 5000 10000]);

    case 'libre'    % ce que les lignes du haut ont posé, sans balayage
        par.balayage = struct('champ', '', 'valeurs', 0);

    otherwise
        error('campagne_chapitre_4:campagne', 'Campagne inconnue : %s.', par.campagne);
end

par = surcharge(par);                % la surcharge a le dernier mot sur les deux switch

%% Carte et vérité
c   = load(par.map.file);
map = struct('h', double(c.mesures).', 'step', double(c.pas(:))' .* [1 1]);

d_x = numel(par.x0);
par.traj.extent = [size(map.h, 2), size(map.h, 1)] .* map.step;

scen.F  = [eye(2), par.dt * eye(2); zeros(2), eye(2)];
scen.Q  = diag(par.sigma_q).^2;
scen.P0 = diag(par.sigma_0).^2;
[scen.x_true, scen.Bu, traj] = trajectoire(par);
n_iter = traj.n_pas;                     % la forme peut en donner moins que demandé

krig    = krigeage();          % construit ici : sa lecture de grille sert dès la vérité
nav     = navigation();
z_true  = krig.lire(map.h, map.step, scen.x_true(1:2, 2:end));
extent = [size(map.h, 2), size(map.h, 1)] .* map.step;
margin = min([min(scen.x_true(1:2, :), [], 2)', extent - max(scen.x_true(1:2, :), [], 2)']);
if margin < 0
    error('campagne_chapitre_4:trajectoire', 'La trajectoire sort de la carte de %.0f km.', ...
          -margin / 1e3);
end

% Le corridor sur lequel l'estimateur statique est tabulé : le vol, plus la place du
% nuage. La marge doit couvrir celui-ci, sinon les particules les plus lointaines
% interrogent la tabulation hors de son domaine et lisent son bord.
pad = max(15e3, 3.5 * max(par.sigma_0(1:2)));
par.krig.corridor = [max(min(scen.x_true(1, :)) - pad, map.step(1)), ...
                     min(max(scen.x_true(1, :)) + pad, extent(1) - map.step(1)), ...
                     max(min(scen.x_true(2, :)) - pad, map.step(2)), ...
                     min(max(scen.x_true(2, :)) + pad, extent(2) - map.step(2))];

fprintf('Carte %d x %d noeuds à %.0f m, champ %.0f %s rms\n', ...
        size(map.h, 2), size(map.h, 1), map.step(1), std(map.h(:)), par.unit);
% Parcouru et non traj.chemin, qui est la longueur de la polyligne entière : quand
% par.mc.n_steps arrête le vol avant sa fin, les deux diffèrent et diviser le second par
% le nombre de pas donne un pas de temps fictif.
parcouru = n_iter * norm(par.x0(3:4)) * par.dt;
fprintf('Vol %s : %d pas, parcouru %.0f km sur %.0f, durée %.0f s, marge %.0f km\n', ...
        traj.forme, n_iter, parcouru / 1e3, traj.chemin / 1e3, n_iter * par.dt, ...
        margin / 1e3);
fprintf('   %.2f km par pas, inclinaison max %.0f deg, vitesse de %.0f à %.0f m/s\n', ...
        parcouru / 1e3 / n_iter, traj.bank, traj.v_min, traj.v_max);
fprintf('   le modèle laisse %.0f m et %.2f m/s par pas, contre %.0f m et %.2f m/s de bruit\n\n', ...
        traj.resid_pos, traj.resid_vit, par.sigma_q(1), par.sigma_q(3));

%% Relevé, modèles et campagne, une passe par valeur balayée
res = struct('name', {}, 'filter', {}, 'model', {}, 'balayage', {}, 'rmse', {}, ...
             'nees', {}, 'track', {}, 'n_used', {}, 'n_modes', {}, 'n_res', {}, ...
             'wall', {}, 'err', {});
valeurs = par.balayage.valeurs;

for b = 1:numel(valeurs)
if ~isempty(par.balayage.champ)
    par = poser(par, par.balayage.champ, valeurs(b));
    fprintf('\n===== %s = %g =====\n', par.balayage.champ, valeurs(b));
end

survey  = krig.releve(map, par);
models  = cellfun(@(k) krig.modele(k, map, survey, par), par.compare.models, ...
                  'UniformOutput', false);
fprintf('\n');

% n_th est dérivé de n_part, que la campagne D balaie : il se recalcule ici.
opt = struct('n_part', par.mc.n_part, 'n_th', 0.5 * par.mc.n_part, ...
             'alpha_reg', par.filt.alpha_reg);

for f = 1:numel(par.compare.filters)
    for m = 1:numel(models)
        kind  = par.compare.filters{f};
        model = models{m};
        err   = zeros(n_iter, par.mc.n_runs);
        d2    = zeros(n_iter, par.mc.n_runs);
        nu    = zeros(n_iter, par.mc.n_runs);
        nm    = zeros(n_iter, par.mc.n_runs);
        co    = zeros(3, n_iter, par.mc.n_runs);
        nres  = zeros(1, par.mc.n_runs);
        track = zeros(2, n_iter, par.mc.n_runs);
        t0    = tic;

        for r = 1:par.mc.n_runs
            rng(par.run.seed + r);
            z_obs = z_true + par.sigma_obs * randn(n_iter, 1);
            out   = nav.run(kind, model, z_obs, scen, opt);
            err(:, r)      = out.err;
            d2(:, r)       = out.d2;
            nu(:, r)       = out.n_used;
            nm(:, r)       = out.n_modes;
            co(:, :, r)    = out.cout;
            nres(r)        = out.n_res;
            track(:, :, r) = out.track;
        end

        e = numel(res) + 1;
        res(e).name   = sprintf('%s / %s', kind, model.name);
        res(e).filter = kind;
        res(e).model  = model.name;
        res(e).balayage = valeurs(b);      % la valeur qui a produit cette ligne
        res(e).rmse   = sqrt(mean(err.^2, 2));
        res(e).nees   = mean(d2, 2, 'omitnan') / d_x;
        res(e).track  = track(:, :, 1);
        res(e).n_used = nu;                % par pas et par essai, comme les modes : c'est la
                                           % QUEUE de cette distribution que l'article oppose
                                           % entre AK et CA-OK, et une moyenne l'efface
        res(e).cout   = mean(co, 3);       % somme des fenêtres, somme de leurs cubes,
                                           % nombre d'ajustements — moyennés sur les essais
        res(e).n_modes = nm;               % par pas et par essai : c'est ce qui sépare
                                           % le CA-OK de l'AK, et sa moyenne le cache
        res(e).n_res  = mean(nres) / n_iter;
        res(e).wall   = toc(t0);
        res(e).err    = err;               % les erreurs brutes, pour les médianes et les
                                           % taux de convergence qu'une RMSE ne donne pas

        fprintf(['%-26s ARMSE %6.0f m, méd. %5.0f m, finale %5.0f m, NEES %6.2f, ' ...
                 'rééch. %3.0f %%, ñ %5.0f, modes %4.1f, %6.1f s\n'], ...
                res(e).name, mean(res(e).rmse), median(err(:)), res(e).rmse(end), ...
                mean(res(e).nees(end-9:end)), 100 * res(e).n_res, mean(nu(:)), ...
                mean(nm(:)), res(e).wall);
    end
end
end                                        % fin du balayage

%% Résultats
% Une étiquette par ligne : le modèle seul ne suffit plus dès qu'un balayage produit
% plusieurs passes, ni dès que plusieurs filtres sont croisés.
if numel(par.compare.filters) > 1
    etiq = arrayfun(@(r) sprintf('%s/%s', r.filter, r.model(1:min(8, end))), res, ...
                    'UniformOutput', false);
else
    etiq = {res.model};
end
if ~isempty(par.balayage.champ)
    etiq = arrayfun(@(i) sprintf('%s %g', etiq{i}, res(i).balayage), 1:numel(res), ...
                    'UniformOutput', false);
end

fprintf('\n%-6s', 'k');  fprintf('%14s', etiq{:});  fprintf('   |');
fprintf('%10s', etiq{:});  fprintf('\n%-6s%s RMSE (m)%s NEES\n', '', ...
        repmat(' ', 1, 14 * numel(res) - 9), repmat(' ', 1, 10 * numel(res) - 5));
for k = unique([1, 5, 10:10:n_iter])
    fprintf('%-6d', k);
    fprintf('%14.0f', arrayfun(@(r) r.rmse(k), res));
    fprintf('   |');
    fprintf('%10.2f', arrayfun(@(r) r.nees(k), res));
    fprintf('\n');
end

if ~exist(par.path.out, 'dir')
    mkdir(par.path.out);
end
stamp = sprintf('%s_%s_%s_%druns', char(datetime('now', 'Format', 'yyyyMMdd_HHmmss')), ...
                par.scenario, par.campagne, ...
                par.mc.n_runs);
save(fullfile(par.path.out, [stamp '.mat']), 'par', 'res', 'scen', 'survey');
fprintf('\nSauvegardé dans resultats/%s.mat\n', stamp);

%% Figures
colours = [0.20 0.30 0.65; 0.85 0.15 0.15; 0.10 0.65 0.25; 0.55 0.35 0.70
           0.90 0.55 0.10; 0.20 0.60 0.70];
t       = (1:n_iter) * par.dt;

figure('Color', 'w', 'Name', 'RMSE'); hold on;
for i = 1:numel(res)
    plot(t, res(i).rmse, 'Color', colours(mod(i-1, size(colours,1)) + 1, :), 'LineWidth', 1.2);
end
% Échelle linéaire : la log écrase le régime établi, qui est ce qu'on compare, pour
% faire de la place au transitoire initial, qui ne dépend que de P0.
hold off; grid on; box on;
xlabel('t (s)'); ylabel('RMSE (m)'); legend(res.name, 'Location', 'northeast');
title('Erreur de position selon la reconstruction de carte');

figure('Color', 'w', 'Name', 'NEES'); hold on;
lo = chi2inv(0.025, par.mc.n_runs * d_x) / (par.mc.n_runs * d_x);
hi = chi2inv(0.975, par.mc.n_runs * d_x) / (par.mc.n_runs * d_x);
fill([t, fliplr(t)], [lo * ones(size(t)), hi * ones(size(t))], [0.9 0.9 0.9], ...
     'EdgeColor', 'none');
for i = 1:numel(res)
    plot(t, res(i).nees, 'Color', colours(mod(i-1, size(colours,1)) + 1, :), 'LineWidth', 1.2);
end
hold off; grid on; box on; set(gca, 'YScale', 'log');
xlabel('t (s)'); ylabel('NEES / d_x');
legend([{'bande de cohérence 95 %'}, res.name], 'Location', 'southwest');
title('Cohérence : l''erreur dans les unités que le filtre déclare');

figure('Color', 'w', 'Name', 'Trajectoires');
ix = max(floor(par.krig.corridor(1) / map.step(1)), 1):floor(par.krig.corridor(2) / map.step(1));
iy = max(floor(par.krig.corridor(3) / map.step(2)), 1):floor(par.krig.corridor(4) / map.step(2));
imagesc(ix * map.step(1) / 1e3, iy * map.step(2) / 1e3, map.h(iy, ix));
set(gca, 'YDir', 'normal'); colormap(gca, gray); hold on;
cb = colorbar; cb.Label.String = sprintf('champ (%s)', par.unit);
plot(survey.X(1, :) / 1e3, survey.X(2, :) / 1e3, 'w.', 'MarkerSize', 2);
plot(scen.x_true(1, :) / 1e3, scen.x_true(2, :) / 1e3, 'w-', 'LineWidth', 3.5);
plot(scen.x_true(1, :) / 1e3, scen.x_true(2, :) / 1e3, 'k-', 'LineWidth', 1.8);
for i = 1:numel(res)
    plot(res(i).track(1, :) / 1e3, res(i).track(2, :) / 1e3, '--', ...
         'Color', colours(mod(i-1, size(colours,1)) + 1, :), 'LineWidth', 1.2);
end
hold off; axis image; box on;
xlabel('x (km)'); ylabel('y (km)');
legend([{'relevé'}, {''}, {'vérité'}, res.name], 'Location', 'northeast');
title('Un essai de chaque modèle sur le champ');

%% -------------------------------------------------------------------- surcharge

function par = surcharge(par)
% Applique la structure `over` du workspace appelant, si elle existe, par-dessus par.
% C'est ce qui permet de balayer une campagne sans éditer ce fichier :
%
%   matlab -batch "cd('...'); over.krig.n_side = 40; over.mc.n_runs = 100; campagne_chapitre_4"
%
% Deux niveaux suffisent, la structure des paramètres n'en a pas d'autres. par.filt.n_th
% est dérivé de par.mc.n_part après la surcharge, donc le surcharger n'a pas d'effet.
    if evalin('caller', 'exist(''over'', ''var'')') ~= 1
        return
    end
    over = evalin('caller', 'over');
    if ~isstruct(over)
        return
    end
    for c1 = fieldnames(over)'
        v = over.(c1{1});
        if isstruct(v) && isfield(par, c1{1}) && isstruct(par.(c1{1}))
            for c2 = fieldnames(v)'
                par.(c1{1}).(c2{1}) = v.(c2{1});
            end
        else
            par.(c1{1}) = v;
        end
    end
end

%% ------------------------------------------------------------------- trajectoire

function [x_true, Bu, info] = trajectoire(par)
% TRAJECTOIRE  La vérité, et la commande qui la rend compatible avec le modèle du filtre.
%
%   [x_true, Bu, info] = trajectoire(par)
%
% x_true est 4 x (n+1), [x, y, vx, vy]', et Bu est 4 x n, l'incrément que la commande
% ajoute à chaque pas. Le filtre propage x_{k+1} = F x_k + Bu(:,k) + bruit.
%
% DEUX FORMES :
%
%   'ligne'     vitesse constante depuis par.x0, aucune commande. C'est le cas des
%               sauvegardes, et le moins informatif qu'une navigation par carte puisse
%               prendre : tout ce que le filtre apprend de l'erreur travers, il
%               l'apprend de la seule forme du champ.
%   'tondeuse'  n branches parallèles reliées par des demi-tours. Les branches sont
%               espacées d'un diamètre de virage, ce qui les rend indépendantes tant que
%               cet espacement dépasse la longueur de corrélation du champ.
%
% POURQUOI UNE COMMANDE. Un virage casse le modèle à vitesse constante. Deux façons de
% le porter : les accéléromètres lisent la manœuvre et la commande la transporte, ce qui
% est la navigation inertielle que le chapitre décrit ; ou personne ne la lit et le bruit
% de modèle doit la couvrir, ce qui demande un Q bien plus large que ce qu'une mesure
% scalaire par pas peut rattraper. La première est écrite ici.
%
% u_k n'est pas l'accélération analytique mais celle qui rend le pas exact :
%
%   u_k = 2 (p_{k+1} - p_k - v_k dt) / dt^2,     v_{k+1} = v_k + u_k dt
%
% La vérité satisfait alors le modèle du filtre exactement, et tout ce qui reste à
% expliquer vient de l'observation, qui est ce qu'on mesure. L'accélération analytique
% laisserait, elle, un résidu de l'ordre de l'accélération centripète fois dt^2 sur deux,
% soit 120 m par pas à 1 g et dt = 5 s — comparable au bruit de modèle lui-même.
    dt = par.dt;
    v  = norm(par.x0(3:4));

    switch par.traj.forme
        case 'ligne'
            n      = par.mc.n_steps;
            x_true = zeros(4, n + 1);
            x_true(:, 1) = par.x0;
            F = [eye(2), dt * eye(2); zeros(2), eye(2)];
            for k = 1:n
                x_true(:, k + 1) = F * x_true(:, k);
            end
            Bu   = zeros(4, n);
            % Les memes champs que la tondeuse, pour qu'un appelant n'ait pas a savoir
            % quelle forme il a demandee. Une ligne droite ne manoeuvre pas, donc le
            % modele est exact et la vitesse constante.
            info = struct('forme', 'ligne', 'chemin', v * n * dt, 'bank', 0, ...
                          'n_pas', n, 'R', Inf, 'v_min', v, 'v_max', v, ...
                          'resid_pos', 0, 'resid_vit', 0);
            return

        case 'tondeuse'
            P = polyligne_tondeuse(par, v);

        otherwise
            error('trajectoire:forme', 'Forme inconnue : %s.', par.traj.forme);
    end

    % Échantillonnage à abscisse curviligne constante : l'avion vole à vitesse constante,
    % donc un pas de temps est une distance.
    s = [0, cumsum(vecnorm(diff(P, 1, 2)))];
    [s, iu] = unique(s, 'stable');      % les jonctions de segments dupliquent des points
    P = P(:, iu);

    n  = min(par.mc.n_steps, floor(s(end) / (v * dt)));
    sk = (0:n) * v * dt;
    p  = [interp1(s, P(1, :), sk); interp1(s, P(2, :), sk)];

    % La vitesse est la tangente au chemin, prise par différences centrées sur la
    % polyligne fine, et la commande est l'accélération que cette vitesse implique. Le
    % bloc vitesse du modèle est alors exact.
    ds = v * dt / 100;
    tangente = @(q) ([interp1(s, P(1, :), min(q + ds, s(end)));
                      interp1(s, P(2, :), min(q + ds, s(end)))] - ...
                     [interp1(s, P(1, :), max(q - ds, 0));
                      interp1(s, P(2, :), max(q - ds, 0))]);
    T   = tangente(sk);
    vel = v * T ./ vecnorm(T, 2, 1);

    x_true = [p; vel];
    u  = diff(vel, 1, 2) / dt;
    Bu = [dt^2 / 2 * eye(2); dt * eye(2)] * u;

    % Ce que le modèle laisse sur la position. Un pas à accélération constante ne suit
    % pas un arc : l'écart vaut l'accélération centripète fois dt^2 sur deux, soit une
    % centaine de mètres à 1 g et dt = 5 s. C'est du même ordre que le bruit de modèle,
    % qui doit donc l'absorber — et c'est le prix d'un pas de temps long sur une
    % trajectoire qui tourne.
    F      = [eye(2), dt * eye(2); zeros(2), eye(2)];
    resid  = x_true(:, 2:end) - (F * x_true(:, 1:end-1) + Bu);

    info = struct('forme', par.traj.forme, 'chemin', s(end), ...
                  'bank', atand(max(vecnorm(u, 2, 1)) / 9.81), 'n_pas', n, ...
                  'R', v^2 / (9.81 * tand(par.traj.bank)), ...
                  'v_min', min(vecnorm(x_true(3:4, :), 2, 1)), ...
                  'v_max', max(vecnorm(x_true(3:4, :), 2, 1)), ...
                  'resid_pos', max(vecnorm(resid(1:2, :), 2, 1)), ...
                  'resid_vit', max(vecnorm(resid(3:4, :), 2, 1)));
end

function P = polyligne_tondeuse(par, v)
% La tondeuse comme polyligne fine : branche, demi-tour, branche. Le rayon vient de
% l'inclinaison demandée, R = v^2 / (g tan phi), et les demi-tours débordent des branches
% de ce rayon, ce dont il faut tenir compte pour rester dans la carte.
    R  = v^2 / (9.81 * tand(par.traj.bank));
    L  = par.traj.extent;
    lo = par.traj.marge + R;
    hi = L(2) - par.traj.marge - R;
    x1 = (L(1) - (par.traj.n_branches - 1) * 2 * R) / 2;

    if hi <= lo || x1 <= par.traj.marge
        error('trajectoire:tondeuse', ...
              ['%d branches à %.0f deg demandent %.0f km de large et %.0f km de haut, ' ...
               'la carte fait %.0f x %.0f km.'], par.traj.n_branches, par.traj.bank, ...
              ((par.traj.n_branches - 1) * 2 * R + 2 * par.traj.marge) / 1e3, ...
              (2 * R + 2 * par.traj.marge) / 1e3, L(1) / 1e3, L(2) / 1e3);
    end

    P = [];
    for i = 1:par.traj.n_branches
        x = x1 + (i - 1) * 2 * R;
        y = linspace(lo, hi, 2000);
        if mod(i, 2) == 0
            y = fliplr(y);
        end
        P = [P, [x * ones(1, 2000); y]];                             %#ok<AGROW>

        if i < par.traj.n_branches
            if mod(i, 2) == 1
                th = linspace(pi, 0, 800);     yc = hi;   % demi-tour par le haut
            else
                th = linspace(pi, 2 * pi, 800); yc = lo;  % et par le bas
            end
            P = [P, [x + R + R * cos(th); yc + R * sin(th)]];        %#ok<AGROW>
        end
    end
end

function par = poser(par, chemin, valeur)
% Écrit `valeur` dans par au bout d'un chemin pointé, par exemple 'krig.n_side'. C'est ce
% qui permet à une campagne de nommer le paramètre qu'elle balaie plutôt que de le coder.
    p = split(string(chemin), '.');
    switch numel(p)
        case 1
            par.(p{1}) = valeur;
        case 2
            par.(p{1}).(p{2}) = valeur;
        otherwise
            error('campagne_chapitre_4:poser', 'Chemin non géré : %s.', chemin);
    end
end
