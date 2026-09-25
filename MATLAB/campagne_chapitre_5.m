% CAMPAGNE_CHAPITRE_5
%
% Navigation magnétique quand il n'y a PAS de relevé : la carte est estimée en même temps
% que l'état. Une campagne croise une liste de méthodes sur le vol et le champ du chapitre
% 4, pour que la comparaison porte sur les algorithmes et non sur l'environnement.
%
% QUATRE FICHIERS : slam.m tient les deux algorithmes du chapitre 5, navigation.m et
% krigeage.m ceux du chapitre 4 dont on tire les bornes, celui-ci les paramètres, la
% vérité, la campagne, les tableaux et les figures.
%
% LES BORNES SONT DU CHAPITRE 4, ET C'EST VOULU. Un SLAM ne se juge pas dans l'absolu :
%
%   'carte'     le filtre du chapitre 4 sur le champ vrai. Ce qu'on obtiendrait si la
%               carte était parfaitement connue, donc la borne haute.
%   'ok'        le même sur un krigeage du relevé. Ce que coûte de reconstruire la carte
%               AVANT la mission, à comparer à ce que coûte de la reconstruire pendant.
%   'inertiel'  le même nuage sans aucune correction. La borne basse, sans laquelle un
%               chiffre de SLAM ne se lit pas.
%   'rbpf'      une carte par particule, corrigée analytiquement le long de sa trajectoire.
%   'akkf'      une seule carte, apprise par moindres carrés récursifs sur les estimées.
%
% Les quatre tournent au MÊME nombre de particules, ce qui n'est pas gratuit : l'AKKF
% inverse trois systèmes N x N par pas, donc N est celui qu'il supporte, et non les 5000
% du chapitre 4. Une comparaison à N différents ne dirait rien sur les méthodes.
%
% L'INCERTITUDE INITIALE EST CELLE DU SLAM, PAS CELLE DU CHAPITRE 4. Huit kilomètres
% là-bas, un ici, et c'est le paramètre le plus structurant du chapitre.
%
% Un filtre qui dispose d'une carte RECONNAÎT un motif : un nuage large lui est favorable,
% puisqu'il couvre les positions candidates et que la carte tranche entre elles — c'est la
% phase d'acquisition que le chapitre 4 mesure sur ses trente premiers pas. Un SLAM n'a
% rien à reconnaître. Il construit la carte depuis ce qu'il voit, donc toute trajectoire
% supposée engendre une carte qui explique les observations, et l'ambiguïté initiale n'est
% pas levée : elle est recopiée dans la carte. Ce qu'un SLAM peut faire est empêcher
% l'erreur de CROÎTRE, pas la réduire, et c'est ce que la campagne doit lire.
%
% La valeur vient du rapport qu'observent les deux codes de référence, sigma_0 valant
% environ un cinquième de la longueur de corrélation dans les deux — 100 m pour un noyau
% de 500 m dans le SLAM par RBPF, 0,3 pour un noyau de 1,5 dans celui par RKHS. Ici la
% longueur de corrélation vaut 5,5 km, d'où le kilomètre.
%
% C'est aussi pourquoi 'inertiel' est dans la liste des méthodes : sans carte du tout,
% l'erreur croît, et la mesure du chapitre est l'écart entre cette croissance et ce que
% les deux SLAM en retiennent.
%
% DEUX OBJETS MESURÉS, contre un seul au chapitre 4 : l'erreur de trajectoire, et l'erreur
% de reconstruction de la carte sur le corridor du vol. La seconde est la moitié du
% problème, et c'est elle qui montre où la carte apprise est bonne — le long de ce qui a
% été survolé — et où elle ne l'est pas.

clearvars -except over; close all; clc;   % `over` survit : voir surcharge() en fin de fichier
here = fileparts(mfilename('fullpath'));
addpath(here);

maxNumCompThreads(feature('numcores') + 2);   % voir campagne_chapitre_4, même raison

%% Paramètres
par.compare.methodes = {'carte', 'ok', 'inertiel', 'rbpf', 'akkf'};

par.scenario   = 'boucle';          % 'boucle', 'reference' ou 'court'
par.campagne   = 'A';               % 'A', 'B', 'C', 'D' ou 'libre'
par.run.seed   = 123456789;
par.mc.n_runs  = 100;

par.unit    = 'nT';
par.filt.alpha_reg = 0.3;
par.conv.seuil     = 3500;       % m, au-delà la piste est perdue
par.conv.confiance = 0.99;
par.conv.n_fin     = 10;
par.conv.k_acq     = 30;

% LA CARTE QUE LE SLAM SE DONNE. Il ne connaît pas le champ, mais il connaît le noyau qui
% le décrit : c'est l'hypothèse de stationnarité du chapitre, et elle remplace le relevé.
% Les valeurs sont celles que le chapitre 4 ajuste par maximum de vraisemblance sur son
% relevé — les donner ici revient à supposer qu'une mission antérieure a mesuré la
% STRUCTURE du champ sans en mesurer les valeurs, ce qui est l'hypothèse la plus favorable
% qu'on puisse faire au SLAM. Les dégrader est un balayage à part (campagne D).
par.slam.noyau   = 'matern32';   % le même qu'au chapitre 4
par.slam.sigma_f = 124;          % nT, écart-type a priori du champ
par.slam.ell     = 4e3;          % m, longueur de corrélation. Pour le champ analytique,
                                 % dont l'autocorrélation tombe à un demi vers lambda/6 ;
                                 % sur l'anomalie du chapitre 4, mettre les 5,5 km que
                                 % fitrgp y ajuste. Entre 3 et 6 km la base retenue ne
                                 % change pratiquement pas, le tri par densité spectrale
                                 % rendant presque le même disque de fréquences.
par.slam.r_kl    = 225;          % rang de la base de Laplace du RBPF. Mesuré : le long du
                                 % vol, la projection du champ laisse 11 nT à r = 225,
                                 % 6 nT à 400 et 5 nT à 625, contre 5 nT de bruit capteur.
                                 % Le résidu de troncature reste donc ici du même ordre que
                                 % le bruit du capteur, ce dont la covariance d'innovation
                                 % tient compte — voir innovation() dans slam.m.
                                 %
                                 % C'EST LA MÉMOIRE QUI TRANCHE, et elle est serrée. Les N
                                 % covariances r x r pèsent 8 r^2 N octets, et le
                                 % rééchantillonnage en fait une copie puisque les indices
                                 % se répètent : 1,3 Go au rang 400 sur cinq cents
                                 % particules, 670 au rang 289, 400 au rang 225. La machine
                                 % a quatre gigaoctets de libres, mais ils fluctuent avec
                                 % ce qui tourne à côté, et un tableau de cette taille
                                 % demande de surcroît de la mémoire CONTIGUË qu'une
                                 % session fragmentée n'a plus. Les rangs 400 et 289 ont
                                 % tous deux échoué en cours de campagne, à des endroits
                                 % différents ; 225 laisse la marge qu'il faut.
                                 %
                                 % Deux façons de remonter en rang sans plus de mémoire,
                                 % aucune essayée : stocker les covariances en simple
                                 % précision, ou n'en porter qu'une seule pour tout le
                                 % nuage. La seconde se défend ici, un nuage d'un kilomètre
                                 % étant petit devant une longueur de corrélation de
                                 % quatre : deux particules ont alors vu les mêmes lieux à
                                 % la résolution de la carte près, donc presque la même
                                 % covariance. Ce ne serait plus le RBPF exact.
par.slam.r_rff   = 500;          % features du RFF de l'AKKF, celui de l'article
par.slam.marge   = 25e3;         % m, marge entre le vol et le bord du domaine de Dirichlet,
                                 % où toute la base s'annule. Trois longueurs de
                                 % corrélation : en deçà le champ reconstruit est tiré vers
                                 % zéro là où l'avion vole, au delà le spectre se raffine
                                 % et il faut plus de modes pour la même résolution.
par.slam.mem_max    = 2.5 * 2^30;   % octets, garde-fou sur les N covariances r x r. La
                                    % machine a 16 Go dont quatre libres, et un dépassement
                                    % part en swap plutôt qu'en erreur.
par.slam.bloc_octets = 64 * 2^20;  % taille du bloc de particules dans la correction des
                                    % cartes. Assez grand pour que pagemtimes travaille,
                                    % assez petit pour que le temporaire ne compte pas.

% L'AKKF. Son noyau d'état porte une distance sur un vecteur qui mélange des mètres et des
% mètres par seconde : les composantes sont donc divisées par leurs échelles avant, et
% c'est sigma_0 qui les donne — l'ordre de grandeur du nuage au moment où il est le plus
% large. L'échelle du noyau d'observation, elle, ne peut pas être fixée d'avance et suit
% la dispersion des observations simulées.
par.akkf.gamma_x   = 1.0;            % échelle du noyau d'état, en unités de sigma_0
par.akkf.echelle_z = 'mediane';      % 'mediane' ou une valeur en nT
par.akkf.lambda_x  = 0.1;            % régularisations du code de référence
par.akkf.lambda_z  = 0.1;
par.akkf.oubli     = 0.995;          % facteur d'oubli du RLS
par.akkf.echelle_p = 1e2;            % covariance initiale des coefficients de carte

% Le relevé n'est lu que par la borne 'ok'. Les valeurs sont celles du chapitre 4.
par.krig.jitter       = 0.30;
par.krig.phase_alea   = true;
par.krig.n_side       = 80;
par.krig.ell_init     = 10e3;
par.krig.n_ml         = 500;
par.krig.n_ml_win     = 250;
par.krig.n_chunk      = 2000;
par.krig.refit        = 'fenetre';
par.krig.window_factor = 2;
par.krig.n_min        = 25;
par.krig.n_surveys    = 1;
par.krig.pas          = 3500;
par.krig.sigma_map    = 10;
par.krig.tab_pitch    = 250;
par.krig.n_max        = 2000;
par.krig.bandwidth    = 1000;
par.krig.noyau        = 'matern32';

% LE CHAMP. Deux sources, et le chapitre commence par la seconde.
%
%   'fichier'     l'anomalie magnétique du chapitre 4. Elle a de l'énergie à toutes les
%                 échelles, ce qui est sa qualité pour une navigation sur carte connue et
%                 son défaut pour un SLAM : une base tronquée en laisse toujours une part,
%                 et sur le domaine que le vol impose il faudrait 1600 modes pour en
%                 représenter 92 %, soit 38 Go de covariances.
%   'analytique'  la somme de trois sinusoïdes du code RKHS, remise à l'échelle de la
%                 trajectoire. Son spectre tient en trois raies, donc une base de rang
%                 modeste le représente presque exactement et ce qui reste mesure les
%                 ALGORITHMES plutôt que la troncature. C'est le champ sur lequel les deux
%                 méthodes ont été mises au point, et le point de départ raisonnable.
%
% La remise à l'échelle est le seul ajustement. Le champ de l'article vit sur dix unités
% avec des longueurs d'onde de cinq à onze ; il est ici dilaté pour que sa longueur d'onde
% principale vaille par.map.lambda, et son amplitude pour que son écart-type vaille celui
% de l'anomalie. Tout le reste du pilote le lit comme une carte tabulée ordinaire.
par.map.source = 'analytique';   % 'analytique' ou 'fichier'
par.map.lambda = 25e3;           % m, longueur d'onde principale du champ analytique
par.map.sigma  = 124;            % nT, son écart-type
par.map.name   = 'carte_magnetometrie_anomalie_mexique.mat';
par.map.file   = fullfile(here, '..', '..', 'MATLAB 2A - sauvegarde', 'Cartes', par.map.name);
par.path.cache = fullfile(here, 'cache');
par.path.out   = fullfile(here, 'resultats');

par = surcharge(par);               % pour que `over.scenario` soit lu avant le switch

switch par.scenario
    % Le vol du chapitre 4, à ceci près que le nuage est plus petit. Tout le reste — la
    % carte, la tondeuse, les bruits, la géométrie — est identique, pour que la
    % comparaison avec les résultats du chapitre 4 tienne.
    case 'reference'
        par.mc.n_steps = 150;
        par.mc.n_part  = 1000;           % PAS 5000 : l'AKKF est cubique en N. À 5000 un
                                         % essai demande des heures, à 1000 il demande des
                                         % secondes, et le chapitre 4 a mesuré que le RPF
                                         % tient encore à 1250 particules.
        par.dt         = 5;
        par.sigma_obs  = 5;              % nT
        par.x0         = [0; 0; 0; 500];
        par.traj.forme      = 'tondeuse';
        par.traj.n_branches = 3;
        par.traj.bank       = 45;
        par.traj.marge      = 35e3;
        par.sigma_0    = [1000, 1000, 2, 2];   % PAS les 8 km du chapitre 4, voir plus haut
        par.sigma_q    = [100, 100, 1, 1];

    % LE VOL QUI DONNE AU SLAM DE QUOI TRAVAILLER. Même carte, même véhicule, mêmes
    % bruits que 'reference', mais un circuit fermé reparcouru : voir polyligne_boucle en
    % fin de fichier pour ce que le recoupement change. La tondeuse de 'reference' reste
    % là comme témoin, et la comparaison des deux mesure exactement ce que la fermeture de
    % boucle apporte.
    case 'boucle'
        par.mc.n_steps = 300;            % quatre tours du circuit
        par.mc.n_part  = 500;            % ce que la mémoire du RBPF autorise au rang 400,
                                         % voir slam.r_kl. Le nuage fait ici un kilomètre
                                         % et non huit, donc cinq cents particules y sont
                                         % aussi denses que les cinq mille du chapitre 4.
        par.dt         = 5;
        par.sigma_obs  = 5;
        par.x0         = [0; 0; 0; 500];
        par.traj.forme   = 'boucle';
        par.traj.bank    = 45;
        par.traj.marge   = 20e3;
        par.traj.largeur = 60e3;
        par.traj.hauteur = 60e3;
        par.sigma_0    = [1000, 1000, 2, 2];
        par.sigma_q    = [100, 100, 1, 1];

    % Une branche unique et soixante pas, pour mettre au point sans payer la tondeuse
    % entière. La carte n'y est apprise que le long d'une ligne, donc l'erreur de
    % reconstruction hors trajectoire y est de l'extrapolation pure, et aucun recoupement
    % ne vient contraindre la position : c'est un cas de mise au point, pas de mesure.
    case 'court'
        par.mc.n_steps = 60;
        par.mc.n_part  = 500;
        par.dt         = 5;
        par.sigma_obs  = 5;
        par.x0         = [80e3; 25e3; 50; 500];
        par.traj.forme = 'ligne';
        par.sigma_0    = [1000, 1000, 2, 2];
        par.sigma_q    = [100, 100, 1, 1];

    otherwise
        error('campagne_chapitre_5:scenario', 'Scénario inconnu : %s.', par.scenario);
end

switch par.campagne
    case 'A'    % LA comparaison : les deux SLAM entre les deux bornes.
        par.compare.methodes = {'carte', 'ok', 'inertiel', 'rbpf', 'akkf'};
        par.mc.n_runs = 100;
        par.balayage  = struct('champ', '', 'valeurs', 0);

    case 'B'    % le rang de la base du RBPF. C'est le paramètre qui décide à la fois de
                % la résolution de la carte et de la mémoire, et les deux se contredisent.
        par.compare.methodes = {'rbpf'};
        par.mc.n_runs = 100;
        par.balayage  = struct('champ', 'slam.r_kl', 'valeurs', [49 100 169 256 400]);

    case 'C'    % les particules, pour comparer les deux SLAM à coût égal plutôt qu'à
                % effectif égal — le coût de l'AKKF croît en N^3 et celui du RBPF en N.
        par.compare.methodes = {'rbpf', 'akkf'};
        par.mc.n_runs = 100;
        par.balayage  = struct('champ', 'mc.n_part', 'valeurs', [250 500 1000 2000]);

    case 'D'    % la longueur de corrélation donnée au SLAM, la vraie valant par.slam.ell.
                % Mesure ce que coûte de se tromper sur la structure du champ, qui est la
                % seule chose que le SLAM reçoive à la place d'un relevé.
        par.compare.methodes = {'rbpf', 'akkf'};
        par.mc.n_runs = 100;
        par.balayage  = struct('champ', 'slam.ell', 'valeurs', [4e3 6e3 8e3 12e3 16e3]);

    case 'libre'
        par.balayage = struct('champ', '', 'valeurs', 0);

    otherwise
        error('campagne_chapitre_5:campagne', 'Campagne inconnue : %s.', par.campagne);
end

par = surcharge(par);                % la surcharge a le dernier mot sur les deux switch

%% Carte et vérité
map = carte(par);

d_x = numel(par.x0);
par.traj.extent = [size(map.h, 2), size(map.h, 1)] .* map.step;

scen.F  = [eye(2), par.dt * eye(2); zeros(2), eye(2)];
scen.Q  = diag(par.sigma_q).^2;
scen.P0 = diag(par.sigma_0).^2;
[scen.x_true, scen.Bu, traj] = trajectoire(par);
n_iter = traj.n_pas;

krig   = krigeage();
nav    = navigation();
sl     = slam();
z_true = krig.lire(map.h, map.step, scen.x_true(1:2, 2:end));

extent = [size(map.h, 2), size(map.h, 1)] .* map.step;
margin = min([min(scen.x_true(1:2, :), [], 2)', extent - max(scen.x_true(1:2, :), [], 2)']);
if margin < 0
    error('campagne_chapitre_5:trajectoire', 'La trajectoire sort de la carte de %.0f km.', ...
          -margin / 1e3);
end

pad = max(15e3, 3.5 * max(par.sigma_0(1:2)));
par.krig.corridor = [max(min(scen.x_true(1, :)) - pad, map.step(1)), ...
                     min(max(scen.x_true(1, :)) + pad, extent(1) - map.step(1)), ...
                     max(min(scen.x_true(2, :)) - pad, map.step(2)), ...
                     min(max(scen.x_true(2, :)) + pad, extent(2) - map.step(2))];

% LA MOYENNE DU CHAMP EST RETIRÉE. Les deux algorithmes écrivent la carte comme un
% processus centré sur une base finie ; un champ de moyenne non nulle devrait dépenser des
% coefficients à la représenter, et la base de Dirichlet ne le peut pas du tout puisqu'elle
% s'annule au bord. On retire donc la moyenne du corridor à la fois du champ vrai et des
% observations, ce qui revient à supposer cette moyenne connue — une hypothèse de plus
% donnée au SLAM, du même ordre que celle du noyau, et à déclarer comme telle.
[ix_c, iy_c] = corridor_indices(map, par.krig.corridor);
par.slam.moyenne = mean(map.h(iy_c, ix_c), 'all');
fprintf('Champ : moyenne %.1f %s retirée, écart-type %.1f %s sur le corridor\n', ...
        par.slam.moyenne, par.unit, std(map.h(iy_c, ix_c), 0, 'all'), par.unit);

fprintf('Carte %d x %d noeuds à %.0f m\n', size(map.h, 2), size(map.h, 1), map.step(1));
parcouru = n_iter * norm(par.x0(3:4)) * par.dt;
fprintf('Vol %s : %d pas, parcouru %.0f km, marge %.0f km\n', ...
        traj.forme, n_iter, parcouru / 1e3, margin / 1e3);

% LE DOMAINE DE DIRICHLET, et la grille sur laquelle les cartes reconstruites sont jugées.
% Le premier est le vol dilaté de par.slam.marge ; la seconde est ce même rectangle
% échantillonné assez fin pour que l'erreur de reconstruction s'y lise, mais pas plus : ce
% n'est qu'une moyenne spatiale.
domaine = [min(scen.x_true(1, :)) - par.slam.marge, max(scen.x_true(1, :)) + par.slam.marge, ...
           min(scen.x_true(2, :)) - par.slam.marge, max(scen.x_true(2, :)) + par.slam.marge];
[Ex, Ey] = meshgrid(linspace(domaine(1), domaine(2), 60), ...
                    linspace(domaine(3), domaine(4), 60));
P_eval   = [Ex(:), Ey(:)]';
h_eval   = krig.lire(map.h, map.step, P_eval) - par.slam.moyenne;

% Le long de la trajectoire, où la carte a réellement été apprise, contre partout, où elle
% est extrapolée. Les deux chiffres ne disent pas la même chose et le second seul ne dirait
% rien sur la méthode.
pres_traj = min(pdist2(P_eval', scen.x_true(1:2, :)'), [], 2) <= 2 * par.slam.ell;

fprintf('Domaine SLAM %.0f x %.0f km, %d points d''évaluation dont %d le long du vol\n\n', ...
        (domaine(2) - domaine(1)) / 1e3, (domaine(4) - domaine(3)) / 1e3, ...
        numel(pres_traj), sum(pres_traj));

%% Campagne, une passe par valeur balayée
res = struct('name', {}, 'methode', {}, 'balayage', {}, 'rmse', {}, 'nees', {}, ...
             'track', {}, 'wall', {}, 'err', {}, 'd2', {}, 'rmse_map', {}, ...
             'armse_map', {}, 'armse_map_traj', {}, 'cout', {}, 'n_res', {});
valeurs = par.balayage.valeurs;
nM      = numel(par.compare.methodes);

for b = 1:numel(valeurs)
if ~isempty(par.balayage.champ)
    par = poser(par, par.balayage.champ, valeurs(b));
    fprintf('\n===== %s = %g =====\n', par.balayage.champ, valeurs(b));
end

opt = struct('n_part', par.mc.n_part, 'n_th', 0.5 * par.mc.n_part, ...
             'alpha_reg', par.filt.alpha_reg, 'R_obs', par.sigma_obs^2, ...
             'P_eval', P_eval, 'mem_max', par.slam.mem_max, 'akkf', par.akkf, ...
             'sigma_f2', par.slam.sigma_f^2, 'bloc_octets', par.slam.bloc_octets);
opt.akkf.echelle_x = par.sigma_0(:);

% Le relevé et le modèle krigé ne servent qu'aux bornes du chapitre 4, et ne sont
% construits que si l'une d'elles est demandée — la tabulation de 'ok' coûte plus cher que
% toute la campagne SLAM.
par.krig.survey_seed = par.run.seed + 1000;
survey = krig.releve(map, par);
modeles = containers.Map();
for m = 1:nM
    k = par.compare.methodes{m};
    if any(strcmp(k, {'carte', 'ok', 'bilineaire'}))
        modeles(k) = krig.modele(k, map, survey, par);
    end
end

err  = zeros(n_iter, par.mc.n_runs, nM);
d2   = zeros(n_iter, par.mc.n_runs, nM);
emap = zeros(numel(pres_traj), par.mc.n_runs, nM);
hmap = zeros(numel(pres_traj), nM);        % la carte moyenne, pour la figure
co   = zeros(4, n_iter, par.mc.n_runs, nM);
nres = zeros(par.mc.n_runs, nM);
track = zeros(2, n_iter, nM);
wall  = zeros(1, nM);

for m = 1:nM
    kind = par.compare.methodes{m};
    t0   = tic;

    % La base est tirée une fois par méthode et non par essai : les fréquences du RFF sont
    % un choix d'algorithme, pas un aléa de mission, et les refaire à chaque essai
    % mélangerait leur variance à celle du filtre. La campagne B, qui balaie le rang, est
    % l'endroit où cette variance se mesure.
    switch kind
        case 'rbpf',     base = sl.base('laplace', par, domaine);
        case 'akkf',     base = sl.base('rff', par, domaine);
        case 'inertiel', base = struct('name', 'aucune', 'r', 0);
        otherwise,       base = [];
    end
    if ~isempty(base) && isfield(base, 'variance_tronquee')
        fprintf('%-6s base %s, rang %d\n', kind, base.name, base.r);
    end

    for r = 1:par.mc.n_runs
        rng(par.run.seed + r);
        z_obs = z_true - par.slam.moyenne + par.sigma_obs * randn(n_iter, 1);

        if isempty(base)
            % Les bornes du chapitre 4 lisent le champ non centré : leur modèle porte la
            % carte telle quelle, donc on leur rend la moyenne qu'on a retirée.
            out = nav.run('RPF', modeles(kind), z_obs + par.slam.moyenne, scen, opt);
            out.h_eval = h_eval;                 % par construction, leur carte est exacte
        else
            out = sl.run(kind, base, z_obs, scen, opt);
        end

        err(:, r, m)   = out.err;
        d2(:, r, m)    = out.d2;
        emap(:, r, m)  = out.h_eval - h_eval;
        hmap(:, m)     = hmap(:, m) + out.h_eval / par.mc.n_runs;
        co(:, :, r, m) = out.cout;
        nres(r, m)     = out.n_res;
        if r == 1, track(:, :, m) = out.track; end
    end
    wall(m) = toc(t0);
end
fprintf('\n');

for m = 1:nM
    E = err(:, :, m);  D = d2(:, :, m);  M = emap(:, :, m);
    e = numel(res) + 1;
    res(e).name     = par.compare.methodes{m};
    res(e).methode  = par.compare.methodes{m};
    res(e).balayage = valeurs(b);
    res(e).rmse     = sqrt(mean(E.^2, 2));
    res(e).nees     = mean(D, 2, 'omitnan') / d_x;
    res(e).track    = track(:, :, m);
    res(e).cout     = mean(co(:, :, :, m), 3);
    res(e).n_res    = mean(nres(:, m)) / n_iter;
    res(e).wall     = wall(m);
    res(e).err      = E;
    res(e).d2       = D;

    % L'erreur de carte, fonction du POINT et non du pas : la racine est prise sur les
    % essais à point fixe, comme le RMSE du chapitre 4 l'est à pas fixe, et sa moyenne
    % spatiale est l'ARMSE.
    res(e).rmse_map       = sqrt(mean(M.^2, 2));
    res(e).armse_map      = mean(res(e).rmse_map);
    res(e).armse_map_traj = mean(res(e).rmse_map(pres_traj));
    res(e).h_map          = hmap(:, m);

    [c_s, c_m] = converges(E, D, d_x, par);
    fprintf(['%-8s ARMSE %6.0f m | seuil %6.0f (%3.0f %%) | méd. acq %5.0f, fin %5.0f' ...
             ' | carte %5.1f %s, le long du vol %5.1f | %6.0f s\n'], ...
            res(e).name, mean(res(e).rmse), ...
            mean(sqrt(mean(E(:, c_s).^2, 2))), 100 * mean(c_s), ...
            median(E(min(par.conv.k_acq, end), :)), median(E(end, :)), ...
            res(e).armse_map, par.unit, res(e).armse_map_traj, res(e).wall);
end
end                                        % fin du balayage

%% Résultats
etiq = {res.name};
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
stamp = sprintf('%s_ch5_%s_%s_%druns', char(datetime('now', 'Format', 'yyyyMMdd_HHmmss')), ...
                par.scenario, par.campagne, par.mc.n_runs);
save(fullfile(par.path.out, [stamp '.mat']), 'par', 'res', 'scen', 'P_eval', 'h_eval', ...
     'pres_traj', '-v7.3');
fprintf('\nSauvegardé dans resultats/%s.mat\n', stamp);

%% Figures
colours = [0.00 0.00 0.00; 0.15 0.35 0.70; 0.95 0.60 0.00; 0.85 0.00 0.75
           0.75 0.02 0.15; 0.00 0.52 0.28];
t = (1:n_iter) * par.dt;

figure('Color', 'w', 'Name', 'RMSE'); hold on;
for i = 1:numel(res)
    plot(t, res(i).rmse, 'Color', colours(mod(i-1, size(colours,1)) + 1, :), 'LineWidth', 1.2);
end
hold off; grid on; box on; set(gca, 'YScale', 'log');
xlabel('t (s)'); ylabel('RMSE (m)'); legend(etiq, 'Location', 'northeast');
title('Erreur de position');

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
legend([{'bande de cohérence 95 %'}, etiq], 'Location', 'southwest');
title('Cohérence');

% LA CARTE RECONSTRUITE, qui n'a pas d'équivalent au chapitre 4. Une colonne par méthode
% qui en apprend une, la vérité en premier, et la même échelle de couleur partout sans
% quoi deux cartes également fausses se ressemblent.
appris = find(~cellfun(@isempty, {res.rmse_map}) & ...
              ismember({res.methode}, {'rbpf', 'akkf'}));
if ~isempty(appris)
    figure('Color', 'w', 'Name', 'Cartes');
    n_col = numel(appris) + 1;
    clim_  = [-3 3] * std(h_eval);
    subplot(2, n_col, 1);
    afficher_carte(Ex, Ey, h_eval, clim_, 'vérité', scen, par);
    for j = 1:numel(appris)
        i = appris(j);
        subplot(2, n_col, 1 + j);
        afficher_carte(Ex, Ey, res(i).h_map, clim_, etiq{i}, scen, par);
        subplot(2, n_col, n_col + 1 + j);
        afficher_carte(Ex, Ey, res(i).rmse_map, [0 max(clim_)], ...
                       sprintf('%s : RMSE carte', etiq{i}), scen, par);
    end
end

figure('Color', 'w', 'Name', 'Trajectoires');
[ix, iy] = corridor_indices(map, par.krig.corridor);
imagesc(ix * map.step(1) / 1e3, iy * map.step(2) / 1e3, map.h(iy, ix));
set(gca, 'YDir', 'normal'); colormap(gca, gray); hold on;
cb = colorbar; cb.Label.String = sprintf('champ (%s)', par.unit);
plot(scen.x_true(1, :) / 1e3, scen.x_true(2, :) / 1e3, 'w-', 'LineWidth', 3.5);
plot(scen.x_true(1, :) / 1e3, scen.x_true(2, :) / 1e3, 'k-', 'LineWidth', 1.8);
for i = 1:numel(res)
    plot(res(i).track(1, :) / 1e3, res(i).track(2, :) / 1e3, '--', ...
         'Color', colours(mod(i-1, size(colours,1)) + 1, :), 'LineWidth', 1.2);
end
rectangle('Position', [domaine(1), domaine(3), domaine(2) - domaine(1), ...
                       domaine(4) - domaine(3)] / 1e3, 'EdgeColor', [0.8 0 0], ...
          'LineStyle', ':', 'LineWidth', 1.2);
hold off; axis image; box on;
xlabel('x (km)'); ylabel('y (km)');
legend([{'vérité'}, {''}, etiq], 'Location', 'northeast');
title('Un essai de chaque méthode, et le domaine de Dirichlet en pointillés');

%% -------------------------------------------------------------- utilitaires
%
% CE QUI SUIT EST REPRIS DE campagne_chapitre_4.m À L'IDENTIQUE — trajectoire, tondeuse,
% surcharge, poser, converges. Les deux pilotes doivent voler exactement le même vol, et
% une copie qui dérive serait pire qu'une duplication visible. Le jour où l'un des deux
% change, sortir ces cinq fonctions dans un fichier commun.

function map = carte(par)
% LE CHAMP, tabulé sur la même grille dans les deux cas pour que rien en aval n'ait à
% savoir d'où il vient.
%
% Le champ analytique est celui du code RKHS, dilaté d'un facteur `s` mètres par unité :
%
%   h(x,y) = A [ sin(k1 x) + 0.5 cos(k2 y) + 0.3 sin(k3 (x+y)) ],   ki = ci / s
%
% avec c = [0.7, 1.2, 0.4]. `s` est fixé par la longueur d'onde principale demandée,
% 2 pi / c1 unités valant par.map.lambda, et A par l'écart-type demandé — celui de la
% somme vaut sqrt(1/2 + 0.5^2/2 + 0.3^2/2) = 0.8185 fois A.
%
% LES TROIS RAIES SONT CE QUI COMPTE. Le champ n'a d'énergie qu'aux fréquences
% ||omega|| = c1/s, c2/s et c3 sqrt(2)/s, soit ici des longueurs d'onde de 25, 15 et 31
% km. Une base dont le spectre les couvre le représente à la précision de ses conditions
% au bord près, là où l'anomalie magnétique en demande dix fois plus.
    c = [0.7, 1.2, 0.4];

    switch par.map.source
        case 'fichier'
            f   = load(par.map.file);
            map = struct('h', double(f.mesures).', 'step', double(f.pas(:))' .* [1 1]);

        case 'analytique'
            % La même grille que la carte du chapitre 4, pour que le corridor, la
            % trajectoire et les marges se comparent sans conversion.
            step = [1000, 1000];
            n    = 201;
            [Gx, Gy] = meshgrid((1:n) * step(1), (1:n) * step(2));

            s = par.map.lambda * c(1) / (2 * pi);
            A = par.map.sigma / sqrt(0.5 * (1 + 0.5^2 + 0.3^2));
            h = A * (sin(c(1) * Gx / s) + 0.5 * cos(c(2) * Gy / s) + ...
                     0.3 * sin(c(3) * (Gx + Gy) / s));
            map = struct('h', h, 'step', step);

        otherwise
            error('campagne_chapitre_5:map', 'Source inconnue : %s.', par.map.source);
    end
end

function afficher_carte(Ex, Ey, v, clim_, titre, scen, ~)
    imagesc(Ex(1, :) / 1e3, Ey(:, 1) / 1e3, reshape(v, size(Ex)));
    set(gca, 'YDir', 'normal'); clim(clim_); axis image; hold on;
    plot(scen.x_true(1, :) / 1e3, scen.x_true(2, :) / 1e3, 'k-', 'LineWidth', 1);
    hold off; title(titre); xlabel('x (km)'); ylabel('y (km)'); colorbar;
end

function [ix, iy] = corridor_indices(map, corridor)
    ix = max(floor(corridor(1) / map.step(1)), 1):floor(corridor(2) / map.step(1));
    iy = max(floor(corridor(3) / map.step(2)), 1):floor(corridor(4) / map.step(2));
end

function par = surcharge(par)
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

function [x_true, Bu, info] = trajectoire(par)
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
            info = struct('forme', 'ligne', 'chemin', v * n * dt, 'bank', 0, ...
                          'n_pas', n, 'R', Inf, 'v_min', v, 'v_max', v, ...
                          'resid_pos', 0, 'resid_vit', 0);
            return

        case 'tondeuse'
            P = polyligne_tondeuse(par, v);

        case 'boucle'
            P = polyligne_boucle(par, v);

        otherwise
            error('trajectoire:forme', 'Forme inconnue : %s.', par.traj.forme);
    end

    s = [0, cumsum(vecnorm(diff(P, 1, 2)))];
    [s, iu] = unique(s, 'stable');
    P = P(:, iu);

    n  = min(par.mc.n_steps, floor(s(end) / (v * dt)));
    sk = (0:n) * v * dt;
    p  = [interp1(s, P(1, :), sk); interp1(s, P(2, :), sk)];

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
    R  = v^2 / (9.81 * tand(par.traj.bank));
    L  = par.traj.extent;
    lo = par.traj.marge + R;
    hi = L(2) - par.traj.marge - R;
    x1 = (L(1) - (par.traj.n_branches - 1) * 2 * R) / 2;

    if isfield(par.traj, 'decalage')
        x1 = x1 + par.traj.decalage(1);
        lo = lo + par.traj.decalage(2);
        hi = hi + par.traj.decalage(2);
    end

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
                th = linspace(pi, 0, 800);     yc = hi;
            else
                th = linspace(pi, 2 * pi, 800); yc = lo;
            end
            P = [P, [x + R + R * cos(th); yc + R * sin(th)]];        %#ok<AGROW>
        end
    end
end

function P = polyligne_boucle(par, v)
% UN CIRCUIT FERMÉ, PARCOURU PLUSIEURS FOIS. C'est la seule forme de ce fichier qui ne
% vienne pas du chapitre 4, et elle y est parce que le chapitre 5 ne mesure rien sans
% elle.
%
% POURQUOI. Un capteur qui lit le champ à sa position n'apprend, à lui seul, rien sur
% cette position : la carte s'ajuste à ce qu'il mesure quelle que soit la trajectoire
% supposée, et une trajectoire fausse portant une carte fausse explique les observations
% aussi bien que le couple vrai. Ce qui brise cette indétermination est le RECOUPEMENT :
% quand le véhicule repasse à portée d'un endroit déjà mesuré, une trajectoire fausse
% prédit une valeur que la carte qu'elle s'est construite contredit. La tondeuse du
% chapitre 4 n'en offre aucun, ses branches étant espacées d'un diamètre de virage, soit
% dix fois la longueur de corrélation du champ : elle est excellente pour couvrir une
% zone, et aveugle pour un SLAM.
%
% Le circuit est un rectangle à coins arrondis du rayon de virage, donc ses côtés font au
% moins deux rayons. Il est reparcouru autant de fois que les pas demandés le permettent,
% chaque tour repassant exactement sur le précédent.
    R  = v^2 / (9.81 * tand(par.traj.bank));
    L  = par.traj.extent;
    W  = max(par.traj.largeur, 2 * R);
    H  = max(par.traj.hauteur, 2 * R);
    c  = L / 2 + [0, 0];                       % centré sur la carte
    if isfield(par.traj, 'decalage')
        c = c + par.traj.decalage;
    end

    marge = min([c - [W, H] / 2, L - c - [W, H] / 2]);
    if marge < par.traj.marge
        error('trajectoire:boucle', ...
              ['un circuit %.0f x %.0f km laisse %.0f km de marge, il en faut %.0f. ' ...
               'Réduire traj.largeur et traj.hauteur, ou traj.marge.'], ...
              W / 1e3, H / 1e3, marge / 1e3, par.traj.marge / 1e3);
    end

    % Un tour, dans le sens direct : côté bas, coin, côté droit, coin, etc.
    x0 = c(1) - W / 2;  x1 = c(1) + W / 2;
    y0 = c(2) - H / 2;  y1 = c(2) + H / 2;
    seg = @(a, b) [linspace(a(1), b(1), 500); linspace(a(2), b(2), 500)];
    arc = @(ct, th0) [ct(1) + R * cos(linspace(th0, th0 + pi/2, 300));
                      ct(2) + R * sin(linspace(th0, th0 + pi/2, 300))];

    tour = [seg([x0 + R, y0], [x1 - R, y0]), arc([x1 - R, y0 + R], -pi/2), ...
            seg([x1, y0 + R], [x1, y1 - R]), arc([x1 - R, y1 - R], 0), ...
            seg([x1 - R, y1], [x0 + R, y1]), arc([x0 + R, y1 - R], pi/2), ...
            seg([x0, y1 - R], [x0, y0 + R]), arc([x0 + R, y0 + R], pi)];

    % Assez de tours pour les pas demandés, un de plus pour ne pas s'arrêter court.
    perimetre = 2 * (W - 2 * R) + 2 * (H - 2 * R) + 2 * pi * R;
    n_tours   = ceil(par.mc.n_steps * v * par.dt / perimetre) + 1;
    P = repmat(tour, 1, n_tours);
end

function par = poser(par, chemin, valeur)
    p = split(string(chemin), '.');
    switch numel(p)
        case 1
            par.(p{1}) = valeur;
        case 2
            par.(p{1}).(p{2}) = valeur;
        otherwise
            error('campagne_chapitre_5:poser', 'Chemin non géré : %s.', chemin);
    end
end

function [c_seuil, c_mahal] = converges(err, d2, d_x, par)
    k       = size(err, 1) - par.conv.n_fin + 1 : size(err, 1);
    c_seuil = median(err(k, :), 1) <= par.conv.seuil;
    c_mahal = median(d2(k, :),  1) <= chi2inv(par.conv.confiance, d_x);
end
