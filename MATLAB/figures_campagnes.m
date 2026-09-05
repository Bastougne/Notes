% FIGURES_CAMPAGNES
%
% Les figures du chapitre 4 tirees des campagnes enregistrees. Un seul script, pour que
% les courbes du manuscrit aient une source unique et suivent d'elles-memes le passage a
% 1000 essais : il relit la campagne la plus recente, sans parametre a retoucher.
%
%   RMSE par maille   grille 3 x 2, une maille par case, carte connue / krigeage /
%                     bilineaire sur le meme axe
%   NEES par maille   la meme grille, echelle log et bande du chi2
%   Tableau           taux de convergence sous les deux criteres, sous include/
%
% Les six mailles tiennent en une grille plutot qu'en six flottants : la degradation se
% lit d'une case a l'autre, ce qu'une figure isolee par maille ne montre pas. La page est
% pleine largeur et haute, donc \includegraphics[width=\textwidth] sur une page entiere.
%
% Ce script ne porte que les courbes VALIDEES. figures_temporelles.m garde pour l'instant
% celles qui ne le sont pas encore (modes, fenetre, coherence) ; elles viendront ici une
% a une, et figures_maille.m disparaitra quand ses deux figures y seront passees.

clear; close all; clc;
here = fileparts(mfilename('fullpath'));

%% Ce que les figures declarent de leur insertion
img_dir   = fullfile(fileparts(here), 'Images');
tab_dir   = fullfile(fileparts(here), 'Latex', 'Manuscrit', 'include');
textwidth = 425;                 % pt, mesure sur le document
tex_frac  = 1.0;                 % \includegraphics[width=\textwidth]
pdf_w     = textwidth * tex_frac;
pdf_h     = 580;                 % haute : la grille occupe une page entiere
fontsize  = 12 * pdf_w / (tex_frac * textwidth);   % donc 12 pt, la figure n'est pas mise a l'echelle
linewidth = 1;
pas_trace = 3;                   % un point sur trois sur les courbes log : a 150 pas le
                                 % NEES est trop bruite pour qu'on distingue les modeles
seuil_conv = 2000;               % m, ou [] pour reprendre celui de la campagne.
                                 % 2000 et non les 3500 de campagne_chapitre_4.m : a 3500
                                 % le krigeage est plat a 99 % sur quatre mailles, donc le
                                 % critere ne separe plus rien. A 2000 la degradation
                                 % apparait (99 -> 71) et l'ecart au bilineaire tient a
                                 % chaque maille. Plus bas, le critere cesse de mesurer la
                                 % convergence pour mesurer la precision finale : a 1000 m
                                 % le krigeage a 4,5 km tombe a 49 % alors que son erreur
                                 % mediane y est bonne, et a 500 m la carte connue elle-meme
                                 % perd des essais.
                                 % LE SEUIL EST UN CHOIX DE DEPOUILLEMENT, PAS DE
                                 % SIMULATION : navigation.m n'enregistre que l'erreur, et
                                 % la convergence s'en deduit ici. Le baisser ne demande
                                 % donc pas de rejouer quoi que ce soit. Attention tout de
                                 % meme, il sert deux fois — au taux du tableau, et a
                                 % selectionner les essais que moyennent les courbes.

%% La campagne
d = dir(fullfile(here, 'resultats', '*_reference_A_*runs.mat'));
if isempty(d)
    error('figures_campagnes:campagne', 'Aucune campagne A sous resultats/.');
end
[~, i] = max([d.datenum]);
S = load(fullfile(here, 'resultats', d(i).name));
fprintf('Campagne A : %s\n', d(i).name);
par = S.par;  res = S.res;
if ~isempty(seuil_conv), par.conv.seuil = seuil_conv; end
fprintf('Seuil de convergence : %g m\n', par.conv.seuil);

% Le RPF seul : ces figures comparent des modeles d'observation, pas des filtres.
res  = res(strcmp({res.filter}, 'RPF'));
noms = cellfun(@char, {res.model}, 'UniformOutput', false);
pas  = unique([res.balayage]) / 1e3;

% 'connue' et non 'carte' : le modele krige s'appelle « OK carte totale », donc 'carte'
% l'attrape aussi et le rangerait dans la colonne de la carte connue, qu'il ecraserait.
% Les trois motifs doivent etre exclusifs.
modeles = {'connue', 'bilin', 'totale'};
etiq    = {'true map', 'bilinear', 'kriging'};

d_x  = numel(par.x0);
n_it = size(res(1).err, 1);
kk   = (1:n_it)';
fin  = n_it - par.conv.n_fin + 1 : n_it;
maha = chi2inv(par.conv.confiance, d_x);

RMSE = nan(n_it, numel(pas), numel(modeles));
NEES = nan(n_it, numel(pas), numel(modeles));
CONV = nan(numel(pas), numel(modeles), 2);        % (:,:,1) seuil, (:,:,2) Mahalanobis
for e = 1:numel(res)
    ip = find(abs(pas - res(e).balayage/1e3) < 1e-9, 1);
    im = find(cellfun(@(m) contains(noms{e}, m), modeles), 1);
    if isempty(ip) || isempty(im), continue, end
    E = res(e).err;
    D = res(e).d2;  D(~isfinite(D)) = NaN;
    CONV(ip, im, 1) = 100 * mean(median(E(fin, :), 1) <= par.conv.seuil);
    CONV(ip, im, 2) = 100 * mean(median(D(fin, :), 1) <= maha, 'omitnan');

    % Les courbes portent les ESSAIS CONVERGES seuls. Sur tous les essais, les quelques
    % pistes perdues valent des dizaines de kilometres et ecrasent l'echelle : les courbes
    % se rejoignent vers 4 km, y compris celle de la carte connue, qui finit pourtant a
    % 102 m. Combien d'essais survivent est la question du tableau, a quelle precision
    % celle de la figure.
    ok = median(E(fin, :), 1) <= par.conv.seuil;
    RMSE(:, ip, im) = sqrt(mean(E(:, ok).^2, 2));
    NEES(:, ip, im) = median(D(:, ok), 2, 'omitnan') / d_x;
end

% La carte connue ne depend pas de la maille : une seule courbe, reprise dans chaque case.
ref_rmse = RMSE(:, 1, 1);
ref_nees = NEES(:, 1, 1);

%% Grille RMSE
y_rmse = max(RMSE(:), [], 'omitnan');
fig = figure('Color', 'w', 'Units', 'points', 'Position', [40 40 pdf_w pdf_h]);
t = tiledlayout(fig, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for ip = 1:numel(pas)
    ax = nexttile(t); hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on');
    plot(ax, kk, ref_rmse,        'k:', 'LineWidth', linewidth);
    plot(ax, kk, RMSE(:, ip, 3),  '-',  'LineWidth', linewidth);
    plot(ax, kk, RMSE(:, ip, 2),  '--', 'LineWidth', linewidth);
    hold(ax, 'off');
    title(ax, sprintf('$%.1f$ km', pas(ip)), 'Interpreter', 'latex', 'FontSize', fontsize);
    set(ax, 'TickLabelInterpreter', 'latex', 'FontSize', fontsize);
    xlim(ax, [1 n_it]); ylim(ax, [0 1.02 * y_rmse]);
    if ip == 1
        legend(ax, etiq([1 3 2]), 'Interpreter', 'latex', 'FontSize', fontsize, ...
               'Location', 'northeast');
    end
    if ip < numel(pas) - 1, set(ax, 'XTickLabel', []); end
    if mod(ip, 2) == 0,     set(ax, 'YTickLabel', []); end
end
xlabel(t, 'time step $k$', 'Interpreter', 'latex', 'FontSize', fontsize);
ylabel(t, '$\mathrm{RMSE}_{\hat{x}}(k)$ (m)', 'Interpreter', 'latex', 'FontSize', fontsize);
exporter(fig, 'RMSE par maille', pdf_w, pdf_h, img_dir);

%% Grille NEES
bas  = chi2inv((1 - par.conv.confiance) / 2, d_x) / d_x;
haut = chi2inv(1 - (1 - par.conv.confiance) / 2, d_x) / d_x;
y_nees = [0.5 * bas, 2 * max([NEES(:); haut], [], 'omitnan')];
kd = 1 : pas_trace : n_it;                        % les pas effectivement traces
fig = figure('Color', 'w', 'Units', 'points', 'Position', [40 40 pdf_w pdf_h]);
t = tiledlayout(fig, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for ip = 1:numel(pas)
    ax = nexttile(t); hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on');
    fill(ax, [1 n_it n_it 1], [bas bas haut haut], [1 0.85 0.88], ...
         'EdgeColor', 'none', 'FaceAlpha', 0.55, 'HandleVisibility', 'off');
    plot(ax, kk(kd), ref_nees(kd),       'k:', 'LineWidth', linewidth);
    plot(ax, kk(kd), NEES(kd, ip, 3),    '-',  'LineWidth', linewidth);
    plot(ax, kk(kd), NEES(kd, ip, 2),    '--', 'LineWidth', linewidth);
    hold(ax, 'off'); set(ax, 'YScale', 'log');
    title(ax, sprintf('$%.1f$ km', pas(ip)), 'Interpreter', 'latex', 'FontSize', fontsize);
    set(ax, 'TickLabelInterpreter', 'latex', 'FontSize', fontsize);
    % Trois graduations, deux decades d'ecart : la log par defaut en pose une par decade,
    % ce qui remplit la case d'etiquettes sans rien apprendre de plus.
    set(ax, 'YTick', 10.^(0:2:4), 'YMinorTick', 'off', 'YMinorGrid', 'off');
    xlim(ax, [1 n_it]); ylim(ax, y_nees);
    if ip == 1
        legend(ax, etiq([1 3 2]), 'Interpreter', 'latex', 'FontSize', fontsize, ...
               'Location', 'northwest');
    end
    if ip < numel(pas) - 1, set(ax, 'XTickLabel', []); end
    if mod(ip, 2) == 0,     set(ax, 'YTickLabel', []); end
end
xlabel(t, 'time step $k$', 'Interpreter', 'latex', 'FontSize', fontsize);
% NEES et non NEES / d_x : la definition du manuscrit porte deja la division par d_x.
ylabel(t, '$\mathrm{NEES}_{\hat{x}}(k)$', 'Interpreter', 'latex', 'FontSize', fontsize);
exporter(fig, 'NEES par maille', pdf_w, pdf_h, img_dir);

%% ARMSE de la reconstruction de la carte
%
% Ce que coute la reconstruction AVANT tout filtrage : l'ecart entre le champ vrai et ce
% que chaque modele en rend, sur le corridor survole. C'est la grandeur qui explique les
% deux figures precedentes, et elle ne depend d'aucun filtre.
%
% ARMSE AU SENS DU MANUSCRIT : la racine est prise sur les tirages de releve A POINT FIXE,
% puis moyennee sur l'espace. Un simple sqrt(mean(.^2)) sur un tirage donnerait une RMSE
% spatiale, ce que le cache stocke deja sous err_carte, et non une ARMSE.
%
% Le krigeage est relu du cache, ou la campagne a laisse ses estimateurs tabules. Le
% bilineaire, lui, ne se cache pas : on refait le releve du meme tirage et on l'interpole.
n_draw = 10;                     % tirages de releve, sur les 100 que la campagne a produits
krig   = krigeage();
cm     = load(par.map.file);
map    = struct('h', double(cm.mesures).', 'step', double(cm.pas(:))' .* [1 1]);

ARMSE = nan(numel(pas), 2);      % (:,1) bilineaire, (:,2) krigeage
for ip = 1:numel(pas)
    f = dir(fullfile(here, 'cache', sprintf('ok_*_p%d_*_matern32_ph_j30.mat', pas(ip) * 1e3)));
    if isempty(f), continue, end
    f = f(1:min(n_draw, numel(f)));
    S2 = 0; S2b = 0; vrai = [];
    for j = 1:numel(f)
        C  = load(fullfile(here, 'cache', f(j).name));
        st = C.tab.step;
        ax = st(2, 1) + (0:size(C.tab.z, 2) - 1) * st(1, 1);
        ay = st(2, 2) + (0:size(C.tab.z, 1) - 1) * st(1, 2);
        [GX, GY] = meshgrid(ax, ay);
        P = [GX(:), GY(:)]';
        if isempty(vrai), vrai = krig.lire(map.h, map.step, P); end

        S2 = S2 + (C.tab.z(:) - vrai).^2;

        % Le meme releve, rejoue depuis sa graine, lu par interpolation lineaire.
        p2 = par;
        p2.krig.survey_seed = C.key.seed;
        p2.krig.pas         = pas(ip) * 1e3;
        p2.krig.jitter      = C.key.jitter;
        p2.krig.phase_alea  = true;
        sv = krig.releve(map, p2);
        F  = scatteredInterpolant(sv.X(1, :)', sv.X(2, :)', sv.z, 'linear', 'nearest');
        S2b = S2b + (F(P(1, :)', P(2, :)') - vrai).^2;
    end
    ARMSE(ip, 2) = mean(sqrt(S2  / numel(f)));
    ARMSE(ip, 1) = mean(sqrt(S2b / numel(f)));
    fprintf('  %.1f km : ARMSE krigeage %.1f %s, bilineaire %.1f %s (%d tirages)\n', ...
            pas(ip), ARMSE(ip, 2), par.unit, ARMSE(ip, 1), par.unit, numel(f));
end

% UN EQUIVALENT SIMPLE. Sur la plage balayee, l'ARMSE suit une loi LOGARITHMIQUE de la
% maille, mieux qu'une puissance, une racine ou une droite : R2 = 0.9996 et 1.7 % d'ecart
% maximal pour le krigeage, contre 0.979 et 6.5 % pour la meilleure loi de puissance. On
% l'ecrit a * ln(p / p0), p0 etant la maille sous laquelle la reconstruction serait
% exacte. Prudence : trois octaves de maille seulement, donc l'equivalent decrit la plage
% et ne s'extrapole pas — a grande maille l'erreur sature vers sigma_f et ne peut pas
% croitre indefiniment.
lois     = zeros(2, 2);
etiq_rec = {'bilineaire', 'krigeage'};
for m = 1:2
    c = polyfit(log(pas(:)), ARMSE(:, m), 1);
    lois(m, :) = [c(1), exp(-c(2) / c(1))];
    fprintf('  equivalent %s : %.1f ln(p / %.2f) nT\n', etiq_rec{m}, c(1), lois(m, 2));
end

% La figure ne porte QUE les deux courbes. L'equivalent et le nombre de tirages sont du
% texte : ils vont dans la caption et dans l'analyse, comme pour les autres figures du
% chapitre, dont aucune n'a de titre.
fig = figure('Color', 'w', 'Units', 'points', 'Position', [80 80 348 280]);
hold on; box on; grid on;
plot(pas, ARMSE(:, 1), 'o--', 'LineWidth', linewidth);
plot(pas, ARMSE(:, 2), 'o-',  'LineWidth', linewidth);
hold off;
fs = 12 * 348 / (0.7 * textwidth);
xlabel('map resolution (km)', 'Interpreter', 'latex', 'FontSize', fs);
ylabel('$\mathrm{ARMSE}_{\hat{z}}$ (nT)', 'Interpreter', 'latex', 'FontSize', fs);
legend({'bilinear', 'kriging'}, 'Interpreter', 'latex', 'FontSize', fs, 'Location', 'best');
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fs);
xlim([min(pas)-0.3 max(pas)+0.3]);
exporter(fig, 'ARMSE de reconstruction', 348, 280, img_dir);

%% Le cout arithmetique par modele
%
% ANALYTIQUE ET NON MESURE. res.cout est nul pour le krigeage ordinaire de cette campagne,
% et le temps mur ne dit rien non plus, le bilineaire y sortant plus lent que le krigeage.
% Les complexites de th:kriging_mean_cost, instanciees a chaque maille, sont donc la seule
% mesure honnete.
%
% ON COMPTE COMME SI TOUT ETAIT RECALCULE A CHAQUE PAS. La tabulation de l'estimateur est
% une optimisation d'implementation : elle ne change pas la complexite de la methode, et
% le manuscrit n'a pas a l'exposer. Deux postes par essai, pour un champ scalaire :
%   setup     n^3/3, la resolution du systeme d'entrainement, une fois
%   queries   n^2 par particule et par pas, soit n_part * n_steps * n^2
% Le bilineaire ne paie pas de setup et lit en huit flops par requete.
n_tilde = nan(1, numel(pas));
for ip = 1:numel(pas)
    e = find([res.balayage] == pas(ip) * 1e3 & contains(noms, 'totale'), 1);
    n_tilde(ip) = max(res(e).n_used(:));
end
n_query = par.mc.n_part * par.mc.n_steps;

cout_setup = n_tilde.^3 / 3;
cout_req   = n_query * n_tilde.^2;
cout_bil   = n_query * 8;                         % huit flops par lecture bilineaire

if ~exist(tab_dir, 'dir'), mkdir(tab_dir); end
fid = fopen(fullfile(tab_dir, 'cout_maille.tex'), 'w', 'n', 'UTF-8');
fprintf(fid, '%% Generated by figures_campagnes.m, do not edit by hand.\n');
fprintf(fid, ['%% Analytic flop counts per run, in Gflop, not timings, for %d particles over\n' ...
              '%% %d steps. Setup solves the training system once; queries cost n^2 per\n' ...
              '%% particle and per step. The true map and the bilinear model are read in eight\n' ...
              '%% flops per query and pay no setup, hence one row each: neither depends on the\n' ...
              '%% map resolution.\n\n'], par.mc.n_part, par.mc.n_steps);
fprintf(fid, '\\begin{tabular}{|l|cc|}\\hline\n');
% L'UNITE EST FACTORISEE DANS L'EN-TETE, LE PREFIXE RESTE EN CELLULE. Une colonne
% entierement en Gflop forcerait 29213 d'un cote et 0.01 de l'autre, donc cinq chiffres
% pour l'un et deux decimales perdues pour l'autre ; un prefixe par valeur garde trois
% chiffres significatifs partout, sans allonger les cellules.
fprintf(fid, 'Observation model & Setup (flops) & Queries (flops)\\\\\\hline\n');
fprintf(fid, 'True map & -- & %s\\\\\n', flops_tex(cout_bil));
fprintf(fid, 'Bilinear & -- & %s\\\\\\hline\n', flops_tex(cout_bil));
for ip = 1:numel(pas)
    fprintf(fid, 'Kriging, $%.1f$km & %s & %s\\\\\n', ...
            pas(ip), flops_tex(cout_setup(ip)), flops_tex(cout_req(ip)));
end
fprintf(fid, '\\hline\n\\end{tabular}\n');
fclose(fid);

fprintf('Cout par essai (Gflop) :\n');
for ip = 1:numel(pas)
    fprintf('  %.1f km, n=%5d : setup %7.1f  requetes %9.1f\n', ...
            pas(ip), n_tilde(ip), cout_setup(ip)/1e9, cout_req(ip)/1e9);
end
%% Convergence et erreur finale contre la maille
%
% Deux figures a deux panneaux, pleine largeur. Elles remplacent le tableau des taux de
% convergence : treize lignes de pourcentages se lisent moins bien que deux courbes, et
% la comparaison entre criteres se fait alors d'un panneau a l'autre.
%
% « Finale » veut dire mediane sur les dix derniers pas, la fenetre qui sert deja aux deux
% criteres de convergence.
%
% ON MEDIANE LES COURBES DEJA CALCULEES, et non les erreurs brutes. La RMSE du manuscrit
% est une racine prise sur les essais A PAS FIXE ; medianer directement le nuage des
% couples (pas, essai) donnerait une mediane d'erreurs, qui n'est pas une RMSE et que
% l'axe ne pourrait pas appeler ainsi. Meme raison pour le NEES, moyenne sur les essais a
% pas fixe. Les deux figures portent donc exactement les grandeurs des grilles.
FIN = cat(3, squeeze(median(RMSE(fin, :, :), 1)), squeeze(median(NEES(fin, :, :), 1)));

sty = {':', '--', '-'};                           % carte connue, bilineaire, krigeage
duo = @(titres, Y, ylab, logy) deal_duo(pas, Y, titres, ylab, logy, sty, etiq, ...
                                        textwidth, fontsize, linewidth);

fig = duo({'Threshold', 'Mahalanobis'}, CONV, 'converged runs (\%)', false);
exporter(fig, 'Convergence par maille', textwidth, 230, img_dir);

fig = duo({'$\mathrm{RMSE}_{\hat{x}}$ (m)', '$\mathrm{NEES}_{\hat{x}}$'}, FIN, ...
          'median over the last ten steps', true);
exporter(fig, 'Erreur finale par maille', textwidth, 230, img_dir);

%% Le tableau des taux de convergence, sous les deux criteres
if false        % remplace par la figure Convergence par maille, garde le temps de valider
if ~exist(tab_dir, 'dir'), mkdir(tab_dir); end
fid = fopen(fullfile(tab_dir, 'convergence_maille.tex'), 'w', 'n', 'UTF-8');
fprintf(fid, '%% Generated by figures_campagnes.m, do not edit by hand.\n');
fprintf(fid, ['%% Converged runs out of %d, regularised filter. Threshold: median error over\n' ...
              '%% the last %d steps below %g m. Mahalanobis: median over the same steps below\n' ...
              '%% the %g quantile of a chi2 with %d degrees of freedom. True map: %.0f%% and %.0f%%.\n\n'], ...
        par.mc.n_runs, par.conv.n_fin, par.conv.seuil, par.conv.confiance, d_x, ...
        CONV(1, 1, 1), CONV(1, 1, 2));
% Une ligne par fonction d'observation, soit les treize configurations comparees : la
% carte connue, puis le krigeage et le bilineaire a chacune des six mailles. Groupe par
% modele et non par maille, pour que les deux series se lisent en colonne.
fprintf(fid, '\\begin{tabular}{|l|cc|}\\hline\n');
fprintf(fid, 'Observation model & Threshold & Mahalanobis\\\\\\hline\n');
fprintf(fid, 'True map & $%.0f\\%%$ & $%.0f\\%%$\\\\\\hline\n', CONV(1, 1, 1), CONV(1, 1, 2));
for im = [3 2]
    nom = [upper(etiq{im}(1)) etiq{im}(2:end)];
    for ip = 1:numel(pas)
        fprintf(fid, '%s, $%.1f$km & $%.0f\\%%$ & $%.0f\\%%$\\\\\n', ...
                nom, pas(ip), CONV(ip, im, 1), CONV(ip, im, 2));
    end
    fprintf(fid, '\\hline\n');
end
fprintf(fid, '\\end{tabular}\n');
fclose(fid);
end

fprintf('Deux grilles et un tableau ecrits\n');

%% ------------------------------------------------------------ fonctions locales

function fig = deal_duo(pas, Y, titres, ylab, logy, sty, etiq, textwidth, fontsize, linewidth)
% Deux panneaux cote a cote sur une meme abscisse de maille, trois modeles par panneau.
% Y est (maille, modele, panneau). ylab est commun aux deux panneaux, ou l'un par panneau
% quand les deux grandeurs n'ont pas la meme unite.
    fig = figure('Color', 'w', 'Units', 'points', 'Position', [40 40 textwidth 230]);
    t = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    for j = 1:2
        ax = nexttile(t); hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on');
        for im = 1:size(Y, 2)
            plot(ax, pas, Y(:, im, j), sty{im}, 'LineWidth', linewidth, 'Marker', 'o');
        end
        hold(ax, 'off');
        if logy, set(ax, 'YScale', 'log'); else, ylim(ax, [0 105]); end
        title(ax, titres{j}, 'Interpreter', 'latex', 'FontSize', fontsize);
        set(ax, 'TickLabelInterpreter', 'latex', 'FontSize', fontsize);
        xlim(ax, [min(pas) - 0.3, max(pas) + 0.3]);
        % Les graduations restent sur LES DEUX panneaux : ils ne partagent pas d'echelle,
        % mieux, pas meme d'unite — des metres a gauche, un rapport a droite. Les masquer
        % a droite ferait lire les valeurs de gauche.
        if iscell(ylab)
            ylabel(ax, ylab{j}, 'Interpreter', 'latex', 'FontSize', fontsize);
        elseif j == 1
            ylabel(ax, ylab, 'Interpreter', 'latex', 'FontSize', fontsize);
        end
        if j == 1
            legend(ax, etiq, 'Interpreter', 'latex', 'FontSize', fontsize, ...
                   'Location', 'best');
        end
    end
    xlabel(t, 'map resolution (km)', 'Interpreter', 'latex', 'FontSize', fontsize);
end

function t = flops_tex(v)
% Trois chiffres significatifs et le prefixe decimal, colle a la valeur comme le reste du
% document ecrit ses unites : $200\text{km}$, $10^2\text{nT}^2$.
    units = {'', 'k', 'M', 'G', 'T', 'P'};
    i     = max(1, min(numel(units), 1 + floor(log10(v) / 3)));
    t     = sprintf('$%.3g$%s', v / 1000^(i - 1), units{i});
end

function exporter(fig, nom, pdf_w, pdf_h, img_dir)
    if ~exist(img_dir, 'dir'), mkdir(img_dir); end
    set(fig, 'PaperUnits', 'points', 'PaperSize', [pdf_w pdf_h], ...
             'PaperPosition', [0 0 pdf_w pdf_h]);
    print(fig, fullfile(img_dir, [nom '.pdf']), '-dpdf', '-vector');
end
