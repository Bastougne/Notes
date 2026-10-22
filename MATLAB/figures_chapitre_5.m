% FIGURES_CHAPITRE_5
%
% Les figures et le tableau de la section de resultats du chapitre 5, tires des campagnes
% enregistrees par campagne_chapitre_5.m sur le scenario 'boucle'.
%
%   Champ et circuit du SLAM   le champ analytique, le circuit et le domaine de la base
%                              de Laplace
%   RMSE et NEES du SLAM       l'erreur de position des six filtres au cours du vol, et
%                              leur NEES avec la bande du chi2 du chapitre 4
%   resultats_slam.tex         les erreurs moyennes, sous Latex/Manuscrit/include
%
% CHAQUE METHODE EST LUE DANS UN FICHIER NOMME, et non plus dans le plus recent qui la
% contient : un balayage range dans resultats/ aurait remplace la campagne du chapitre. La
% campagne A du 27 septembre porte les deux bornes, l'inertiel et le RBPF ; son AKKF lisait
% ses moments en coordonnees brutes et a perdu ses cent essais. AKKF-SLAM vient de la
% valeur r_RFF = 500 du balayage du 13 octobre, dont la base RFF est retiree a chaque essai,
% et l'AKKF sur carte connue du lot du 20 octobre (akkf.obs_simulees = 1). Les essais sont
% apparies d'un fichier a l'autre (meme graine, meme vol, memes bruits d'observation) ; le
% script verifie qu'ils partagent le vol, le nombre de particules et le nombre d'essais.
%
% Les conventions sont celles de figures_campagnes.m : page figee, police calculee pour
% sortir a 12 pt dans le manuscrit, une couleur et un trait par methode.

clear; close all; clc;
here = fileparts(mfilename('fullpath'));

img_dir   = fullfile(fileparts(here), 'Images');
tab_dir   = fullfile(fileparts(here), 'Latex', 'Manuscrit', 'include');
res_dir   = fullfile(here, 'resultats');
textwidth = 425;                          % pt, mesure sur le document
linewidth = 1;
pas_trace = 3;                            % un point sur trois sur la NEES, comme au chapitre 4

% Cle de style, methode dans le fichier, fichier sous resultats/, valeur balayee (NaN : aucune).
% L'ordre est celui de la legende, qui se remplit colonne par colonne : les deux cartes
% connues, puis les deux bornes de navigation, puis les deux SLAM.
camp_A  = '20260927_032929_ch5_boucle_A_100runs.mat';
sources = {
    'carte',      'carte',    camp_A,                                                       NaN
    'akkf_carte', 'akkf',     fullfile('ch5_lots', '20261020_121849_ch5_boucle_libre_100runs.mat'), 1
    'ok',         'ok',       camp_A,                                                       NaN
    'inertiel',   'inertiel', camp_A,                                                       NaN
    'rbpf',       'rbpf',     camp_A,                                                       NaN
    'akkf',       'akkf',     fullfile('ch5_lots', '20261013_194404_ch5_boucle_libre_100runs.mat'), 500};

%% Lecture
R = struct('cle', {}, 'res', {}, 'fichier', {});
ref = [];
for ic = 1:size(sources, 1)
    S = load(fullfile(res_dir, sources{ic, 3}), 'par', 'res', 'scen', 'P_eval', ...
             'h_eval', 'pres_traj');
    sel = strcmp({S.res.methode}, sources{ic, 2});
    if ~isnan(sources{ic, 4}), sel = sel & [S.res.balayage] == sources{ic, 4}; end
    i = find(sel, 1);
    if isempty(i)
        error('figures_chapitre_5:lots', '%s ne contient pas %s.', sources{ic, 3}, sources{ic, 1});
    end
    if isempty(ref)
        ref = S;
    elseif ~isequal(S.scen.x_true, ref.scen.x_true) || ...
           S.par.mc.n_part ~= ref.par.mc.n_part || S.par.mc.n_runs ~= ref.par.mc.n_runs
        error('figures_chapitre_5:lots', '%s ne partage pas le scenario de %s.', ...
              sources{ic, 3}, R(1).fichier);
    end
    R(end + 1) = struct('cle', sources{ic, 1}, 'res', S.res(i), 'fichier', sources{ic, 3}); %#ok<SAGROW>
    if strcmp(sources{ic, 1}, 'rbpf'), par_rbpf = S.par; end
    if strcmp(sources{ic, 1}, 'akkf')
        fprintf(['AKKF-SLAM : lambda_x %g, lambda_z %g, oubli %g, echelle_p %g, r_rff %d, ' ...
                 'base par essai %d\n'], S.par.akkf.lambda_x, S.par.akkf.lambda_z, ...
                S.par.akkf.oubli, S.par.akkf.echelle_p, S.par.slam.r_rff, S.par.slam.rff_par_essai);
    end
end
par    = ref.par;
x_true = ref.scen.x_true;
n_it   = size(R(1).res.err, 1);
d_x    = size(x_true, 1);
lire   = @(c) R(strcmp({R.cle}, c)).res;     % la sortie d'une methode, par sa cle
for ic = 1:numel(R)
    fprintf('%-10s lu dans %s\n', R(ic).cle, R(ic).fichier);
end

%% Le champ, le circuit et le domaine de la base de Laplace
% Le champ analytique est recalcule comme dans campagne_chapitre_5.m, sur la grille de la
% carte du chapitre 4. La figure est inseree a 0.7\textwidth, donc sa police est grossie
% d'autant pour sortir a 12 pt.
pdf_w = 348;  pdf_h = 280;
fs    = 12 * pdf_w / (0.7 * textwidth);

c = [0.7, 1.2, 0.4];
s = par.map.lambda * c(1) / (2 * pi);
A = par.map.sigma / sqrt(0.5 * (1 + 0.5^2 + 0.3^2));
g = (1:201) * 1000;
[Gx, Gy] = meshgrid(g, g);
h = A * (sin(c(1) * Gx / s) + 0.5 * cos(c(2) * Gy / s) + 0.3 * sin(c(3) * (Gx + Gy) / s));
fprintf(['Champ : amplitude %.1f nT, longueurs d''onde %.2f, %.2f et %.2f km, ' ...
         'ecart-type %.1f nT\n'], A, 2 * pi * s ./ c / 1e3, std(h(:)));

dom = [min(x_true(1, :)) - par.slam.marge, max(x_true(1, :)) + par.slam.marge, ...
       min(x_true(2, :)) - par.slam.marge, max(x_true(2, :)) + par.slam.marge] / 1e3;
vue = dom + 8 * [-1 1 -1 1];              % km, un peu plus large que le domaine

fig = figure('Color', 'w', 'Units', 'points', 'Position', [80 80 pdf_w pdf_h]);
t   = tiledlayout(fig, 1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ax  = nexttile(t);
imagesc(ax, g / 1e3, g / 1e3, h); hold(ax, 'on');
set(ax, 'YDir', 'normal'); axis(ax, 'image');
rectangle(ax, 'Position', [dom(1), dom(3), dom(2) - dom(1), dom(4) - dom(3)], ...
          'EdgeColor', 'k', 'LineStyle', '--', 'LineWidth', linewidth);
plot(ax, x_true(1, :) / 1e3, x_true(2, :) / 1e3, 'k-', 'LineWidth', linewidth);
plot(ax, x_true(1, 1) / 1e3, x_true(2, 1) / 1e3, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
hold(ax, 'off');
xlim(ax, vue(1:2)); ylim(ax, vue(3:4));
colormap(ax, parula);
cb = colorbar(ax);
cb.Layout.Tile = 'east';                  % un element de la grille, pour que son label reste dans la page
cb.TickLabelInterpreter = 'latex';  cb.FontSize = fs;
ylabel(cb, '$h(x)$ (nT)', 'Interpreter', 'latex', 'FontSize', fs);
xlabel(ax, '$p^{(x)}$ (km)', 'Interpreter', 'latex', 'FontSize', fs);
ylabel(ax, '$p^{(y)}$ (km)', 'Interpreter', 'latex', 'FontSize', fs);
set(ax, 'TickLabelInterpreter', 'latex', 'FontSize', fs, 'Layer', 'top');
box(ax, 'on');
exporter(fig, 'Champ et circuit du SLAM', pdf_w, pdf_h, img_dir);

%% RMSE et NEES au cours du vol
kk = (1:n_it)';
kd = (1:pas_trace:n_it)';
bas  = chi2inv((1 - par.conv.confiance) / 2, d_x) / d_x;
haut = chi2inv(1 - (1 - par.conv.confiance) / 2, d_x) / d_x;
courbes(kk, kd, R, [bas haut], textwidth, linewidth, img_dir, 'RMSE et NEES du SLAM');

%% Le tableau des erreurs moyennes
% Trois tiers du vol pour l'erreur de position, le dernier pour la NEES, et l'erreur de
% reconstruction du champ le long du vol pour les methodes qui en estiment un. L'inertiel
% n'en estime pas : sa ligne donne l'erreur du champ nul, a laquelle les deux autres se
% comparent. L'AKKF sur carte connue apprend bien une carte a cote, mais elle ne lui sert
% pas : sa case reste vide comme celles des deux bornes, et le texte la cite.
tiers = {1:100, 101:200, 201:n_it};
fid = fopen(fullfile(tab_dir, 'resultats_slam.tex'), 'w', 'n', 'UTF-8');
fprintf(fid, '%% Generated by figures_chapitre_5.m, do not edit by hand.\n');
fprintf(fid, ['%% Scenario %s, %d runs, %d particles. Position errors are averaged over\n' ...
              '%% three thirds of the mission, the NEES over the last one, and the field\n' ...
              '%% error over the evaluation points within %g km of the trajectory.\n'], ...
        par.scenario, par.mc.n_runs, par.mc.n_part, 2 * par.slam.ell / 1e3);
for ic = 1:numel(R)
    fprintf(fid, '%% %s read from %s\n', R(ic).cle, R(ic).fichier);
end
fprintf(fid, '\n\\begin{tabular}{|l|ccc|c|c|}\\hline\n');
% L'en-tete de la derniere colonne tient sur deux lignes : sur une seule, le tableau
% debordait de la page de neuf points.
fprintf(fid, ['\\multirow{2}{*}{Filter} & \\multicolumn{3}{c|}{$\\mathrm{ARMSE}_{\\hat{x}}$ (m)}' ...
              ' & \\multirow{2}{*}{NEES} & Field ARMSE\\\\\\cline{2-4}\n']);
fprintf(fid, ' & $[\\![1,100]\\!]$ & $[\\![101,200]\\!]$ & $[\\![201,%d]\\!]$ & & (nT)\\\\\\hline\n', n_it);
for ic = 1:numel(R)
    [~, ~, nom] = style_methode(R(ic).cle);
    r = R(ic).res;
    fprintf(fid, '%s', nom);
    for j = 1:3
        fprintf(fid, ' & $%.0f$', mean(r.rmse(tiers{j})));
    end
    fprintf(fid, ' & $%s$', nombre(mean(r.nees(tiers{3}))));
    if any(strcmp(R(ic).cle, {'carte', 'akkf_carte', 'ok'}))
        fprintf(fid, ' & --');
    else
        fprintf(fid, ' & $%.1f$', r.armse_map_traj);
    end
    fprintf(fid, '\\\\\\hline\n');
end
fprintf(fid, '\\end{tabular}\n');
fclose(fid);

% Ce que le texte cite et que le tableau ne porte pas.
fprintf('\n%-10s %8s %8s %8s %8s %8s | NEES k=100 %8s fin | carte vol / hors vol\n', 'methode', ...
        'k=1', 'k=30', 'k=78', 'k=157', 'k=300', '');
for ic = 1:numel(R)
    r = R(ic).res;
    p = ref.pres_traj;
    fprintf('%-10s %8.0f %8.0f %8.0f %8.0f %8.0f | %8.2f %8.2f | %6.1f / %6.1f\n', R(ic).cle, ...
            r.rmse([1 30 78 157 n_it]), r.nees(100), mean(r.nees(n_it - 9:n_it)), ...
            mean(r.rmse_map(p)), mean(r.rmse_map(~p)));
end
for ic = 1:numel(R)
    r = R(ic).res.rmse;
    [m1, k1] = max(r(1:100));  [m2, k2] = min(r(60:120));
    fprintf(['%-10s max sur [1,100] : %5.0f m au pas %3d ; min sur [60,120] : %5.0f m au pas %3d ;' ...
             ' min et max apres 100 : %5.0f et %5.0f m\n'], R(ic).cle, m1, k1, m2, k2 + 59, ...
            min(r(101:end)), max(r(101:end)));
end
nees_rbpf = lire('rbpf').nees;
fprintf('La NEES du RBPF depasse 1 au pas %d, et le haut de la bande (%.2f) au pas %d\n', ...
        find(nees_rbpf > 1, 1), haut, find(nees_rbpf > haut, 1));
k_cr = find(lire('akkf').rmse < lire('inertiel').rmse & kk > 100, 1);
fprintf('L''AKKF passe sous l''inertiel au pas %d ; apres le pas 100, son RMSE y est inferieur %d pas sur %d\n', ...
        k_cr, sum(lire('akkf').rmse(101:end) < lire('inertiel').rmse(101:end)), n_it - 100);
for c = {'akkf', 'akkf_carte'}
    y = lire(c{1}).nees;
    fprintf('NEES %-10s : mediane %.2f, dernier tiers %.2f, min et max apres 100 : %.1f et %.1f\n', ...
            c{1}, median(y), mean(y(tiers{3})), min(y(101:end)), max(y(101:end)));
end
fprintf('Temps par essai : %s\n', strjoin(arrayfun(@(i) sprintf('%s %.1f s', R(i).cle, ...
        R(i).res.wall / par.mc.n_runs), 1:numel(R), 'UniformOutput', false), ', '));
fprintf('Rang de la base de Laplace du RBPF : %d ; domaine %.1f x %.1f km\n', ...
        par_rbpf.slam.r_kl, dom(2) - dom(1), dom(4) - dom(3));

%% -------------------------------------------------------------- utilitaires

function courbes(kk, kd, R, bande, textwidth, linewidth, img_dir, nom_fig)
% Les six filtres sur deux panneaux empiles, l'erreur de position en haut et la NEES en
% bas, avec une seule legende au-dessus : dans un panneau elle couvrirait une courbe, les
% filtres sur carte connue longeant le bas de l'axe et les trois autres le haut. La figure est inseree
% a \textwidth, donc sa police n'est pas mise a l'echelle.
    pdf_w = textwidth;  pdf_h = 400;  fs = 12;
    n_it  = numel(kk);
    fig   = figure('Color', 'w', 'Units', 'points', 'Position', [40 40 pdf_w pdf_h]);
    t     = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    ax = nexttile(t); hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on');
    leg = cell(1, numel(R));
    ymax = 0;
    for ic = 1:numel(R)
        [c, l, leg{ic}] = style_methode(R(ic).cle);
        plot(ax, kk, R(ic).res.rmse, l, 'Color', c, 'LineWidth', linewidth);
        ymax = max(ymax, max(R(ic).res.rmse));
    end
    hold(ax, 'off');
    xlim(ax, [1 n_it]);  ylim(ax, [0, 1.05 * ymax]);
    set(ax, 'TickLabelInterpreter', 'latex', 'FontSize', fs, 'XTickLabel', []);
    ylabel(ax, '$\mathrm{RMSE}_{\hat{x}}(k)$ (m)', 'Interpreter', 'latex', 'FontSize', fs);
    lg = legend(ax, leg, 'Interpreter', 'latex', 'FontSize', fs, 'NumColumns', 3);
    lg.ItemTokenSize = [18 18];
    lg.Layout.Tile   = 'north';

    ax = nexttile(t); hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on');
    fill(ax, [1 n_it n_it 1], [bande(1) bande(1) bande(2) bande(2)], ...
         [1 0.85 0.88], 'EdgeColor', 'none', 'FaceAlpha', 0.55, 'HandleVisibility', 'off');
    lim = [Inf, 0];
    for ic = 1:numel(R)
        [c, l] = style_methode(R(ic).cle);
        y = R(ic).res.nees;
        plot(ax, kk(kd), y(kd), l, 'Color', c, 'LineWidth', linewidth);
        lim = [min(lim(1), min(y(kd))), max(lim(2), max(y(kd)))];
    end
    hold(ax, 'off');
    lim = [min(0.5 * lim(1), 0.7 * bande(1)), max(2 * lim(2), 1.4 * bande(2))];
    dec = floor(log10(lim(1))) : ceil(log10(lim(2)));
    set(ax, 'YScale', 'log', 'YTick', 10.^dec, 'YMinorGrid', 'off');
    xlim(ax, [1 n_it]);  ylim(ax, lim);
    set(ax, 'TickLabelInterpreter', 'latex', 'FontSize', fs);
    ylabel(ax, '$\mathrm{NEES}_{\hat{x}}(k)$', 'Interpreter', 'latex', 'FontSize', fs);
    xlabel(ax, 'Time step $k$', 'Interpreter', 'latex', 'FontSize', fs);
    exporter(fig, nom_fig, pdf_w, pdf_h, img_dir);
end

function [c, l, nom] = style_methode(cle)
% COULEUR, TRAIT ET NOM D'UNE METHODE. Les deux bornes gardent ceux de figures_campagnes.m,
% pour qu'un lecteur du chapitre 4 les reconnaisse ; la methode proposee prend l'ambre que
% l'AK y porte, celle de l'etat de l'art le rouge, et l'inertiel un gris tirete.
    switch cle
        case 'carte',    c = [0.00 0.00 0.00];  l = ':';   nom = 'Known map';
        case 'akkf_carte', c = [0.95 0.60 0.00]; l = ':';  nom = 'Known map (AKKF)';
        case 'ok',       c = [0.15 0.35 0.70];  l = '-';   nom = 'Static kriging';
        case 'inertiel', c = [0.45 0.45 0.45];  l = '--';  nom = 'Inertial navigation';
        case 'rbpf',     c = [0.75 0.02 0.15];  l = '-';   nom = 'RBPF-SLAM';
        case 'akkf',     c = [0.95 0.60 0.00];  l = '-';   nom = 'AKKF-SLAM';
        otherwise
            error('figures_chapitre_5:style', 'Cle de methode inconnue : %s.', cle);
    end
end

function s = nombre(v)
% Deux decimales sous dix, un entier au-dela.
    if v < 10, s = sprintf('%.2f', v); else, s = sprintf('%.0f', v); end
end

function exporter(fig, nom, pdf_w, pdf_h, img_dir)
% Repris de figures_campagnes.m : le PDF n'est remplace que si le dessin a change, la date
% et l'identifiant que MATLAB inscrit a chaque export mis a part.
    if ~exist(img_dir, 'dir'), mkdir(img_dir); end
    set(fig, 'PaperUnits', 'points', 'PaperSize', [pdf_w pdf_h], ...
             'PaperPosition', [0 0 pdf_w pdf_h]);
    cible = fullfile(img_dir, [nom '.pdf']);
    tmp   = [tempname '.pdf'];
    print(fig, tmp, '-dpdf', '-vector');
    if isfile(cible) && strcmp(sans_horodatage(tmp), sans_horodatage(cible))
        delete(tmp);
    else
        movefile(tmp, cible, 'f');
    end
end

function s = sans_horodatage(f)
    fid = fopen(f, 'r');
    s   = char(fread(fid, Inf, '*uint8')');
    fclose(fid);
    s = regexprep(s, '/(CreationDate|ModDate)\s*\([^)]*\)', '');
    s = regexprep(s, '/ID\s*\[[^\]]*\]', '');
end
