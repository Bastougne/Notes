% FIGURES_CHAPITRE_5
%
% Les figures et les tableaux de la section de resultats du chapitre 5, tires des lots de la
% file version 2 (2026-10-22, MATLAB/resultats/ch5_v2/lots), reperes par leur etiquette :
%
%   principal            RPF sur carte connue, inertiel, RBPF-SLAM (r_KL = 300) et
%                        AKKF-SLAM (r_RFF = 500)
%   akkf_variantes       AKKF sur carte connue (akkf.obs_simulees = 1) et AKKF-SLAM avec la
%                        variance predictive des observations simulees (2)
%   akkf_lz1             lambda_z = 1, sans (0) et avec (2) la variance predictive
%   akkf_carte_lz1       AKKF sur carte connue avec lambda_z = 1 (2026-10-26)
%   rangs_akkf           AKKF-SLAM, r_RFF de 50 a 2000
%   rangs_akkf_combine   AKKF-SLAM avec lambda_z = 1 et la variance predictive, r_RFF de 50 a 1000
%   rangs_rbpf_*         RBPF-SLAM, r_KL de 50 a 500 (la valeur 300 vient du lot principal)
%
% LA VERSION PRINCIPALE D'AKKF-SLAM (2026-10-26, demande de Bastien) est celle a variance
% predictive et lambda_z = 1, lue dans akkf_lz1, et son temoin sur carte connue a lui aussi
% lambda_z = 1. L'algorithme de l'article (lot principal) et les deux modifications prises
% separement ne sont plus que des variantes, citees par le texte sans etre tracees.
%
% UN LOT ABSENT EST SAUTE, pour que le script tourne aussi pendant que la file avance. Quand
% une etiquette apparait dans plusieurs fichiers, le plus recent l'emporte.
%
% Sorties, dans les dossiers de la structure `sortie` si elle existe (la file y met les
% siens), sinon dans ceux du manuscrit :
%
%   Champ et circuit du SLAM.pdf       le champ analytique, le circuit et le domaine du RBPF
%   RMSE et NEES du SLAM.pdf           l'erreur de position et la NEES au cours du vol
%   Cartes reconstruites du SLAM.pdf   une carte reconstruite par chaque SLAM sur le premier
%                                      essai, et la RMSE de reconstruction sur tous les essais
%   Rangs du SLAM.pdf                  precision, NEES, carte et flops contre le rang
%   resultats_slam.tex, rangs_slam.tex les tableaux
%
% chiffres_ch5.txt, ce que le texte cite, va dans resultats/ch5_v2/tableaux et jamais dans le
% manuscrit (sortie.chiffres).
%
% AUCUN TEMPS MESURE (2026-10-26, demande de Bastien). Le cout est compte en flops
% analytiques par essai, comme au chapitre 4 (slam.m, flops_pas), et la memoire en octets.
%
% Les conventions sont celles de figures_campagnes.m : page figee, police calculee pour
% sortir a 12 pt dans le manuscrit, une couleur et un trait par methode.

clearvars -except sortie; close all; clc;
here = fileparts(mfilename('fullpath'));
if ~exist('sortie', 'var'), sortie = struct(); end
if ~isfield(sortie, 'lots'), sortie.lots = fullfile(here, 'resultats', 'ch5_v2', 'lots'); end
if ~isfield(sortie, 'img'),  sortie.img  = fullfile(fileparts(here), 'Images'); end
if ~isfield(sortie, 'tab'),  sortie.tab  = fullfile(fileparts(here), 'Latex', 'Manuscrit', 'include'); end
if ~isfield(sortie, 'chiffres'), sortie.chiffres = fullfile(here, 'resultats', 'ch5_v2', 'tableaux'); end
for d = {sortie.img, sortie.tab, sortie.chiffres}
    if ~exist(d{1}, 'dir'), mkdir(d{1}); end
end
textwidth = 425;                          % pt, mesure sur le document
linewidth = 1;
pas_trace = 3;                            % un point sur trois sur la NEES, comme au chapitre 4

%% Catalogue des lots
L = dir(fullfile(sortie.lots, '*.mat'));
[~, o] = sort({L.name});  L = L(o);       % horodates : le plus recent en dernier
fichiers = containers.Map();
for i = 1:numel(L)
    S = load(fullfile(L(i).folder, L(i).name), 'par');
    if isfield(S.par, 'etiquette') && ~isempty(S.par.etiquette)
        fichiers(S.par.etiquette) = fullfile(L(i).folder, L(i).name);
    end
end
fprintf('Lots trouves : %s\n', strjoin(fichiers.keys, ', '));
lots = containers.Map();
for c = fichiers.keys
    lots(c{1}) = load(fichiers(c{1}), 'par', 'res', 'scen', 'P_eval', 'h_eval', 'pres_traj');
end
% LA NEES DU MANUSCRIT (2026-10-26) : mediane sur tous les essais, normalisee par la mediane
% de chi2_{d_x}, pour qu'un filtre coherent vaille 1, comme au chapitre 4. Les lots v2
% portent la moyenne divisee par d_x de l'ancien campagne_chapitre_5.m : on la recalcule
% depuis les distances de Mahalanobis, enregistrees essai par essai.
for c = lots.keys
    S = lots(c{1});
    for e = 1:numel(S.res)
        D = S.res(e).d2;  D(~isfinite(D)) = NaN;
        S.res(e).nees = median(D, 2, 'omitnan') / chi2inv(0.5, size(S.scen.x_true, 1));
    end
    lots(c{1}) = S;
end

% Chaque configuration : cle de style, etiquette du lot, methode, valeur balayee (NaN : aucune).
configs = {
    'carte',      'principal',      'carte',    NaN
    'akkf_carte', 'akkf_carte_lz1', 'akkf',     1
    'inertiel',   'principal',      'inertiel', NaN
    'rbpf',       'principal',      'rbpf',     NaN
    'akkf',       'akkf_lz1',       'akkf',     2};
variantes = {
    'akkf_orig',   'principal',      'akkf',     NaN
    'akkf_pv',     'akkf_variantes', 'akkf',     2
    'akkf_lz1',    'akkf_lz1',       'akkf',     0
    'akkf_carte01', 'akkf_variantes', 'akkf',    1};

R = struct('cle', {}, 'res', {}, 'lot', {});
ref = [];
for ic = 1:size(configs, 1)
    [r, S] = choisir(lots, configs(ic, :));
    if isempty(r), continue; end
    if isempty(ref)
        ref = S;
    elseif ~isequal(S.scen.x_true, ref.scen.x_true) || S.par.mc.n_runs ~= ref.par.mc.n_runs
        warning('figures_chapitre_5:lots', '%s ne partage pas le scenario de reference : ignore.', ...
                configs{ic, 2});
        continue
    end
    R(end + 1) = struct('cle', configs{ic, 1}, 'res', r, 'lot', configs{ic, 2}); %#ok<SAGROW>
end
V = struct('cle', {}, 'res', {}, 'lot', {});
for ic = 1:size(variantes, 1)
    r = choisir(lots, variantes(ic, :));
    if ~isempty(r)
        V(end + 1) = struct('cle', variantes{ic, 1}, 'res', r, 'lot', variantes{ic, 2}); %#ok<SAGROW>
    end
end

ch = fopen(fullfile(sortie.chiffres, 'chiffres_ch5.txt'), 'w', 'n', 'UTF-8');
ecrire = @(varargin) fprintf(ch, varargin{:});
ecrire('Chiffres du chapitre 5, %s\nLots : %s\n\n', char(datetime('now')), ...
       strjoin(cellfun(@(c) sprintf('%s = %s', c, fichiers(c)), fichiers.keys, ...
       'UniformOutput', false), sprintf('\n       ')));

if ~isempty(ref)
    par    = ref.par;
    x_true = ref.scen.x_true;
    n_it   = size(R(1).res.err, 1);
    d_x    = size(x_true, 1);
    kk     = (1:n_it)';
    tiers  = {1:round(n_it / 3), round(n_it / 3) + 1:round(2 * n_it / 3), round(2 * n_it / 3) + 1:n_it};

    %% Le champ, le circuit et le domaine de la base de Laplace
    pdf_w = 348;  pdf_h = 280;
    fs    = 12 * pdf_w / (0.7 * textwidth);
    g = (1:201) * 1000;
    [Gx, Gy] = meshgrid(g, g);
    h = champ(par, Gx, Gy);
    dom = [min(x_true(1, :)) - par.slam.marge, max(x_true(1, :)) + par.slam.marge, ...
           min(x_true(2, :)) - par.slam.marge, max(x_true(2, :)) + par.slam.marge] / 1e3;
    vue = dom + 8 * [-1 1 -1 1];
    ecrire(['Champ : A = %g nT, poids %s, longueurs d''onde %s km, ecart-type %.1f nT sur la ' ...
            'carte, moyenne retiree %.2f nT\nDomaine U : %.1f x %.1f km\n\n'], ...
           par.map.amplitude, mat2str(par.map.poids), mat2str(par.map.longueurs / 1e3), ...
           std(h(:)), par.slam.moyenne, dom(2) - dom(1), dom(4) - dom(3));

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
    cb.Layout.Tile = 'east';
    cb.TickLabelInterpreter = 'latex';  cb.FontSize = fs;
    ylabel(cb, '$h(x)$ (nT)', 'Interpreter', 'latex', 'FontSize', fs);
    xlabel(ax, '$p^{(x)}$ (km)', 'Interpreter', 'latex', 'FontSize', fs);
    ylabel(ax, '$p^{(y)}$ (km)', 'Interpreter', 'latex', 'FontSize', fs);
    set(ax, 'TickLabelInterpreter', 'latex', 'FontSize', fs, 'Layer', 'top');
    box(ax, 'on');
    exporter(fig, 'Champ et circuit du SLAM', pdf_w, pdf_h, sortie.img);

    %% RMSE et NEES au cours du vol
    kd = (1:pas_trace:n_it)';
    bas  = chi2inv((1 - par.conv.confiance) / 2, d_x) / chi2inv(0.5, d_x);
    haut = chi2inv(1 - (1 - par.conv.confiance) / 2, d_x) / chi2inv(0.5, d_x);
    courbes(kk, kd, R, [bas haut], textwidth, linewidth, sortie.img, 'RMSE et NEES du SLAM');

    %% Les cartes reconstruites
    i_rbpf = find(strcmp({R.cle}, 'rbpf'));  i_akkf = find(strcmp({R.cle}, 'akkf'));
    if ~isempty(i_rbpf) && ~isempty(i_akkf)
        cartes(ref, R([i_rbpf i_akkf]), textwidth, linewidth, sortie.img, 'Cartes reconstruites du SLAM');
    end

    %% Le tableau des erreurs moyennes
    fid = fopen(fullfile(sortie.tab, 'resultats_slam.tex'), 'w', 'n', 'UTF-8');
    fprintf(fid, '%% Generated by figures_chapitre_5.m, do not edit by hand.\n');
    fprintf(fid, ['%% Scenario %s, %d runs, %d particles. Position errors are averaged over\n' ...
                  '%% three thirds of the mission, the NEES over the last one, and the field\n' ...
                  '%% error over the evaluation points within %g km of the trajectory.\n'], ...
            par.scenario, par.mc.n_runs, par.mc.n_part, 2 * par.slam.ell / 1e3);
    for ic = 1:numel(R)
        fprintf(fid, '%% %s read from lot %s\n', R(ic).cle, R(ic).lot);
    end
    fprintf(fid, '\n\\begin{tabular}{|l|ccc|c|c|}\\hline\n');
    fprintf(fid, ['\\multirow{2}{*}{Filter} & \\multicolumn{3}{c|}{$\\mathrm{ARMSE}_{\\hat{x}}$ (m)}' ...
                  ' & \\multirow{2}{*}{NEES} & $\\mathrm{ARMSE}_{\\hat{h}}$\\\\\\cline{2-4}\n']);
    fprintf(fid, ' & $[\\![%d,%d]\\!]$ & $[\\![%d,%d]\\!]$ & $[\\![%d,%d]\\!]$ & & (nT)\\\\\\hline\n', ...
            tiers{1}([1 end]), tiers{2}([1 end]), tiers{3}([1 end]));
    for ic = 1:numel(R)
        [~, ~, nom] = style_methode(R(ic).cle);
        r = R(ic).res;
        fprintf(fid, '%s', nom);
        for j = 1:3
            fprintf(fid, ' & $%.0f$', mean(r.rmse(tiers{j})));
        end
        fprintf(fid, ' & $%s$', nombre(mean(r.nees(tiers{3}))));
        if any(strcmp(R(ic).cle, {'carte', 'akkf_carte'}))
            fprintf(fid, ' & --');
        else
            fprintf(fid, ' & $%.1f$', r.armse_map_traj);
        end
        fprintf(fid, '\\\\\\hline\n');
    end
    fprintf(fid, '\\end{tabular}\n');
    fclose(fid);

    %% Ce que le texte cite
    % conv : sous le seuil de la campagne ; conv2 : sous les 2000 m du chapitre 4, ceux que
    % le texte cite pour justifier une RMSE sur tous les essais.
    ecrire('%-12s %7s | %6s %6s %6s | %6s %6s %6s %6s %6s | %7s %7s | %5s %5s | %6s %6s\n', ...
           'config', 'ARMSE', 'T1', 'T2', 'T3', 'k=1', 'k=30', 'k=47', 'k=78', 'k=fin', ...
           'NEESmed', 'NEES_T3', 'conv', 'conv2', 'carte', 'horsC');
    for ic = 1:numel(R) + numel(V)
        if ic <= numel(R), x = R(ic); else, x = V(ic - numel(R)); end
        r = x.res;  p = ref.pres_traj;
        ks = min([1 30 47 78 n_it], n_it);
        fin = median(r.err(end - min(9, n_it - 1):end, :), 1);
        ecrire('%-12s %7.0f | %6.0f %6.0f %6.0f | %6.0f %6.0f %6.0f %6.0f %6.0f | %7.2f %7.2f | %4.0f%% %4.0f%% | %6.1f %6.1f\n', ...
               x.cle, mean(r.rmse), mean(r.rmse(tiers{1})), mean(r.rmse(tiers{2})), ...
               mean(r.rmse(tiers{3})), r.rmse(ks), median(r.nees), mean(r.nees(tiers{3})), ...
               100 * mean(fin <= par.conv.seuil), 100 * mean(fin <= 2000), ...
               mean(r.rmse_map(p)), mean(r.rmse_map(~p)));
    end
    ecrire('\n');
    for ic = 1:numel(R) + numel(V)
        if ic <= numel(R), x = R(ic); else, x = V(ic - numel(R)); end
        r = x.res.rmse;  y = x.res.nees;
        n1 = min(100, n_it);
        [m1, k1] = max(r(1:n1));
        apres = min(101, n_it):n_it;
        ecrire(['%-12s max sur [1,%d] : %5.0f m au pas %3d ; apres le pas 100 : RMSE entre %5.0f ' ...
                'et %5.0f m, NEES entre %.2f et %.2f\n'], x.cle, n1, m1, k1, min(r(apres)), ...
               max(r(apres)), min(y(apres)), max(y(apres)));
    end
    i_in = find(strcmp({R.cle}, 'inertiel'));
    if ~isempty(i_rbpf)
        y = R(i_rbpf).res.nees;
        ecrire('\nNEES de RBPF-SLAM : au-dessus de 1 au pas %d, au-dessus de la bande (%.2f) au pas %d\n', ...
               premier(y > 1), haut, premier(y > haut));
    end
    if ~isempty(i_akkf) && ~isempty(i_in)
        a = R(i_akkf).res.rmse;  b = R(i_in).res.rmse;
        jusqua_fin = flipud(cumprod(flipud(double(a < b))));   % 1 si a < b jusqu'a la fin
        ecrire(['AKKF-SLAM sous l''inertiel : a partir du pas %d, definitivement a partir du pas %d ' ...
                '(et %d pas sur %d apres le pas 100)\n'], premier(a < b), ...
               premier(jusqua_fin == 1), sum(a(101:end) < b(101:end)), max(n_it - 100, 0));
        if ~isempty(i_rbpf)
            ecrire('AKKF-SLAM sous RBPF-SLAM : %d pas sur %d\n', sum(a < R(i_rbpf).res.rmse), n_it);
        end
    end
    i_orig = find(strcmp({V.cle}, 'akkf_orig'));
    if ~isempty(i_orig) && ~isempty(i_in)
        a = V(i_orig).res.rmse;  b = R(i_in).res.rmse;
        ecrire('Algorithme de l''article sous l''inertiel : a partir du pas %d\n', premier(a < b & kk > 100));
    end
    ecrire('\n');
end

%% Les rangs
% La figure et le tableau ne portent que les deux methodes du chapitre ; l'algorithme de
% l'article (akkf_orig) n'est que dans chiffres_ch5.txt.
tous = rangs(lots, ref);
pts  = tous(ismember({tous.cle}, {'rbpf', 'akkf'}));
if ~isempty(pts)
    figure_rangs(pts, R, textwidth, linewidth, sortie.img, 'Rangs du SLAM');
    fid = fopen(fullfile(sortie.tab, 'rangs_slam.tex'), 'w', 'n', 'UTF-8');
    % LE NOM DE LA METHODE EN LIGNE D'EN-TETE DE SON BLOC (2026-10-26), et non dans une colonne :
    % avec « Modified AKKF-SLAM », la colonne des filtres rendait le tableau trop large de 109 pt.
    fprintf(fid, '%% Generated by figures_chapitre_5.m, do not edit by hand.\n');
    fprintf(fid, ['%% Analytic flop counts per run, not timings, leading terms only (slam.m,\n' ...
                  '%% flops_pas). Memory counts eight bytes per stored double: the N coefficient\n' ...
                  '%% vectors and covariance matrices of RBPF-SLAM, the single set of AKKF-SLAM\n' ...
                  '%% and three N x N matrices.\n\n']);
    fprintf(fid, '\\begin{tabular}{|c|c|c|c|c|c|}\\hline\n');
    % $\mathrm{ARMSE}_{\hat{h}}$ ET NON « Field ARMSE » (2026-10-26, Bastien) : le tableau
    % debordait, et la notation est celle de l'etat, $\mathrm{ARMSE}_{\hat{x}}$.
    fprintf(fid, ['Rank & $\\mathrm{ARMSE}_{\\hat{x}}$ (m) & NEES & $\\mathrm{ARMSE}_{\\hat{h}}$ (nT) & ' ...
                  'Cost (flops) & Memory (bytes)\\\\\\hline\n']);
    for i = 1:numel(pts)
        if i == 1 || ~strcmp(pts(i).cle, pts(i - 1).cle)
            [~, ~, nom] = style_methode(pts(i).cle);
            fprintf(fid, '\\multicolumn{6}{|c|}{%s}\\\\\\hline\n', nom);
        end
        fprintf(fid, '$%d$ & $%.0f$ & $%s$ & $%.1f$ & %s & %s\\\\\\hline\n', ...
                pts(i).rang, pts(i).armse, nombre(pts(i).nees), pts(i).carte, ...
                flops_tex(pts(i).flops), flops_tex(pts(i).memoire));
    end
    fprintf(fid, '\\end{tabular}\n');
    fclose(fid);
    ecrire('%-11s %6s | %7s %8s %7s | %6s %10s %10s\n', 'filtre', 'rang', 'ARMSE', 'NEES_T3', ...
           'carte', 'reech', 'flops', 'octets');
    for i = 1:numel(tous)
        ecrire('%-11s %6d | %7.0f %8.2f %7.1f | %6.3f %10.3g %10.3g\n', tous(i).cle, tous(i).rang, ...
               tous(i).armse, tous(i).nees, tous(i).carte, tous(i).reech, tous(i).flops, ...
               tous(i).memoire);
    end
end
fclose(ch);
fprintf('FIGURES OK\n');

%% -------------------------------------------------------------- utilitaires

function [r, S] = choisir(lots, cfg)
% La sortie d'une configuration {cle, etiquette, methode, valeur}, ou [] si son lot manque.
    r = [];  S = [];
    if ~lots.isKey(cfg{2}), return; end
    S   = lots(cfg{2});
    sel = strcmp({S.res.methode}, cfg{3});
    if ~isnan(cfg{4}), sel = sel & [S.res.balayage] == cfg{4}; end
    i = find(sel, 1);
    if ~isempty(i), r = S.res(i); end
end

function h = champ(par, X, Y)
% Le champ analytique du pilote (campagne_chapitre_5, fonction carte), sur une grille.
    A = par.map.amplitude;  w = par.map.poids;  L = par.map.longueurs;
    h = A * (w(1) * sin(2 * pi * X / L(1)) + w(2) * cos(2 * pi * Y / L(2)) + ...
             w(3) * sin(2 * pi * (X + Y) / L(3)));
end

function k = premier(masque)
    k = find(masque, 1);
    if isempty(k), k = NaN; end
end

function pts = rangs(lots, ref)
% Un point par (filtre, rang) : RBPF-SLAM depuis le lot principal et les lots rangs_rbpf_*,
% AKKF-SLAM et sa variante depuis leurs balayages.
%
% LE COUT EST ANALYTIQUE (slam.m, flops_pas), par essai : le pas de RBPF-SLAM au taux de
% reechantillonnage moyen du lot, celui d'AKKF-SLAM avec la variance predictive si le lot
% l'emploie. La memoire compte huit octets par double stocke : N vecteurs de coefficients
% et N covariances pour RBPF-SLAM ; un seul jeu de coefficients et sa covariance pour
% AKKF-SLAM, plus trois matrices N x N au plus en meme temps (Sigma, la matrice de Gram et
% le gain ou la matrice de reechantillonnage).
    sl  = slam();
    pts = struct('cle', {}, 'rang', {}, 'armse', {}, 'nees', {}, 'carte', {}, 'reech', {}, ...
                 'flops', {}, 'memoire', {});
    sources = {'rbpf', 'principal', 'slam.r_kl'; 'rbpf', 'rangs_rbpf_bas', 'slam.r_kl'
               'rbpf', 'rangs_rbpf_400', 'slam.r_kl'; 'rbpf', 'rangs_rbpf_500', 'slam.r_kl'
               'akkf', 'rangs_akkf_combine', 'slam.r_rff'; 'akkf_orig', 'rangs_akkf', 'slam.r_rff'};
    for s = 1:size(sources, 1)
        if ~lots.isKey(sources{s, 2}), continue; end
        S = lots(sources{s, 2});
        if ~isempty(ref) && ~isequal(S.scen.x_true, ref.scen.x_true), continue; end
        methode = strtok(sources{s, 1}, '_');
        for i = find(strcmp({S.res.methode}, methode))
            r = S.res(i);
            if isempty(S.par.balayage.champ)
                p = split(sources{s, 3}, '.');  rang = S.par.(p{1}).(p{2});
            else
                rang = r.balayage;
            end
            n_it = numel(r.rmse);  t3 = round(2 * n_it / 3) + 1:n_it;
            N = S.par.mc.n_part;
            if strcmp(methode, 'rbpf')
                fl  = n_it * sl.flops('rbpf', N, rang, r.n_res);
                mem = 8 * N * (rang^2 + rang);
            else
                pv  = isfield(S.par.akkf, 'obs_simulees') && S.par.akkf.obs_simulees == 2;
                fl  = n_it * sl.flops('akkf', N, rang, true, pv);
                mem = 8 * (rang^2 + rang + 3 * N^2);
            end
            pts(end + 1) = struct('cle', sources{s, 1}, 'rang', rang, 'armse', mean(r.rmse), ...
                                  'nees', mean(r.nees(t3)), 'carte', r.armse_map_traj, ...
                                  'reech', r.n_res, 'flops', fl, 'memoire', mem); %#ok<AGROW>
        end
    end
    if isempty(pts), return; end
    [~, o] = sortrows([cellfun(@(c) find(strcmp(c, {'rbpf', 'akkf', 'akkf_orig'})), {pts.cle})', ...
                       [pts.rang]']);
    pts = pts(o);
end

function courbes(kk, kd, R, bande, textwidth, linewidth, img_dir, nom_fig)
% Les filtres sur deux panneaux empiles, l'erreur de position en haut et la NEES en bas, avec
% une seule legende au-dessus. La figure est inseree a \textwidth.
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

function cartes(ref, R, textwidth, linewidth, img_dir, nom_fig)
% Comme dans l'article FUSION 2026 : en haut, la carte reconstruite par chaque SLAM sur le
% premier essai, avec la trajectoire vraie en noir et celle qu'il a estimee en pointille
% bleu ; en bas, la RMSE de reconstruction sur tous les essais, point par point. Chaque
% ligne a son echelle de couleur, commune aux deux methodes, et celle du haut est l'etendue
% du champ vrai sur le domaine.
    pdf_w = textwidth;  pdf_h = 400;  fs = 12;
    P = ref.P_eval;  n = round(sqrt(size(P, 2)));
    xs = reshape(P(1, :), n, n);  ys = reshape(P(2, :), n, n);
    xs = xs(1, :) / 1e3;  ys = ys(:, 1)' / 1e3;
    xt = ref.scen.x_true(1:2, :) / 1e3;
    lim_h = [min(ref.h_eval), max(ref.h_eval)];
    lim_e = [0, max(arrayfun(@(x) prctile(x.res.rmse_map, 99), R))];

    fig = figure('Color', 'w', 'Units', 'points', 'Position', [40 40 pdf_w pdf_h]);
    t   = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    for ligne = 1:2
        % Une grille par ligne, pour que chacune porte sa barre de couleur sans retrecir
        % la carte de droite.
        tl = tiledlayout(t, 1, numel(R), 'TileSpacing', 'compact', 'Padding', 'none');
        tl.Layout.Tile = ligne;
        for j = 1:numel(R)
            ax = nexttile(tl);
            if ligne == 1
                imagesc(ax, xs, ys, reshape(R(j).res.h_exemple, n, n), lim_h);
            else
                imagesc(ax, xs, ys, reshape(R(j).res.rmse_map, n, n), lim_e);
            end
            hold(ax, 'on');
            set(ax, 'YDir', 'normal'); axis(ax, 'image'); colormap(ax, parula);
            plot(ax, xt(1, :), xt(2, :), 'k-', 'LineWidth', linewidth);
            if ligne == 1
                tr = R(j).res.track / 1e3;
                plot(ax, tr(1, :), tr(2, :), ':', 'Color', [0 0 1], 'LineWidth', linewidth);
                [~, ~, nom] = style_methode(R(j).cle);
                title(ax, nom, 'Interpreter', 'latex', 'FontSize', fs, 'FontWeight', 'normal');
            end
            hold(ax, 'off');
            xlim(ax, xs([1 end])); ylim(ax, ys([1 end]));
            set(ax, 'TickLabelInterpreter', 'latex', 'FontSize', fs, 'Layer', 'top');
            if ligne == 2, xlabel(ax, '$p^{(x)}$ (km)', 'Interpreter', 'latex', 'FontSize', fs);
            else, set(ax, 'XTickLabel', []); end
            if j == 1, ylabel(ax, '$p^{(y)}$ (km)', 'Interpreter', 'latex', 'FontSize', fs);
            else, set(ax, 'YTickLabel', []); end
            box(ax, 'on');
        end
        cb = colorbar(ax);
        cb.Layout.Tile = 'east';
        cb.TickLabelInterpreter = 'latex';  cb.FontSize = fs;
        if ligne == 1, et = '$\hat{h}(x)$ (nT)'; else, et = '$\mathrm{RMSE}_{\hat{h}}$ (nT)'; end
        ylabel(cb, et, 'Interpreter', 'latex', 'FontSize', fs);
    end
    exporter(fig, nom_fig, pdf_w, pdf_h, img_dir);
end

function figure_rangs(pts, R, textwidth, linewidth, img_dir, nom_fig)
% Quatre panneaux contre le rang de la base, en echelle logarithmique : erreur de position
% sur tout le vol, NEES du dernier tiers, erreur de carte le long du vol, et flops par essai.
% L'inertiel et la carte nulle donnent les references en gris.
    pdf_w = textwidth;  pdf_h = 400;  fs = 12;
    cles = unique({pts.cle}, 'stable');
    i_in = find(strcmp({R.cle}, 'inertiel'), 1);
    champs = {'armse', 'nees', 'carte', 'flops'};
    etiq   = {'$\mathrm{ARMSE}_{\hat{x}}$ (m)', 'NEES', '$\mathrm{ARMSE}_{\hat{h}}$ (nT)', ...
              'Cost per run (flops)'};
    fig = figure('Color', 'w', 'Units', 'points', 'Position', [40 40 pdf_w pdf_h]);
    t   = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    h_leg = gobjects(0);  leg = {};
    for p = 1:4
        ax = nexttile(t); hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on');
        for c = 1:numel(cles)
            q = pts(strcmp({pts.cle}, cles{c}));
            [col, l, nom] = style_methode(cles{c});
            hp = plot(ax, [q.rang], [q.(champs{p})], l, 'Color', col, 'LineWidth', linewidth, ...
                      'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', col);
            if p == 1, h_leg(end + 1) = hp; leg{end + 1} = nom; end %#ok<AGROW>
        end
        if ~isempty(i_in) && any(p == [1 3])
            if p == 1, v = mean(R(i_in).res.rmse); else, v = R(i_in).res.armse_map_traj; end
            hp = yline(ax, v, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', linewidth);
            if p == 1, h_leg(end + 1) = hp; leg{end + 1} = 'Inertial navigation'; end %#ok<AGROW>
        end
        hold(ax, 'off');
        set(ax, 'XScale', 'log', 'TickLabelInterpreter', 'latex', 'FontSize', fs);
        if any(p == [2 4]), set(ax, 'YScale', 'log'); end
        ylabel(ax, etiq{p}, 'Interpreter', 'latex', 'FontSize', fs);
        if p >= 3, xlabel(ax, 'Rank', 'Interpreter', 'latex', 'FontSize', fs);
        else, set(ax, 'XTickLabel', []); end
    end
    lg = legend(h_leg, leg, 'Interpreter', 'latex', 'FontSize', fs, 'NumColumns', 2);
    lg.ItemTokenSize = [18 18];
    lg.Layout.Tile = 'north';
    exporter(fig, nom_fig, pdf_w, pdf_h, img_dir);
end

function [c, l, nom] = style_methode(cle)
% COULEUR, TRAIT ET NOM D'UNE METHODE. La carte connue garde le noir pointille de
% figures_campagnes.m, et l'AKKF qui la lit l'ambre de la methode proposee ; l'etat de l'art
% prend le rouge, l'inertiel un gris tirete. Les variantes d'AKKF-SLAM gardent l'ambre.
    switch cle
        case 'carte',      c = [0.00 0.00 0.00];  l = ':';   nom = 'Known map';
        case 'akkf_carte', c = [0.95 0.60 0.00];  l = ':';   nom = 'Known map (AKKF)';
        case 'ok',         c = [0.15 0.35 0.70];  l = '-';   nom = 'Static kriging';
        case 'inertiel',   c = [0.45 0.45 0.45];  l = '--';  nom = 'Inertial navigation';
        case 'rbpf',       c = [0.75 0.02 0.15];  l = '-';   nom = 'RBPF-SLAM';
        case 'akkf',       c = [0.95 0.60 0.00];  l = '-';   nom = 'AKKF-SLAM';
        case 'akkf_orig',  c = [0.60 0.35 0.00];  l = '-';   nom = 'AKKF-SLAM of the original article';
        case 'akkf_pv',    c = [0.95 0.60 0.00];  l = '--';  nom = 'AKKF-SLAM ($\alpha_{\mathrm{RKHS},z}=0.1$)';
        case 'akkf_lz1',   c = [0.95 0.60 0.00];  l = '-.';  nom = 'AKKF-SLAM (noise variance only)';
        case 'akkf_carte01', c = [0.95 0.60 0.00]; l = ':';  nom = 'Known map (AKKF, $\alpha_{\mathrm{RKHS},z}=0.1$)';
        otherwise
            error('figures_chapitre_5:style', 'Cle de methode inconnue : %s.', cle);
    end
end

function s = nombre(v)
% Deux decimales sous dix, un entier au-dela.
    if v < 10, s = sprintf('%.2f', v); else, s = sprintf('%.0f', v); end
end

function t = flops_tex(v)
% Repris de figures_campagnes.m : trois chiffres significatifs et le prefixe decimal, separe
% de la valeur par une espace ($333$ k, pas $333$k).
    units = {'', 'k', 'M', 'G', 'T', 'P'};
    i     = max(1, min(numel(units), 1 + floor(log10(v) / 3)));
    t     = sprintf('$%.3g$ %s', v / 1000^(i - 1), units{i});
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
