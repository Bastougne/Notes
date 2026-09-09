% FIGURES_TEMPORELLES
%
% Les métriques du chapitre 4 au cours du vol, tracées depuis les campagnes enregistrées.
% Deux figures exportées séparément sous Images/, en PDF vectoriel.
%
%   Modes            ce que trouve le mean-shift, pas par pas
%   Fenêtre          échantillons retenus par requête, avec le plancher n_min
%
% Ce qui reste ici n'est que ce dont les sections AK et CAK ont besoin. Les RMSE, NEES et
% cohérence au cours du vol qu'il portait ont été retirés le 2026-09-04 : leurs figures
% ne sont plus dans le manuscrit, remplacées par les grilles par maille de
% figures_campagnes.m, qui est le script des courbes validées.
%
% La page est figée comme les autres figures du manuscrit et la police en découle. Voir
% figures_chapitre_3.m pour la chaîne complète.

clear; close all; clc;
here = fileparts(mfilename('fullpath'));

img_dir   = fullfile(fileparts(here), 'Images');
textwidth = 425;                 % pt, mesuré sur le document
tex_frac  = 0.7;                 % \includegraphics[width=0.7\textwidth]
pdf_w     = 348;                 % page figée, en points
pdf_h     = 280;
fontsize  = 12 * pdf_w / (tex_frac * textwidth);   % 14,0, donc 12 pt dans le manuscrit
linewidth = 1;

%% Les données
B = recent(here, '*_reference_B_*runs.mat');
A = recent(here, '*_reference_A_*runs.mat');
fprintf('Campagne B : %s\n', B.nom);
fprintf('Campagne A : %s\n', A.nom);

par   = B.S.par;
res   = B.S.res;
n_it  = numel(res(1).rmse);
k     = (1:n_it)';
seuil = par.conv.seuil;
noms  = cellfun(@char, {res.model}, 'UniformOutput', false);   % char et non string
etiq  = etiquettes(noms);                                       % legendes du manuscrit
styles = {'-', '--', ':', '-.', '-'};

% Le pas où commence le premier demi-tour, relevé sur la trajectoire vraie plutôt que
% recalculé : c'est lui qui lève l'ambiguïté, en échantillonnant le champ dans une
% seconde direction, et les courbes ne se lisent pas sans lui.
k_vir = premier_virage(B.S.scen.x_true);
fprintf('Premier virage au pas %d sur %d\n', k_vir, n_it);

%% 4. Modes trouvés par le mean-shift
j = find(strcmp(noms, 'CA-OK'), 1);
fig = figure('Color', 'w', 'Units', 'points', 'Position', [80 80 pdf_w pdf_h]);
plot(k, mean(res(j).n_modes, 2), '-', 'LineWidth', linewidth);
box on; grid on; hold on; marqueur(k_vir, fontsize); hold off;
set(gca, 'YScale', 'log');
xlabel('$k$', 'Interpreter', 'latex', 'FontSize', fontsize);
ylabel('Number of modes', 'Interpreter', 'latex', 'FontSize', fontsize);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fontsize);
xlim([1 n_it]);
exporter(fig, 'Modes au cours du vol', pdf_w, pdf_h, img_dir);

%% 5. Fenêtre : échantillons retenus par requête, et son plancher
sel = find(ismember(noms, {'OK adaptatif', 'CA-OK'}));
fig = figure('Color', 'w', 'Units', 'points', 'Position', [80 80 pdf_w pdf_h]);
hold on; box on; grid on;
for j = sel
    plot(k, res(j).cout(1, :), styles{j}, 'LineWidth', linewidth);
end
yline(par.krig.n_min, ':', 'LineWidth', linewidth, 'HandleVisibility', 'off');
marqueur(k_vir, fontsize);
hold off; set(gca, 'YScale', 'log');
xlabel('$k$', 'Interpreter', 'latex', 'FontSize', fontsize);
ylabel('$\tilde{N}_\mathrm{w}(k)$', 'Interpreter', 'latex', 'FontSize', fontsize);
legend(etiq(sel), 'Interpreter', 'latex', 'FontSize', fontsize, 'Location', 'northeast');
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fontsize);
xlim([1 n_it]);
exporter(fig, 'Fenetre au cours du vol', pdf_w, pdf_h, img_dir);

%% ------------------------------------------------------------ fonctions locales

function out = recent(here, motif)
% La campagne la plus récente qui corresponde au motif. Les figures suivent ainsi les
% dernières mesures sans qu'on ait à les nommer.
    d = dir(fullfile(here, 'resultats', motif));
    if isempty(d)
        error('figures_temporelles:donnees', 'Aucune campagne ne correspond a %s.', motif);
    end
    [~, i] = max([d.datenum]);
    out.nom = d(i).name;
    out.S   = load(fullfile(here, 'resultats', d(i).name));
end

function k = premier_virage(x_true)
% Le premier pas ou le cap change. Mesure sur la vitesse vraie, donc valable quelle que
% soit la forme du vol : une ligne droite rend simplement un marqueur vide.
    v   = x_true(3:4, :);
    n   = vecnorm(v, 2, 1);
    cos_a = sum(v(:, 1:end-1) .* v(:, 2:end), 1) ./ max(n(1:end-1) .* n(2:end), eps);
    ang = acosd(min(max(cos_a, -1), 1));
    k   = find(ang > 0.5, 1);
    if isempty(k), k = NaN; end
end

function marqueur(k_vir, fontsize)
% Le debut du premier demi-tour. C'est lui qui leve l'ambiguite en echantillonnant le
% champ dans une seconde direction, donc les courbes ne se lisent pas sans lui.
    if isnan(k_vir), return, end
    xl = xline(k_vir, '--', 'first turn', 'LineWidth', 0.75, 'Color', [0.4 0.4 0.4], ...
               'HandleVisibility', 'off', 'Interpreter', 'latex', 'FontSize', fontsize, ...
               'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'left');
    uistack(xl, 'bottom');
end

function c = etiquettes(noms)
% Les legendes du manuscrit, en anglais et sans accent : l'interpreteur LaTeX de MATLAB
% ne rend ni les caracteres accentues ni le souligne nu.
% La correspondance se fait sur un prefixe purement ASCII, les noms de modeles portant
% des accents que l'on ne veut pas avoir a reecrire ici.
    dico = { 'carte',    'True map'
             'bilin',    'Bilinear'
             'OK carte', 'Static OK'
             'OK adapt', 'Adaptive OK'
             'CA-OK',    'CA-OK' };
    c = noms;
    for i = 1:numel(noms)
        j = find(cellfun(@(p) strncmp(noms{i}, p, numel(p)), dico(:, 1)), 1);
        if isempty(j)
            c{i} = strrep(noms{i}, '_', '\_');
        else
            c{i} = dico{j, 2};
        end
    end
end

function exporter(fig, nom, pdf_w, pdf_h, img_dir)
% Page figee, donc meme echelle et meme police apparente pour toutes les figures.
    if ~exist(img_dir, 'dir')
        mkdir(img_dir);
    end
    set(fig, 'PaperUnits', 'points', 'PaperSize', [pdf_w pdf_h], ...
             'PaperPosition', [0 0 pdf_w pdf_h]);
    print(fig, fullfile(img_dir, [nom '.pdf']), '-dpdf', '-vector');
end
