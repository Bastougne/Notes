% TEMPO_MATERN
%
% Radial profile of the Matérn kernels of half-integer smoothness and of their
% squared exponential limit, drawn with unit variance and unit correlation
% length so that only the shape near the origin and the tail weight differ.
%
% Closed forms are those of the manuscript, chapter 3, th:Matérn and the
% half-integer expression that follows it.

clear; close all; clc;

%% Paramètres
ell       = 1;                   % longueur caractéristique
r_max     = 5;                   % distance maximale tracée
n_plot    = 1000;                % résolution
img_dir   = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'Images');
                                 % le manuscrit et les notes y lisent leurs figures,
                                 % résolu depuis le script et non depuis pwd
img_name  = 'Noyaux de Matérn';
linewidth = 1;                   % partagé avec les scripts de krigeage
fontsize  = 14;                  % 14 x 298/348 = 12,0 pt à l'inclusion
pdf_w     = 348;                 % page figée en points, partagée par les 3 scripts
pdf_h     = 280;

r = linspace(0, r_max, n_plot);

%% Noyaux de Matérn, formes fermées pour nu demi-entier
k_12  = exp(-r ./ ell);
k_32  = (1 + sqrt(3) .* r ./ ell) .* exp(-sqrt(3) .* r ./ ell);
k_52  = (1 + sqrt(5) .* r ./ ell + 5 .* r.^2 ./ (3 .* ell.^2)) ...
        .* exp(-sqrt(5) .* r ./ ell);
k_inf = exp(-r.^2 ./ (2 .* ell.^2));   % nu -> inf : exponentiel quadratique

%% Tracé
fig = figure('Color', 'w', 'Units', 'points', 'Position', [80 80 pdf_w pdf_h]);
hold on; box on; grid on;

plot(r, k_12,  ':',  'LineWidth', linewidth);
plot(r, k_32,  '-.', 'LineWidth', linewidth);
plot(r, k_52,  '--', 'LineWidth', linewidth);
plot(r, k_inf, '-',  'LineWidth', linewidth);

xlabel('$r=\|x-x''\|$', 'Interpreter', 'latex', 'FontSize', fontsize);
ylabel('$k_\nu(r)$',    'Interpreter', 'latex', 'FontSize', fontsize);

legend({'$\nu=1/2$', '$\nu=3/2$', '$\nu=5/2$', '$\nu=+\infty$ (SE)'}, ...
       'Interpreter', 'latex', 'FontSize', fontsize, 'Location', 'northeast');

set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fontsize);
xlim([0 r_max]); ylim([-0.02 1.02]);
hold off;

%% Export
% print plutôt que exportgraphics : ce dernier recadre au contenu, donc la page
% varierait d'une figure à l'autre. Ici la page vaut exactement pdf_w x pdf_h,
% ce qui garantit une même échelle, donc une même police apparente, partout.
if ~exist(img_dir, 'dir')
    mkdir(img_dir);
end
set(fig, 'PaperUnits', 'points', 'PaperSize', [pdf_w pdf_h], ...
         'PaperPosition', [0 0 pdf_w pdf_h]);
print(fig, fullfile(img_dir, [img_name '.pdf']), '-dpdf', '-vector');
