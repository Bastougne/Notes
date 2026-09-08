% FIGURES_MULTIMODAL
%
% Le regime multimodal en fonction de la maille du releve, a sigma_0 = 12 km et capteur
% inchange a 5 nT. Trois figures sous Images/, en PDF vectoriel.
%
%   Fraction multimodale  part des pas ou le mean-shift trouve plus d'un mode, avec la
%                         bande de la longueur de correlation du champ
%   Convergence           les trois estimateurs contre la maille
%   Erreur finale         idem, mediane sur les essais convergés
%
% Page figee comme les autres figures du manuscrit, police derivee de la fraction
% d'insertion declaree.

clear; close all; clc;
here = fileparts(mfilename('fullpath'));

img_dir   = fullfile(fileparts(here), 'Images');
textwidth = 425;
tex_frac  = 0.7;
pdf_w     = 348;
pdf_h     = 280;
fontsize  = 12 * pdf_w / (tex_frac * textwidth);
linewidth = 1;

ell_bas = 4.5;                   % km, longueur de correlation ajustee, min sur les tirages
ell_haut = 5.7;                  % km, max

%% Rassembler les lots complets du regime multimodal
% TOUT LOT MULTIMODAL, QUEL QUE SOIT LE NOMBRE D'ESSAIS. Le motif etait fige sur
% '2026090*_reference_libre_100runs.mat' : le nombre d'essais ET la date etaient en dur,
% donc le lot a 1000 essais ne serait jamais apparu et la figure serait restee muette
% sans rien signaler. Ce qui identifie le regime multimodal est l'incertitude initiale
% portee a 12 km, testee plus bas — ni la date, ni le nombre d'essais.
d = dir(fullfile(here, 'resultats', '*_reference_libre_*runs.mat'));
d = d([d.datenum] > datenum([2026 9 1 23 30 0]));   % avant, le releve n'etait pas secoue

pas = []; frac = []; conv = []; fin = [];
for f = 1:numel(d)
    S = load(fullfile(here, 'resultats', d(f).name));
    par = S.par; res = S.res;
    if par.sigma_obs ~= 5 || par.dt ~= 5 || par.sigma_0(1) < 10e3, continue, end
    noms = cellfun(@char, {res.model}, 'UniformOutput', false);
    io = trouve(noms, 'totale');  ia = trouve(noms, 'adaptatif');  ic = trouve(noms, 'CA-OK');
    if isempty(io) || isempty(ia) || isempty(ic), continue, end

    n_it = numel(res(1).rmse);
    m    = mean(res(ic).n_modes, 2);
    kv   = find(m < 1.1, 1);  if isempty(kv), kv = n_it; end

    c = zeros(1, 3); e = zeros(1, 3);
    for t = 1:3
        j  = [io ia ic];  j = j(t);
        E  = res(j).err;
        ok = median(E(n_it-par.conv.n_fin+1:n_it, :), 1) <= par.conv.seuil;
        c(t) = 100 * mean(ok);
        e(t) = median(E(end, ok));
    end
    pas(end+1)  = par.krig.pas / 1e3; %#ok<AGROW>
    frac(end+1) = 100 * kv / n_it;    %#ok<AGROW>
    conv(end+1, :) = c;               %#ok<AGROW>
    fin(end+1, :)  = e;               %#ok<AGROW>
end
[pas, o] = sort(pas);  frac = frac(o);  conv = conv(o, :);  fin = fin(o, :);
fprintf('%d mailles : %s\n', numel(pas), mat2str(pas));

leg = {'static OK', 'adaptive OK', 'CA-OK'};
sty = {'-o', '--s', ':^'};

%% 1. Fraction multimodale, avec la bande de la longueur de correlation
fig = figure('Color', 'w', 'Units', 'points', 'Position', [80 80 pdf_w pdf_h]);
hold on; box on; grid on;
fill([ell_bas ell_haut ell_haut ell_bas], [0 0 105 105], [1 0.85 0.88], ...
     'EdgeColor', 'none', 'FaceAlpha', 0.55, 'HandleVisibility', 'off');
plot(pas, frac, '-o', 'LineWidth', linewidth, 'MarkerFaceColor', 'auto');
hold off;
xlabel('survey pitch (km)', 'Interpreter', 'latex', 'FontSize', fontsize);
ylabel('multimodal steps (\%)', 'Interpreter', 'latex', 'FontSize', fontsize);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fontsize);
xlim([min(pas)-0.2 max(pas)+0.2]); ylim([0 105]);
exporter(fig, 'Fraction multimodale contre maille', pdf_w, pdf_h, img_dir);

%% 2. Convergence
fig = figure('Color', 'w', 'Units', 'points', 'Position', [80 80 pdf_w pdf_h]);
hold on; box on; grid on;
for t = 1:3
    plot(pas, conv(:, t), sty{t}, 'LineWidth', linewidth);
end
hold off;
xlabel('survey pitch (km)', 'Interpreter', 'latex', 'FontSize', fontsize);
ylabel('converged runs (\%)', 'Interpreter', 'latex', 'FontSize', fontsize);
legend(leg, 'Interpreter', 'latex', 'FontSize', fontsize, 'Location', 'southwest');
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fontsize);
xlim([min(pas)-0.2 max(pas)+0.2]);
exporter(fig, 'Convergence contre maille multimodale', pdf_w, pdf_h, img_dir);

%% 3. Erreur finale
fig = figure('Color', 'w', 'Units', 'points', 'Position', [80 80 pdf_w pdf_h]);
hold on; box on; grid on;
for t = 1:3
    plot(pas, fin(:, t), sty{t}, 'LineWidth', linewidth);
end
hold off;
xlabel('survey pitch (km)', 'Interpreter', 'latex', 'FontSize', fontsize);
ylabel('median final error (m)', 'Interpreter', 'latex', 'FontSize', fontsize);
legend(leg, 'Interpreter', 'latex', 'FontSize', fontsize, 'Location', 'northwest');
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fontsize);
xlim([min(pas)-0.2 max(pas)+0.2]);
exporter(fig, 'Erreur finale contre maille multimodale', pdf_w, pdf_h, img_dir);

fprintf('Trois figures ecrites dans %s\n', img_dir);

%% ------------------------------------------------------------ fonctions locales

function i = trouve(noms, motif)
    i = find(~cellfun(@isempty, strfind(noms, motif)), 1); %#ok<STRIFCND>
end

function exporter(fig, nom, pdf_w, pdf_h, img_dir)
    if ~exist(img_dir, 'dir'), mkdir(img_dir); end
    set(fig, 'PaperUnits', 'points', 'PaperSize', [pdf_w pdf_h], ...
             'PaperPosition', [0 0 pdf_w pdf_h]);
    print(fig, fullfile(img_dir, [nom '.pdf']), '-dpdf', '-vector');
end
