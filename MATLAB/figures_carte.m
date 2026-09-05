% FIGURES_CARTE
%
% La carte d'anomalie magnétique du chapitre 4, avec la trajectoire de référence par-dessus,
% exportée en PDF vectoriel sous Images/.
%
% Le champ est celui que lisent tous les essais du chapitre. La trajectoire est relue d'une
% campagne enregistrée plutôt que recalculée, pour que la figure montre exactement le vol
% qui a produit les résultats et non une reconstruction qui pourrait en diverger.
%
% La page est figée comme les autres figures du manuscrit, et la police en découle. Voir
% figures_chapitre_3.m pour la chaîne complète.

clear; close all; clc;
here = fileparts(mfilename('fullpath'));

%% Ce que la figure déclare de son insertion
img_dir   = fullfile(fileparts(here), 'Images');
textwidth = 425;                 % pt, mesuré sur le document
tex_frac  = 0.7;                 % \includegraphics[width=0.7\textwidth]
pdf_w     = 348;                 % page figée, en points
pdf_h     = 280;
fontsize  = 12 * pdf_w / (tex_frac * textwidth);   % 14,0, donc 12 pt dans le manuscrit
linewidth = 1;

%% La carte
carte = load(fullfile(here, '..', '..', 'MATLAB 2A - sauvegarde', 'Cartes', ...
                      'carte_magnetometrie_anomalie_mexique.mat'));

% TRANSPOSEE, comme campagne_chapitre_4.m. Dans la sauvegarde, mesures indexe (x, y) ;
% grid_read veut des lignes en y et des colonnes en x, ce qu'attend aussi imagesc. Sans
% la transposition la figure montre le champ retourné sur sa diagonale, avec une
% trajectoire correcte par-dessus — invisible ici, la carte étant carrée.
h    = double(carte.mesures).';
pas  = double(carte.pas(:))' .* [1 1];                 % m

% Le noeud (1,1) est en (pas, pas) et non en (0,0) : c'est l'origine que grid_read
% suppose pour les grilles des sauvegardes, et le relevé est place dessus.
x_km = (1 : size(h, 2)) * pas(1) / 1e3;
y_km = (1 : size(h, 1)) * pas(2) / 1e3;

fprintf('Carte %d x %d noeuds a %.0f m, champ de %.0f a %.0f nT, ecart-type %.0f nT\n', ...
        size(h, 2), size(h, 1), pas(1), min(h(:)), max(h(:)), std(h(:)));

%% Le relevé, construit par krigeage.m et non redérivé ici
%
% Passer par krig.releve garantit que les croix tombent exactement où le filtre lit ses
% échantillons : centres de cellules, meme origine, et surtout meme secouage. Un réseau
% exact rendrait le biais d'interpolation périodique, donc répété d'une branche à l'autre
% du vol au lieu de se moyenner ; le manuscrit annonce un secouage de 30 % de la maille,
% c'est celui-ci.
par.unit               = 'nT';
par.run.seed           = 123456789;
par.krig.pas           = 7.5e3;      % m, la maille. Compromis de lisibilité : à 5 km les
                                     % croix saturent la carte, à 10 km elles la désertent
par.krig.jitter        = 0.30;       % fraction de maille, comme annoncé au manuscrit
par.krig.sigma_map     = 10;         % nT, scénario de référence
par.krig.n_ml          = 500;
par.krig.ell_init      = 10e3;
par.krig.noyau         = 'matern32';

krig   = krigeage();
survey = krig.releve(struct('h', h, 'step', pas), par);
srv_km = survey.X / 1e3;

%% La trajectoire, relue d'une campagne
d = dir(fullfile(here, 'resultats', '*_A_*runs.mat'));
if isempty(d)
    error('figures_carte:trajectoire', ...
          'Aucune campagne A enregistree sous resultats/ : la trajectoire ne peut pas etre tracee.');
end
[~, i] = max([d.datenum]);
S = load(fullfile(here, 'resultats', d(i).name), 'scen');
vol = S.scen.x_true(1:2, :) / 1e3;                     % km
fprintf('Trajectoire relue de %s, %d pas\n', d(i).name, size(vol, 2) - 1);

%% La figure
fig = figure('Color', 'w', 'Units', 'points', 'Position', [80 80 pdf_w pdf_h]);
imagesc(x_km, y_km, h); hold on;
set(gca, 'YDir', 'normal');
axis image;
% Traits fins mais pas sous-pixel : a 0.25 pt les croix se fondent dans le fond a
% l'impression et virent au sombre au lieu de rester rouges.
plot(srv_km(1, :), srv_km(2, :), 'x', 'Color', [0.85 0 0], ...
     'MarkerSize', 2.5, 'LineWidth', 0.6);
plot(vol(1, :), vol(2, :), 'k-', 'LineWidth', linewidth);
plot(vol(1, 1), vol(2, 1), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
hold off;

colormap(gca, parula);
c = colorbar;
c.Label.String = '$h(x)$ (nT)';
c.Label.Interpreter = 'latex';
c.Label.FontSize = fontsize;
c.TickLabelInterpreter = 'latex';
c.FontSize = fontsize;

% Les axes portent les composantes de position de l'etat du scenario, pas des indices.
xlabel('$p^{(x)}$ (km)', 'Interpreter', 'latex', 'FontSize', fontsize);
ylabel('$p^{(y)}$ (km)', 'Interpreter', 'latex', 'FontSize', fontsize);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fontsize, 'Layer', 'top');
box on;

% La barre de couleur et son label vivent HORS des axes, et MATLAB dimensionne les axes
% avant de savoir ce que le label mesure : sur une page figee, il sort donc du PDF. Les
% deux positions sont imposees plutot que laissees automatiques, la largeur reservee a
% droite couvrant les graduations et le label. axis image ne s'y oppose pas, il ne fixe
% que les rapports d'aspect et laisse la boite de trace se centrer dans Position.
ax = gca;
drawnow;
ax.Units = 'normalized';
c.Units  = 'normalized';

% Les marges de gauche et du bas sont MESUREES et non devinees : TightInset donne la place
% que reclament le label et les graduations de chaque cote, donc une valeur en dur finit
% par rogner des que la police ou les chiffres d'axe changent. Seule la bande de droite
% reste imposee, TightInset ignorant la barre de couleur, qui vit hors des axes.
ti     = ax.TightInset;
marg   = 0.015;
cb_w   = 0.035;                  % la barre elle-meme
cb_pad = 0.21;                   % ce qui vit A SA DROITE : graduations puis label
gap    = 0.02;                   % entre les axes et la barre
left   = ti(1) + marg;
bot    = ti(2) + marg;
top    = ti(4) + marg;
cb_x   = 1 - cb_pad - cb_w;      % la barre est placee depuis le bord droit, pas depuis
                                 % les axes : c'est son label qui touche le bord de page
ax.Position = [left, bot, cb_x - gap - left, 1 - bot - top];
c.Position  = [cb_x, bot, cb_w, 1 - bot - top];
fprintf('Marges mesurees : gauche %.3f, bas %.3f, haut %.3f\n', ti(1), ti(2), ti(4));

%% Export
if ~exist(img_dir, 'dir'), mkdir(img_dir); end
set(fig, 'PaperUnits', 'points', 'PaperSize', [pdf_w pdf_h], ...
         'PaperPosition', [0 0 pdf_w pdf_h]);
print(fig, fullfile(img_dir, 'Carte anomalie magnetique.pdf'), '-dpdf', '-vector');
fprintf('Ecrit dans %s\n', fullfile(img_dir, 'Carte anomalie magnetique.pdf'));
