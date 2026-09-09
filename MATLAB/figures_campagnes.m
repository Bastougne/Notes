% FIGURES_CAMPAGNES
%
% Les figures du chapitre 4 tirees des campagnes enregistrees. Un seul script, pour que
% les courbes du manuscrit aient une source unique et suivent d'elles-memes le passage a
% 1000 essais : il relit la campagne la plus recente, sans parametre a retoucher.
%
%   Carte anomalie magnetique  le champ, le releve secoue et la trajectoire de reference
%   ARMSE de reconstruction    l'erreur de carte contre la maille, et sa loi de puissance
%   RMSE par maille            grille 3 x 2, une maille par case, carte connue / krigeage /
%                              bilineaire sur le meme axe
%   NEES par maille            la meme grille, echelle log et bande du chi2
%   Convergence par maille     les deux criteres de convergence, deux panneaux
%   Erreur finale par maille   RMSE et NEES medianes sur les dix derniers pas
%
% et cinq tableaux sous Latex/Manuscrit/include/ :
%
%   cout_maille.tex   le cout analytique par maille          campagne A
%   ak_vs_ok.tex      adaptatif contre statique              campagne B
%   cak_vs_ak.tex     clusterise contre adaptatif            campagne B
%   cout_cak.tex      les trois postes de cout du CA-OK      campagne B
%   fenetrage.tex     les parametres de fenetre, un par ligne  balayages CA-OK
%
% Les six mailles tiennent en une grille plutot qu'en six flottants : la degradation se
% lit d'une case a l'autre, ce qu'une figure isolee par maille ne montre pas. La page est
% pleine largeur et haute, donc \includegraphics[width=\textwidth] sur une page entiere.
%
% figures_carte.m a ete absorbe ici le 7 septembre : une figure isolee dans son propre
% fichier finit par diverger des autres sur la page, la police ou la campagne qu'elle
% relit. Il reste figures_temporelles.m, qui porte les deux figures au cours du vol.
%
% LES TROIS TABLEAUX DE CAMPAGNE ETAIENT ECRITS EN DUR DANS LE MANUSCRIT. Ils sont
% desormais generes, comme cout_maille.tex l'etait deja : douze lots passent a 1000
% essais, et une transcription a la main par lot est une erreur qui attend son tour.

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
etiq    = {'True map', 'Bilinear', 'Kriging'};

d_x  = numel(par.x0);
n_it = size(res(1).err, 1);
kk   = (1:n_it)';
fin  = n_it - par.conv.n_fin + 1 : n_it;
maha = chi2inv(par.conv.confiance, d_x);

RMSE = nan(n_it, numel(pas), numel(modeles));
NEES = nan(n_it, numel(pas), numel(modeles));
CONV = nan(numel(pas), numel(modeles), 2);        % (:,:,1) seuil, (:,:,2) Mahalanobis
FIN  = nan(numel(pas), numel(modeles), 2);        % moyenne sur les essais des criteres
for e = 1:numel(res)
    ip = find(abs(pas - res(e).balayage/1e3) < 1e-9, 1);
    im = find(cellfun(@(m) contains(noms{e}, m), modeles), 1);
    if isempty(ip) || isempty(im), continue, end
    E = res(e).err;
    D = res(e).d2;  D(~isfinite(D)) = NaN;
    % Les deux criteres, essai par essai : ce sont eux que CONV seuille, et dont FIN
    % moyenne les valeurs. Medianer les courbes agregees donnerait autre chose.
    crit_e = median(E(fin, :), 1);
    crit_d = median(D(fin, :), 1) / d_x;
    CONV(ip, im, 1) = 100 * mean(crit_e <= par.conv.seuil);
    CONV(ip, im, 2) = 100 * mean(crit_d * d_x <= maha, 'omitnan');
    % MEDIANE SUR TOUS LES ESSAIS, divergés compris. Une moyenne aurait exige de les
    % ecarter — un nuage effondre rend une distance de Mahalanobis astronomique, et
    % quelques essais suffisent alors a emporter la moyenne, jusqu'a 7,6e15 sur la carte
    % connue. La mediane y est insensible, donc la figure porte les memes essais que les
    % taux de convergence et n'a aucune restriction a justifier.
    FIN(ip, im, 1)  = median(crit_e);
    FIN(ip, im, 2)  = median(crit_d, 'omitnan');

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

%% La carte d'anomalie, avec le releve et la trajectoire de reference
%
% Vient de figures_carte.m, supprime : les figures du chapitre sont produites par ce
% script et par figures_temporelles.m, et une figure isolee dans son propre fichier
% finit par diverger des autres sur la page, la police ou la campagne qu'elle relit.
%
% ELLE A SA PROPRE PAGE. Les grilles par maille occupent \textwidth sur une page
% entiere ; celle-ci est inseree a 0.7\textwidth, donc sa police doit etre grossie
% d'autant pour sortir a 12 pt dans le manuscrit. D'ou un jeu de dimensions local
% plutot que celui du haut de fichier.
carte_frac = 1.0;                                     % \includegraphics[width=\textwidth]
carte_w    = textwidth * carte_frac;                  % page figee, en points
carte_h    = 340;                                     % l'ancienne hauteur, mise a la nouvelle
                                                      % largeur : 280 * 425 / 348
carte_font = 12 * carte_w / (carte_frac * textwidth); % 12,0 : a pleine largeur la figure
                                                      % n'est pas mise a l'echelle

% Le champ est celui que lisent tous les essais du chapitre. TRANSPOSEE, comme
% campagne_chapitre_4.m : dans la sauvegarde, mesures indexe (x, y) ; grid_read veut des
% lignes en y et des colonnes en x, ce qu'attend aussi imagesc. Sans la transposition la
% figure montre le champ retourne sur sa diagonale, avec une trajectoire correcte
% par-dessus — invisible ici, la carte etant carree.
carte = load(fullfile(here, '..', '..', 'MATLAB 2A - sauvegarde', 'Cartes', ...
                      'carte_magnetometrie_anomalie_mexique.mat'));
h_map = double(carte.mesures).';
pas_m = double(carte.pas(:))' .* [1 1];               % m

% Le noeud (1,1) est en (pas, pas) et non en (0,0) : c'est l'origine que grid_read
% suppose pour les grilles des sauvegardes, et le releve est place dessus.
x_km = (1 : size(h_map, 2)) * pas_m(1) / 1e3;
y_km = (1 : size(h_map, 1)) * pas_m(2) / 1e3;

fprintf('Carte %d x %d noeuds a %.0f m, champ de %.0f a %.0f nT, ecart-type %.0f nT\n', ...
        size(h_map, 2), size(h_map, 1), pas_m(1), ...
        min(h_map(:)), max(h_map(:)), std(h_map(:)));

% LE RELEVE PASSE PAR krigeage.m ET N'EST PAS REDERIVE ICI, pour que les croix tombent
% exactement ou le filtre lit ses echantillons : centres de cellules, meme origine, et
% surtout meme secouage. Un reseau exact rendrait le biais d'interpolation periodique,
% donc repete d'une branche a l'autre du vol au lieu de se moyenner.
par_srv.unit           = par.unit;
par_srv.run.seed       = par.run.seed;
par_srv.krig.pas       = 7.5e3;      % m. Compromis de lisibilite, et non la maille de
                                     % reference : a 5 km les croix saturent la carte,
                                     % a 10 km elles la desertent
par_srv.krig.jitter    = par.krig.jitter;
par_srv.krig.sigma_map = par.krig.sigma_map;
par_srv.krig.n_ml      = 500;
par_srv.krig.ell_init  = 10e3;
par_srv.krig.noyau     = par.krig.noyau;

krig    = krigeage();
survey  = krig.releve(struct('h', h_map, 'step', pas_m), par_srv);
srv_km  = survey.X / 1e3;

% La trajectoire est relue de la campagne deja chargee plutot que recalculee, pour que la
% figure montre exactement le vol qui a produit les resultats.
vol = S.scen.x_true(1:2, :) / 1e3;                    % km

fig = figure('Color', 'w', 'Units', 'points', 'Position', [80 80 carte_w carte_h]);
imagesc(x_km, y_km, h_map); hold on;
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
c.Label.String        = '$h(x)$ (nT)';
c.Label.Interpreter   = 'latex';
c.Label.FontSize      = carte_font;
c.TickLabelInterpreter = 'latex';
c.FontSize            = carte_font;

% Les axes portent les composantes de position de l'etat du scenario, pas des indices.
xlabel('$p^{(x)}$ (km)', 'Interpreter', 'latex', 'FontSize', carte_font);
ylabel('$p^{(y)}$ (km)', 'Interpreter', 'latex', 'FontSize', carte_font);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', carte_font, 'Layer', 'top');
box on;

% La barre de couleur et son label vivent HORS des axes, et MATLAB dimensionne les axes
% avant de savoir ce que le label mesure : sur une page figee, il sort donc du PDF. Les
% marges de gauche et du bas sont MESUREES par TightInset et non devinees — une valeur en
% dur finit par rogner des que la police ou les chiffres d'axe changent. Seule la bande de
% droite reste imposee, TightInset ignorant la barre de couleur.
ax = gca;
drawnow;
ax.Units = 'normalized';
c.Units  = 'normalized';
ti     = ax.TightInset;
marg   = 0.015;
cb_w   = 0.035;                  % la barre elle-meme
% Ce qui vit A DROITE de la barre : les graduations, larges de cinq caracteres pour
% « -1000 », puis le label tourne, plus une marge. On l'exprime en POINTS et on divise
% par la largeur de page : en fraction, cette bande retrecit quand la page s'elargit ou
% que la police diminue, et une valeur en dur laisse une bande blanche des que l'un des
% deux change — ce qui est arrive au passage de 348 pt et 14 pt a 425 pt et 12 pt.
cb_pad = (3.4 * carte_font + 10) / carte_w;
gap    = 0.02;                   % entre les axes et la barre
left   = ti(1) + marg;
bot    = ti(2) + marg;
top    = ti(4) + marg;
cb_x   = 1 - cb_pad - cb_w;      % la barre est placee depuis le bord droit, pas depuis
                                 % les axes : c'est son label qui touche le bord de page
ax.Position = [left, bot, cb_x - gap - left, 1 - bot - top];
c.Position  = [cb_x, bot, cb_w, 1 - bot - top];

exporter(fig, 'Carte anomalie magnetique', carte_w, carte_h, img_dir);

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
xlabel(t, 'Time step $k$', 'Interpreter', 'latex', 'FontSize', fontsize);
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
xlabel(t, 'Time step $k$', 'Interpreter', 'latex', 'FontSize', fontsize);
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
xlabel('Map resolution (km)', 'Interpreter', 'latex', 'FontSize', fs);
ylabel('$\mathrm{ARMSE}_{\hat{z}}$ (nT)', 'Interpreter', 'latex', 'FontSize', fs);
legend({'Bilinear', 'Kriging'}, 'Interpreter', 'latex', 'FontSize', fs, 'Location', 'best');
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fs);
xlim([min(pas)-0.3 max(pas)+0.3]);

% Les marges par defaut de MATLAB laissent 13 % a gauche et 9,5 % a droite : le trace
% n'est pas centre dans la page, et il reste une bande blanche a droite qu'aucune
% etiquette ne remplit. On les mesure et on remplit la page.
ax = gca; drawnow; ti = ax.TightInset; marg = 0.02;
ax.Position = [ti(1) + marg, ti(2) + marg, ...
               1 - ti(1) - ti(3) - 2 * marg, 1 - ti(2) - ti(4) - 2 * marg];
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

% LA MEMOIRE, en octets. Il n'existe pas d'unite analogue au flop pour la taille : on
% compte des doubles, huit octets chacun, et les prefixes decimaux sont les memes.
%   krigeage    la matrice de Gram, 8 n^2, qui domine tout le reste
%   bilineaire  les echantillons seuls, trois doubles par point avec leurs coordonnees
%   carte vraie la grille entiere du champ
mem_carte = 8 * numel(map.h);
mem_bil   = 24 * n_tilde;
mem_krig  = 8 * n_tilde.^2 + mem_bil;

if ~exist(tab_dir, 'dir'), mkdir(tab_dir); end
fid = fopen(fullfile(tab_dir, 'cout_maille.tex'), 'w', 'n', 'UTF-8');
fprintf(fid, '%% Generated by figures_campagnes.m, do not edit by hand.\n');
fprintf(fid, ['%% Analytic flop counts per run, in Gflop, not timings, for %d particles over\n' ...
              '%% %d steps. Setup solves the training system once; queries cost n^2 per\n' ...
              '%% particle and per step. The true map and the bilinear model are read in eight\n' ...
              '%% flops per query and pay no setup, hence one row each: neither depends on the\n' ...
              '%% map resolution.\n\n'], par.mc.n_part, par.mc.n_steps);
fprintf(fid, '\\begin{tabular}{|l|cc|c|}\\hline\n');
% L'UNITE EST FACTORISEE DANS L'EN-TETE, LE PREFIXE RESTE EN CELLULE. Une colonne
% entierement en Gflop forcerait 29213 d'un cote et 0.01 de l'autre, donc cinq chiffres
% pour l'un et deux decimales perdues pour l'autre ; un prefixe par valeur garde trois
% chiffres significatifs partout, sans allonger les cellules.
fprintf(fid, ['Observation model & Setup (flops) & Queries (flops)' ...
              ' & Memory (bytes)\\\\\\hline\n']);
fprintf(fid, 'True map & -- & %s & %s\\\\\n', flops_tex(cout_bil), flops_tex(mem_carte));
fprintf(fid, 'Bilinear & -- & %s & %s\\\\\\hline\n', flops_tex(cout_bil), flops_tex(mem_bil(1)));
for ip = 1:numel(pas)
    fprintf(fid, 'Kriging, $%.1f$km & %s & %s & %s\\\\\n', ...
            pas(ip), flops_tex(cout_setup(ip)), flops_tex(cout_req(ip)), ...
            flops_tex(mem_krig(ip)));
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

% Les rapports bilineaire / krigeage, que l'analyse cite : c'est l'ecart entre les deux
% qui porte l'argument, et il n'est pas le meme sur la position et sur la consistance.
lab = {'RMSE', 'NEES'};   % un litteral cellule ne s'indexe pas directement en MATLAB
nom_metrique = {'RMSE', 'NEES'};
for j = 1:2
    rc = FIN(:, 3, j) ./ FIN(:, 1, j);
    fprintf('  rapport krig/carte %s : de %.1f a %.1f (', nom_metrique{j}, min(rc), max(rc));
    fprintf('%.1f ', rc); fprintf(')\n');
    r = FIN(:, 2, j) ./ FIN(:, 3, j);
    fprintf('  rapport bilin/krig %s : de %.1f a %.1f (', ...
            nom_metrique{j}, min(r), max(r));
    fprintf('%.1f ', r); fprintf(')\n');
end

sty = {':', '--', '-'};                           % carte connue, bilineaire, krigeage
duo = @(titres, Y, ylab, logy) deal_duo(pas, Y, titres, ylab, logy, sty, etiq, ...
                                        textwidth, fontsize, linewidth);

fig = duo({'Threshold', 'Mahalanobis'}, CONV, 'Converged runs (\%)', false);
exporter(fig, 'Convergence par maille', textwidth, 230, img_dir);

% NI RMSE NI NEES : les deux panneaux portent, par essai, la mediane des dix derniers pas,
% puis la mediane de ces valeurs sur tous les essais.
% moyennee ensuite sur les essais convergés. Ce sont exactement les deux quantites que les
% criteres de convergence seuillent, et aucune des deux ne passe par la racine des carres
% moyens sur les essais qui definit la RMSE du manuscrit. Les titres le disent donc en
% clair, et le label porte la double agregation, sur deux lignes plutot que dans une
% parenthese illisible.
fig = duo({'Position error (m)', 'Normalised error squared'}, FIN, ...
          {'Median over the last ten steps,', 'then over the $n_\mathrm{MC}$ runs'}, true);
exporter(fig, 'Erreur finale par maille', textwidth, 230, img_dir);


%% Les trois tableaux de la campagne B, sous include/
%
% LES CHIFFRES DE CAMPAGNE NE SE RECOPIENT PLUS A LA MAIN. Les tableaux tab:AK, tab:CAK et
% tab:cout_CAK etaient ecrits en dur dans 4_adaptive_kriging.tex, donc a retaper a chaque
% lot qui atterrit — ils l'ont deja ete deux fois, et douze lots de transcription sont une
% erreur qui attend. Ils sortent maintenant d'ici, comme cout_maille.tex.
%
% La campagne B est LA comparaison : les cinq modeles a maille et filtre fixes. On la relit
% separement de la campagne A, qui balaie la maille et ne porte ni AK ni CA-OK.
dB = dir(fullfile(here, 'resultats', '*_reference_B_*runs.mat'));
if isempty(dB)
    warning('figures_campagnes:campagneB', ...
            'Aucune campagne B sous resultats/ : les trois tableaux ne sont pas ecrits.');
else
[~, iB] = max([dB.datenum]);
SB   = load(fullfile(here, 'resultats', dB(iB).name));
parB = SB.par;  resB = SB.res;
nomB = cellfun(@char, {resB.model}, 'UniformOutput', false);
fprintf('Campagne B : %s\n', dB(iB).name);

% LE SEUIL EST CELUI DU HAUT DE CE FICHIER, comme pour les grilles par maille. Le chapitre
% l'annonce explicitement — « a fixed threshold of $2000$m, under which the median over the
% last ten steps of a filter must be » — et une figure calculee a 2000 m voisinant un
% tableau calcule a 3500 dirait deux choses differentes sous la meme definition.
%
% C'est un CHOIX DE DEPOUILLEMENT et non de simulation : navigation.m n'enregistre que
% l'erreur, la convergence s'en deduit ici, et le changer ne demande de rejouer aucun
% essai. La campagne porte le sien dans son fichier, qui sert de valeur de repli quand
% seuil_conv est vide.
seuilB = seuil_conv;
if isempty(seuilB), seuilB = parB.conv.seuil; end
n_itB = size(resB(1).err, 1);
finB  = n_itB - parB.conv.n_fin + 1 : n_itB;
k_acq = min(parB.conv.k_acq, n_itB);

M_ok  = metriques(resB, nomB, 'OK carte totale', finB, k_acq, seuilB);
M_ak  = metriques(resB, nomB, 'OK adaptatif',    finB, k_acq, seuilB);
M_cak = metriques(resB, nomB, 'CA-OK',           finB, k_acq, seuilB);

if ~exist(tab_dir, 'dir'), mkdir(tab_dir); end

% Les deux tableaux d'exactitude partagent leurs colonnes et ne different que par la paire
% de lignes comparee, d'ou les trois anonymes plutot que deux blocs recopies.
entete = @(fid) fprintf(fid, ['\\begin{tabular}{lrrrr}\n\t\\hline\n' ...
    '\t& converged & acquisition & end & $\\tilde{N}_\\mathrm{w}$ \\\\\n\t\\hline\n']);
ligne  = @(fid, nom, M) fprintf(fid, ...
    '\t%s & $%.0f\\%%$ & $%.0f$ m & $%.0f$ m & $%.0f$ \\\\\n', nom, M(1), M(2), M(3), M(4));
pied   = @(fid) fprintf(fid, '\t\\hline\n\\end{tabular}\n');
tete   = @(fid) fprintf(fid, ...
    ['%% Generated by figures_campagnes.m, do not edit by hand.\n' ...
     '%% Survey pitch %.2f km, %s filter, n_MC = %d over %d survey draws.\n' ...
     '%% Converged: median error over the last %d steps below %g m. Acquisition and end\n' ...
     '%% are median errors at steps %d and %d.\n\n'], ...
    parB.krig.pas / 1e3, resB(1).filter, parB.mc.n_runs, parB.krig.n_surveys(1), ...
    parB.conv.n_fin, seuilB, k_acq, n_itB);

fid = fopen(fullfile(tab_dir, 'ak_vs_ok.tex'), 'w', 'n', 'UTF-8');
tete(fid);  entete(fid);
ligne(fid, 'Ordinary kriging', M_ok);
ligne(fid, 'Adaptive kriging', M_ak);
pied(fid);  fclose(fid);

fid = fopen(fullfile(tab_dir, 'cak_vs_ak.tex'), 'w', 'n', 'UTF-8');
tete(fid);  entete(fid);
ligne(fid, 'Adaptive kriging',           M_ak);
ligne(fid, 'Clustered adaptive kriging', M_cak);
pied(fid);  fclose(fid);

% LE COUT, ANALYTIQUE ET PAR ESSAI. res.cout enregistre quatre lignes par pas : la somme
% des fenetres, la somme de leurs cubes, le nombre d'ajustements d'hyperparametres, et les
% evaluations de distance du mean-shift. On en tire trois postes :
%
%   factorisations  somme des n_c^3 / 3, enregistree telle quelle, PLUS l'ajustement local
%                   des hyperparametres, cubique lui aussi et plafonne a n_ml_win
%   requetes        N (somme n_c / n_modes)^2 : chaque mode ne resout QUE ses propres
%                   particules, donc le terme dominant est divise deux fois — une fois par
%                   la taille de fenetre, une fois par le partage des particules. C'est
%                   l'erreur que j'avais faite en moyennant le nombre de modes sur la
%                   mission au lieu de le prendre pas par pas.
%   clusterisation  cinq flops par distance en dimension deux, COMPTES A TOUS LES PAS, y
%                   compris ceux ou le mean-shift ne rend qu'un mode et ou la methode
%                   degenere en krigeage adaptatif.
n_opt = 30;                      % iterations de l'optimiseur du maximum de vraisemblance
C_ak  = cout_modele(resB, nomB, 'OK adaptatif', parB, n_opt);
C_cak = cout_modele(resB, nomB, 'CA-OK',        parB, n_opt);
G     = 1e9;

fid = fopen(fullfile(tab_dir, 'cout_cak.tex'), 'w', 'n', 'UTF-8');
fprintf(fid, ['%% Generated by figures_campagnes.m, do not edit by hand.\n' ...
    '%% Analytic flop counts per run, in Gflop, not timings, for %d particles over %d\n' ...
    '%% steps at a survey pitch of %.2f km. Factorisations include the local hyperparameter\n' ...
    '%% fit, cubic as well. Queries are N (sum n_c / n_modes)^2. Clustering is five flops\n' ...
    '%% per two-dimensional distance, counted at every step.\n\n'], ...
    parB.mc.n_part, parB.mc.n_steps, parB.krig.pas / 1e3);
fprintf(fid, ['\\begin{tabular}{lrrrr}\n\t\\hline\n' ...
    '\t& factorisations & queries & clustering & total \\\\\n\t\\hline\n']);
fprintf(fid, '\tAdaptive kriging & $%.2f$ & $%.2f$ & --- & $%.2f$ \\\\\n', ...
        C_ak(1) / G, C_ak(2) / G, sum(C_ak) / G);
fprintf(fid, '\tClustered adaptive kriging & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ \\\\\n', ...
        C_cak(1) / G, C_cak(2) / G, C_cak(3) / G, sum(C_cak) / G);
fprintf(fid, '\t\\hline\n\\end{tabular}\n');
fclose(fid);

% CONTROLE A LA CONSOLE. Le seuil alternatif dit si la conclusion « indiscernables » tient
% a la valeur choisie ; le rapport de couts est celui qu'annonce la prose du chapitre.
alt   = parB.conv.seuil;         % celui que portait la campagne, pour comparaison
A_ok  = metriques(resB, nomB, 'OK carte totale', finB, k_acq, alt);
A_ak  = metriques(resB, nomB, 'OK adaptatif',    finB, k_acq, alt);
A_cak = metriques(resB, nomB, 'CA-OK',           finB, k_acq, alt);
fprintf(['Trois tableaux ecrits. Convergence a %g m : OK %.0f %%, AK %.0f %%, CA-OK %.0f %% ' ...
         '| a %g m : %.0f, %.0f, %.0f\n'], seuilB, M_ok(1), M_ak(1), M_cak(1), ...
        alt, A_ok(1), A_ak(1), A_cak(1));
fprintf('Cout par essai : AK %.2f Gflop, CA-OK %.2f, rapport %.1f\n', ...
        sum(C_ak) / G, sum(C_cak) / G, sum(C_ak) / sum(C_cak));
end


%% Le tableau des parametres de fenetre, sous include/
%
% Les deux — bientot trois — parametres qui definissent la fenetre du krigeage adaptatif :
% le facteur de dilatation alpha_W, qui la fait suivre la covariance declaree, le plancher
% n_min, qui l'empeche de se vider quand la piste est acquise, et le plafond n_max, qui la
% borne. Le chapitre affirmait leur inertie sans la montrer ; elle se montre ici.
%
% CHAQUE LOT EST RECONNU A SON CHAMP BALAYE, jamais a son nom de fichier : les lots 'libre'
% produisent tous le meme '*_reference_libre_*runs.mat'. Un parametre dont le lot n'a pas
% encore tourne est signale et saute, donc le tableau grandit tout seul quand n_max
% atterrira, sans rien a retoucher ici.
fen = {'krig.window_factor', '\alpha_W'
       'krig.n_min',         'n_{\min}'
       'krig.n_max',         'n_{\max}'};

lots = cell(size(fen, 1), 1);
for k = 1:size(fen, 1)
    lots{k} = lot_balaye(here, fen{k, 1}, par.mc.n_runs);
    if isempty(lots{k})
        fprintf('  parametre de fenetre sans lot a %d essais : %s\n', par.mc.n_runs, fen{k, 1});
    end
end

if all(cellfun(@isempty, lots))
    warning('figures_campagnes:fenetrage', ...
            'Aucun balayage de fenetre sous resultats/ : fenetrage.tex n''est pas ecrit.');
else
if ~exist(tab_dir, 'dir'), mkdir(tab_dir); end
fid = fopen(fullfile(tab_dir, 'fenetrage.tex'), 'w', 'n', 'UTF-8');
fprintf(fid, ['%% Generated by figures_campagnes.m, do not edit by hand.\n' ...
    '%% Window parameters of the adaptive estimator, swept one at a time around the\n' ...
    '%% operating point, clustered adaptive kriging under the regularised filter.\n' ...
    '%% Converged: median error over the last ten steps below %g m.\n\n'], seuil_conv);
fprintf(fid, ['\\begin{tabular}{lrrrr}\n\t\\hline\n' ...
    '\t& $\\tilde{N}_\\mathrm{w}$ & converged & acquisition & end \\\\\n\t\\hline\n']);

for k = 1:size(fen, 1)
    S = lots{k};
    if isempty(S), continue, end
    n_i = size(S.res(1).err, 1);
    fi  = n_i - S.par.conv.n_fin + 1 : n_i;
    ka  = min(S.par.conv.k_acq, n_i);
    for e = 1:numel(S.res)
        M = metriques_e(S.res, e, fi, ka, seuil_conv);
        fprintf(fid, '\t$%s = %g$ & $%.0f$ & $%.0f\\%%$ & $%.0f$ m & $%.0f$ m \\\\\n', ...
                fen{k, 2}, S.res(e).balayage, M(4), M(1), M(2), M(3));
    end
    fprintf(fid, '\t\\hline\n');
end
fprintf(fid, '\\end{tabular}\n');
fclose(fid);

% CONTROLE : l'amplitude relative de chaque colonne sur chaque balayage. C'est elle que la
% prose du chapitre cite, et elle seule dit si « inerte » est encore le mot juste.
for k = 1:size(fen, 1)
    S = lots{k};
    if isempty(S), continue, end
    n_i = size(S.res(1).err, 1);
    fi  = n_i - S.par.conv.n_fin + 1 : n_i;
    ka  = min(S.par.conv.k_acq, n_i);
    T   = cell2mat(arrayfun(@(e) metriques_e(S.res, e, fi, ka, seuil_conv), ...
                            (1:numel(S.res))', 'UniformOutput', false));
    fprintf(['  %-20s n~ x%.1f  |  convergence %.0f a %.0f %%  |  acquisition %.0f a %.0f m ' ...
             '(%+.0f %%)  |  fin %.0f a %.0f m (%+.0f %%)\n'], fen{k, 1}, ...
            max(T(:, 4)) / min(T(:, 4)), min(T(:, 1)), max(T(:, 1)), ...
            min(T(:, 2)), max(T(:, 2)), 100 * (max(T(:, 2)) / min(T(:, 2)) - 1), ...
            min(T(:, 3)), max(T(:, 3)), 100 * (max(T(:, 3)) / min(T(:, 3)) - 1));
end
end

fprintf('figures_campagnes : six figures sous Images/ et cinq tableaux sous include/\n');
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
        if j == 1
            legend(ax, etiq, 'Interpreter', 'latex', 'FontSize', fontsize, ...
                   'Location', 'best');
        end
    end
    % Le label commun est porte par la disposition et non par le panneau de gauche, ce qui
    % le centre sur les deux et accepte un tableau de lignes tel quel.
    ylabel(t, ylab, 'Interpreter', 'latex', 'FontSize', fontsize);
    xlabel(t, 'Map resolution (km)', 'Interpreter', 'latex', 'FontSize', fontsize);
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

function M = metriques(res, noms, motif, fin, k_acq, seuil)
% Les quatre colonnes des tableaux d'exactitude, dans l'ordre ou elles y figurent :
% convergence en pourcentage, erreur mediane a l'acquisition, erreur mediane en fin de
% vol, et taille moyenne de la fenetre. Ce sont exactement les grandeurs que
% campagne_chapitre_4.m imprime en fin de campagne, calculees de la meme facon, pour que
% le tableau du manuscrit et le journal de la campagne ne puissent pas diverger.
    e = find(strcmp(noms, motif), 1);
    if isempty(e)
        error('figures_campagnes:modele', 'Modele absent de la campagne : %s.', motif);
    end
    M = metriques_e(res, e, fin, k_acq, seuil);
end

function C = cout_modele(res, noms, motif, par, n_opt)
% Les trois postes de cout analytique par essai, en flops : factorisations, requetes,
% clusterisation. Voir le commentaire du bloc appelant pour ce que chacun compte.
%
% LE NOMBRE DE MODES EST PRIS PAS PAR PAS. Le moyenner sur la mission fait apparaitre le
% CA-OK 23 % plus cher que l'AK au lieu de treize fois moins : le mean-shift rend
% cinquante-quatre modes au premier pas et un seul des le seizieme, donc une moyenne de
% mission ne decrit aucun pas reel.
    e = find(strcmp(noms, motif), 1);
    if isempty(e)
        error('figures_campagnes:modele', 'Modele absent de la campagne : %s.', motif);
    end
    c   = res(e).cout;
    sn  = c(1, :);  sn3 = c(2, :);  nf = c(3, :);  nd = c(4, :);
    nmk = max(mean(res(e).n_modes, 2)', 1);      % modes par pas, moyennes sur les essais
    nc  = sn ./ nmk;                             % taille moyenne d'une fenetre de mode
    C   = [sum(sn3) / 3 + sum(nf) * n_opt * min(mean(nc), par.krig.n_ml_win)^3 / 3, ...
           par.mc.n_part * sum(nc.^2), ...
           5 * sum(nd)];
end

function S = lot_balaye(here, champ, n_runs)
% Le lot le plus recent dont le balayage porte sur ce champ, AU NOMBRE D'ESSAIS DEMANDE,
% ou [] s'il n'a pas tourne.
%
% ON IDENTIFIE PAR LE CONTENU ET NON PAR LE NOM. Les lots 'libre' produisent tous le meme
% '*_reference_libre_*runs.mat' : seul par.balayage.champ dit ce qu'un fichier contient.
% On ne charge d'abord que par, qui pese quelques kilo-octets, et le reste seulement si
% c'est le bon lot — les fichiers font une centaine de mega-octets chacun.
%
% LE NOMBRE D'ESSAIS EST UNE CONDITION, PAS UNE PREFERENCE. Sans lui, le premier tirage de
% ce tableau est alle chercher un balayage de n_max a 100 essais, d'une campagne anterieure
% et d'un autre scenario : il annoncait 250 m d'acquisition la ou les autres blocs en
% donnent 475, dans le meme tableau et sans rien signaler. Un lot d'un autre nombre
% d'essais est desormais refuse et nomme.
    S = [];
    d = dir(fullfile(here, 'resultats', '*_reference_libre_*runs.mat'));
    if isempty(d), return, end
    [~, ordre] = sort([d.datenum], 'descend');
    autre = '';
    for k = ordre(:)'
        f = fullfile(here, 'resultats', d(k).name);
        p = load(f, 'par');
        if ~isfield(p.par, 'balayage') || ~strcmp(p.par.balayage.champ, champ), continue, end
        if p.par.mc.n_runs ~= n_runs
            if isempty(autre)
                autre = sprintf('%s (%d essais)', d(k).name, p.par.mc.n_runs);
            end
            continue
        end
        fprintf('  %-20s <- %s\n', champ, d(k).name);
        S = load(f);
        return
    end
    if ~isempty(autre)
        fprintf('  %-20s ECARTE : %s, il en faut %d\n', champ, autre, n_runs);
    end
end

function M = metriques_e(res, e, fin, k_acq, seuil)
% Les quatre grandeurs d'une ligne de campagne, pour l'entree e : convergence en
% pourcentage, erreur mediane a l'acquisition, erreur mediane en fin de vol, taille
% moyenne de la fenetre. Memes definitions que celles qu'imprime campagne_chapitre_4.m.
    E = res(e).err;
    M = [100 * mean(median(E(fin, :), 1) <= seuil), ...
         median(E(k_acq, :)), ...
         median(E(end, :)), ...
         mean(res(e).n_used, 'all')];
end
