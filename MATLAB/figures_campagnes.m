% FIGURES_CAMPAGNES
%
% Les figures du chapitre 4 tirees des campagnes enregistrees. Un seul script, pour que
% les courbes du manuscrit aient une source unique et suivent d'elles-memes le passage a
% 1000 essais : il relit la campagne la plus recente, sans parametre a retoucher.
%
%   Carte anomalie magnetique  le champ, le releve secoue et la trajectoire de reference
%   ARMSE de reconstruction    l'erreur de carte contre la maille, et sa loi de puissance
%   Erreur finale par maille   RMSE et NEES medianes sur les dix derniers pas
%   Echantillons retenus par maille  la fenetre de l'AK a chaque pas, une courbe par maille
%   Modes au cours du vol      ce que trouve le mean-shift, pas par pas, campagne B
%   Fenetre au cours du vol    echantillons retenus par requete et plancher n_min, campagne B
%
% puis cinq comparaisons par paires, chacune avec la carte connue en borne, en trois
% figures : 'RMSE <etape>', grille 3 x 2 d'une maille par case ; 'NEES <etape>', la meme en
% echelle log avec la bande du chi2 ; et la convergence contre la maille, deux criteres en
% deux panneaux. Les etapes, dans l'ordre du chapitre :
%
%   par maille              bilineaire contre statique, sous les noms historiques
%   eclairci par maille     statique contre releve eclairci
%   adaptatif par maille    statique contre AK
%   reunion par maille      AK contre reunion des voisins
%   clusterise par maille   AK contre CA-OK
%
% COULEUR, TRAIT ET NOM DE CHAQUE METHODE sortent d'une seule fonction, style_modele, en fin
% de fichier : figures, legendes et en-tetes de tableaux la consultent tous, et renommer une
% methode ne touche qu'une ligne.
%
% et neuf tableaux sous Latex/Manuscrit/include/ :
%
%   cout_maille.tex   le cout analytique par maille          campagne A
%   cout_<etape>.tex  le cout de chaque comparaison, sauf la premiere
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
% relit. figures_temporelles.m l'a ete a son tour le 15 septembre, avec les deux figures au
% cours du vol : le chapitre n'a plus qu'un script de figures, donc qu'une palette.
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
etiq    = {'connue', 'bilin', 'statique'};    % les cles de style_modele, dans le meme ordre

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
% Vient de figures_carte.m, supprime : les figures du chapitre sont produites par ce seul
% script, et une figure isolee dans son propre fichier finit par diverger des autres sur la
% page, la police ou la campagne qu'elle relit.
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

%% Les bornes communes aux grilles
%
% La bande du chi2 ou une NEES honnete doit tomber, et le sous-echantillonnage des pas
% traces : a 150 pas la NEES est trop bruitee pour qu'on distingue les modeles.
bas  = chi2inv((1 - par.conv.confiance) / 2, d_x) / d_x;
haut = chi2inv(1 - (1 - par.conv.confiance) / 2, d_x) / d_x;
kd   = 1 : pas_trace : n_it;

% Le premier virage, trace sur les grilles de RMSE : c'est lui qui leve l'ambiguite le long
% de la premiere branche, et les retards d'acquisition entre methodes se lisent contre lui.
k_vir = premier_virage(S.scen.x_true);

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
% La palette commune et non l'ordre par defaut de MATLAB, qui mettait le bilineaire en bleu
% et le krigeage en orange, a rebours de toutes les autres figures.
[cb, lb, eb] = style_modele('bilin');
[ck, lk, ek] = style_modele('statique');
plot(pas, ARMSE(:, 1), lb, 'Color', cb, 'Marker', 'o', 'LineWidth', linewidth);
plot(pas, ARMSE(:, 2), lk, 'Color', ck, 'Marker', 'o', 'LineWidth', linewidth);
hold off;
fs = 12 * 348 / (0.7 * textwidth);
xlabel('Map resolution (km)', 'Interpreter', 'latex', 'FontSize', fs);
ylabel('$\mathrm{ARMSE}_{\hat{z}}$ (nT)', 'Interpreter', 'latex', 'FontSize', fs);
legend({eb, ek}, 'Interpreter', 'latex', 'FontSize', fs, 'Location', 'southeast');
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
n_opt   = 30;                    % iterations de l'optimiseur du maximum de vraisemblance

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
fprintf(fid, ['%% Analytic flop counts per run, not timings, for %d particles over %d steps.\n' ...
              '%% Setup solves the training system once; queries cost n^2 per particle and per\n' ...
              '%% step; total is their sum, the figure every comparison table of the chapter\n' ...
              '%% reports. The known map and bilinear interpolation are read in eight flops per\n' ...
              '%% query and pay no setup, hence one row each: neither depends on the map\n' ...
              '%% resolution.\n\n'], par.mc.n_part, par.mc.n_steps);
fprintf(fid, '\\begin{tabular}{|l|ccc|c|}\\hline\n');
% L'UNITE EST FACTORISEE DANS L'EN-TETE, LE PREFIXE RESTE EN CELLULE. Une colonne
% entierement en Gflop forcerait 29213 d'un cote et 0.01 de l'autre, donc cinq chiffres
% pour l'un et deux decimales perdues pour l'autre ; un prefixe par valeur garde trois
% chiffres significatifs partout, sans allonger les cellules.
[~, ~, e_connue] = style_modele('connue');
[~, ~, e_bilin]  = style_modele('bilin');
[~, ~, e_stat]   = style_modele('statique');
% LE TOTAL EST ECRIT, et pas seulement ses deux termes : c'est lui que reportent les tableaux
% de comparaison du chapitre, et un lecteur qui lit 29.2T ici et 29.3T la-bas croit a une
% erreur, alors que l'ecart n'est que la factorisation. L'unite passe sur une ligne d'en-tete
% commune, sans quoi la cinquieme colonne ne tient plus dans \textwidth.
fprintf(fid, ' & \\multicolumn{3}{c|}{Cost (flops)} & Memory\\\\\n');
fprintf(fid, 'Observation model & Setup & Queries & Total & (bytes)\\\\\\hline\n');
fprintf(fid, '%s & -- & %s & %s & %s\\\\\n', e_connue, flops_tex(cout_bil), ...
        flops_tex(cout_bil), flops_tex(mem_carte));
fprintf(fid, '%s & -- & %s & %s & %s\\\\\\hline\n', e_bilin, flops_tex(cout_bil), ...
        flops_tex(cout_bil), flops_tex(mem_bil(1)));
for ip = 1:numel(pas)
    fprintf(fid, '%s, $%.1f$km & %s & %s & %s & %s\\\\\n', e_stat, ...
            pas(ip), flops_tex(cout_setup(ip)), flops_tex(cout_req(ip)), ...
            flops_tex(cout_setup(ip) + cout_req(ip)), flops_tex(mem_krig(ip)));
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

duo = @(titres, Y, ylab, logy) deal_duo(pas, Y, titres, ylab, logy, etiq, ...
                                        textwidth, fontsize, linewidth);

% La convergence contre la maille est desormais produite par le bloc des comparaisons,
% pour que les cinq figures de convergence du chapitre sortent du meme code.

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


%% Les comparaisons par paires, sur les six mailles
%
% LE CHAPITRE PROGRESSE PAR ETAPES, et chaque etape a sa figure. Le bilineaire contre le
% krigeage statique, puis le statique contre l'adaptatif, puis l'adaptatif contre le
% clusterise. Chacune porte la carte connue en pointille, comme borne que personne ne
% depasse, et ne montre qu'UN pas de la progression : trois courbes par case, jamais sept.
%
% Les controles negatifs se lisent de la meme facon, contre l'etape qu'ils mettent en
% cause — le releve eclairci contre le statique, puisqu'il pretend en resoudre le cout.
%
% TOUTES LES COURBES SONT ASSEMBLEES UNE FOIS ICI. La campagne A porte la carte connue, le
% bilineaire et le statique ; le lot des quatre modeles porte l'AK, le CA-OK et les deux
% controles. Les deux ont tourne sur les memes cent tirages — les releves sont des
% fonctions deterministes de par.run.seed + 1000*s — donc les courbes se comparent d'un
% lot a l'autre sans reserve.
%
% ON IDENTIFIE LE SECOND LOT PAR SON CONTENU. Le balayage multimodal porte lui aussi
% krig.pas sur un lot 'libre' a 1000 essais ; seul sigma_0 les separe. Sans ce test, les
% figures iraient chercher le mauvais des deux des que le multimodal aura fini.
dE = dir(fullfile(here, 'resultats', '*_reference_libre_*runs.mat'));
[~, ordE] = sort([dE.datenum], 'descend');
SE = [];
for k = ordE(:)'
    f = fullfile(here, 'resultats', dE(k).name);
    p = load(f, 'par');
    if ~isfield(p.par, 'balayage') || ~strcmp(p.par.balayage.champ, 'krig.pas'), continue, end
    if p.par.mc.n_runs ~= par.mc.n_runs || p.par.sigma_0(1) > 10e3, continue, end
    fprintf('Lot des quatre modeles : %s\n', dE(k).name);
    SE = load(f);
    break
end

if isempty(SE)
    warning('figures_campagnes:comparaisons', ...
            'Pas de balayage de maille au scenario de reference : comparaisons sautees.');
else
% Le vivier : les deux lots mis bout a bout, identifies par motif exclusif. On n'apparie
% jamais sur 'OK' ni sur 'carte' seuls, quatre modeles contenant le premier et deux le
% second.
%
% LE RPF SEUL DANS LES DEUX LOTS, comme plus haut pour la campagne A : l'appariement prend la
% premiere entree qui convient, et un lot portant plusieurs filtres lui en offrirait un
% autre. Et les deux lots n'ont pas tourne avec la meme version de navigation.m : chacun
% recoit, vides, les champs qui ne sont que dans l'autre, sans quoi la concatenation echoue.
rA = reshape(res, 1, []);
rE = reshape(SE.res(strcmp({SE.res.filter}, 'RPF')), 1, []);
for ch = setdiff(fieldnames(rE), fieldnames(rA))', [rA.(ch{1})] = deal([]); end
for ch = setdiff(fieldnames(rA), fieldnames(rE))', [rE.(ch{1})] = deal([]); end
pool = [orderfields(rA), orderfields(rE)];
nomP = cellfun(@char, {pool.model}, 'UniformOutput', false);
balP = [pool.balayage] / 1e3;

% Premiere colonne, la cle de style_modele, qui donne couleur, trait et nom ; seconde, le
% motif qui retrouve le modele dans le vivier.
cand = {'connue',   'connue'
        'bilin',    'bilin'
        'statique', 'totale'
        'ak',       'adaptatif'
        'cak',      'CA-OK'
        'eclairci', 'clairci'
        'reunion',  'union'};
n_mod = size(cand, 1);

CRB  = nan(n_it, numel(pas), n_mod);       % RMSE par pas
CNE  = nan(n_it, numel(pas), n_mod);       % NEES par pas
CCV  = nan(numel(pas), n_mod, 2);          % convergence, seuil et Mahalanobis
CNU  = nan(numel(pas), n_mod);             % echantillons retenus
CGF  = nan(numel(pas), n_mod);             % cout par essai, en flops
trouve = false(numel(pas), n_mod);         % l'entree existe-t-elle dans le vivier

for im = 1:n_mod
    for ip = 1:numel(pas)
        e = find(contains(nomP, cand{im, 2}) & abs(balP - pas(ip)) < 1e-9, 1);
        if isempty(e), continue, end
        trouve(ip, im) = true;
        E = pool(e).err;
        D = pool(e).d2;  D(~isfinite(D)) = NaN;
        crit_e = median(E(fin, :), 1);
        crit_d = median(D(fin, :), 1);
        CCV(ip, im, 1) = 100 * mean(crit_e <= par.conv.seuil);
        CCV(ip, im, 2) = 100 * mean(crit_d <= maha, 'omitnan');
        ok = crit_e <= par.conv.seuil;     % memes essais converges pour toutes les courbes
        CRB(:, ip, im) = sqrt(mean(E(:, ok).^2, 2));
        CNE(:, ip, im) = median(D(:, ok), 2, 'omitnan') / d_x;
        CNU(ip, im)    = mean(pool(e).n_used, 'all');

        % LE COUT : mesure quand le modele l'enregistre, analytique sinon. L'AK, le CA-OK
        % et la reunion comptent leurs fenetres pas par pas dans res.cout, et cout_e les
        % chiffre comme tab:cout_CAK, ajustements locaux et mean-shift compris. La carte
        % connue, le bilineaire, le statique et l'eclairci n'ont rien a compter, leur
        % fenetre etant fixe. Pour ces derniers on applique la formule du tableau par
        % maille — une factorisation en n^3/3, une requete en n^2 par particule et par pas.
        if any(pool(e).cout(:))
            CGF(ip, im) = sum(cout_e(pool(e), SE.par, n_opt));
        elseif CNU(ip, im) > 0 && ~strcmp(cand{im, 1}, 'bilin')
            % Le bilineaire enregistre la taille du releve dans n_used mais ne resout rien :
            % il lit en huit flops, comme dans cout_maille.tex, et non en n^2.
            CGF(ip, im) = CNU(ip, im)^3 / 3 + n_query * CNU(ip, im)^2;
        else
            CGF(ip, im) = n_query * 8;     % huit flops par lecture, sans entrainement
        end
    end
end

% CE QUI A ETE TROUVE, maille par maille. Un modele absent d'une maille laisse une case sans
% sa courbe, et une case incomplete se remarque moins qu'une ligne a la console.
for im = 1:n_mod
    fprintf('  %-18s trouve a %d mailles sur %d\n', cand{im, 1}, sum(trouve(:, im)), numel(pas));
end

% La carte connue ne depend pas de la maille : si elle n'a pas tourne a chacune, la courbe
% d'une maille est reprise dans les cases vides, comme le faisaient les grilles d'origine.
ic = find(strcmp(cand(:, 1), 'connue'));
iv = find(trouve(:, ic), 1);
if ~isempty(iv)
    for ip = find(~trouve(:, ic))'
        CRB(:, ip, ic) = CRB(:, iv, ic);   CNE(:, ip, ic) = CNE(:, iv, ic);
        CCV(ip, ic, :) = CCV(iv, ic, :);   CNU(ip, ic)    = CNU(iv, ic);
        CGF(ip, ic)    = CGF(iv, ic);
    end
end

% LES CINQ COMPARAISONS DU CHAPITRE. Chaque ligne donne les trois modeles a superposer, le
% nom des deux grilles, celui de la figure de convergence et celui du tableau de cout. Les
% legendes viennent de style_modele. La premiere sort sous les noms de fichier historiques,
% que le manuscrit cite deja.
%
% L'ordre suit la progression : le bilineaire d'abord, que le krigeage statique remplace ;
% puis les deux controles negatifs du sous-echantillonnage et du fenetrage, chacun contre
% l'estimateur qu'il pretend egaler ; puis les deux contributions, fenetrage et clusters.
etapes = {[1 3 2], 'par maille',            'Convergence par maille', ''
          [1 3 6], 'eclairci par maille',   'Convergence eclairci',   'cout_eclairci'
          [1 3 4], 'adaptatif par maille',  'Convergence adaptatif',  'cout_adaptatif'
          [1 4 7], 'reunion par maille',    'Convergence reunion',    'cout_reunion'
          [1 4 5], 'clusterise par maille', 'Convergence clusterise', 'cout_clusterise'};
for s = 1:size(etapes, 1)
    sel  = etapes{s, 1};
    cles = cand(sel, 1)';
    leg  = cell(size(cles));
    for im = 1:numel(cles), [~, ~, leg{im}] = style_modele(cles{im}); end

    grille(kk, kk, CRB(:, :, sel), pas, cles, '$\mathrm{RMSE}_{\hat{x}}(k)$ (m)', false, [], k_vir, ...
           pdf_w, pdf_h, fontsize, linewidth, img_dir, ['RMSE ' etapes{s, 2}]);
    grille(kk, kd, CNE(:, :, sel), pas, cles, '$\mathrm{NEES}_{\hat{x}}(k)$', true, [bas haut], NaN, ...
           pdf_w, pdf_h, fontsize, linewidth, img_dir, ['NEES ' etapes{s, 2}]);

    f2 = figure('Color', 'w', 'Units', 'points', 'Position', [40 40 textwidth 230]);
    t2 = tiledlayout(f2, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    titres = {'Threshold', 'Mahalanobis'};
    for j = 1:2
        ax = nexttile(t2); hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on');
        for im = 1:numel(sel)
            [c, l] = style_modele(cles{im});
            plot(ax, pas, CCV(:, sel(im), j), l, 'Color', c, ...
                 'LineWidth', linewidth, 'Marker', 'o');
        end
        hold(ax, 'off');
        ylim(ax, [0 105]);
        xlim(ax, [min(pas) - 0.3, max(pas) + 0.3]);
        title(ax, titres{j}, 'Interpreter', 'latex', 'FontSize', fontsize);
        set(ax, 'TickLabelInterpreter', 'latex', 'FontSize', fontsize);
        if j == 1
            lg = legend(ax, leg, 'Interpreter', 'latex', 'FontSize', fontsize, ...
                        'Orientation', 'horizontal');
            lg.ItemTokenSize = [18 18];      % propriete cachee, que legend() refuse en argument
            lg.Layout.Tile   = 'north';
        end
    end
    ylabel(t2, 'Converged runs (\%)', 'Interpreter', 'latex', 'FontSize', fontsize);
    xlabel(t2, 'Map resolution (km)', 'Interpreter', 'latex', 'FontSize', fontsize);
    exporter(f2, etapes{s, 3}, textwidth, 230, img_dir);

    % Le tableau de cout de la premiere etape existe deja sous cout_maille.tex, avec ses
    % colonnes de setup et de requetes : on ne le refait pas.
    if isempty(etapes{s, 4}), continue, end
    fid = fopen(fullfile(tab_dir, [etapes{s, 4} '.tex']), 'w', 'n', 'UTF-8');
    fprintf(fid, ['%% Generated by figures_campagnes.m, do not edit by hand.\n' ...
        '%% %s against %s.\n' ...
        '%% Analytic flop counts per run, not timings, for %d particles over %d steps.\n' ...
        '%% Where the model records its windows, they are counted as in cout_cak.tex:\n' ...
        '%% factorisations with the local fits, queries, clustering. Where the training set\n' ...
        '%% is fixed, n^3/3 for the setup plus N n^2 per query. The known map reads in eight\n' ...
        '%% flops per query.\n\n'], ...
        leg{2}, lower(leg{3}), par.mc.n_part, par.mc.n_steps);
    fprintf(fid, '\\begin{tabular}{|c|ccc|}\\hline\n');
    % Le meme nombre que la colonne Total de cout_maille.tex, et l'en-tete le dit.
    fprintf(fid, ' & \\multicolumn{%d}{c|}{Total cost (flops)}\\\\\n', numel(sel));
    fprintf(fid, 'Map resolution');
    for im = 1:numel(sel), fprintf(fid, ' & %s', leg{im}); end
    fprintf(fid, '\\\\\\hline\n');
    for ip = 1:numel(pas)
        fprintf(fid, '$%.1f$km', pas(ip));
        for im = 1:numel(sel)
            fprintf(fid, ' & %s', flops_tex(CGF(ip, sel(im))));
        end
        fprintf(fid, '\\\\\n');
    end
    fprintf(fid, '\\hline\n\\end{tabular}\n');
    fclose(fid);
end

% LA FENETRE DE L'AK AU COURS DU VOL, UNE COURBE PAR MAILLE, sur un seul axe : c'est le
% nombre d'echantillons retenus que le texte introduit pour expliquer les ecarts de cout
% d'une maille a l'autre. Superposees et non en grille, parce que les six courbes ne se
% separent qu'a l'acquisition et se rejoignent sur le plancher n_min des que la piste est
% acquise : c'est cet ecart-la que la figure doit montrer.
%
% UNE SEULE METHODE, DONC UNE SEULE TEINTE : l'orange de l'AK, du plus sombre pour la maille
% la plus fine au plus clair pour la plus lache. La couleur reste celle de la methode, et la
% maille se lit dans sa valeur.
fen_w  = 348;  fen_h = 280;                  % insere a 0.7\textwidth, comme l'ARMSE
fen_fs = 12 * fen_w / (0.7 * textwidth);
teinte = interp1([0 0.5 1], [0.45 0.18 0.00; 0.93 0.50 0.05; 1.00 0.76 0.38], ...
                 linspace(0, 1, numel(pas)));
[~, ~, e_ak] = style_modele('ak');
fig = figure('Color', 'w', 'Units', 'points', 'Position', [80 80 fen_w fen_h]);
hold on; box on; grid on;
fen_leg = {};
fen_max = 0;
for ip = 1:numel(pas)
    e = find(contains(nomP, 'adaptatif') & abs(balP - pas(ip)) < 1e-9, 1);
    if isempty(e), continue, end
    plot(kk, pool(e).cout(1, :), '-', 'Color', teinte(ip, :), 'LineWidth', linewidth);
    fen_leg{end + 1} = sprintf('$%.1f$ km', pas(ip));   %#ok<SAGROW>
    fen_max = max(fen_max, max(pool(e).cout(1, :)));
end
% Le plancher en pointille et le premier virage en tirete, comme sur la fenetre du
% chapitre au cours du vol.
yline(SE.par.krig.n_min, ':', 'LineWidth', linewidth, 'HandleVisibility', 'off');
marqueur(premier_virage(S.scen.x_true), fen_fs);
hold off; set(gca, 'YScale', 'log');
% De la marge au-dessus du pic de la maille la plus fine : sans elle, MATLAB cale l'axe sur
% lui et la courbe touche le cadre.
ylim([10^floor(log10(SE.par.krig.n_min)), 2 * fen_max]);
xlabel('$k$', 'Interpreter', 'latex', 'FontSize', fen_fs);
ylabel('$\tilde{N}^r(k)$', 'Interpreter', 'latex', 'FontSize', fen_fs);
lg = legend(fen_leg, 'Interpreter', 'latex', 'FontSize', fen_fs, 'Location', 'northeast');
lg.Title.String      = e_ak;
lg.Title.Interpreter = 'latex';
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fen_fs);
xlim([1 n_it]);
exporter(fig, 'Echantillons retenus par maille', fen_w, fen_h, img_dir);

fprintf('%d comparaisons ecrites : %s\n', size(etapes, 1), strjoin(etapes(:, 2)', ', '));
for im = 1:n_mod
    if all(isnan(CGF(:, im))), continue, end
    fprintf('  %-18s n~ %5.0f a %5.0f, cout %6.2f a %8.2f Gflop\n', cand{im, 1}, ...
            min(CNU(:, im)), max(CNU(:, im)), min(CGF(:, im))/1e9, max(CGF(:, im))/1e9);
end
end

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
    '\t& converged & acquisition & end & $\\tilde{N}^r$ \\\\\n\t\\hline\n']);
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
[~, ~, e_stat] = style_modele('statique');
[~, ~, e_ak]   = style_modele('ak');
[~, ~, e_cak]  = style_modele('cak');
ligne(fid, e_stat, M_ok);
ligne(fid, e_ak,   M_ak);
pied(fid);  fclose(fid);

fid = fopen(fullfile(tab_dir, 'cak_vs_ak.tex'), 'w', 'n', 'UTF-8');
tete(fid);  entete(fid);
ligne(fid, e_ak,  M_ak);
ligne(fid, e_cak, M_cak);
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
fprintf(fid, '\t%s & $%.2f$ & $%.2f$ & --- & $%.2f$ \\\\\n', e_ak, ...
        C_ak(1) / G, C_ak(2) / G, sum(C_ak) / G);
fprintf(fid, '\t%s & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ \\\\\n', e_cak, ...
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


%% Les deux figures au cours du vol, sous Images/
%
% Venues de figures_temporelles.m, absorbe le 15 septembre : toutes les figures et tous les
% tableaux du chapitre sortent de ce seul script, et la palette avec eux.
%
%   Modes au cours du vol     ce que trouve le mean-shift, pas par pas
%   Fenetre au cours du vol   echantillons retenus par requete, avec le plancher n_min
%
% Toutes deux relisent la campagne B et sont inserees a 0.7\textwidth, sur la page figee de
% 348 x 280 pt : une police de 14 y sort a 12 pt, comme pour l'ARMSE.
if ~isempty(dB)
vol_w  = 348;  vol_h = 280;
vol_fs = 12 * vol_w / (0.7 * textwidth);
k_vol  = (1:n_itB)';

% Le pas ou commence le premier demi-tour, releve sur la trajectoire vraie plutot que
% recalcule : c'est lui qui leve l'ambiguite, en echantillonnant le champ dans une seconde
% direction, et les courbes ne se lisent pas sans lui.
k_vir = premier_virage(SB.scen.x_true);
fprintf('Premier virage au pas %d sur %d\n', k_vir, n_itB);

j      = find(strcmp(nomB, 'CA-OK'), 1);
[c, l] = style_modele('cak');
fig = figure('Color', 'w', 'Units', 'points', 'Position', [80 80 vol_w vol_h]);
plot(k_vol, mean(resB(j).n_modes, 2), l, 'Color', c, 'LineWidth', linewidth);
box on; grid on; hold on; marqueur(k_vir, vol_fs); hold off;
set(gca, 'YScale', 'log');
xlabel('$k$', 'Interpreter', 'latex', 'FontSize', vol_fs);
ylabel('Number of modes', 'Interpreter', 'latex', 'FontSize', vol_fs);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', vol_fs);
xlim([1 n_itB]);
exporter(fig, 'Modes au cours du vol', vol_w, vol_h, img_dir);

fig = figure('Color', 'w', 'Units', 'points', 'Position', [80 80 vol_w vol_h]);
hold on; box on; grid on;
vol_cle   = {'ak', 'cak'};
vol_motif = {'OK adaptatif', 'CA-OK'};
vol_leg   = cell(1, 2);
for m = 1:2
    j = find(strcmp(nomB, vol_motif{m}), 1);
    [c, l, vol_leg{m}] = style_modele(vol_cle{m});
    plot(k_vol, resB(j).cout(1, :), l, 'Color', c, 'LineWidth', linewidth);
end
% Le plancher en pointille : c'est ainsi que la legende du manuscrit le designe.
yline(parB.krig.n_min, ':', 'LineWidth', linewidth, 'HandleVisibility', 'off');
marqueur(k_vir, vol_fs);
hold off; set(gca, 'YScale', 'log');
xlabel('$k$', 'Interpreter', 'latex', 'FontSize', vol_fs);
ylabel('$\tilde{N}^r(k)$', 'Interpreter', 'latex', 'FontSize', vol_fs);
legend(vol_leg, 'Interpreter', 'latex', 'FontSize', vol_fs, 'Location', 'northeast');
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', vol_fs);
xlim([1 n_itB]);
exporter(fig, 'Fenetre au cours du vol', vol_w, vol_h, img_dir);
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
fen = {'krig.window_factor', '\alpha_{\mathrm{cov}}'
       'krig.n_min',         '\tilde{N}_{\min}'
       'krig.n_max',         '\tilde{N}_{\max}'};

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
    '\t& $\\tilde{N}^r$ & converged & acquisition & end \\\\\n\t\\hline\n']);

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

fprintf('figures_campagnes : vingt et une figures sous Images/ et neuf tableaux sous include/\n');
%% ------------------------------------------------------------ fonctions locales

function fig = deal_duo(pas, Y, titres, ylab, logy, etiq, textwidth, fontsize, linewidth)
% Deux panneaux cote a cote sur une meme abscisse de maille, trois modeles par panneau.
% Y est (maille, modele, panneau). ylab est commun aux deux panneaux, ou l'un par panneau
% quand les deux grandeurs n'ont pas la meme unite. etiq porte les cles de style_modele.
    fig = figure('Color', 'w', 'Units', 'points', 'Position', [40 40 textwidth 230]);
    t = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    leg = cell(1, size(Y, 2));
    for j = 1:2
        ax = nexttile(t); hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on');
        for im = 1:size(Y, 2)
            [c, l, leg{im}] = style_modele(etiq{im});
            plot(ax, pas, Y(:, im, j), l, 'Color', c, 'LineWidth', linewidth, 'Marker', 'o');
        end
        hold(ax, 'off');
        if logy
            % Une marge d'un facteur 1,5 de part et d'autre : sans elle, MATLAB cale l'axe
            % sur les valeurs extremes et coupe en deux les marqueurs de la carte connue.
            y   = Y(:, :, j);  y = y(isfinite(y) & y > 0);
            lim = [min(y) / 1.5, 1.5 * max(y)];
            % Et une graduation par decade, que MATLAB espace d'office des qu'on fixe l'axe.
            set(ax, 'YScale', 'log', 'YLim', lim, ...
                    'YTick', 10.^(ceil(log10(lim(1))) : floor(log10(lim(2)))));
        else
            ylim(ax, [0 105]);
        end
        title(ax, titres{j}, 'Interpreter', 'latex', 'FontSize', fontsize);
        set(ax, 'TickLabelInterpreter', 'latex', 'FontSize', fontsize);
        xlim(ax, [min(pas) - 0.3, max(pas) + 0.3]);
        % Les graduations restent sur LES DEUX panneaux : ils ne partagent pas d'echelle,
        % mieux, pas meme d'unite — des metres a gauche, un rapport a droite. Les masquer
        % a droite ferait lire les valeurs de gauche.
        if j == 1
            lg = legend(ax, leg, 'Interpreter', 'latex', 'FontSize', fontsize, ...
                        'Orientation', 'horizontal');
            lg.ItemTokenSize = [18 18];      % propriete cachee, que legend() refuse en argument
            lg.Layout.Tile   = 'north';
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
    C = cout_e(res(e), par, n_opt);
end

function C = cout_e(r, par, n_opt)
% Les trois postes pour une entree deja trouvee. A part de cout_modele pour que les
% tableaux de cout par maille chiffrent exactement comme tab:cout_CAK.
    c   = r.cout;
    sn  = c(1, :);  sn3 = c(2, :);  nf = c(3, :);  nd = c(4, :);
    nmk = mean(r.n_modes, 2)';                   % modes par pas, moyennes sur les essais
    if isempty(nmk), nmk = ones(size(sn)); end   % un modele qui n'enregistre pas ses modes
    nmk = max(nmk, 1);
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

function fig = grille(kk, kd, Y, pas, cles, ylab, logy, bande, k_vir, pdf_w, pdf_h, fontsize, linewidth, img_dir, nom_fig)
% Une grille 3 x 2, une maille par case, les modeles superposes dans chacune. Y est
% (pas de temps, maille, modele) et cles donne la cle style_modele de chaque modele, d'ou
% se deduisent couleur, trait et nom. kd donne les pas effectivement traces : a 150 pas la
% NEES est trop bruitee pour qu'on distingue les modeles, la RMSE ne l'est pas. k_vir est le
% pas du premier virage, marque d'un tirete dans chaque case, ou NaN pour s'en passer.
%
% Les six mailles tiennent en une grille plutot qu'en six flottants : la degradation se lit
% d'une case a l'autre, ce qu'une figure isolee par maille ne montre pas.
    n_it = size(Y, 1);
    % LES BORNES SUIVENT CE QUI EST TRACE, et non tous les pas : la NEES du bilineaire a un
    % pic hors des pas retenus par kd, et l'axe montait a 1e6 pour des courbes qui ne
    % depassent pas 1e3. La bande du chi2 y entre aussi, pour que ses deux bords se voient.
    Yk   = Y(kd, :, :);
    y_ok = Yk(isfinite(Yk) & Yk > 0);
    if logy
        lim = [0.5 * min(y_ok), 2 * max(y_ok)];
        if ~isempty(bande)
            lim = [min(lim(1), 0.7 * bande(1)), max(lim(2), 1.4 * bande(2))];
        end
    else
        lim = [0, 1.02 * max(y_ok)];
    end
    fig  = figure('Color', 'w', 'Units', 'points', 'Position', [40 40 pdf_w pdf_h]);
    t    = tiledlayout(fig, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    for ip = 1:numel(pas)
        ax = nexttile(t); hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on');
        if ~isempty(bande)
            % La bande du chi2 : ou une NEES honnete doit tomber. Hors legende.
            fill(ax, [1 n_it n_it 1], [bande(1) bande(1) bande(2) bande(2)], ...
                 [1 0.85 0.88], 'EdgeColor', 'none', 'FaceAlpha', 0.55, ...
                 'HandleVisibility', 'off');
        end
        leg = cell(1, numel(cles));
        for im = 1:numel(cles)
            [c, l, leg{im}] = style_modele(cles{im});
            plot(ax, kk(kd), Y(kd, ip, im), l, 'Color', c, 'LineWidth', linewidth);
        end
        % Le premier virage, etiquete dans la premiere case seulement : six fois la meme
        % etiquette chargerait la grille sans rien apprendre de plus.
        marqueur(k_vir, fontsize, ip == 1, ax);
        hold(ax, 'off');
        title(ax, sprintf('$%.1f$ km', pas(ip)), 'Interpreter', 'latex', 'FontSize', fontsize);
        set(ax, 'TickLabelInterpreter', 'latex', 'FontSize', fontsize);
        xlim(ax, [1 n_it]);
        ylim(ax, lim);
        if logy
            % UNE GRADUATION PAR DECADE QUAND L'AXE EN COUVRE PEU, UNE SUR DEUX AU-DELA, et
            % alors sans le quadrillage intermediaire : sur les cinq decades de la NEES du
            % bilineaire, cinq etiquettes et huit pointilles par decade chargeaient les cases
            % sans rien apprendre. C'est l'aspect qu'avait cette figure avant le 15 septembre.
            dec  = floor(log10(lim(1))) : ceil(log10(lim(2)));
            saut = 1 + (numel(dec) > 5);
            set(ax, 'YScale', 'log', 'YTick', 10.^dec(mod(dec, saut) == 0));
            if saut > 1, set(ax, 'YMinorTick', 'off', 'YMinorGrid', 'off'); end
        end
        if ip == 1
            % Au-dessus de la grille et non dans une case : dans la premiere, elle cachait
            % le pic de NEES du krigeage vers le pas 70. Echantillons de trait de 18 pt et non
            % 30, ici comme dans les autres legendes horizontales : avec le nom complet du
            % CA-OK, la legende debordait de la page a droite.
            lg = legend(ax, leg, 'Interpreter', 'latex', 'FontSize', fontsize, ...
                        'Orientation', 'horizontal');
            lg.ItemTokenSize = [18 18];      % propriete cachee, que legend() refuse en argument
            lg.Layout.Tile   = 'north';
        end
        if ip < numel(pas) - 1, set(ax, 'XTickLabel', []); end
        if mod(ip, 2) == 0,     set(ax, 'YTickLabel', []); end
    end
    xlabel(t, 'Time step $k$', 'Interpreter', 'latex', 'FontSize', fontsize);
    ylabel(t, ylab, 'Interpreter', 'latex', 'FontSize', fontsize);
    exporter(fig, nom_fig, pdf_w, pdf_h, img_dir);
end

function [c, l, nom] = style_modele(cle)
% COULEUR, TRAIT ET NOM D'UNE METHODE, les memes dans toutes les figures et tous les
% tableaux du chapitre. Un lecteur qui a appris la couleur du krigeage statique sur une
% figure la retrouve sur les autres, et la FORME du trait dit a quelle famille la methode
% appartient :
%
%   pointille  la borne superieure : la carte connue, qu'aucun estimateur ne depasse
%   plein      les trois estimateurs que le chapitre defend : statique, AK, CA-OK
%   brise      les alternatives auxquelles il se compare : bilineaire, eclairci, reunion
%
% La forme porte donc l'argument et la couleur l'identite, ce qui laisse une figure lisible
% en noir et blanc. LES COULEURS SEPARENT D'ABORD LES PAIRES QU'UNE MEME FIGURE SUPERPOSE :
% statique et AK, bleu et orange ; AK et CA-OK, orange et vert ; AK et reunion, orange et
% violet ; statique et eclairci, bleu et magenta. L'AK etait sarcelle, trop proche du bleu
% statique pour que la figure qui les compare se lise.
%
% ON APPARIE SUR UNE CLE COURTE ET EXACTE, jamais sur un nom de modele ni sur une
% etiquette : 'OK' et 'carte' figurent dans plusieurs noms de modeles, et l'appariement par
% motif a deja fausse deux figures en silence. Une cle inconnue est une erreur, pas un gris.
%
% Les noms sont ceux du manuscrit, arretes le 15 septembre. Celui de la reunion ne l'est
% pas encore.
    switch cle
        case 'connue',   c = [0.00 0.00 0.00];  l = ':';   nom = 'Known map';
        case 'bilin',    c = [0.85 0.15 0.15];  l = '--';  nom = 'Bilinear interpolation';
        case 'statique', c = [0.20 0.30 0.65];  l = '-';   nom = 'Static kriging';
        case 'eclairci', c = [0.80 0.30 0.60];  l = '-.';  nom = 'Random subsampling';
        case 'ak',       c = [0.93 0.50 0.05];  l = '-';   nom = 'Adaptive kriging';
        case 'cak',      c = [0.10 0.62 0.25];  l = '-';   nom = 'Clustered adaptive kriging';
        case 'reunion',  c = [0.50 0.35 0.70];  l = '--';  nom = 'Neighbour union';  % provisoire
        otherwise
            error('figures_campagnes:style', 'Cle de methode inconnue : %s.', cle);
    end
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

function marqueur(k_vir, fontsize, etiquette, ax)
% Le debut du premier demi-tour. C'est lui qui leve l'ambiguite en echantillonnant le
% champ dans une seconde direction, donc les courbes ne se lisent pas sans lui. etiquette
% dit s'il faut ecrire « first turn » le long du trait, ax sur quels axes le tracer.
    if nargin < 3, etiquette = true; end
    if nargin < 4, ax = gca; end
    if isnan(k_vir), return, end
    txt = '';
    if etiquette, txt = 'first turn'; end
    xl = xline(ax, k_vir, '--', txt, 'LineWidth', 0.75, 'Color', [0.4 0.4 0.4], ...
               'HandleVisibility', 'off', 'Interpreter', 'latex', 'FontSize', fontsize, ...
               'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'left');
    % En haut de l'axe : en bas, l'etiquette traversait les courbes posees sur le plancher
    % n_min, qu'elles atteignent toutes avant le premier virage.
    uistack(xl, 'bottom');
end
