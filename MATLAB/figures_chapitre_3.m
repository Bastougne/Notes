% FIGURES_CHAPITRE_3
%
% Les deux illustrations du chapitre 3, exportées en PDF vectoriel sous Images/ :
%
%   "Noyaux de Matérn"        profil radial des noyaux de Matérn de régularité
%                             demi-entière et de leur limite exponentielle quadratique,
%                             à variance et longueur de corrélation unitaires, de sorte
%                             que seules la forme près de l'origine et la masse de queue
%                             les distinguent. Formes fermées de th:Matérn et de
%                             l'expression demi-entière qui le suit. Incluse dans
%                             3_Gaussian_processes.tex.
%   "Exemple de krigeage 3D"  krigeage ordinaire d'un champ scalaire sur un ensemble
%                             d'indices bidimensionnel, tracé en surfaces — la troisième
%                             dimension est la valeur du champ et non une entrée. Le
%                             champ est l'analogue bidimensionnel de celui du chapitre 4,
%                             f(x1, x2) = (x1 sin(x1) + x2 cos(x2)) / 2, dont l'amplitude
%                             croît avec chaque coordonnée : le comportement en
%                             extrapolation reste visible aux coins du domaine. Équations
%                             eq:Gram_matrix, eq:ordinary_kriging_equations et
%                             eq:kriged_ordinary_mean. Incluse dans 3_kriging.tex.
%
% LA PAGE EST FIGÉE à pdf_w x pdf_h points pour les deux, et l'export passe par print
% plutôt que par exportgraphics : ce dernier recadre au contenu, donc la page varierait
% d'une figure à l'autre et la police apparente avec elle. À 348 x 280 points inclus à
% 0.7\textwidth, une police de 14 pt en donne 12,0 dans le manuscrit.

clear; close all; clc;
rng(123456789);                  % reproductibilité de la partie krigée

%% Ce que les deux figures partagent
img_dir   = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'Images');
                                 % le manuscrit et les notes y lisent leurs figures,
                                 % résolu depuis le script et non depuis pwd
linewidth = 1;
fontsize  = 14;
pdf_w     = 348;                 % page figée, en points
pdf_h     = 280;

%% ------------------------------------------------------- Noyaux de Matérn
ell_k  = 1;                      % longueur caractéristique
r_max  = 5;                      % distance maximale tracée
n_r    = 1000;                   % résolution

r = linspace(0, r_max, n_r);

k_12  = exp(-r ./ ell_k);
k_32  = (1 + sqrt(3) .* r ./ ell_k) .* exp(-sqrt(3) .* r ./ ell_k);
k_52  = (1 + sqrt(5) .* r ./ ell_k + 5 .* r.^2 ./ (3 .* ell_k.^2)) ...
        .* exp(-sqrt(5) .* r ./ ell_k);
k_inf = exp(-r.^2 ./ (2 .* ell_k.^2));   % nu -> inf : exponentiel quadratique

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

exporter(fig, 'Noyaux de Matérn', pdf_w, pdf_h, img_dir, true);

%% -------------------------------------------------- Exemple de krigeage 3D
x_min     = 0;                   % domaine carré [x_min, x_max]^2
x_max     = 10;
n_side    = 6;                   % plan d'expérience : n_side^2 points de grille secoués
jitter    = 0.35;                % amplitude du secouage, en pas de grille
sigma_n   = 0.5;                 % écart-type du bruit d'observation
sigma_f   = 5;                   % écart-type du noyau
ell       = 2;                   % longueur de corrélation
R_query   = 0;                   % bruit en requête : 0 estime le champ lui-même
n_plot    = 60;                  % grille de prédiction, n_plot^2 points
z_lim     = [-15 15];            % axe vertical, élargi pour la bande à 3 sigma
alpha_band  = 0.55;              % opacité de l'enveloppe +/- 3 sigma, voir plus bas
alpha_field = 0.35;              % opacité de la surface du champ vrai
vector_pdf  = true;              % false rastérise, voir la note dans exporter

R_tilde   = sigma_n^2;           % variance du bruit d'apprentissage

f = @(X) (X(:,1) .* sin(X(:,1)) + X(:,2) .* cos(X(:,2))) / 2;

step    = (x_max - x_min) / n_side;
centres = x_min + step * (0.5:1:n_side-0.5);
[T1, T2] = meshgrid(centres, centres);
X_train  = [T1(:), T2(:)] + jitter * step * (2 * rand(n_side^2, 2) - 1);
n_train  = size(X_train, 1);
z_train  = f(X_train) + sigma_n * randn(n_train, 1);

side    = linspace(x_min, x_max, n_plot);
[P1, P2] = meshgrid(side, side);
X_star   = [P1(:), P2(:)];
f_star   = reshape(f(X_star), n_plot, n_plot);

% Matrice de Gram, eq:Gram_matrix
G = se_kernel(X_train, X_train, sigma_f, ell) + R_tilde * eye(n_train);
L = chol(G, 'lower');            % résolution par Cholesky plutôt que inv(G)

K_star = se_kernel(X_star, X_train, sigma_f, ell);     % k(x, x_tilde)
k_xx   = sigma_f^2;                                    % k(x, x), constant ici

V      = L \ K_star';                                  % L^{-1} k(x_tilde, x)
R_sk   = k_xx + R_query - sum(V.^2, 1)';               % covariance du krigeage simple

Ginv_z = L' \ (L \ z_train);

% Krigeage ordinaire, eq:ordinary_kriging_equations et eq:kriged_ordinary_mean
C       = ones(n_train, 1);                            % C = 1_n kron I_1
Ginv_C  = L' \ (L \ C);
Minv_ok = C' * Ginv_C;                                 % M^{-1}, gardé implicite

m_ok    = Minv_ok \ (C' * Ginv_z);                     % moyenne constante estimée
z_ok    = m_ok + K_star * (L' \ (L \ (z_train - C * m_ok)));

U_ok    = 1 - K_star * Ginv_C;                         % U' = I - k(x,x_tilde) G^{-1} C
R_mok   = sum((U_ok / Minv_ok) .* U_ok, 2);            % U' M U, le prix d'estimer m
R_ok    = R_mok + R_sk;

Z_ok    = reshape(z_ok,  n_plot, n_plot);
S_ok    = reshape(sqrt(max(R_ok, 0)), n_plot, n_plot); % contre les -0 d'arrondi

report('Krigeage ordinaire', f_star(:), z_ok, R_ok);
fprintf('\nMoyenne estimée m_ok = %.2f sur %d échantillons\n', m_ok, n_train);

figure_krigeage_3d(P1, P2, f_star, Z_ok, S_ok, X_train, z_train, z_lim, ...
                   linewidth, fontsize, alpha_band, alpha_field, pdf_w, pdf_h, ...
                   vector_pdf, img_dir);

%% ------------------------------------------------------------ fonctions locales

function exporter(fig, nom, pdf_w, pdf_h, img_dir, vectoriel)
% L'export partagé par les deux figures : page figée, donc même échelle et même police
% apparente partout. La transparence peut forcer MATLAB à rastériser une partie de la
% page ; si les surfaces sortent plates ou moirées, passer vectoriel à false.
%
% Le PDF n'est remplacé que si le dessin a changé, comme dans figures_campagnes.m. MATLAB
% inscrit dans chaque export sa date de création et un identifiant /ID recalculé à chaque
% fois, et git signalait sinon chaque figure à chaque passage, même inchangée.
    if ~exist(img_dir, 'dir')
        mkdir(img_dir);
    end
    set(fig, 'PaperUnits', 'points', 'PaperSize', [pdf_w pdf_h], ...
             'PaperPosition', [0 0 pdf_w pdf_h]);
    cible = fullfile(img_dir, [nom '.pdf']);
    tmp   = [tempname '.pdf'];
    if vectoriel
        print(fig, tmp, '-dpdf', '-vector');
    else
        print(fig, tmp, '-dpdf', '-image', '-r300');
    end
    if isfile(cible) && strcmp(sans_horodatage(tmp), sans_horodatage(cible))
        delete(tmp);
    else
        movefile(tmp, cible, 'f');
    end
end

function s = sans_horodatage(f)
% Les octets d'un PDF sans sa date ni son identifiant. Lu en octets et non en caractères :
% un décodage confondrait des flux compressés différents.
    fid = fopen(f, 'r');
    s   = char(fread(fid, Inf, '*uint8')');
    fclose(fid);
    s = regexprep(s, '/(CreationDate|ModDate)\s*\([^)]*\)', '');
    s = regexprep(s, '/ID\s*\[[^\]]*\]', '');
end

function K = se_kernel(Xa, Xb, sigma_f, ell)
% Noyau exponentiel quadratique en dimension d'entrée quelconque, un point par ligne.
% k(x, x') = sigma_f^2 exp(-||x - x'||^2 / (2 ell^2)).
    D2 = sum(Xa.^2, 2) - 2 * Xa * Xb' + sum(Xb.^2, 2)';
    K  = sigma_f^2 * exp(-max(D2, 0) / (2 * ell^2));   % max contre l'arrondi
end

function report(nom, f_true, z_hat, R_hat)
% Erreur quadratique moyenne, couverture à +/- 3 sigma et largeur moyenne de bande.
    s        = sqrt(max(R_hat, 0));
    rmse     = sqrt(mean((f_true - z_hat).^2));
    coverage = mean(abs(f_true - z_hat) <= 3 * s);
    fprintf('%s : RMSE = %6.3f, couverture = %5.1f %%, bande moyenne = %5.2f\n', ...
            nom, rmse, 100 * coverage, mean(6 * s));
end

function figure_krigeage_3d(P1, P2, F, Z, S, X_train, z_train, z_lim, ...
                            linewidth, fontsize, alpha_band, alpha_field, ...
                            pdf_w, pdf_h, vector_pdf, img_dir)
% Surfaces du champ vrai, de la moyenne krigée et de son enveloppe à +/- 3 sigma, avec
% les échantillons. Pas de titre : la légende LaTeX porte le nommage.
%
% Les figures unidimensionnelles remplissent la bande d'un aplat opaque, les courbes
% étant tracées par-dessus dans le plan. Ce n'est pas possible ici : l'enveloppe
% supérieure est entre la caméra et le champ, donc une opacité totale masquerait tout ce
% qu'elle est censée encadrer. alpha_band est le plus proche de l'aplat rose de ces
% figures qui reste lisible.
    pink = [1 0.85 0.88];
    fig  = figure('Color', 'w', 'Units', 'points', 'Position', [80 80 pdf_w pdf_h]);

    surf(P1, P2, Z + 3 * S, 'FaceColor', pink, 'FaceAlpha', alpha_band, ...
         'EdgeColor', 'none'); hold on;
    lower = surf(P1, P2, Z - 3 * S, 'FaceColor', pink, 'FaceAlpha', alpha_band, ...
                 'EdgeColor', 'none');
    lower.Annotation.LegendInformation.IconDisplayStyle = 'off';  % une seule entrée

    surf(P1, P2, F, 'FaceColor', 'b', 'FaceAlpha', alpha_field, 'EdgeColor', 'none');
    mesh(P1, P2, Z, 'EdgeColor', 'r', 'FaceColor', 'none', ...
         'LineWidth', 0.5 * linewidth);
    plot3(X_train(:,1), X_train(:,2), z_train, 'kx', 'MarkerSize', 8, ...
          'LineWidth', linewidth);
    hold off;

    box on; grid on; view(-37.5, 42);
    xlim([P1(1,1) P1(1,end)]); ylim([P2(1,1) P2(end,1)]); zlim(z_lim);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fontsize);
    % Pas de légende : la légende LaTeX nomme déjà chaque élément par sa couleur.

    exporter(fig, 'Exemple de krigeage 3D', pdf_w, pdf_h, img_dir, vector_pdf);
end
