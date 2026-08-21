% KRIGEAGE_3D
%
% Ordinary kriging of a scalar field over a two-dimensional index set, drawn as
% surfaces. The name follows the figure it produces, "Exemple de krigeage 3D",
% the third dimension being the field value rather than an input.
%
% The field is the two-dimensional analogue of the one used by
% krigeage_simple_ordinaire.m, so that the two figures read as a family:
%
%   f(x1, x2) = (x1 sin(x1) + x2 cos(x2)) / 2
%
% Its amplitude grows with each coordinate, which keeps the extrapolation
% behaviour visible at the far corners of the domain.
%
% Ordinary kriging is the estimator the thesis retains, so it is the one drawn
% here. Equations are those of the manuscript, chapter 3:
%   - Gram matrix        : eq:Gram_matrix
%   - simple kriging     : eq:simple_kriging_equations    (th:simple_kriging)
%   - ordinary kriging   : eq:ordinary_kriging_equations  (th:ordinary_kriging)
%                          eq:kriged_ordinary_mean
%
% Scalar case d_z = 1 over d_x = 2, so that C = ones(n_tilde, 1) and m_ok is a
% scalar. The kernel routine below is written for any d_x.

clear; close all; clc;
rng(123456789);                  % reproducibility

%% Parameters
x_min     = 0;                   % square domain [x_min, x_max]^2
x_max     = 10;
n_side    = 6;                   % training design: n_side^2 jittered grid points
jitter    = 0.35;                % jitter amplitude, in units of the grid step
sigma_n   = 0.5;                 % observation noise standard deviation
sigma_f   = 5;                   % kernel standard deviation
ell       = 2;                   % correlation length
R_query   = 0;                   % query noise: 0 estimates the field itself
n_plot    = 60;                  % prediction grid, n_plot^2 points
z_lim     = [-15 15];            % vertical axis, widened for the 3 sigma band
img_dir   = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'Images');
                                 % the manuscript and the notes read figures here,
                                 % resolved from the script rather than from pwd
img_name  = 'Exemple de krigeage 3D';
linewidth = 1;
fontsize  = 14;                  % 14 x 298/348 = 12.0 pt once included at 0.7\textwidth
pdf_w     = 348;                 % frozen page, in points, shared by the three scripts
pdf_h     = 280;                 % same page as the 1D family, thanks to the view angle
alpha_band  = 0.55;              % opacity of the +/- 3 sigma envelope, see below
alpha_field = 0.35;              % opacity of the true field surface
vector_pdf = true;               % false rasterises, see the note in save_kriging_3d

R_tilde   = sigma_n^2;           % training noise variance

%% Field and data
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

%% Gram matrix, eq:Gram_matrix
G = se_kernel(X_train, X_train, sigma_f, ell) + R_tilde * eye(n_train);
L = chol(G, 'lower');            % solve via Cholesky rather than inv(G)

K_star = se_kernel(X_star, X_train, sigma_f, ell);     % k(x, x_tilde)
k_xx   = sigma_f^2;                                    % k(x, x), constant here

V      = L \ K_star';                                  % L^{-1} k(x_tilde, x)
R_sk   = k_xx + R_query - sum(V.^2, 1)';               % simple kriging covariance

Ginv_z = L' \ (L \ z_train);

%% Ordinary kriging, eq:ordinary_kriging_equations and eq:kriged_ordinary_mean
C       = ones(n_train, 1);                            % C = 1_n kron I_1
Ginv_C  = L' \ (L \ C);
Minv_ok = C' * Ginv_C;                                 % M^{-1}, kept implicit

m_ok    = Minv_ok \ (C' * Ginv_z);                     % estimated constant mean
z_ok    = m_ok + K_star * (L' \ (L \ (z_train - C * m_ok)));

U_ok    = 1 - K_star * Ginv_C;                         % U' = I - k(x,x_tilde) G^{-1} C
R_mok   = sum((U_ok / Minv_ok) .* U_ok, 2);            % U' M U, cost of estimating m
R_ok    = R_mok + R_sk;

Z_ok    = reshape(z_ok,  n_plot, n_plot);
S_ok    = reshape(sqrt(max(R_ok, 0)), n_plot, n_plot); % guard against -0 round-off

%% Diagnostics
report('Ordinary kriging', f_star(:), z_ok, R_ok);
fprintf('\nEstimated mean m_ok = %.2f over %d samples\n', m_ok, n_train);

%% Figure
save_kriging_3d(img_name, P1, P2, f_star, Z_ok, S_ok, X_train, z_train, ...
                z_lim, linewidth, fontsize, alpha_band, alpha_field, ...
                pdf_w, pdf_h, vector_pdf, img_dir);

%% Local functions
function K = se_kernel(Xa, Xb, sigma_f, ell)
% Squared exponential kernel over any input dimension, one point per row.
% k(x, x') = sigma_f^2 exp(-||x - x'||^2 / (2 ell^2)).
    D2 = sum(Xa.^2, 2) - 2 * Xa * Xb' + sum(Xb.^2, 2)';
    K  = sigma_f^2 * exp(-max(D2, 0) / (2 * ell^2));   % max guards round-off
end

function report(name, f_true, z_hat, R_hat)
% Root mean squared error, +/- 3 sigma coverage and mean band width.
    s        = sqrt(max(R_hat, 0));
    rmse     = sqrt(mean((f_true - z_hat).^2));
    coverage = mean(abs(f_true - z_hat) <= 3 * s);
    fprintf('%s : RMSE = %6.3f, coverage = %5.1f %%, mean band = %5.2f\n', ...
            name, rmse, 100 * coverage, mean(6 * s));
end

function save_kriging_3d(name, P1, P2, F, Z, S, X_train, z_train, z_lim, ...
                         linewidth, fontsize, alpha_band, alpha_field, ...
                         pdf_w, pdf_h, vector_pdf, img_dir)
% Surfaces of the true field, of the kriging mean and of its +/- 3 sigma
% envelope, with the samples. No title: the caption carries the naming.
%
% The one-dimensional figures fill the band with an opaque patch, the curves
% being drawn on top of it in the plane. That is not available here: the upper
% envelope lies between the camera and the field, so full opacity would hide
% everything it is meant to bracket. alpha_band is therefore the closest the
% surfaces get to the flat pink of those figures while keeping them readable.
    pink = [1 0.85 0.88];
    fig  = figure('Color', 'w', 'Units', 'points', 'Position', [80 80 pdf_w pdf_h]);

    surf(P1, P2, Z + 3 * S, 'FaceColor', pink, 'FaceAlpha', alpha_band, ...
         'EdgeColor', 'none'); hold on;
    lower = surf(P1, P2, Z - 3 * S, 'FaceColor', pink, 'FaceAlpha', alpha_band, ...
                 'EdgeColor', 'none');
    lower.Annotation.LegendInformation.IconDisplayStyle = 'off';  % one entry only

    surf(P1, P2, F, 'FaceColor', 'b', 'FaceAlpha', alpha_field, 'EdgeColor', 'none');
    mesh(P1, P2, Z, 'EdgeColor', 'r', 'FaceColor', 'none', ...
         'LineWidth', 0.5 * linewidth);
    plot3(X_train(:,1), X_train(:,2), z_train, 'kx', 'MarkerSize', 8, ...
          'LineWidth', linewidth);
    hold off;

    box on; grid on; view(-37.5, 42);
    xlim([P1(1,1) P1(1,end)]); ylim([P2(1,1) P2(end,1)]); zlim(z_lim);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fontsize);
    % No legend: the LaTeX caption already names every element by colour.

    if ~exist(img_dir, 'dir')
        mkdir(img_dir);
    end
    % print rather than exportgraphics: the latter crops to the content, so the
    % page would not match the one-dimensional family. Here it is exactly
    % pdf_w x pdf_h, hence the same scale and the same apparent font size.
    target = fullfile(img_dir, [name '.pdf']);
    set(fig, 'PaperUnits', 'points', 'PaperSize', [pdf_w pdf_h], ...
             'PaperPosition', [0 0 pdf_w pdf_h]);
    if vector_pdf
        % Transparency may force MATLAB to rasterise part of the page. If the
        % surfaces come out flat or banded, set vector_pdf = false.
        print(fig, target, '-dpdf', '-vector');
    else
        print(fig, target, '-dpdf', '-image', '-r300');
    end
end
