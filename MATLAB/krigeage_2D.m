% KRIGEAGE_2D
%
% Simple, ordinary and universal kriging of the field f(x) = x*sin(x) + m(x), on a
% one-dimensional index set. Figures are written to Images/ under the names
% "Exemple de krigeage ... 2D.pdf", to be assembled in LaTeX.
%
% Five estimators are compared on the very same samples, one figure each so
% that they can be assembled in LaTeX rather than as MATLAB subplots:
%   1. simple kriging given the true mean m,
%   2. simple kriging given a wrong mean,
%   3. ordinary kriging, which estimates the constant mean from the data,
%   4. universal kriging over an affine trend basis,
%   5. universal kriging over a high-order polynomial basis.
%
% Case 2 is the first point of the script. Its error covariance R_sk does not
% depend on m at all, so a wrong mean shifts the estimate without widening the
% band by a single unit: the estimator is not merely inaccurate, it is
% over-confident. Ordinary and universal kriging pay for the unknown mean
% through the extra terms R_mok and R_muk, which widen the band precisely where
% the data no longer constrain the mean. That price grows with the number of
% estimated coefficients, and shrinks as samples accumulate.
%
% That price is tabulated rather than drawn. Nesting a second band per figure was
% tried and dropped: the contribution of the mean is a couple of percent of the
% half-width for ordinary kriging, so its two boundaries read as one thick line
% precisely in the case the manuscript argues for. The figures therefore keep a
% single band, and the decomposition goes to the table written by latex_table.
%
% Case 5 is the second point. A basis rich enough to chase the samples absorbs
% structure that the GP should have carried, and extrapolates on its own terms
% once the data run out. More coefficients therefore do not buy a better field,
% they buy a trend that misbehaves off-support, which is the objection raised
% against fitting a smooth global basis to a field whose information lies in
% its local anomaly.
%
% Equations are those of the manuscript, chapter 3:
%   - Gram matrix        : eq:Gram_matrix
%   - simple kriging     : eq:simple_kriging_equations    (th:simple_kriging)
%   - ordinary kriging   : eq:ordinary_kriging_equations  (th:ordinary_kriging)
%                          eq:kriged_ordinary_mean
%   - universal kriging  : eq:universal_kriging_equations (th:universal_kriging)
%                          eq:kriged_universal_mean
%
% Scalar case d_z = 1, so that C = ones(n_tilde, 1) and m_ok is a scalar.

clear; close all; clc;
rng(123456789);                  % reproducibility

%% Parameters (matched to the figure)
m_true    = @(x) 10 * ones(size(x));   % true mean of the field, and the one handed to
                                 % the well-informed simple kriging. Any function of x
                                 % is admissible, a constant being what the figures use
m_wrong   = @(x) zeros(size(x)); % mean handed to the misinformed simple kriging
x_min     = -3;                  % plotting domain
x_max     = 13;
x_train   = (0.5:1:9.5)';        % training locations, unit spacing
x_gap     = [2 5];               % interval cleared of training points, [] for none.
                                 % A hole inside the sampled range is the realistic
                                 % case for a map, and the only place where the
                                 % estimators can be told apart without leaving it.
                                 % The gap region reported below is not this interval
                                 % but the one the surviving points actually bound.
sigma_n   = 1;                   % observation noise standard deviation
sigma_f   = 5;                   % kernel standard deviation
ell       = 1.2;                   % correlation length
R_query   = 0;                   % query noise: 0 estimates the field f(x),
                                 % sigma_n^2 would predict a new observation z
n_plot    = 400;                 % resolution of the prediction grid
n_mc      = 1000;                % noise realisations behind the tabulated errors.
                                 % n_mc = 1 falls back to the single-draw behaviour,
                                 % and the figures always show the first draw
y_lim     = [-15 30];            % vertical range shared by the five panels, to be
                                 % widened if a non-constant mean moves the field
deg_poly  = 5;                   % degree of the over-rich trend basis, case 5
root      = fileparts(fileparts(mfilename('fullpath')));  % repository root, one level up
img_dir   = fullfile(root, 'Images');                     % the manuscript and the notes read their figures here
tab_dir   = fullfile(root, 'Latex', 'Manuscrit', 'include');   % read by \input, next to the chapters
linewidth = 1;
pdf_w     = 348;                 % frozen page, in points, shared by the three scripts
pdf_h     = 280;
textwidth = 425;                 % pt of text on the page, measured on the class
tex_frac  = 0.49;                % the five panels are placed two per row in chapter 4,
                                 % not one per row like the 3D and Matérn figures
fontsize  = 12 * pdf_w / (tex_frac * textwidth);   % 20.1, so that the labels render at
                                 % 12.0 pt, the body size. At 14 they came out at 8.4:
                                 % the convention had been measured on the single-panel
                                 % figures and never rechecked once these went in pairs
alpha_band = 0.55;               % opacity of the band, shared with krigeage_3D.m

if ~isempty(x_gap)
    x_train = x_train(x_train < x_gap(1) | x_train > x_gap(2));
end
n_train   = numel(x_train);
R_tilde   = sigma_n^2;           % training noise variance

%% Data
% One column of observations per Monte Carlo realisation. The estimators are all
% linear in the observations, so the whole experiment goes through the same
% triangular solves as a single draw, and every quantity below carries its
% realisations along its second dimension. Column 1 is drawn first, hence is the
% draw a single-realisation run would have produced, and it is the one the
% figures show.
f        = @(x) x .* sin(x) + 10;
z_mc     = f(x_train) + sigma_n * randn(n_train, n_mc);
z_train  = z_mc(:, 1);
x_star   = linspace(x_min, x_max, n_plot)';
f_star   = f(x_star);

%% Gram matrix, eq:Gram_matrix
G = se_kernel(x_train, x_train, sigma_f, ell) + R_tilde * eye(n_train);
L = chol(G, 'lower');            % solve via Cholesky rather than inv(G)

K_star   = se_kernel(x_star, x_train, sigma_f, ell);   % k(x, x_tilde)
k_xx     = sigma_f^2;                                  % k(x, x), constant here

V        = L \ K_star';                                % L^{-1} k(x_tilde, x)
R_sk     = k_xx + R_query - sum(V.^2, 1)';             % simple kriging covariance

Ginv_z   = L' \ (L \ z_mc);

%% Simple kriging, eq:simple_kriging_equations
% The mean is known, so it is subtracted before kriging and added back after. It
% is evaluated at the training inputs on one side and at the query points on the
% other, which is what leaves a non-constant mean admissible without touching the
% weights. R_sk is computed once and reused for both means: it does not involve m,
% and no more does it involve the observations, hence it is shared by every draw.
z_right  = simple_kriging(m_true,  x_train, x_star, L, K_star, z_mc);
z_wrong  = simple_kriging(m_wrong, x_train, x_star, L, K_star, z_mc);

%% Ordinary kriging, eq:ordinary_kriging_equations and eq:kriged_ordinary_mean
C        = ones(n_train, 1);                           % C = 1_n kron I_1
Ginv_C   = L' \ (L \ C);
Minv_ok  = C' * Ginv_C;                                % M^{-1}, kept implicit

m_ok     = Minv_ok \ (C' * Ginv_z);                    % estimated mean, one per draw
z_ok     = m_ok + K_star * (L' \ (L \ (z_mc - C * m_ok)));

U_ok     = 1 - K_star * Ginv_C;                        % U' = I - k(x,x_tilde) G^{-1} C
R_mok    = sum((U_ok / Minv_ok) .* U_ok, 2);           % U' M U, cost of estimating m
R_ok     = R_mok + R_sk;

%% Universal kriging, eq:universal_kriging_equations and eq:kriged_universal_mean
% Affine trend basis g = [1, x], i.e. r_UK = 2, which the GP then corrects.
g_aff_train = [ones(n_train, 1), x_train];
g_aff_star  = [ones(n_plot, 1),  x_star];
[z_uk, R_uk, m_uk] = universal_kriging(g_aff_train, g_aff_star, L, K_star, ...
                                       z_mc, R_sk);

% Polynomial basis of degree deg_poly, i.e. r_UK = deg_poly + 1 coefficients for
% n_train samples. The variable is centred and scaled on the training range, so
% that what the figure shows is the effect of the basis and not of conditioning.
g_pol_train = poly_basis(x_train, deg_poly, x_train);
g_pol_star  = poly_basis(x_star,  deg_poly, x_train);
[z_pol, R_pol, m_pol] = universal_kriging(g_pol_train, g_pol_star, L, K_star, ...
                                          z_mc, R_sk);

%% Diagnostics
% Every quantity of the table is a function of x collapsed over the grid, so the
% grid is split at the training range first. A single average over [x_min, x_max]
% would bury the whole point, since the cost of an estimated mean is nil where
% the samples constrain it and only opens up beyond them; it would also depend on
% how much extrapolation the plotting domain happens to contain, which is a free
% parameter of the figure and not a property of the estimator.
hull = x_star >= min(x_train) & x_star <= max(x_train);
if isempty(x_gap)
    regions = struct('name', {'in', 'out'}, ...
                     'mask', {hull, ~hull});
else
    % The gap is bounded by the two training points that survive on either side of
    % the cleared interval, not by the interval itself: what an estimator has to
    % bridge is the distance between the samples it still has, and taking the
    % widest such hole keeps the definition independent of how x_gap was written.
    d         = diff(x_train);
    [~, w]    = max(d);
    gap_edges = [x_train(w), x_train(w + 1)];
    gap       = x_star >= gap_edges(1) & x_star <= gap_edges(2);
    regions   = struct('name', {'in', 'gap', 'out'}, ...
                       'mask', {hull & ~gap, gap, ~hull});
    fprintf('Gap bounded by the samples at %.1f and %.1f\n\n', gap_edges);
end

% Simple kriging is handed R_sk as its own reference, which makes its mean share
% exactly zero rather than undefined. That is the honest value: its covariance
% carries no mean estimation term at all, and a column of zeros states the point
% the section is making better than a column of dashes would.
stats = [kriging_stats('Simple, correct mean',  f_star, z_right, R_sk,  R_sk, regions)
         kriging_stats('Simple, wrong mean',    f_star, z_wrong, R_sk,  R_sk, regions)
         kriging_stats('Ordinary',              f_star, z_ok,    R_ok,  R_sk, regions)
         kriging_stats('Universal, affine',     f_star, z_uk,    R_uk,  R_sk, regions)
         kriging_stats('Universal, degree $5$', f_star, z_pol,   R_pol, R_sk, regions)];

print_stats(stats, regions, n_mc);
print_dispersion(stats, regions);
print_head_to_head(stats, regions, 'Ordinary');
latex_tables(stats, regions, tab_dir, n_mc);

fprintf('\nField mean over the samples %.2f, estimated at %.2f on the figure draw\n', ...
        mean(m_true(x_train)), m_ok(1));
fprintf('   and at %.2f on average over the %d realisations, standard deviation %.2f\n', ...
        mean(m_ok), n_mc, std(m_ok));
fprintf('Estimated affine trend m_uk = [%.2f, %.2f] on the figure draw\n', m_uk(1, 1), m_uk(2, 1));
% poly_basis works in a reduced variable, where the trend Gram matrix stays well
% conditioned, so its coefficients are not those of the monomial basis announced
% by the scenario. The estimator is left as it is and only the reported
% coefficients are converted, which is exact and costs nothing.
x0_pol = mean(x_train);
s_pol  = max(x_train - x0_pol);
a_pol  = monomial_coefficients(m_pol(:, 1), x0_pol, s_pol);
fprintf('Estimated degree %d trend, %d coefficients for %d samples,\n', ...
        deg_poly, size(m_pol, 1), n_train);
fprintf('   reduced variable t = (x - %.2f) / %.2f : [%s]\n', ...
        x0_pol, s_pol, strtrim(sprintf('%.2f ', m_pol(:, 1))));
% Printed to six significant figures rather than two: in the monomial basis the
% terms cancel over three orders of magnitude, around 800 against a result near
% 5 at the middle of the training range, so a coefficient rounded to two decimals
% no longer reproduces the trend it describes.
fprintf('   monomial basis (1, x, ..., x^%d)        : [%s]\n', ...
        deg_poly, strtrim(sprintf('%.6g ', a_pol)));
fprintf('   check at x = 5, reduced %.4f, monomial %.4f\n', ...
        polyval(flipud(m_pol(:, 1))', (5 - x0_pol) / s_pol), polyval(fliplr(a_pol), 5));

%% Figures, one file each
% Identical axes throughout, so that the five can be placed side by side, and no
% title, the caption of the assembled figure naming the panels instead.
save_kriging('Exemple de krigeage simple 2D', x_star, f_star, z_right(:, 1), R_sk, ...
             x_train, z_train, linewidth, fontsize, alpha_band, pdf_w, pdf_h, y_lim, img_dir);

save_kriging('Exemple de krigeage simple à mauvaise moyenne 2D', x_star, f_star, z_wrong(:, 1), R_sk, ...
             x_train, z_train, linewidth, fontsize, alpha_band, pdf_w, pdf_h, y_lim, img_dir);

save_kriging('Exemple de krigeage ordinaire 2D', x_star, f_star, z_ok(:, 1), R_ok, ...
             x_train, z_train, linewidth, fontsize, alpha_band, pdf_w, pdf_h, y_lim, img_dir);

save_kriging('Exemple de krigeage universel 2D', x_star, f_star, z_uk(:, 1), R_uk, ...
             x_train, z_train, linewidth, fontsize, alpha_band, pdf_w, pdf_h, y_lim, img_dir);

save_kriging('Exemple de krigeage universel à mauvaise base 2D', x_star, f_star, z_pol(:, 1), R_pol, ...
             x_train, z_train, linewidth, fontsize, alpha_band, pdf_w, pdf_h, y_lim, img_dir);

%% Local functions
function z_hat = simple_kriging(m, x_train, x_star, L, K_star, Z)
% Simple kriging under a known mean function, eq:simple_kriging_equations. The
% mean is a handle rather than a scalar, so the residual z_tilde - m(x_tilde) and
% the recentring m(x) are evaluated on their own grids. A constant handle
% reproduces the scalar case exactly, the weights being unchanged either way.
% Z holds one column of observations per realisation, and the two mean terms
% broadcast over them.
    z_hat = m(x_star) + K_star * (L' \ (L \ (Z - m(x_train))));
end

function [z_hat, R_hat, m_hat] = universal_kriging(g_train, g_star, L, K_star, ...
                                                   Z, R_sk)
% Universal kriging over an arbitrary trend basis, eq:universal_kriging_equations.
% The estimation cost is the sandwich U' M U of eq:kriged_universal_mean, with
% M = (g(x_tilde)' G^{-1} g(x_tilde))^{-1} kept implicit.
%
% Z holds one column of observations per realisation, so m_hat comes out r_UK by
% n_mc and z_hat n_plot by n_mc. R_hat stays a single column: neither M nor U
% involves the observations, so the announced covariance is the same for all of
% them, which is why the covariance table cannot be a Monte Carlo average.
    Ginv_g = L' \ (L \ g_train);
    Minv   = g_train' * Ginv_g;                        % M^{-1}, r_UK x r_UK

    m_hat  = Minv \ (g_train' * (L' \ (L \ Z)));       % estimated coefficients
    z_hat  = g_star * m_hat + K_star * (L' \ (L \ (Z - g_train * m_hat)));

    U      = g_star - K_star * Ginv_g;                 % U' = g(x) - k(x,x_tilde) G^{-1} g(x_tilde)
    R_muk  = sum((U / Minv) .* U, 2);                  % U' M U
    R_hat  = R_muk + R_sk;
end

function a_x = monomial_coefficients(a_t, x0, s)
% Rewrites a polynomial given by its coefficients in the reduced variable
% t = (x - x0)/s as its coefficients in the monomial basis (1, x, ..., x^deg).
% Both are ascending in degree. The substitution is carried out by Horner in t,
% each step multiplying by t seen as the degree one polynomial (x - x0)/s.
    deg = numel(a_t) - 1;
    lin = [1 / s, -x0 / s];                            % t, descending as polyval wants
    p   = 0;
    for j = deg:-1:0
        p      = conv(p, lin);
        p(end) = p(end) + a_t(j + 1);
    end
    a_x = fliplr(p(end - deg:end));                    % drop the leading zero, ascend
end

function g = poly_basis(x, deg, x_ref)
% Monomial basis [1, t, ..., t^deg] in the variable t centred and scaled on the
% range of x_ref. Rescaling spans the same functions, so the estimator is
% unchanged, but the Gram matrix of the basis stays well conditioned.
    x0 = mean(x_ref);
    s  = max(x_ref - x0);
    t  = (x - x0) / s;
    g  = t .^ (0:deg);                                 % implicit expansion
end

function K = se_kernel(xa, xb, sigma_f, ell)
% Squared exponential kernel, k(r) = sigma_f^2 exp(-r^2 / (2 ell^2)).
    D = xa - xb';                                      % implicit expansion
    K = sigma_f^2 * exp(-D.^2 / (2 * ell^2));
end

function s = kriging_stats(name, f_true, z_hat, R_hat, R_inner, regions)
% Collapses the x-dependent quantities into scalars, once per region of the grid.
% Each field is a row vector indexed like regions. Three different collapses are
% involved, and they are not interchangeable:
%
%   rmse  : mean over the region of the pointwise RMSE, the root being taken over
%           the realisations at fixed x. The RMSE is a function of the index,
%           space here and time wherever a trajectory is estimated, and any
%           aggregate is an arithmetic mean of that function. Rooting a single
%           pooled average instead would sit above this by Jensen, and would no
%           longer be the mean of a curve that can be plotted on its own,
%   hw    : arithmetic mean of the confidence half-width, that is the area of the
%           half-band divided by the length of the region. Free of the draw,
%           since no kriging covariance involves the observations. The band is
%           still drawn and measured at three standard deviations, only the
%           manuscript no longer names the multiple,
%   share : arithmetic mean of the pointwise 1 - R_inner / R_hat, the fraction of
%           the kriging covariance contributed by the estimation of the mean.
%
% The share is taken on covariances and not on half-widths, so that it is a true
% proportion: R_hat = R_mok + R_sk is an additive decomposition, which standard
% deviations would not preserve. It also averages ratios rather than taking a
% ratio of averages, so that every x weighs the same, a region where the band is
% widest being otherwise left to dictate the figure on its own.
    sd  = sqrt(max(R_hat, 0));                         % guard against -0 round-off
    err = f_true - z_hat;
    if isempty(R_inner)
        r = [];                                        % no reference, nothing to split
    else
        r = 1 - max(R_inner, 0) ./ max(R_hat, eps);
    end

    s.name    = name;
    s.rmse    = zeros(1, numel(regions));
    s.rmse_sd = zeros(1, numel(regions));
    s.draw    = zeros(numel(regions), size(err, 2));
    s.hw      = zeros(1, numel(regions));
    s.share   = nan(1, numel(regions));
    for j = 1:numel(regions)
        w = regions(j).mask;
        e = err(w, :);
        s.draw(j, :) = sqrt(mean(e.^2, 1));            % RMSE of each realisation
        s.rmse(j)    = mean(sqrt(mean(e.^2, 2)));      % mean over E of the pointwise RMSE
        s.rmse_sd(j) = std(s.draw(j, :));
        s.hw(j)      = mean(3 * sd(w));
        if ~isempty(r)
            s.share(j) = mean(r(w));
        end
    end
end

function print_stats(stats, regions, n_mc)
% The numbers of the table, in a form that survives a console.
    n = numel(regions);
    for j = 1:n
        fprintf('%-14s %4d grid points\n', regions(j).name, sum(regions(j).mask));
    end
    fprintf('%-14s %4d noise realisations\n', 'Monte Carlo', n_mc);
    fprintf('\n%-22s', '');
    for label = {'ARMSE', 'confidence hw', 'mean share'}
        fprintf(' %*s', 8 * n - 1, label{1});
    end
    fprintf('\n%-22s', '');
    for j = 1:3 * n
        fprintf(' %7s', regions(mod(j - 1, n) + 1).name(1:min(7, end)));
    end
    fprintf('\n');
    for i = 1:numel(stats)
        fprintf('%-22s', stats(i).name);
        fprintf(' %7.2f', stats(i).rmse);
        fprintf(' %7.2f', stats(i).hw);
        for j = 1:n
            fprintf(' %7s', share_str(stats(i).share(j)));
        end
        fprintf('\n');
    end
end

function print_dispersion(stats, regions)
% Spread of the per-realisation RMSE. It says whether a gap between two rows of the
% table is a property of the estimators or an accident of one draw, which no
% averaged column can answer on its own.
    fprintf('\nStandard deviation of the per-realisation RMSE\n');
    print_region_head(regions);
    for i = 1:numel(stats)
        fprintf('%-22s', stats(i).name);
        fprintf(' %7.2f', stats(i).rmse_sd);
        fprintf('\n');
    end
end

function print_head_to_head(stats, regions, ref_name)
% Fraction of the realisations in which each estimator beats the reference one,
% region by region. A ranking read off the averaged table says nothing about how
% often it holds, whereas a figure near 50% means the two estimators are
% indistinguishable however far apart their averages happen to fall.
    ref = find(strcmp({stats.name}, ref_name), 1);
    if isempty(ref)
        return
    end
    fprintf('\nFraction of realisations with a lower RMSE than %s\n', ref_name);
    print_region_head(regions);
    for i = 1:numel(stats)
        if i == ref
            continue
        end
        fprintf('%-22s', stats(i).name);
        for j = 1:numel(regions)
            fprintf(' %6.1f%%', 100 * mean(stats(i).draw(j, :) < stats(ref).draw(j, :)));
        end
        fprintf('\n');
    end
end

function print_region_head(regions)
% The region names, once, over a single group of columns.
    fprintf('%-22s', '');
    for j = 1:numel(regions)
        fprintf(' %7s', regions(j).name(1:min(7, end)));
    end
    fprintf('\n');
end

function t = share_str(v)
% Percentage, or a dash for the two estimators that do not estimate a mean.
    if isnan(v)
        t = '-';
    else
        t = sprintf('%.1f%%', 100 * v);
    end
end

function latex_tables(stats, regions, tab_dir, n_mc)
% Writes the comparison as two LaTeX tabulars under tab_dir, so that the numbers
% of the manuscript cannot drift away from those of the figures. Bare tabulars,
% with neither float nor caption: those belong to the manuscript, which alone
% knows where each table sits and what it is being used to argue.
%
% One file each, since the two carry different claims and therefore need a
% caption apiece. The cut separates what the estimator got wrong from what it
% says it knows: kriging_rmse holds the error alone, kriging_covariance the
% declared uncertainty and the part of it contributed by the estimated mean.
    if ~exist(tab_dir, 'dir')
        mkdir(tab_dir);
    end
    n    = numel(regions);
    cols = repmat('c', 1, n);

    fid = fopen(fullfile(tab_dir, 'kriging_rmse.tex'), 'w', 'n', 'UTF-8');
    write_header(fid, n, n_mc);
    fprintf(fid, '\\begin{tabular}{|l|%s|}\\hline\n', cols);
    fprintf(fid, 'Estimator & \\multicolumn{%d}{c|}{ARMSE}\\\\\\hline\n', n);
    write_region_row(fid, regions, 1);
    for i = 1:numel(stats)
        fprintf(fid, '%s', stats(i).name);
        fprintf(fid, ' & %.2f', stats(i).rmse);
        fprintf(fid, '\\\\\\hline\n');
    end
    fprintf(fid, '\\end{tabular}\n');
    fclose(fid);

    fid = fopen(fullfile(tab_dir, 'kriging_covariance.tex'), 'w', 'n', 'UTF-8');
    write_header(fid, n, n_mc);
    fprintf(fid, '\\begin{tabular}{|l|%s|%s|}\\hline\n', cols, cols);
    fprintf(fid, ['Estimator & \\multicolumn{%d}{c|}{Confidence half-width}' ...
                  ' & \\multicolumn{%d}{c|}{$R_{\\hat{m}}/R_{\\mathrm{K}}$}\\\\\\hline\n'], n, n);
    write_region_row(fid, regions, 2);
    for i = 1:numel(stats)
        fprintf(fid, '%s', stats(i).name);
        fprintf(fid, ' & %.2f', stats(i).hw);
        for j = 1:n
            fprintf(fid, ' & %s', share_tex(stats(i).share(j)));
        end
        fprintf(fid, '\\\\\\hline\n');
    end
    fprintf(fid, '\\end{tabular}\n');
    fclose(fid);
end

function write_header(fid, n, n_mc)
% Warns off hand edits, which the next run of the script would overwrite.
    fprintf(fid, '%% Generated by krigeage_2D.m, do not edit by hand.\n');
    fprintf(fid, '%% Each metric is averaged separately over %d regions of the grid.\n', n);
    fprintf(fid, ['%% Errors are pooled over %d noise realisations. The declared covariances' ...
                  ' are the same for all of them,\n%% since no kriging covariance involves' ...
                  ' the observations.\n\n'], n_mc);
end

function write_region_row(fid, regions, groups)
% The region names, repeated once per column group.
    for j = repmat(1:numel(regions), 1, groups)
        fprintf(fid, '& %s ', regions(j).name);
    end
    fprintf(fid, '\\\\\\hline\n');
end

function t = share_tex(v)
% Same as share_str, with the percent sign escaped for LaTeX.
    if isnan(v)
        t = '--';
    else
        t = sprintf('%.1f\\%%', 100 * v);
    end
end

function save_kriging(name, x, f_true, z_hat, R_hat, x_train, z_train, ...
                      linewidth, fontsize, alpha_band, pdf_w, pdf_h, y_lim, img_dir)
% One figure per estimator, exported as vector PDF under img_dir. No title: the
% caption of the assembled figure carries the naming. The band uses the same
% pink and the same opacity as the envelope of krigeage_3D.m, so that the
% one- and two-dimensional figures read as a single family.
%
% A single band throughout, carrying the full covariance. Nesting a second one to
% expose the contribution of the estimated mean was tried and reverted: for
% ordinary kriging the two boundaries sit a couple of percent of the half-width
% apart, so they read as one thick line rather than as two. That comparison is
% carried by the generated table instead.
    pink      = [1 0.85 0.88];
    pink_edge = [0.85 0.55 0.60];                      % boundary of the band
    s         = sqrt(max(R_hat, 0));                   % guard against -0 round-off
    hi        = z_hat + 3 * s;
    lo        = z_hat - 3 * s;

    fig = figure('Color', 'w', 'Units', 'points', 'Position', [80 80 pdf_w pdf_h]);
    fill([x; flipud(x)], [hi; flipud(lo)], pink, ...
         'FaceAlpha', alpha_band, 'EdgeColor', 'none'); hold on;

    % Outlining the band would detach it from the page instead of letting it fade
    % into the white, but it adds a third pink to a figure that already carries a
    % fill and two curves. Kept here, commented, rather than deleted: it is worth
    % re-enabling if the band ever has to be read against a tinted background.
    % pink_edge is left defined for the same reason, hence unused as it stands.
    %plot(x, hi, '-', 'Color', pink_edge, 'LineWidth', 0.35 * linewidth);
    %plot(x, lo, '-', 'Color', pink_edge, 'LineWidth', 0.35 * linewidth);

    plot(x, f_true, 'b-',  'LineWidth', linewidth);
    plot(x, z_hat,  'r--', 'LineWidth', linewidth);
    plot(x_train, z_train, 'kx', 'MarkerSize', 8, 'LineWidth', linewidth);
    hold off; box on; xlim([x(1) x(end)]); ylim(y_lim);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fontsize);
    % No legend: the LaTeX caption already names every element by colour.

    % print rather than exportgraphics: the latter crops to the content, so the
    % page would vary from one panel to the next. Here every page is exactly
    % pdf_w x pdf_h, which keeps one scale, one apparent font size, and panels
    % that align when assembled.
    if ~exist(img_dir, 'dir')
        mkdir(img_dir);
    end
    set(fig, 'PaperUnits', 'points', 'PaperSize', [pdf_w pdf_h], ...
             'PaperPosition', [0 0 pdf_w pdf_h]);
    print(fig, fullfile(img_dir, [name '.pdf']), '-dpdf', '-vector');
end
