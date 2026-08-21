% FIGURES_CHAPITRE_2
%
% The four illustrations of chapter 2, exported as vector PDF under Images/ so
% that they join the figure conventions of the kriging family:
%   2.1 "Étapes du filtre de Kalman avec équations" -- the prediction/correction
%       cycle, redrawn with the manuscript's own notation,
%   2.2 "Rééchantillonnage CDF"                     -- the inverse cdf of the
%       particle weights read by systematic resampling,
%   2.3 "Appauvrissement particulaire"              -- degeneracy, resampling and
%       the impoverishment it causes,
%   2.4 "Régularisation par noyau"                  -- a weighted empirical
%       measure beside its kernel regularisation.
%
% SIZING. The convention was measured on the seven kriging PDFs: a 348 pt canvas
% included at 0.7\textwidth = 298 pt renders MATLAB's fontsize 14 at exactly the
% 12 pt of the body text. What matters there is the ratio 348/298 = 1.168, not
% the number 348, so any canvas 1.168 times its included width keeps fontsize 14
% and lands at 12 pt. Each figure therefore declares the fraction of \textwidth
% it is included at, and its page follows; changing that fraction in the LaTeX
% means changing it here too, or the apparent font size drifts. The shape is
% carried separately by an aspect ratio, so that a page recomputed for a new
% fraction reproduces exactly the same layout at a different scale.
%
% Figures 2.1, 2.3 and 2.4 are schematics rather than plots, so their coordinates
% are laid out by hand in a unit square with the axes turned off. Only figure 2.2
% carries data.

clear; close all; clc;

%% Shared parameters
textwidth = 425;                 % pt, measured on the class: 12pt A4, margin 2.5cm,
                                 % bindingoffset 1cm, hence 150 mm of text
oversize  = 348 / 298;           % canvas width over included width, the ratio the
                                 % seven kriging PDFs were built on
fontsize  = 14;                  % renders at 12.0 pt under that ratio, whatever the width
linewidth = 1;
root      = fileparts(fileparts(mfilename('fullpath')));  % repository root, one level up
img_dir   = fullfile(root, 'Images');

canvas    = @(frac) oversize * frac * textwidth;   % the page an include at frac needs

% Fraction of \textwidth each figure is included at, and the shape it is drawn
% in. These four numbers are the only ones to touch when a caption is resized.
w_21 = canvas(0.9);  h_21 = 0.605 * w_21;   % the Kalman cycle
w_22 = canvas(0.7);  h_22 = 0.805 * w_22;   % the cdf, on the kriging family page
w_23 = canvas(0.9);  h_23 = 0.468 * w_23;   % the two state axes
w_24 = canvas(0.8);  h_24 = 0.468 * w_24;   % the two panels

blue      = [0.20 0.30 0.65];    % Gaussian blobs and kernels
green     = [0.10 0.65 0.25];    % resampled particles, regularised density
red       = [0.85 0.15 0.15];    % weighted particles before resampling
grey      = [0.35 0.35 0.35];
alpha_disc = 0.80;               % opacity of the particles of figure 2.3, low enough
                                 % that overlapping copies all stay visible

%% Figure 2.1, the Kalman cycle
% Laid out on a unit square, the boxes at fixed anchors and every arrow drawn
% between two of them. Two departures from the sketch it replaces: the incoming
% measurement is written z_k rather than y_k, and the second half-step is named
% the correction rather than the update, so that the figure speaks the notation
% and the vocabulary of the chapter.
fig = new_figure(w_21, h_21);
layout_axes();
asp = h_21 / w_21;               % points per y-unit over points per x-unit

txt(0.140, 0.900, {'\textbf{Prior knowledge}', '\textbf{of state}'}, fontsize);
txt(0.400, 0.905, {'$P_{k-1|k-1}$', '$\hat{x}_{k-1|k-1}$'}, fontsize);
txt(0.628, 0.905, {'\textbf{Prediction}'}, fontsize);
txt(0.628, 0.640, {'$P_{k|k-1}$', '$\hat{x}_{k|k-1}$'}, fontsize);
txt(0.628, 0.355, {'\textbf{Correction}'}, fontsize);
txt(0.400, 0.355, {'$P_{k|k}$', '$\hat{x}_{k|k}$'}, fontsize);
txt(0.400, 0.635, {'\textbf{Next time step}', '$k\leftarrow k+1$'}, fontsize);
txt(0.880, 0.412, {'$z_k$'}, fontsize);          % above, so that the word below it
txt(0.880, 0.355, {'\textbf{Observation}'}, fontsize);   % sits on the row it enters
txt(0.400, 0.115, {'\textbf{Output estimate}', '\textbf{of state}'}, fontsize);

arrow([0.275 0.338], [0.905 0.905], linewidth, asp);   % prior       -> P_{k-1|k-1}
arrow([0.470 0.520], [0.905 0.905], linewidth, asp);   % P_{k-1|k-1} -> prediction
arrow([0.628 0.628], [0.860 0.720], linewidth, asp);   % prediction  -> P_{k|k-1}
arrow([0.628 0.628], [0.560 0.450], linewidth, asp);   % P_{k|k-1}   -> correction
arrow([0.762 0.715], [0.355 0.355], linewidth, asp);   % observation -> correction
arrow([0.520 0.470], [0.355 0.355], linewidth, asp);   % correction  -> P_{k|k}
arrow([0.400 0.400], [0.445 0.550], linewidth, asp);   % P_{k|k}     -> next step
arrow([0.400 0.400], [0.725 0.825], linewidth, asp);   % next step   -> P_{k-1|k-1}
arrow([0.400 0.400], [0.270 0.185], linewidth, asp);   % P_{k|k}     -> output

% The three densities last: each opens an axes of its own, which would otherwise
% capture the text and the arrows drawn after it.
gaussian_blob(0.045, 0.540, 0.135, 0.215, -0.55, blue, linewidth);
gaussian_blob(0.840, 0.055, 0.125, 0.200,  0.10, blue, linewidth);
gaussian_blob(0.125, 0.010, 0.125, 0.200, -0.30, blue, linewidth);

save_pdf(fig, 'Étapes du filtre de Kalman avec équations', w_21, h_21, img_dir);

%% Figure 2.2, the cdf read by Kitagawa resampling
% Notation of alg:Kitagawa throughout: the particles are indexed by i and the
% resampled ones by j, cdf_i is the cumulative weight of the first i particles,
% and the reading points are the Kitagawa samples u_j = u_1 + (j-1)/N with u_1
% drawn once in [0, 1/N]. The staircase maps u_j to the parent index i^{(j)}, so
% the width of the step of particle i is exactly w_k^{(i)}: a particle whose step
% is wider than the spacing 1/N is certain to be drawn, and drawn as many times
% as its step holds reading points.
w      = [0.14 0.06 0.05 0.14 0.05 0.10 0.04 0.12 0.16 0.14];
w      = w / sum(w);
N      = numel(w);
cdf    = [0, cumsum(w)];         % cdf_0 = 0 and cdf_i, one step per particle
i_star = 4;                      % the particle the brace measures, one of the widest
u_1    = 0.03;                   % the single draw in [0, 1/N] the whole grid rests on
u      = u_1 + (0:N-1) / N;      % the Kitagawa samples u_j
parent = arrayfun(@(uj) find(uj < cdf(2:end), 1), u);   % the parent indices i^{(j)}

y_min  = -0.6;
y_max  = N + 0.4;
ax_pos = [0.150 0.150 0.820 0.820];
marker = 5;                      % points. Half of it puts the tip on the staircase,
gap    = 2;                      % and this much more leaves the curve unbroken
head   = (0.5 * marker + gap) * (y_max - y_min) / (ax_pos(4) * h_22);
brace_inset = 0.010;             % so that the arm of the brace stops short of the riser
                                 % the step it measures ends on

fig = new_figure(w_22, h_22);
axes('Units', 'normalized', 'Position', ax_pos); hold on;

% The reading points first, so that the staircase and the axis line sit above them.
for j = 1:N
    plot([u(j) u(j)], [0, parent(j) - head], '-', 'Color', green, 'LineWidth', linewidth);
    plot(u(j), parent(j) - head, '^', 'Color', green, 'MarkerFaceColor', green, ...
         'MarkerSize', marker);
end

stairs(cdf, [1:N, N], '-', 'Color', blue, 'LineWidth', linewidth);
plot([0 1], [0 0], '-', 'Color', red, 'LineWidth', linewidth);

% The width of one step, named, and its height carried back to the axis it is
% read on. Together they say what the staircase is: a weight in, a parent out.
curly_brace(cdf(i_star) + brace_inset, cdf(i_star+1) - brace_inset, ...
            i_star + 0.28, 0.50, grey, linewidth);
text(mean(cdf(i_star:i_star+1)), i_star + 1.30, '$w_k^{(i)}$', 'Interpreter', 'latex', ...
     'FontSize', fontsize, 'HorizontalAlignment', 'center');
plot([cdf(i_star) 0.035], [i_star i_star], ':', 'Color', grey, 'LineWidth', linewidth);
plot(0.035, i_star, '<', 'Color', grey, 'MarkerFaceColor', grey, 'MarkerSize', 4);
text(0.055, i_star + 0.45, '$i^{(j)}$', 'Interpreter', 'latex', 'FontSize', fontsize);

hold off; box on;
xlim([0 1]); ylim([y_min y_max]);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fontsize, ...
         'XTick', 0:0.2:1, 'YTick', 2:2:N);
% \mathrm rather than the manuscript's \operatorname, which MATLAB's interpreter
% does not know. Both set the name upright, so the figure and the text agree.
xlabel('$\mathrm{cdf}$', 'Interpreter', 'latex', 'FontSize', fontsize);
ylabel('parent index', 'Interpreter', 'latex', 'FontSize', fontsize);

save_pdf(fig, 'Rééchantillonnage CDF', w_22, h_22, img_dir);

%% Figure 2.3, degeneracy and impoverishment
% Two state axes carrying the same particle set, before and after resampling.
% Areas are proportional to weights, so the upper line shows the mass gathered on
% a handful of particles, and the lower one shows what survives once those few
% have each been drawn several times. The copies are spread on a small circle
% rather than superimposed, so that they can be counted: the support has lost
% every location but two, and resampling has created none.
before = struct('x', {0.10 0.16 0.21 0.27 0.30 0.33 0.36 0.41 0.44 0.55 0.62 0.68 0.75 0.80 0.84}, ...
                'r', {0.010 0.018 0.012 0.014 0.022 0.026 0.020 0.016 0.055 0.075 0.017 0.021 0.011 0.008 0.006});
kept    = [9 10];                % the two particles resampling concentrates on
centres = [0.440 0.550];         % where their copies land, at the same state values
r_after = 0.030;                 % one radius for all of them: resampling resets every
                                 % weight to 1/N, so what a heavy parent buys is not a
                                 % bigger particle but more copies of itself
copies  = [3 5];

fig = new_figure(w_23, h_23);
layout_axes(); axis equal; xlim([0 1]); ylim([0 h_23 / w_23]);
y_hi = 0.355; y_lo = 0.075;

plot([0.03 0.90], [y_hi y_hi], 'k-', 'LineWidth', 0.6 * linewidth);
plot([0.03 0.90], [y_lo y_lo], 'k-', 'LineWidth', 0.6 * linewidth);
txt(0.955, y_hi, {'State'}, fontsize);
txt(0.955, y_lo, {'State'}, fontsize);

% Every disc is collected before anything is drawn, so that the fills can go down
% first and the outlines all together on top. Painted one by one, a late disc
% would bury the outline of an early one, and the copies of the lower axis are
% precisely the ones that need to stay countable.
discs = struct('x', {}, 'y', {}, 'r', {}, 'c', {});
for i = 1:numel(before)
    face = red;
    if any(i == kept)
        face = green;
    end
    discs(end+1) = struct('x', before(i).x, 'y', y_hi, 'r', before(i).r, 'c', face); %#ok<SAGROW>
end
for c = 1:numel(centres)
    for i = 1:copies(c)
        a = 2*pi*(i-1) / copies(c);
        discs(end+1) = struct('x', centres(c) + 0.42*r_after*cos(a), ...
                              'y', y_lo + 0.42*r_after*sin(a), ...
                              'r', r_after, 'c', green); %#ok<SAGROW>
    end
end
for d = discs
    disc_fill(d, alpha_disc);
end
for d = discs
    disc_edge(d, linewidth);
end

txt(0.175, y_hi - 0.090, {'\textbf{Degeneracy}'}, fontsize);
txt(0.510, y_hi - 0.125, {'\textbf{Resampling}'}, fontsize);
txt(0.175, y_lo - 0.050, {'\textbf{Impoverishment}'}, fontsize);
arrow([0.290 0.290], [0.215 0.130], linewidth, 1);
arrow([0.800 0.800], [0.215 0.130], linewidth, 1);

save_pdf(fig, 'Appauvrissement particulaire', w_23, h_23, img_dir);

%% Figure 2.4, a weighted empirical measure and its regularisation
% The same five particles on both panels, at the same abscissae, so that the two
% approximations can be read against each other. On the left each particle
% carries its weight as a spike, on the right it carries a kernel of that weight,
% and the sum of the kernels is the density the regularised filter resamples
% from. The bandwidth is deliberately generous, to show the gaps being bridged.
x_p     = [0.06 0.15 0.24 0.68 0.82];   % two clusters, the hole five times the spacing
w_p     = [0.16 0.28 0.20 0.24 0.12];   % inside a cluster, so that it reads as a hole
h_band  = 0.12;                          % wide enough that the kernels bridge the hole
                                         % rather than leaving it empty, which is the
                                         % whole difference between the two panels
x_grid  = linspace(-0.30, 1.20, 500);
kernels = w_p .* exp(-0.5 * ((x_grid' - x_p) / h_band).^2) / (h_band * sqrt(2*pi));
density = sum(kernels, 2);
y_max   = 1.15 * max(density);
x_lim   = [-0.22 1.12];
y_lim   = [-0.12 * y_max, y_max];   % just enough under the baseline for the samples
                                    % label, the other two having moved to the corridor
spikes  = w_p / max(w_p) * 0.88 * y_max;   % the tallest weight fills the panel

% The two panels are pushed apart to open a corridor between them, in which the
% three names common to both are stacked once instead of twice.
pos_l = [0.010 0.080 0.390 0.890];
pos_r = [0.600 0.080 0.390 0.890];

fig = new_figure(w_24, h_24);

axes('Units', 'normalized', 'Position', pos_l); hold on;
for i = 1:numel(x_p)
    plot([x_p(i) x_p(i)], [0, spikes(i)], '-', 'Color', blue, 'LineWidth', linewidth);
end
plot(x_p, zeros(size(x_p)), 'o', 'Color', blue, 'MarkerFaceColor', 'w', ...
     'MarkerSize', 4, 'LineWidth', linewidth);
plot(x_lim, [0 0], 'k-', 'LineWidth', 0.6 * linewidth);
hold off; axis off; xlim(x_lim); ylim(y_lim);

axes('Units', 'normalized', 'Position', pos_r); hold on;
plot(x_grid, kernels, '-', 'Color', blue, 'LineWidth', linewidth);
plot(x_grid, density, '-', 'Color', green, 'LineWidth', 1.4 * linewidth);
plot(x_p, zeros(size(x_p)), 'o', 'Color', blue, 'MarkerFaceColor', 'w', ...
     'MarkerSize', 4, 'LineWidth', linewidth);
plot(x_lim, [0 0], 'k-', 'LineWidth', 0.6 * linewidth);

hold off; axis off; xlim(x_lim); ylim(y_lim);

% The shared column, drawn in an axes spanning the whole page so that a leader
% may start between the panels and end inside either of them. Every leader is
% horizontal, so each label sits at the height of what it names and the four of
% them space themselves out. "samples" is the exception twice over: it stands at
% the height of the two baselines, and it carries no leader, which at that height
% would be taken for the axis itself.
axes('Units', 'normalized', 'Position', [0 0 1 1]); hold on;
axis off; xlim([0 1]); ylim([0 1]);
x_col = mean([pos_l(1) + pos_l(3), pos_r(1)]);
lead  = 0.078;                   % where a leader starts, clear of the widest label

p_base = to_fig(pos_l, x_lim, y_lim, 0, 0);
text(x_col, p_base(2), 'samples', 'Interpreter', 'latex', ...
     'FontSize', fontsize, 'HorizontalAlignment', 'center');

% The kernels are named on the left flank of the tallest one, high enough up that
% the flank is unmistakably its own and not the foot of its neighbour.
x_kern = x_p(2) - 0.9 * h_band;
p_kern = to_fig(pos_r, x_lim, y_lim, x_kern, interp1(x_grid, kernels(:, 2), x_kern));
label_with_leader(x_col, p_kern, lead, {'kernels'}, fontsize, grey, linewidth);

% The weights are named three quarters up the tallest spike of the right cluster,
% where the leader clears the shorter one beside it.
p_wght = to_fig(pos_l, x_lim, y_lim, x_p(4), 0.75 * spikes(4));
label_with_leader(x_col, p_wght, -lead, {'weights'}, fontsize, grey, linewidth);

% The sum, near the top of it, and on the flank facing its label. The level is
% above anything the right hump reaches, so the last crossing is the descending
% side of the left one -- the first thing a leader coming from the right meets,
% which is what keeps it from cutting across the curve to reach the other flank.
y_dens = 0.88 * max(density);
x_dens = x_grid(find(density >= y_dens, 1, 'last'));
p_dens = to_fig(pos_r, x_lim, y_lim, x_dens, y_dens);
label_with_leader(x_col+0.38, p_dens, -lead, {'density', 'estimate'}, fontsize, grey, linewidth);
hold off;

save_pdf(fig, 'Régularisation par noyau', w_24, h_24, img_dir);

fprintf('Four PDFs written to %s\n', img_dir);

%% Local functions
function label_with_leader(x_col, target, lead, lines, fontsize, colour, linewidth)
% One entry of the shared column: the words centred in the corridor at the height
% of what they name, and a horizontal leader running out to it. lead is signed,
% positive to reach into the right panel and negative into the left one.
    text(x_col, target(2), lines, 'Interpreter', 'latex', 'FontSize', fontsize, ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    plot([x_col + lead, target(1)], [target(2), target(2)], '-', ...
         'Color', colour, 'LineWidth', 0.6 * linewidth);
end

function p = to_fig(pos, x_lim, y_lim, x, y)
% Where a point of an axes lands on the page, in the normalised coordinates of a
% full-page overlay. Lets a leader drawn between two panels end on a curve inside
% one of them, which no axes of its own could reach.
    p = [pos(1) + pos(3) * (x - x_lim(1)) / (x_lim(2) - x_lim(1)), ...
         pos(2) + pos(4) * (y - y_lim(1)) / (y_lim(2) - y_lim(1))];
end

function fig = new_figure(pdf_w, pdf_h)
% A frozen page in points, as in the kriging scripts: print writes exactly this
% size, where exportgraphics would crop to the content and let the page drift.
    fig = figure('Color', 'w', 'Units', 'points', 'Position', [80 80 pdf_w pdf_h]);
end

function ax = layout_axes()
% A single invisible axes spanning the whole page, in which the schematics are
% laid out by hand. Coordinates are then read directly as fractions of the width.
    ax = axes('Units', 'normalized', 'Position', [0 0 1 1]); hold on;
    axis off; xlim([0 1]); ylim([0 1]);
end

function txt(x, y, lines, fontsize)
% A block of centred lines. Bold goes through \textbf rather than FontWeight,
% which the latex interpreter ignores.
    text(x, y, lines, 'Interpreter', 'latex', 'FontSize', fontsize, ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end

function arrow(x, y, linewidth, aspect)
% A straight arrow in the coordinates of the current axes, drawn as a shaft plus
% a filled head rather than as an annotation, which would live in figure
% coordinates instead. aspect is the height of one y-unit over the width of one
% x-unit, in points, so that the head stays isotropic on an axes whose two
% directions do not share a scale.
    plot(x, y, 'k-', 'LineWidth', linewidth);
    dx = x(2) - x(1);
    dy = (y(2) - y(1)) * aspect;
    d  = hypot(dx, dy);
    dx = dx / d; dy = dy / d;
    nx = -dy;    ny = dx;
    L  = 0.020; W = 0.008;                             % head length and half-width
    bx = x(2) - L * dx;         by = y(2) - L * dy / aspect;
    patch('XData', [x(2), bx + W*nx, bx - W*nx], ...
          'YData', [y(2), by + W*ny/aspect, by - W*ny/aspect], ...
          'FaceColor', 'k', 'EdgeColor', 'none');
end

function disc_fill(d, alpha)
% The body of one particle, its area carrying its weight. Translucent, so that a
% disc drawn over another still lets it be seen rather than erasing it.
    t = linspace(0, 2*pi, 120);
    patch('XData', d.x + d.r * cos(t), 'YData', d.y + d.r * sin(t), ...
          'FaceColor', d.c, 'FaceAlpha', alpha, 'EdgeColor', 'none');
end

function disc_edge(d, linewidth)
% Its outline, drawn in a second pass over every fill, so that each circle is
% closed and countable however many others overlap it.
    t = linspace(0, 2*pi, 120);
    plot(d.x + d.r * cos(t), d.y + d.r * sin(t), 'k-', 'LineWidth', 0.4 * linewidth);
end

function gaussian_blob(left, bottom, w, h, rho, blue, linewidth)
% A small bivariate Gaussian density in its own axes, with two arrows for axes.
% Filled contours rather than an image, so that the PDF stays vector. The caller's
% axes is restored on the way out, since creating an axes also makes it current
% and would silently capture everything drawn afterwards.
    parent = gcf;
    prev   = get(parent, 'CurrentAxes');
    ax     = axes('Parent', parent, 'Units', 'normalized', 'Position', [left bottom w h]);
    hold on;
    t = linspace(-3, 3, 60);
    [X, Y] = meshgrid(t, t);
    Q = (X.^2 - 2*rho*X.*Y + Y.^2) / (1 - rho^2);
    contourf(X, Y, exp(-0.5 * Q), 10, 'LineStyle', 'none');
    colormap(ax, white_to(blue));
    plot([-3.4 -3.4], [-3.4 3.2], 'k-', 'LineWidth', linewidth);
    plot([-3.4 3.2], [-3.4 -3.4], 'k-', 'LineWidth', linewidth);
    hold off; axis off equal; xlim([-3.7 3.7]); ylim([-3.7 3.7]);
    if ~isempty(prev)
        set(parent, 'CurrentAxes', prev);
    end
end

function map = white_to(colour)
% A colormap running from the page white to the given colour, so that the tail of
% the density fades into the page instead of ending on a visible square.
    s   = linspace(0, 1, 64)';
    map = 1 - s * (1 - colour);
end

function curly_brace(x1, x2, y, depth, colour, linewidth)
% An overbrace spanning [x1, x2], its two arms at y and its tip at y + depth,
% MATLAB having no such shape. Built from four quarter ellipses joined by a
% straight run at mid-height, the two radii taken separately because the axes it
% is drawn on has one unit of weight on one side and one unit of index on the
% other.
    rx = (x2 - x1) / 4;
    ry = depth / 2;
    xm = (x1 + x2) / 2;
    t  = linspace(0, pi/2, 24);

    arm_x  = x1 + rx - rx * cos(t);      arm_y  = y + ry * sin(t);
    rise_x = xm - rx + rx * sin(t);      rise_y = y + depth - ry * cos(t);

    half_x = [arm_x, rise_x];
    half_y = [arm_y, rise_y];
    plot([half_x, 2*xm - fliplr(half_x)], [half_y, fliplr(half_y)], ...
         '-', 'Color', colour, 'LineWidth', linewidth);
end

function save_pdf(fig, name, pdf_w, pdf_h, img_dir)
% Vector PDF on a frozen page, the export path shared with the kriging scripts.
% print truncates its target before writing, so never preview the same file in
% the same command: a concurrent read destroys it and breaks the LaTeX build.
    if ~exist(img_dir, 'dir')
        mkdir(img_dir);
    end
    set(fig, 'PaperUnits', 'points', 'PaperSize', [pdf_w pdf_h], ...
             'PaperPosition', [0 0 pdf_w pdf_h]);
    print(fig, fullfile(img_dir, [name '.pdf']), '-dpdf', '-vector');
end
