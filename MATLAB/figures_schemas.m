% FIGURES_SCHEMAS
%
% The hand-laid diagrams of the manuscript, exported as vector PDF under Images/.
% Unlike the per-chapter scripts, nothing here is read from a campaign: every
% coordinate is placed by hand so that the drawing says one thing and says it once.
%
%   "Construction de la fenêtre AK"     the AK window in two regimes, sec:AK_PF
%   "Trois estimateurs ponctuels"       mean, median and mode on one law, sec:bayesian
%   "Schéma de l'application krigeage"  the data flow of ch. 4, sec:kriging_navigation
%   "Schéma de l'application slam"      the same flow without a prior survey, sec:SLAM
%
% SIZING. Same convention as figures_chapitre_2.m, measured on the kriging PDFs: a
% canvas 348/298 times its included width renders MATLAB's fontsize 14 at the 12 pt
% of the body text. Each figure declares the fraction of \textwidth it is included
% at, and changing that fraction in the LaTeX means changing it here too.
%
% SLIDE VERSIONS (2026-11-22). Called by figures_soutenance.m with the diapo structure
% in its workspace, the script also draws the slide version of the AK window, at the
% size Latex/Soutenance/main.tex includes it and with text at diapo.fs. The manuscript
% figures then go to diapo.controle, where figures_soutenance compares them with Images/.

clearvars -except diapo; close all; clc;

%% Shared parameters
textwidth = 425;                 % pt, measured on the class
oversize  = 348 / 298;
fontsize  = 14;
linewidth = 1;
root      = fileparts(fileparts(mfilename('fullpath')));
img_dir   = fullfile(root, 'Images');
en_diapo  = exist('diapo', 'var') && isstruct(diapo);
if en_diapo
    img_dir = diapo.controle;
end

canvas = @(frac) oversize * frac * textwidth;

red   = [0.85 0.15 0.15];        % the window, as in the article it replaces
blue  = [0.20 0.30 0.65];        % the particles
grey  = [0.45 0.45 0.45];        % the confidence ellipse

%% Fenêtrage du krigeage adaptatif
% Two panels on the same survey and the same centre, differing only by the spread of
% the particle set. On the left the cloud is wide and r_cov governs; on the right it
% has concentrated, r_cov would retain almost nothing, and r_min takes over. The two
% radii are drawn to scale against a survey of unit spacing: everything the reader is
% asked to compare is a length on the page.
%
% THE ELLIPSE IS THE 99 % REGION of the predicted position, not a 3-sigma contour: the
% manuscript sizes the window with chi2_{d_pos}(99 %), and alpha_AK dilates the ball
% that circumscribes it. Drawing the ellipse and the ball in the same picture is the
% whole point, since the text argues for the ball against the ellipse.
w_ak = canvas(0.9);  h_ak = 0.52 * w_ak;

chi2   = 9.21;                   % chi2_2(99 %)
alpha  = 2;                      % alpha_AK, the value used in the chapter
pas    = 1;                      % survey spacing, the unit of the drawing
jitter = 0.30;                   % as in sec:scenario

rng(7);                          % one draw, kept for reproducibility

% Left panel: spread cloud. lam are the eigenvalues of P_pos, tilted by theta.
gauche = struct('lam', [0.244 0.080], 'theta', 25, 'n_min', 12, 'r_min', 0);
% Right panel: the same, concentrated by a factor of eight in each direction.
droite = struct('lam', [0.030 0.012], 'theta', 25, 'n_min', 12, 'r_min', 1.95);

fig = new_figure(w_ak, h_ak);

panneau(0.030, [], gauche, chi2, alpha, pas, jitter, fontsize, linewidth, ...
        red, blue, grey);
panneau(0.520, [], droite, chi2, alpha, pas, jitter, fontsize, linewidth, ...
        red, blue, grey);

save_pdf(fig, 'Construction de la fenêtre AK', w_ak, h_ak, img_dir);

% THE SLIDE VERSION, under the block of the AK slide, at the 0.7\textwidth main.tex
% includes it at, and 95 pt high, what the block leaves above the bottom of the page.
% The same two panels from the same draw, hence the second rng(7). The squares are half
% the size they have in the manuscript, so the crosses and the centre shrink by 0.6 to
% keep the survey readable; the particles keep their dots, already the smallest mark.
if en_diapo
    w_d = 0.7 * diapo.textwidth;  h_d = 95;
    rng(7);
    fig = new_figure(w_d, h_d);
    panneau(0.030, 0.6, gauche, chi2, alpha, pas, jitter, diapo.fs, linewidth, ...
            red, blue, grey);
    panneau(0.520, 0.6, droite, chi2, alpha, pas, jitter, diapo.fs, linewidth, ...
            red, blue, grey);
    save_pdf(fig, 'Construction de la fenêtre AK', w_d, h_d, diapo.img);
end

%% Trois estimateurs ponctuels
% One posterior and three point estimates, drawn to carry the fact that opens
% sec:analytical_estimation: on a multimodal law the MMSE estimate falls where the
% density is almost zero. The likelihood has two modes because a scalar field
% measurement does not determine a position, which is the situation of the navigation
% problem itself, and the prior then breaks the tie, so that the ML and the MAP
% estimates do not coincide either. The three values are computed from the curves.
%
% The prior and the posterior share the scale their densities actually have, so that the
% contraction the measurement brings is read off the drawing. Only the likelihood is
% rescaled, since read as a function of the state it is not a density and carries no
% scale of its own, and it is held just below the posterior to leave it the tallest.
w_es = canvas(0.75);  h_es = 0.56 * w_es;

gauss = @(x, m, s) exp(-0.5*((x - m)/s).^2) / (s*sqrt(2*pi));
u     = linspace(-4.2, 4.2, 4000);
% The prior is broad and slightly off-centre, which is enough to flip which of the two
% likelihood modes carries the posterior: the ML estimate stays on the taller mode of
% the likelihood, the MAP moves to the other, and the MMSE lands between the two.
prior = gauss(u, 0.90, 2.20);
lik   = 0.52*gauss(u, -1.80, 0.45) + 0.48*gauss(u, 1.80, 0.55);
post  = prior .* lik;  post = post / trapz(u, post);
s     = max(post);                   % the reference height of the drawing

[~, i_ml]  = max(lik);
[~, i_map] = max(post);
u_mmse     = trapz(u, u .* post);

fig = new_figure(w_es, h_es);
axes('Units', 'normalized', 'Position', [0.04 0.10 0.93 0.78]); hold on;
xlim([-4.2 4.2]); ylim([0 1.30*s]); axis off;

plot(u, prior,                ':',  'Color', grey, 'LineWidth', linewidth);
plot(u, 0.85*s*lik/max(lik),  '--', 'Color', grey, 'LineWidth', linewidth);
plot(u, post,                 '-',  'Color', blue, 'LineWidth', 1.6*linewidth);
plot([-4.2 4.2], [0 0], 'k-', 'LineWidth', 0.7*linewidth);

estime(u(i_ml),  post, u, s, '$\hat{\psi}_{\mathrm{ML}}$',   fontsize, linewidth, red);
estime(u(i_map), post, u, s, '$\hat{\psi}_{\mathrm{MAP}}$',  fontsize, linewidth, red);
estime(u_mmse,   post, u, s, '$\hat{\psi}_{\mathrm{MMSE}}$', fontsize, linewidth, red);

text(3.9, 0.30*s, 'prior',      'Interpreter', 'latex', 'FontSize', fontsize, ...
     'Color', grey, 'HorizontalAlignment', 'right');
text(-3.9, 0.62*s, 'likelihood', 'Interpreter', 'latex', 'FontSize', fontsize, ...
     'Color', grey);
text(-3.9, 0.95*s, 'posterior',  'Interpreter', 'latex', 'FontSize', fontsize, ...
     'Color', blue);
hold off;
save_pdf(fig, 'Trois estimateurs ponctuels', w_es, h_es, img_dir);

%% Schéma de l'application
% The whole method as a data flow, redrawn after the figure of the FUSION 2024 article
% from the photographs it was built with, under Images/sources. The composition is that
% figure's: the vehicle on the left carrying the two sensors, the field read by the
% magnetometer, the two measurements merging into the estimation block, the filter above
% and the map below exchanging a selection and a reconstruction, and on the right the
% same aircraft seen through the posterior it ends up with.
%
% Three things are adapted rather than copied. The article's figure is two and a half
% times as wide as it is tall, which body text at 12 pt cannot hold here, since a label
% takes twice the room it took on a page twice as wide. The drawing is therefore less
% flat, the estimation block wider than its share, and the label of the inertial unit
% shortened. The palette is the manuscript's rather than the poster's navy on lavender.
%
% In sec:kriging_navigation a survey flown beforehand holds the map, and the filter only
% chooses which of its points to use. In sec:SLAM nothing is known beforehand, and the
% map is estimated from the observations collected along the way, which changes the two
% labels of the lower block and nothing else.
%
% Each photograph is reduced to the pixels its printed size actually uses, at the 300 dpi
% a photograph asks for on paper. The reduction is not what the exported file stores, as
% print pins an embedded image to 576 dpi of the canvas and ignores -r, but a file
% resampled up from a reduced source still weighs close to half of what the whole one
% gives, and nothing of the photograph is lost at the printed size.
w_ap  = canvas(1.00);
H     = 0.590;                       % the height of the page, in units of its width
h_ap  = H * w_ap;
src   = fullfile(root, 'Images', 'sources');   % not img_dir, which a slide run redirects
pouce = textwidth / 72;              % the printed width of the figure, in inches
cible = @(w) round(300 * w * pouce);
lw_fl = 1.6;                         % the arrows carry the drawing, so they are heavier
r_arr = 0.022;                       % the radius the article rounds its bends at

drone = charge(fullfile(src, 'Drone.png'),                     [0 0 0 0.30861],     cible(0.118));
globe = charge(fullfile(src, 'Champ magnétique terrestre.png'),[0 0 0 0],           cible(0.095));
magne = charge(fullfile(src, 'Magnétomètre.png'),              [0 0 0 0],           cible(0.098));
inert = charge(fullfile(src, 'Centrale inertielle.jpg'),       [0.14235 0.05488 0.02708 0.18913], cible(0.097));
cycle = charge(fullfile(src, 'Cycle du filtre.png'),           [0 0 0 0],           cible(0.175));
surfa = charge(fullfile(src, 'Surface krigée.png'),            [0.07804 0 0.06604 0], cible(0.195));
estim = charge(fullfile(src, 'Drone flou.png'),                [0.11371 0 0.10207 0.30861], cible(0.145));

for ch = ["krigeage", "slam"]
    fig = new_figure(w_ap, h_ap);
    axes('Units', 'normalized', 'Position', [0 0 1 1]); hold on;
    axis off; xlim([0 1]); ylim([0 H]);
    yy = @(d) H - d;                 % a depth below the top edge, as an ordinate

    % The field, and the magnetometer that reads it. The two names share the top line,
    % which is why the globe sits further right than the block it stands above.
    txt_s(0.445, yy(0.018), {'Earth''s magnetic field'}, fontsize, blue);
    txt_s(0.170, yy(0.018), {'Magnetometer'},            fontsize, blue);
    vignette(0.445, yy(0.100), 0.095, globe);
    vignette(0.170, yy(0.100), 0.098, magne);
    chemin([0.390 0.232], yy([0.100 0.100]), r_arr, lw_fl, blue, true);

    % The vehicle and what it carries, each sensor reached by a bend of its own.
    vignette(0.068, yy(0.255), 0.118, drone);
    chemin([0.072 0.072 0.130], yy([0.227 0.100 0.100]), r_arr, lw_fl, blue, true);
    chemin([0.072 0.072 0.130], yy([0.283 0.350 0.350]), r_arr, lw_fl, blue, true);
    vignette(0.170, yy(0.350), 0.097, inert);
    txt_s(0.170, yy(0.420), {'Inertial unit'}, fontsize, blue);

    % The two measurements, joined into the single stream the filter consumes.
    chemin([0.165 0.165 0.285], yy([0.145 0.245 0.245]), r_arr, lw_fl, blue, true);
    chemin([0.165 0.165 0.235], yy([0.307 0.245 0.245]), r_arr, lw_fl, blue, false);

    % The estimation block, dashed because it is the contribution and not the hardware.
    rectangle('Position', [0.295 yy(0.568) 0.490 0.416], 'EdgeColor', blue, ...
              'LineStyle', '--', 'LineWidth', linewidth, 'Curvature', 0.09);
    txt_s(0.540, yy(0.185), {'Particle filter'}, fontsize, blue);
    vignette(0.540, yy(0.258), 0.175, cycle);

    if ch == "krigeage"
        bas  = {'Kriging'};
        desc = {'Data', 'selection'};
    else
        bas  = {'Map estimation'};
        desc = {'Estimate and', 'observation'};
    end
    chemin([0.495 0.495], yy([0.320 0.388]), r_arr, lw_fl, blue, true);
    chemin([0.570 0.570], yy([0.388 0.320]), r_arr, lw_fl, blue, true);
    txt_r(0.483, yy(0.3555), desc, fontsize, blue);
    txt_l(0.582, yy(0.3555), {'Reconstructed', 'map'}, fontsize, blue);

    vignette(0.540, yy(0.455), 0.195, surfa);
    txt_s(0.540, yy(0.538), bas, fontsize, blue);

    % The survey flown beforehand, which only sec:kriging_navigation has. Its absence is
    % what separates the two drawings beyond the two labels of the lower block.
    if ch == "krigeage"
        semis(0.130, yy(0.489), 0.130, 0.078, 6, 4, blue);
        txt_s(0.130, yy(0.556), {'Training data'}, fontsize, blue);
        chemin([0.206 0.293], yy([0.489 0.489]), r_arr, lw_fl, blue, true);
    end

    % The output, as the same aircraft seen through the posterior: blurred, inside the
    % ellipse the filter reports, around the point it would hand over if asked for one.
    chemin([0.797 0.845], yy([0.255 0.255]), r_arr, lw_fl, blue, true);
    vignette(0.918, yy(0.255), 0.145, estim);
    phi = linspace(0, 2*pi, 181);
    th  = -28 * pi/180;
    plot(0.915 + 0.059*cos(phi)*cos(th) - 0.0185*sin(phi)*sin(th), ...
         yy(0.255) + 0.059*cos(phi)*sin(th) + 0.0185*sin(phi)*cos(th), ...
         'k-', 'LineWidth', linewidth);
    plot(0.915 + [-0.009 0.009], yy(0.255) * [1 1], 'k-', 'LineWidth', 1.6*linewidth);
    plot(0.915 * [1 1], yy(0.255) + [-0.009 0.009], 'k-', 'LineWidth', 1.6*linewidth);
    txt_s(0.918, yy(0.345), {'Navigation', 'estimation'}, fontsize, blue);

    hold off;
    save_pdf(fig, char("Schéma de l'application " + ch), w_ap, h_ap, img_dir);
end

fprintf('Written to %s\n', img_dir);

%% Local functions
function panneau(left, ech, par, chi2, alpha, pas, jitter, fs, lw, red, blue, grey)
% One panel of the AK window figure, in an axes of its own occupying the left half or
% the right half of the page. Data units are survey spacings throughout. ech scales the
% crosses and the centre, [] leaving them at the manuscript's size.
    if isempty(ech), ech = 1; end
    ax = axes('Units', 'normalized', 'Position', [left 0.02 0.45 0.96]);
    hold on; axis equal off;

    demi = 3.6;                             % half-width of the drawn region
    xlim([-demi demi]); ylim([-demi*0.98 demi*0.98]);

    % The survey: a regular grid shaken by up to jitter * pas, as in sec:scenario.
    [GX, GY] = meshgrid(-4:pas:4, -4:pas:4);
    GX = GX + jitter * pas * (2*rand(size(GX)) - 1);
    GY = GY + jitter * pas * (2*rand(size(GY)) - 1);
    P = [GX(:) GY(:)]';

    % The predicted covariance, its 99 % ellipse, and the radii it induces.
    R     = [cosd(par.theta) -sind(par.theta); sind(par.theta) cosd(par.theta)];
    r_cov = alpha * sqrt(chi2 * max(par.lam));
    r     = max(r_cov, par.r_min);

    t  = linspace(0, 2*pi, 200);
    el = R * [sqrt(chi2*par.lam(1))*cos(t); sqrt(chi2*par.lam(2))*sin(t)];

    % The particles, drawn from the same covariance so that the ellipse is theirs.
    X = R * (sqrt(par.lam(:)) .* randn(2, 160));

    % Retained points first as filled discs, so the crosses read on top of them.
    garde = vecnorm(P) <= r;
    plot(P(1, ~garde), P(2, ~garde), 'x', 'Color', [0.55 0.55 0.55], ...
         'MarkerSize', 5*ech, 'LineWidth', 0.8*lw);
    plot(P(1, garde), P(2, garde), 'x', 'Color', 'k', ...
         'MarkerSize', 6*ech, 'LineWidth', 1.3*lw);

    plot(X(1, :), X(2, :), '.', 'Color', blue, 'MarkerSize', 4);
    plot(el(1, :), el(2, :), '-', 'Color', grey, 'LineWidth', lw);

    % The window itself, and the radius it was not given when r_min takes over.
    if par.r_min > r_cov
        plot(r_cov*cos(t), r_cov*sin(t), ':', 'Color', red, 'LineWidth', lw);
    end
    plot(r*cos(t), r*sin(t), '-', 'Color', red, 'LineWidth', 1.4*lw);
    plot(0, 0, '+', 'Color', red, 'MarkerSize', 8*ech, 'LineWidth', 1.2*lw);

    % Each radius drawn as a segment, so that its name sits on a length rather than
    % beside a curve. When r_min governs, both are named: the reader has to see that
    % the covariance radius was too small, not guess it from a dotted circle. The two
    % directions avoid 25 degrees, which is the major axis of the ellipse.
    rayon(r_cov, 125, '$r_{\mathrm{cov},k}$', fs, lw, red);
    if par.r_min > r_cov
        rayon(par.r_min, 295, '$r_{\min,k}$', fs, lw, red);
    end
    hold off;
end

function estime(u0, courbe, u, h, nom, fs, lw, col)
% One point estimate: a stem up to a common height so that the three names align, a dot
% on the posterior at that abscissa, and the name above. The dot is what carries the
% figure, since on a multimodal law one of the three sits in the trough. The stem and the
% name are placed in units of h, the height of the tallest curve of the drawing.
    y0 = interp1(u, courbe, u0);
    plot([u0 u0], [0 1.08*h], ':', 'Color', col, 'LineWidth', lw);
    plot(u0, y0, 'o', 'Color', col, 'MarkerFaceColor', 'w', 'MarkerSize', 5, ...
         'LineWidth', lw);
    text(u0, 1.19*h, nom, 'Interpreter', 'latex', 'FontSize', fs, 'Color', col, ...
         'HorizontalAlignment', 'center');
end

function txt_s(x, y, lignes, fs, col)
% A block of centred lines, used for the labels that carry no box.
    text(x, y, lignes, 'Interpreter', 'latex', 'FontSize', fs, 'Color', col, ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end

function txt_l(x, y, lignes, fs, col)
% The same, anchored on its left edge, for a label that follows an arrow.
    text(x, y, lignes, 'Interpreter', 'latex', 'FontSize', fs, 'Color', col, ...
         'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
end

function txt_r(x, y, lignes, fs, col)
% The same, anchored on its right edge, for a label that precedes an arrow.
    text(x, y, lignes, 'Interpreter', 'latex', 'FontSize', fs, 'Color', col, ...
         'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
end

function chemin(px, py, r, lw, col, tete)
% An arrow of the flow diagram, straight or bent, as a shaft plus a filled head rather
% than an annotation, which would live in figure coordinates instead. Every bend is a
% right angle and is rounded at radius r, as the article rounds its own. The drawing
% gives its two axes the same scale, so the arc and the head are built in its
% coordinates directly. A branch that joins another one passes tete false, since the
% two then share the single head the article draws where they meet.
    X = px(1);  Y = py(1);
    for i = 2:numel(px)-1
        u = [px(i)-px(i-1), py(i)-py(i-1)];  u = u / norm(u);
        v = [px(i+1)-px(i), py(i+1)-py(i)];  v = v / norm(v);
        a = [px(i), py(i)] - r*u;            % where the shaft enters the arc
        c = a + r*v;                         % its centre, the bend being a right angle
        t = atan2(a(2)-c(2), a(1)-c(1)) + linspace(0, sign(u(1)*v(2)-u(2)*v(1))*pi/2, 12);
        X = [X, a(1), c(1)+r*cos(t)];
        Y = [Y, a(2), c(2)+r*sin(t)];
    end
    plot([X, px(end)], [Y, py(end)], '-', 'Color', col, 'LineWidth', lw);
    if ~tete
        return
    end
    dx = px(end) - px(end-1);
    dy = py(end) - py(end-1);
    d  = hypot(dx, dy);
    dx = dx / d;  dy = dy / d;
    nx = -dy;     ny = dx;
    L  = 0.020;  W = 0.0085;
    bx = px(end) - L * dx;  by = py(end) - L * dy;
    patch('XData', [px(end), bx + W*nx, bx - W*nx], ...
          'YData', [py(end), by + W*ny, by - W*ny], ...
          'FaceColor', col, 'EdgeColor', 'none');
end

function rayon(r, phi, nom, fs, lw, red)
% One radius drawn from the centre of the window, named just outside the circle it
% measures, along its own direction. Placing the name outside keeps it off the
% training points and off the particles, whatever the radius.
    plot([0 r*cosd(phi)], [0 r*sind(phi)], '-', 'Color', red, 'LineWidth', lw);
    text((r + 0.45)*cosd(phi), (r + 0.45)*sind(phi), nom, 'Interpreter', 'latex', ...
         'FontSize', fs, 'Color', red, 'HorizontalAlignment', 'center');
end

function fig = new_figure(pdf_w, pdf_h)
% A frozen page in points, as in the kriging scripts: print writes exactly this size.
    fig = figure('Color', 'w', 'Units', 'points', 'Position', [80 80 pdf_w pdf_h]);
end

function save_pdf(fig, name, pdf_w, pdf_h, img_dir)
% Vector PDF on a frozen page, written through a temporary file and moved over the
% target only if the drawing changed, so that a rerun does not dirty the repository.
    if ~exist(img_dir, 'dir')
        mkdir(img_dir);
    end
    set(fig, 'PaperUnits', 'points', 'PaperSize', [pdf_w pdf_h], ...
             'PaperPosition', [0 0 pdf_w pdf_h]);
    target = fullfile(img_dir, [name '.pdf']);
    tmp    = [tempname '.pdf'];
    print(fig, tmp, '-dpdf', '-vector');
    if isfile(target) && strcmp(without_timestamps(tmp), without_timestamps(target))
        delete(tmp);
    else
        movefile(tmp, target, 'f');
    end
end

function s = without_timestamps(f)
% The bytes of a PDF without its creation date and document ID, the only fields that
% change between two exports of an identical figure.
    fid = fopen(f, 'r');
    s   = char(fread(fid, Inf, '*uint8')');
    fclose(fid);
    s = regexprep(s, '/(CreationDate|ModDate)\s*\([^)]*\)', '');
    s = regexprep(s, '/ID\s*\[[^\]]*\]', '');
end

function I = charge(fichier, crop, cible)
% A source photograph of the application schematic, cropped by the four fractions the
% poster cropped it by, flattened on white, trimmed of its blank border and reduced to
% roughly cible pixels across. Flattening matters beyond appearance: it leaves no
% transparent pixel in the figure, and so no reason for print to fall back to a raster
% page for the whole drawing. The reduction averages over the area each output pixel
% covers rather than picking one input pixel in each, which a wing one pixel wide would
% not survive, and it works at any ratio rather than at whole factors only.
    [A, ~, alpha] = imread(fichier);
    A = double(A);
    if ~isempty(alpha)
        a = double(alpha) / 255;
        A = A .* a + 255 * (1 - a);
    end

    [ny, nx, ~] = size(A);
    A = A(round(crop(2)*ny)+1 : ny-round(crop(4)*ny), ...
          round(crop(1)*nx)+1 : nx-round(crop(3)*nx), :);

    plein = any(A < 250, 3);                 % what is not the white of the page
    lig   = find(any(plein, 2));
    col   = find(any(plein, 1));
    A     = A(lig(1):lig(end), col(1):col(end), :);

    [ny, nx, ~] = size(A);
    if nx > cible
        A = reduit(A, round(ny * cible / nx));
        A = permute(reduit(permute(A, [2 1 3]), cible), [2 1 3]);
    end
    I = uint8(A);
end

function B = reduit(A, m)
% An image reduced to m rows, each output row being the average of the input rows it
% covers, weighted by how much of each it covers. Two calls with a transpose between
% them reduce both directions.
    n     = size(A, 1);
    bords = linspace(1, n+1, m+1);
    W     = zeros(m, n);
    for i = 1:m
        a = bords(i);  b = bords(i+1);
        for j = floor(a) : ceil(b)-1
            W(i, j) = min(b, j+1) - max(a, j);
        end
    end
    W = W ./ sum(W, 2);
    B = reshape(W * A(:, :), [m, size(A, 2), size(A, 3)]);
end

function vignette(cx, cy, w, I)
% A photograph placed by its centre. The two axes of the drawing share a scale, so the
% height follows from the pixels alone and no aspect ratio has to be carried around. The
% ordinates run top row first, since the drawing keeps its vertical axis pointing up.
    h = w * size(I, 1) / size(I, 2);
    image('XData', [cx-w/2, cx+w/2], 'YData', [cy+h/2, cy-h/2], 'CData', I);
end

function semis(cx, cy, w, h, nx, ny, col)
% The prior survey, drawn as the jittered grid of sample points it is. The jitter is a
% fixed hash of the node indices rather than a draw, so that the exported file stays the
% same from one run to the next and save_pdf has nothing to rewrite.
    [I, J] = meshgrid(1:nx, 1:ny);
    r1 = mod(sin(12.9898*I + 78.2330*J) * 43758.5453, 1);
    r2 = mod(sin(39.3468*I + 11.1350*J) * 24634.6345, 1);
    X  = cx - w/2 + (I - 0.5 + 0.6*(r1 - 0.5)) * w / nx;
    Y  = cy - h/2 + (J - 0.5 + 0.6*(r2 - 0.5)) * h / ny;
    plot(X(:), Y(:), '.', 'Color', col, 'MarkerSize', 10);
end
