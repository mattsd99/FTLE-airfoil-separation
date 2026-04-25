% run_bwd_ftle.m
% Backward FTLE at 6 phases with white streamlines.
% t0 = 0, 0.2, 0.4, 0.6, 0.8, 1.0 s
% Saves individual PNGs, an animated GIF, and a 2x3 comparison figure.

clear; close all; clc;

% -------------------------------------------------------------------------
% SETTINGS
% -------------------------------------------------------------------------

DATA_DIR     = './';
FILE_PATTERN = 'naca0012_angle12_%d.csv';

START_FRAMES = [1, 5, 9, 13, 17, 21];   % t0 = 0, 0.2, 0.4, 0.6, 0.8, 1.0 s
T_LENGTH     = 10;                        % integration length in frames
DT           = 0.05;                      % physical time between frames

VMIN = 2.0;
VMAX = 8.0;

PLOT_X1 = -0.1;  PLOT_X2 = 0.6;
PLOT_Y1 = -0.2;  PLOT_Y2 = 0.2;

STREAM_DENSITY = 30;
STREAM_STEP    = 0.002;
STREAM_MAXVERT = 500;

GIF_FILE  = 'ftle_backward_phases.gif';
GIF_DELAY = 0.8;

N_ROWS = 2;
N_COLS = 3;

% -------------------------------------------------------------------------
% COMPUTE FTLE + SAVE INDIVIDUAL PNGS + BUILD GIF
% -------------------------------------------------------------------------

n_runs      = numel(START_FRAMES);
gif_started = false;

% Store results for the comparison figure
ftle_all = cell(1, n_runs);
vel_all  = cell(1, n_runs);
sxcoor   = [];
sycoor   = [];
af_x     = [];
af_y     = [];

for w = 1:n_runs

    t0   = START_FRAMES(w);
    t0_s = (t0 - 1) * DT;
    T_s  = T_LENGTH * DT;

    fprintf('\n=== Phase %d/%d  |  t0 = %.2f s   T = %.2f s ===\n', w, n_runs, t0_s, T_s);

    % Compute FTLE and load velocity (store for reuse in comparison figure)
    [ftle_all{w}, sxcoor, sycoor, af_x, af_y] = ftle_naca(t0, t0 + T_LENGTH);
    [Ug, Vg, xg, yg] = load_velocity(DATA_DIR, FILE_PATTERN, t0);
    vel_all{w} = {Ug, Vg, xg, yg};

    % Streamline seeds
    [sx, sy] = meshgrid(linspace(PLOT_X1, PLOT_X2, STREAM_DENSITY), ...
                        linspace(PLOT_Y1, PLOT_Y2, STREAM_DENSITY));
    sx = sx(:)';
    sy = sy(:)';

    [xmesh, ymesh] = meshgrid(sxcoor, sycoor);

    fig = figure('Color','k', 'Position',[100 100 800 600], 'Visible','off');
    ax  = axes('Color','k', 'XColor','w', 'YColor','w');

    pcolor(ax, xmesh, ymesh, ftle_all{w}');
    shading(ax, 'interp');
    colormap(ax, inferno_cmap());
    caxis(ax, [VMIN VMAX]);

    cb = colorbar(ax);
    cb.Color = 'w';
    cb.Label.String = 'FTLE (1/s)';
    cb.Label.Color  = 'w';

    hold(ax, 'on');
    h = streamline(ax, xg, yg, Ug, Vg, sx, sy, [STREAM_STEP, STREAM_MAXVERT]);
    set(h, 'Color','w', 'LineWidth',0.5);
    fill(ax, af_x, af_y, 'k', 'EdgeColor',[0.6 0.6 0.6]);

    axis(ax, 'equal');
    xlim(ax, [PLOT_X1 PLOT_X2]);
    ylim(ax, [PLOT_Y1 PLOT_Y2]);
    xlabel(ax, 'x/c', 'Color','w');
    ylabel(ax, 'y/c', 'Color','w');
    title(ax, sprintf('Backward FTLE  |  t_0 = %.2f s   T = %.2f s', t0_s, T_s), ...
          'Color','w', 'FontSize',12);

    png_name = sprintf('ftle_backward_t0%.2fs.png', t0_s);
    exportgraphics(fig, png_name, 'Resolution',150, 'BackgroundColor','black');
    fprintf('  Saved %s\n', png_name);

    frame     = getframe(fig);
    [img, cm] = rgb2ind(frame2im(frame), 256);
    if ~gif_started
        imwrite(img, cm, GIF_FILE, 'gif', 'Loopcount',inf, 'DelayTime',GIF_DELAY);
        gif_started = true;
    else
        imwrite(img, cm, GIF_FILE, 'gif', 'WriteMode','append', 'DelayTime',GIF_DELAY);
    end

    close(fig);

end

fprintf('\nSaved %s\n', GIF_FILE);

% -------------------------------------------------------------------------
% STATIC 2x3 COMPARISON FIGURE  (reuses stored results — no recomputation)
% -------------------------------------------------------------------------

fprintf('Building comparison figure...\n');

[xmesh, ymesh] = meshgrid(sxcoor, sycoor);

fig2 = figure('Color','k', 'Position',[100 100 400*N_COLS 350*N_ROWS]);

gap_x  = 0.02;   gap_y  = 0.03;
left   = 0.04;   right  = 0.96;
bottom = 0.10;   top    = 0.93;
pw = (right - left   - (N_COLS+1)*gap_x) / N_COLS;
ph = (top   - bottom - (N_ROWS+1)*gap_y) / N_ROWS;

for w = 1:n_runs
    t0  = START_FRAMES(w);
    row = ceil(w / N_COLS);
    col = mod(w-1, N_COLS) + 1;

    ax = subplot(N_ROWS, N_COLS, w);

    pcolor(ax, xmesh, ymesh, ftle_all{w}');
    shading(ax, 'interp');
    colormap(ax, inferno_cmap());
    caxis(ax, [VMIN VMAX]);
    hold(ax, 'on');

    Ug = vel_all{w}{1};  Vg = vel_all{w}{2};
    xg = vel_all{w}{3};  yg = vel_all{w}{4};

    [sx, sy] = meshgrid(linspace(PLOT_X1, PLOT_X2, STREAM_DENSITY), ...
                        linspace(PLOT_Y1, PLOT_Y2, STREAM_DENSITY));
    sx = sx(:)';  sy = sy(:)';

    h = streamline(ax, xg, yg, Ug, Vg, sx, sy, [STREAM_STEP, STREAM_MAXVERT]);
    set(h, 'Color','w', 'LineWidth',0.4);

    fill(ax, af_x, af_y, 'k', 'EdgeColor',[0.6 0.6 0.6]);

    axis(ax, 'equal');
    xlim(ax, [PLOT_X1 PLOT_X2]);
    ylim(ax, [PLOT_Y1 PLOT_Y2]);

    ax.Position = [left + (col-1)*(pw+gap_x) + gap_x, ...
                   top  - row*(ph+gap_y) + gap_y, ...
                   pw, ph];

    xlabel(ax, 'x/c', 'Color','w');
    ylabel(ax, 'y/c', 'Color','w');
    title(ax, sprintf('t_0 = %.2f s', (t0-1)*DT), 'Color','w', 'FontSize',11);
    ax.Color = 'k';  ax.XColor = 'w';  ax.YColor = 'w';
end

cb2 = colorbar('Location','southoutside', ...
               'Units','normalized', ...
               'Position',[0.15 0.03 0.70 0.02]);
colormap(inferno_cmap());
caxis([VMIN VMAX]);
cb2.Color = 'w';
cb2.Label.String = 'FTLE (1/s)';
cb2.Label.Color  = 'w';
cb2.FontSize     = 11;

sgtitle(sprintf('Backward FTLE + Streamlines  |  T = %.2f s', T_LENGTH*DT), ...
        'Color','w', 'FontSize',13);

exportgraphics(fig2, 'ftle_backward_comparison.png', 'Resolution',200, 'BackgroundColor','black');
fprintf('Saved ftle_backward_comparison.png\n');
close(fig2);

fprintf('\nAll done.\n');

% =========================================================================
function [Ug, Vg, xg, yg] = load_velocity(data_dir, file_pattern, frame_idx)
d = readmatrix(fullfile(data_dir, sprintf(file_pattern, frame_idx)));
[~, uid] = unique(d(:,1:2), 'rows');
d  = d(uid, :);
xg = linspace(min(d(:,1)), max(d(:,1)), 500);
yg = linspace(min(d(:,2)), max(d(:,2)), 375);
[Xg, Yg] = meshgrid(xg, yg);
Fu = scatteredInterpolant(d(:,1), d(:,2), d(:,3), 'linear','none');
Fv = scatteredInterpolant(d(:,1), d(:,2), d(:,4), 'linear','none');
Ug = Fu(Xg, Yg);  Ug(isnan(Ug)) = 0;
Vg = Fv(Xg, Yg);  Vg(isnan(Vg)) = 0;
end

% =========================================================================
function cmap = inferno_cmap()
c = [0.00 0.00 0.00;
     0.18 0.02 0.33;
     0.50 0.05 0.47;
     0.76 0.22 0.26;
     0.96 0.57 0.10;
     0.99 0.99 0.64];
t  = linspace(0, 1, size(c,1));
tq = linspace(0, 1, 256);
cmap = max(0, min(1, interp1(t, c, tq)));
end