function [ftle, sxcoor, sycoor, af_x, af_y] = ftle_naca(start_frame, end_frame)
% ftle_naca  Compute backward FTLE for a single integration window.
%
%   [ftle, sxcoor, sycoor, af_x, af_y] = ftle_naca(start_frame, end_frame)
%
%   Inputs:
%     start_frame : first frame index
%     end_frame   : last frame index  (T = (end-start) * DT)
%
%   Outputs:
%     ftle        : (NX x NY) FTLE field  (NaN inside airfoil)
%     sxcoor      : 1-D seed x coordinates
%     sycoor      : 1-D seed y coordinates
%     af_x, af_y  : airfoil outline coordinates (for plotting)
%
%   Example:
%     ftle = ftle_naca(1, 11);   % backward FTLE, T = 10 frames
%
%   Required files in same folder:
%     velocity_field.m, RungeKutta_ba.m, calculate_FTLE.m

% -------------------------------------------------------------------------
% SETTINGS  — edit here
% -------------------------------------------------------------------------

DATA_DIR     = './';
FILE_PATTERN = 'naca0012_angle12_%d.csv';
DT           = 0.05;

% Seed grid
NX = 1200;   X1 = -1.0;   X2 = 2.0;
NY = 800;    Y1 = -1.0;   Y2 = 1.0;

% -------------------------------------------------------------------------
% AIRFOIL MASK
% -------------------------------------------------------------------------

[af_x, af_y] = naca0012_geometry(12.0);

sxcoor = linspace(X1, X2, NX);
sycoor = linspace(Y1, Y2, NY);
[xmesh, ymesh] = meshgrid(sxcoor, sycoor);

sigma_x0 = xmesh';
sigma_y0 = ymesh';

airfoil_mask = inpolygon(sigma_x0, sigma_y0, af_x, af_y);

% -------------------------------------------------------------------------
% LOAD FRAMES
% -------------------------------------------------------------------------

all_idx = start_frame : end_frame;
n_steps = numel(all_idx) - 1;
T_win   = n_steps * DT;

fprintf('Loading %d frames (backward, T=%.2f)...\n', numel(all_idx), T_win);

d0 = readmatrix(fullfile(DATA_DIR, sprintf(FILE_PATTERN, all_idx(1))));
xi = linspace(min(d0(:,1)), max(d0(:,1)), 500);
yi = linspace(min(d0(:,2)), max(d0(:,2)), 375);
[Xg, Yg] = meshgrid(xi, yi);

x_lo = xi(1);   x_hi = xi(end);
y_lo = yi(1);   y_hi = yi(end);

frames(numel(all_idx)) = struct('U',[],'V',[],'t',[]);

for k = 1:numel(all_idx)
    d  = readmatrix(fullfile(DATA_DIR, sprintf(FILE_PATTERN, all_idx(k))));
    Fu = scatteredInterpolant(d(:,1), d(:,2), d(:,3), 'linear','none');
    Fv = scatteredInterpolant(d(:,1), d(:,2), d(:,4), 'linear','none');
    Ui = Fu(Xg,Yg);  Ui(isnan(Ui)) = 0;
    Vi = Fv(Xg,Yg);  Vi(isnan(Vi)) = 0;
    frames(k).U = Ui;
    frames(k).V = Vi;
    frames(k).t = (start_frame - 1)*DT + (k-1)*DT;
    fprintf('  frame %d loaded\n', all_idx(k));
end

% -------------------------------------------------------------------------
% ADVECT PARTICLES  —  backward integration (dt < 0)
% -------------------------------------------------------------------------

sigma_x    = sigma_x0;
sigma_y    = sigma_y0;
CalcFTLE   = true(NX, NY);
LeftDomain = false(NX, NY);
ftle       = zeros(NX, NY);

for step = 1:n_steps
    fa = frames(step);
    fb = frames(step+1);
    dt = -abs(fb.t - fa.t);   % negative = backward

    [sigma_x_new, sigma_y_new] = RungeKutta_ba(sigma_x, sigma_y, fa, fb, xi, yi, abs(fb.t - fa.t));

    % Particles that just left the domain — compute their FTLE now
    just_left = (sigma_x_new < x_lo | sigma_x_new > x_hi | ...
                 sigma_y_new < y_lo | sigma_y_new > y_hi) & ~LeftDomain;

    [esc_i, esc_j] = find(just_left & CalcFTLE);
    for n = 1:numel(esc_i)
        ix = esc_i(n);   iy = esc_j(n);
        ftle(ix,iy)     = calculate_FTLE(ix, iy, sigma_x, sigma_y, T_win, sxcoor, sycoor);
        CalcFTLE(ix,iy) = false;
    end

    LeftDomain = LeftDomain | just_left;

    sigma_x(~LeftDomain) = sigma_x_new(~LeftDomain);
    sigma_y(~LeftDomain) = sigma_y_new(~LeftDomain);

    fprintf('  step %d/%d\n', step, n_steps);
end

% Remaining particles still inside at the final step
[rem_i, rem_j] = find(CalcFTLE);
for n = 1:numel(rem_i)
    ix = rem_i(n);   iy = rem_j(n);
    ftle(ix,iy) = calculate_FTLE(ix, iy, sigma_x, sigma_y, T_win, sxcoor, sycoor);
end

ftle(airfoil_mask) = NaN;

% =========================================================================
% LOCAL GEOMETRY FUNCTION
% =========================================================================

function [xr, yr] = naca0012_geometry(aoa_deg)
n    = 400;
beta = linspace(0, pi, n);
xc   = 0.5*(1 - cos(beta));
yt   = 5*0.12*(0.2969*sqrt(xc) - 0.1260*xc - 0.3516*xc.^2 ...
              + 0.2843*xc.^3   - 0.1015*xc.^4);
x  = [xc,  xc(end-1:-1:2)];
y  = [yt, -yt(end-1:-1:2)];
a  = deg2rad(aoa_deg);
xr = x*cos(a) + y*sin(a) - (x(1)*cos(a) + y(1)*sin(a));
yr = -x*sin(a) + y*cos(a) - (-x(1)*sin(a) + y(1)*cos(a));
end

end
