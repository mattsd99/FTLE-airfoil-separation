function [ux, uy] = velocity_field(px, py, fa, fb, alpha, xi, yi)
% Interpolate velocity at particle positions (px, py).
% Time-blends between two frames: alpha=0 gives fa, alpha=1 gives fb.
% Points outside the grid return zero velocity.

ua = interp2(xi, yi, fa.U, px, py, 'linear', 0);
va = interp2(xi, yi, fa.V, px, py, 'linear', 0);
ub = interp2(xi, yi, fb.U, px, py, 'linear', 0);
vb = interp2(xi, yi, fb.V, px, py, 'linear', 0);

ux = (1 - alpha)*ua + alpha*ub;
uy = (1 - alpha)*va + alpha*vb;
