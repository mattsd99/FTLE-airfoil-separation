function [x_new, y_new] = RungeKutta_ba(x, y, fa, fb, xi, yi, dt)
% Advance particles one time step backward using RK4.
% Self-contained — does not require RungeKutta_f.m

dt_neg = -dt;   % negative = backward integration

alpha = 0.5;

[k1x, k1y] = velocity_field(x,                   y,                   fa, fb, alpha, xi, yi);
[k2x, k2y] = velocity_field(x + 0.5*dt_neg*k1x,  y + 0.5*dt_neg*k1y,  fa, fb, alpha, xi, yi);
[k3x, k3y] = velocity_field(x + 0.5*dt_neg*k2x,  y + 0.5*dt_neg*k2y,  fa, fb, alpha, xi, yi);
[k4x, k4y] = velocity_field(x +     dt_neg*k3x,  y +     dt_neg*k3y,  fa, fb, alpha, xi, yi);

x_new = x + (dt_neg/6)*(k1x + 2*k2x + 2*k3x + k4x);
y_new = y + (dt_neg/6)*(k1y + 2*k2y + 2*k3y + k4y);
