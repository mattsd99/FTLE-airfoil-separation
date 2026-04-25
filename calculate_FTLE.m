function sigma = calculate_FTLE(ix, iy, flowmap_x, flowmap_y, T, sxcoor, sycoor)
% Compute FTLE at grid point (ix, iy) using central differences.
% Returns 0 at boundary points where neighbours are not available.
%
% Builds the flow map Jacobian A, then computes:
%   sigma = log( sqrt( lambda_max(A'*A) ) ) / T

[nx, ny] = size(flowmap_x);

% Need neighbours on both sides — skip boundary points
if ix == 1 || ix == nx || iy == 1 || iy == ny
    sigma = 0;
    return
end

% Flow map Jacobian via central differences
A11 = (flowmap_x(ix+1,iy) - flowmap_x(ix-1,iy)) / (sxcoor(ix+1) - sxcoor(ix-1));
A12 = (flowmap_x(ix,iy+1) - flowmap_x(ix,iy-1)) / (sycoor(iy+1) - sycoor(iy-1));
A21 = (flowmap_y(ix+1,iy) - flowmap_y(ix-1,iy)) / (sxcoor(ix+1) - sxcoor(ix-1));
A22 = (flowmap_y(ix,iy+1) - flowmap_y(ix,iy-1)) / (sycoor(iy+1) - sycoor(iy-1));

A = [A11 A12; A21 A22];

% Largest eigenvalue of the Cauchy-Green tensor C = A'*A
lambda_max = max(eig(A'*A));
sigma      = log(lambda_max) / (2*T);
