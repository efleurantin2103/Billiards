function orbit = BilliardMapOrbit(coef, f, V0, N)
% BILLIARDMAPORBIT - Compute forward orbit under billiard map
%
% INPUTS:
%   coef : 2xM boundary Fourier coefficients
%   f    : Function handle for map, f(V, coef)
%   V0   : 2x1 initial point
%   N    : Number of iterations
%
% OUTPUTS:
%   orbit : 2x(N+1) orbit points

orbit = zeros(2, N+1);
orbit(:, 1) = V0;
V = V0;

for n = 1:N
    F = f(V, coef);
    V = F(:, 1);
    orbit(:, n+1) = V;
end

end