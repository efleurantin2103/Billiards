function orbit = BilliardMapINVOrbit(coef, V0, N)
% BILLIARDMAPINVORBIT - Compute backward orbit under inverse billiard map
%
% SYNTAX:
%   orbit = BilliardMapINVOrbit(coef, V0, N)
%
% INPUTS:
%   coef  : 2xM matrix of boundary Fourier coefficients
%   V0    : 2x1 initial point [x; y]
%   N     : Number of backward iterations
%
% OUTPUTS:
%   orbit : 2x(N+1) matrix of orbit points, orbit(:,1) = V0
%
% DESCRIPTION:
%   Generates backward trajectory by iterating the inverse billiard map
%   N times starting from V0.
%
% DEPENDENCIES:
%   RealF3INV.m

% Initialize orbit array
orbit = zeros(2, N+1);
orbit(:, 1) = V0;

% Iterate inverse map N times
V = V0;
for n = 1:N
    F = RealF3INV(V, coef);
    V = F(:, 1);
    orbit(:, n+1) = V;
end

end