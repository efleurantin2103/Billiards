function [FX, DFX] = newtonPerK_2dmap_mod(coef, X0, A, f, Df, k)
% NEWTONPERK_2DMAP_MOD - Shooting map for period-k orbit Newton iteration
%
% SYNTAX:
%   [FX, DFX] = newtonPerK_2dmap_mod(coef, X0, A, f, Df, k)
%
% INPUTS:
%   coef : 2xM boundary Fourier coefficients
%   X0   : 2k x 1 vector of orbit points [x1; y1; x2; y2; ...; xk; yk]
%   A    : 2k x 1 perturbation vector (typically zeros)
%   f    : Function handle for map evaluation
%   Df   : Function handle for Jacobian evaluation
%   k    : Period of orbit
%
% OUTPUTS:
%   FX  : 2k x 1 shooting map residual F(X) - X
%   DFX : 2k x 2k Jacobian of shooting map
%
% DESCRIPTION:
%   Constructs the shooting map F^k(X) - X and its derivative for
%   Newton iteration to find periodic orbits of period k.

% Dimension of state vector
d = size(X0, 1);

% Initialize output arrays
FX = zeros(d, 1);
DFX = -eye(2*k);

% First component: f(x_k) - x_1
FF = f(X0(end-1:end), coef);
FX(1:2) = FF(:, 1) - X0(1:2) - A(end-1:end);
DFX(1:2, end-1:end) = Df(X0(end-1:end), coef);

% Remaining components: f(x_{j-1}) - x_j for j = 2, ..., k
for j = 2:k
    % Map residual
    FF1 = f(X0(2*j-3:2*j-2), coef);
    FX(2*j-1:2*j) = FF1(:, 1) - X0(2*j-1:2*j) - A(2*j-3:2*j-2);
    
    % Jacobian block
    DFX(2*j-1:2*j, 2*j-3:2*j-2) = Df(X0(2*j-3:2*j-2), coef);
end

end