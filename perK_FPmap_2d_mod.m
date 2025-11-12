function [TX, DTX] = perK_FPmap_2d_mod(coef, X0, A, f, Df, k)
% PERK_FPMAP_2D_MOD - Evaluate period-k map and derivative
%
% SYNTAX:
%   [TX, DTX] = perK_FPmap_2d_mod(coef, X0, A, f, Df, k)
%
% INPUTS:
%   coef : 2xM boundary Fourier coefficients
%   X0   : 2k x 1 orbit points [x1; y1; ...; xk; yk]
%   A    : 2k x 1 perturbation vector
%   f    : Map function handle
%   Df   : Jacobian function handle
%   k    : Period
%
% OUTPUTS:
%   TX  : 2k x 1 mapped orbit T(X) = [f(xk); f(x1); ...; f(x_{k-1})]
%   DTX : 2k x 2k Jacobian of period-k map

d = size(X0, 1);
TX = zeros(d, 1);
DTX = zeros(d, d);

% First component: f(x_k)
FF = f(X0(end-1:end), coef);
TX(1:2) = FF(:, 1) - A(end-1:end);
DTX(1:2, end-1:end) = Df(X0(end-1:end), coef);

% Remaining components: f(x_{j-1}) for j = 2, ..., k
for j = 2:k
    FF1 = f(X0(2*j-3:2*j-2), coef);
    TX(2*j-1:2*j) = FF1(:, 1) - A(2*j-3:2*j-2);
    DTX(2*j-1:2*j, 2*j-3:2*j-2) = Df(X0(2*j-3:2*j-2), coef);
end

end