function [TX, DTX] = perK_FPmap_billiard(coef, X0, Guess, k)
% PERK_FPMAP_BILLIARD - Period-k map with guess optimization
%
% SYNTAX:
%   [TX, DTX] = perK_FPmap_billiard(coef, X0, Guess, k)
%
% INPUTS:
%   coef  : 2xM boundary Fourier coefficients
%   X0    : 2k x 1 orbit points
%   Guess : 2k x 3 initial guesses for map evaluation
%   k     : Period
%
% OUTPUTS:
%   TX  : 2k x 1 mapped orbit
%   DTX : 2k x 2k Jacobian

d = size(X0, 1);
TX = zeros(d, 1);
DTX = zeros(d, d);

% First component
FF = BilliardMap(X0(end-1:end), coef, Guess(end-1:end, :));
TX(1:2) = FF(:, 1);
DTX(1:2, end-1:end) = FF(:, 2:3);

% Remaining components
for j = 2:k
    FF1 = BilliardMap(X0(2*j-3:2*j-2), coef, Guess(2*j-3:2*j-2, :));
    TX(2*j-1:2*j) = FF1(:, 1);
    DTX(2*j-1:2*j, 2*j-3:2*j-2) = FF1(:, 2:3);
end

end