function P_lambda = P_comp_lambda_1variable(P, lambda, N)
% P_COMP_LAMBDA_1VARIABLE - Compose parameterization with scaling
%
% SYNTAX:
%   P_lambda = P_comp_lambda_1variable(P, lambda, N)
%
% INPUTS:
%   P      : dim x (N+1) Taylor coefficients
%   lambda : Scaling parameter
%   N      : Maximum degree
%
% OUTPUTS:
%   P_lambda : dim x (N+1) scaled coefficients
%
% DESCRIPTION:
%   Computes P(lambda*theta) by scaling each coefficient:
%   P_lambda_n = lambda^n * P_n

dim = size(P, 1);
P_lambda = zeros(dim, N+1);

for n = 0:N
    P_lambda(:, n+1) = lambda^n * P(:, n+1);
end

end