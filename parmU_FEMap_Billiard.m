function [phiP, PhiP, DTP_mu] = parmU_FEMap_Billiard(mu, P, s, coef, NFFT)
% PARMU_FEMAP_BILLIARD - Unstable manifold conjugacy equation residual
%
% SYNTAX:
%   [phiP, PhiP, DTP_mu] = parmU_FEMap_Billiard(mu, P, s, coef, NFFT)
%
% INPUTS:
%   mu   : Scaling parameter (1/lambda for unstable manifold)
%   P    : dim x (NFFT+1) Taylor coefficients of parameterization
%   s    : Target norm for linear coefficient
%   coef : 2xM boundary Fourier coefficients
%   NFFT : Maximum Taylor degree
%
% OUTPUTS:
%   phiP   : Scalar constraint ||P_1||^2 - s = 0
%   PhiP   : dim x (NFFT+1) conjugacy equation residual P(theta) - F(P(mu*theta))
%   DTP_mu : dim x dim x (NFFT+1) Fourier coefficients of DF(P(mu*theta))
%
% DESCRIPTION:
%   Computes residuals for unstable manifold conjugacy equation:
%   P(theta) = F(P(mu*theta)) with normalization constraint.

% Compute P(mu*theta)
P_mu = P_comp_lambda_1variable(P, mu, NFFT);

% Compute F(P(mu*theta)) and DF(P(mu*theta))
[FP_mu, DTP_mu] = get_fourier_coeffs_extended(coef, P_mu, NFFT);

% Normalization constraint
phiP = P(:, 2).' * P(:, 2) - s;

% Conjugacy equation residual
PhiP = P - FP_mu;

end