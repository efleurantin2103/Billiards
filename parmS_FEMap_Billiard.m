function [phiP, PhiP, DTP] = parmS_FEMap_Billiard(mu, P, s, coef, NFFT)
% PARMS_FEMAP_BILLIARD - Stable manifold conjugacy equation residual
%
% SYNTAX:
%   [phiP, PhiP, DTP] = parmS_FEMap_Billiard(mu, P, s, coef, NFFT)
%
% INPUTS:
%   mu   : Scaling parameter (1/lambda for stable manifold)
%   P    : dim x (NFFT+1) Taylor coefficients of parameterization
%   s    : Target norm for linear coefficient
%   coef : 2xM boundary Fourier coefficients
%   NFFT : Maximum Taylor degree
%
% OUTPUTS:
%   phiP : Scalar constraint ||P_1||^2 - s = 0
%   PhiP : dim x (NFFT+1) conjugacy equation residual P(mu*theta) - F(P(theta))
%   DTP  : dim x dim x (NFFT+1) Fourier coefficients of DF(P(theta))
%
% DESCRIPTION:
%   Computes residuals for stable manifold conjugacy equation:
%   P(mu*theta) = F(P(theta)) with normalization constraint.

P_mu = P_comp_lambda_1variable(P, mu, NFFT);
[FP, DTP] = get_fourier_coeffs_extended(coef, P, NFFT);

phiP = P(:, 2).' * P(:, 2) - s;
PhiP = P_mu - FP;

end