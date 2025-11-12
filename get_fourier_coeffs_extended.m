function [FP, DFP] = get_fourier_coeffs_extended(coef, P, NFFT)
% GET_FOURIER_COEFFS_EXTENDED - Compute Fourier coefficients of F(P(z))
%
% SYNTAX:
%   [FP, DFP] = get_fourier_coeffs_extended(coef, P, NFFT)
%
% INPUTS:
%   coef : 2xM boundary Fourier coefficients
%   P    : dim x num_coeffs Taylor coefficients of parameterization
%   NFFT : Maximum degree for output coefficients
%
% OUTPUTS:
%   FP  : dim x (NFFT+1) Fourier coefficients of F(P(z))
%   DFP : dim x dim x (NFFT+1) Fourier coefficients of DF(P(z))
%
% DESCRIPTION:
%   Computes Fourier series representation of the composition F∘P
%   by evaluating at points on the unit circle and solving for coefficients.

[dim, num_coeffs] = size(P);
K = dim / 2;

% Initialize output arrays
FP = zeros(dim, NFFT+1);
DFP = zeros(dim, dim, NFFT+1);

% Process each coordinate pair
for coord = 1:2:dim-1
    
    % Determine polynomial degree and FFT size
    M = min(num_coeffs-1, NFFT);
    M_extended = 2*M + 1;
    M_extended = 2^nextpow2(M_extended);
    numpts = M_extended + 1;
    
    % Generate evaluation points on unit circle
    k_vals = 0:M_extended;
    z_k = exp(2*pi*1i*k_vals/(M_extended+1));
    
    % Evaluate P(z_k) using Horner's method
    P_vals = zeros(2, numpts);
    P_coeffs_1 = P(coord, :);
    P_coeffs_2 = P(coord+1, :);
    
    for j = 1:numpts
        P_vals(1, j) = horner(z_k(j), P_coeffs_1, num_coeffs-1);
        P_vals(2, j) = horner(z_k(j), P_coeffs_2, num_coeffs-1);
    end
    
    % Evaluate F(P(z_k)) and DF(P(z_k))
    F_vals = zeros(2, numpts);
    DF2x2_vals = zeros(2, 2, numpts);
    
    % First point: use real evaluation
    x0 = [P_vals(1,1); P_vals(2,1)];
    [F, DF] = RealF3(real(x0), coef);
    F_vals(:, 1) = F(:, 1);
    DF2x2_vals(:, :, 1) = DF;
    prev_output = F;
    
    % Remaining points: use complex evaluation with previous output
    for j = 2:numpts
        x = [P_vals(1,j); P_vals(2,j)];
        [F, DF] = ComplexF3(x, coef, prev_output);
        F_vals(:, j) = F(:, 1);
        DF2x2_vals(:, :, j) = DF;
        prev_output = F;
    end
    
    % Construct Vandermonde matrix for coefficient recovery
    A = zeros(numpts, numpts);
    for i = 1:numpts
        for j = 1:(M_extended+1)
            A(i, j) = z_k(i)^(j-1);
        end
    end
    
    % Solve for Fourier coefficients of F
    b1 = F_vals(1, :).';
    b2 = F_vals(2, :).';
    coeffs_1_ext = A \ b1;
    coeffs_2_ext = A \ b2;
    
    % Truncate to desired number of coefficients
    coeffs_1 = coeffs_1_ext(1:min(num_coeffs, M+1));
    coeffs_2 = coeffs_2_ext(1:min(num_coeffs, M+1));
    
    for k = 1:min(num_coeffs, M+1)
        FP(coord, k) = coeffs_1(k);
        FP(coord+1, k) = coeffs_2(k);
    end
    
    % Solve for Fourier coefficients of DF components
    B11 = reshape(DF2x2_vals(1, 1, :), [numpts, 1]);
    B12 = reshape(DF2x2_vals(1, 2, :), [numpts, 1]);
    B21 = reshape(DF2x2_vals(2, 1, :), [numpts, 1]);
    B22 = reshape(DF2x2_vals(2, 2, :), [numpts, 1]);
    
    coeffs_11_ext = A \ B11;
    coeffs_12_ext = A \ B12;
    coeffs_21_ext = A \ B21;
    coeffs_22_ext = A \ B22;
    
    coeffs_11 = coeffs_11_ext(1:min(num_coeffs, M+1));
    coeffs_12 = coeffs_12_ext(1:min(num_coeffs, M+1));
    coeffs_21 = coeffs_21_ext(1:min(num_coeffs, M+1));
    coeffs_22 = coeffs_22_ext(1:min(num_coeffs, M+1));
    
    for k = 1:min(num_coeffs, M+1)
        DFP(coord, coord, k) = coeffs_11(k);
        DFP(coord, coord+1, k) = coeffs_12(k);
        DFP(coord+1, coord, k) = coeffs_21(k);
        DFP(coord+1, coord+1, k) = coeffs_22(k);
    end
end

% Reorder coordinates: move last pair to first position
FP_temp = FP;
DFP_temp = DFP;

FP(1:2, :) = FP_temp(end-1:end, :);
DFP(1:2, :, :) = DFP_temp(end-1:end, :, :);

for k = 2:(dim/2)
    FP(2*k-1:2*k, :) = FP_temp(2*k-3:2*k-2, :);
    DFP(2*k-1:2*k, :, :) = DFP_temp(2*k-3:2*k-2, :, :);
end

end

% =========================================================================
% HELPER FUNCTION
% =========================================================================

function Px = horner(x, A, order)
% HORNER - Evaluate polynomial using Horner's method
%
% INPUTS:
%   x     : Evaluation point
%   A     : Polynomial coefficients [a_0, a_1, ..., a_order]
%   order : Polynomial degree
%
% OUTPUT:
%   Px : Polynomial value at x

if order == 0
    p_i = A(1);
else
    p_i = A(order+1) * x + A(order);
    for i = 1:(order-1)
        p_i = p_i * x + A(order-i);
    end
end

Px = p_i;

end