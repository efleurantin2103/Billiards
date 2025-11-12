function [lambda_new, P_new] = stableNewton_maps_1var_Billiard(...
    coef, lambda0, P0, NFFT)
% STABLENEWTON_MAPS_1VAR_BILLIARD_PER2 - Solve stable manifold equation
%
% SYNTAX:
%   [lambda_new, P_new] = stableNewton_maps_1var_Billiard_Per2(...
%       coef, lambda0, P0, NFFT, MP)
%
% INPUTS:
%   coef    : 2xM boundary Fourier coefficients
%   lambda0 : Initial eigenvalue guess (should be < 1)
%   P0      : 2K x (NFFT+1) initial parameterization
%   NFFT    : Maximum Taylor degree
%
% OUTPUTS:
%   lambda_new : Refined eigenvalue
%   P_new      : 2K x (NFFT+1) refined parameterization
%
% DESCRIPTION:
%   Solves P(mu*theta) = F(P(theta)) with normalization constraint
%   using Newton iteration. Includes adaptive stopping to prevent
%   divergence and rollback on error increase.

K = size(P0, 1) / 2;
dim = 2 * K;

% Normalization scale
scale = P0(:, 2).' * P0(:, 2);

if abs(lambda0) > 1
    error('Initial eigenvalue must be < 1 for stable manifold');
end

% Initialize state
mu = lambda0;
P = P0;

% Initial residual
[Phi0, PhiP, DTP_mu] = parmS_FEMap_Billiard(mu, P, scale, coef, NFFT);
ErrorPhi = norm(evaluate_taylor(PhiP, 1), inf);

fprintf('Initial error: %.3e\n', ErrorPhi);

prevErrors = [];

% Newton iteration
for m = 1:10
    % Store previous state for potential rollback
    mu_prev = mu;
    P_prev = P;
    PhiP_prev = PhiP;
    
    % Construct Jacobian blocks
    P_mu = P_comp_lambda_1variable(P, mu, NFFT);
    
    DTP_matrix = zeros(dim*(NFFT+1), dim*(NFFT+1));
    Mu_matrix = zeros(dim*(NFFT+1), dim*(NFFT+1));
    
    for n = 0:NFFT
        thisMatrix = DTP_mu(:, :, n+1);
        rowIndex = n*dim + 1;
        columnIndex = 1;
        
        % Fill DTP block diagonal
        for j = 0:NFFT-n
            DTP_matrix(rowIndex:rowIndex+dim-1, columnIndex:columnIndex+dim-1) = thisMatrix;
            rowIndex = rowIndex + dim;
            columnIndex = columnIndex + dim;
        end
        
        % Fill Mu block diagonal
        Mu_matrix(n*dim+1:(n+1)*dim, n*dim+1:(n+1)*dim) = mu^n * eye(dim);
    end
    
    % Partial derivatives
    partialPhi0 = zeros(1, dim*(NFFT+1));
    partialPhi0(dim+1:2*dim) = 2 * P(:, 2).';
    
    P_prime_mu_theta = [];
    for n = 0:NFFT
        P_prime_mu_theta = [P_prime_mu_theta; n * mu^(n-1) * P_mu(:, n+1)];
    end
    
    % Assemble full Jacobian
    DPhi_matrix = [0, partialPhi0; P_prime_mu_theta, Mu_matrix - DTP_matrix];
    
    % Newton step
    Phi1_vector = reshape(PhiP, [dim*(NFFT+1), 1]);
    Phi_vector = [Phi0; Phi1_vector];
    
    P_vector = reshape(P, [dim*(NFFT+1), 1]);
    X_state = [mu; P_vector];
    
    Delta = -DPhi_matrix \ Phi_vector;
    X_state = X_state + Delta;
    
    % Extract updated state
    mu = X_state(1);
    P_vector_new = X_state(2:end);
    
    P = reshape(P_vector_new, [dim, NFFT+1]);
    
    % Compute new residual
    [Phi0, PhiP, DTP_mu] = parmS_FEMap_Billiard(mu, P, scale, coef, NFFT);
    ErrorPhi = norm(evaluate_taylor(PhiP, 1), inf);
    
    fprintf('Iteration %d: mu = %.6f, Error = %.3e\n', m, mu, ErrorPhi);
    
    % Adaptive stopping
    prevErrors = [prevErrors, ErrorPhi];
    
    if m > 2
        if ErrorPhi > prevErrors(end-1)
            if prevErrors(end-1) < 1e-14
                fprintf('Converged below tolerance (1e-14).\n');
                break;
            else
                fprintf('Error increased, reverting to previous state.\n');
                mu = mu_prev;
                P = P_prev;
                PhiP = PhiP_prev;
                break;
            end
        end
    end
end

lambda_new = mu;
P_new = P;

fprintf('Final error history: ');
fprintf('%.2e ', prevErrors);
fprintf('\n');

end