function MapManifold_Billiard(coef, f, Df, numOrbits, numIterates, bounds, ...
    xplotinit, X0, per, scale_u, scale_s, num_orbit_steps, tablenum, level)
% MAPMANIFOLDBILLIARD - Compute and visualize stable/unstable manifolds
%
% SYNTAX:
%   MapManifold_Billiard(coef, f, Df, numOrbits, numIterates, bounds, ...
%       xplotinit, X0, per, scale_u, scale_s, num_orbit_steps, tablenum, level)
%
% DESCRIPTION:
%   Computes stable and unstable manifolds of periodic orbits for billiard
%   maps using parameterization method. Uses Newton iteration to refine
%   periodic orbits, then solves conjugacy equations for manifolds.
%
% INPUTS:
%   coef            : 2xM boundary Fourier coefficients
%   f, Df           : Map and derivative function handles
%   numOrbits       : Number of sample orbits for background
%   numIterates     : Iterations per sample orbit
%   bounds          : Plot bounds [xmin, xmax, ymin, ymax]
%   xplotinit       : Initial y-coordinate for sample orbits
%   X0              : 2*per x 1 initial guess for periodic orbit
%   per             : Period of orbit
%   scale_u, scale_s: Scaling for unstable/stable eigenvectors
%   num_orbit_steps : Iterations for global manifold computation
%   tablenum        : Table identifier for file saving
%   level           : Orbit level (for symmetric pairs)

% =========================================================================
% SECTION 1: GENERATE BACKGROUND ORBITS
% =========================================================================

x_axis = linspace(bounds(3), bounds(4), numOrbits);
orbits = [];

for k = 1:numOrbits
    x0 = x_axis(k);
    y0 = xplotinit;
    orbit1 = BilliardMapOrbit(coef, f, [x0; y0], numIterates);
    orbits = [orbits, orbit1];
end

% Initial visualization
figure
hold on
plot(mod(orbits(2, :), 1), orbits(1, :), 'g.')
axis(bounds)

% =========================================================================
% SECTION 2: REFINE PERIODIC ORBIT WITH NEWTON ITERATION
% =========================================================================

A = zeros(size(X0));
[FX, DFX] = newtonPerK_2dmap_mod(coef, X0, A, f, Df, per);

X = X0;
for n = 1:8
    Delta = -DFX \ FX;
    X = X + Delta;
    [FX, DFX] = newtonPerK_2dmap_mod(coef, X, A, f, Df, per);
    thisNewtonError = norm(inv(DFX) * FX, inf);
end

fprintf('Refined periodic orbit found. Final Newton error: %.3e\n', thisNewtonError);
plot(X(2:2:2*per), X(1:2:2*per), 'k.', 'MarkerSize', 20)

% =========================================================================
% SECTION 3: STABILITY ANALYSIS
% =========================================================================

% Compute linearization at periodic orbit
[TX, DTX] = perK_FPmap_2d_mod(coef, X, A, f, Df, per);
error = norm(TX - X, inf);

% Eigenvalue decomposition
[V, Sigma] = eig(DTX);
diagsigma = diag(Sigma);

%fprintf('Eigenvalues raised to period: %f\n', real(diagsigma.^per));

% Identify stable and unstable eigenvalues (real, with appropriate magnitude)
stableind = find((abs(real(diagsigma)) < 1) .* (abs(imag(diagsigma)) < 2*eps) == 1);
stableind = stableind(1);

unstableind = find((abs(real(diagsigma)) > 1) .* (abs(imag(diagsigma)) < 2*eps) == 1);
unstableind = unstableind(1);

% Extract eigenvalues and eigenvectors
% Switch Xi_u, Xi_s orientations for boundary cases
lambda_u = real(Sigma(unstableind, unstableind));
Xi_u = -V(:, unstableind);

lambda_s = real(Sigma(stableind, stableind));
Xi_s = -V(:, stableind);

% Verify eigenvector accuracy
eigErrorCheck_u = norm(DTX * Xi_u - lambda_u * Xi_u);
eigErrorCheck_s = norm(DTX * Xi_s - lambda_s * Xi_s);

fprintf('lambda_u = %.6f (should be > 1)\n', lambda_u);
fprintf('lambda_s = %.6f (should be < 1)\n', lambda_s);
fprintf('Eigenvector errors: u = %.2e, s = %.2e\n', eigErrorCheck_u, eigErrorCheck_s);

% =========================================================================
% SECTION 4: MANIFOLD PARAMETERIZATION - SETUP
% =========================================================================

% Initialize Taylor polynomial degree
NFFT = 60;

% Initialize parameterization tensors
% Each has 2*per components (x,y for each point in orbit)
% Each component is a Taylor series of degree NFFT
Pu = zeros(2*per, NFFT+1);  % Unstable manifold
Ps = zeros(2*per, NFFT+1);  % Stable manifold

% Set initial Taylor coefficients
Pu(:, 1) = X;                % Constant term: periodic orbit
Pu(:, 2) = -scale_u * Xi_u;  % Linear term: scaled eigenvector

Ps(:, 1) = X;
Ps(:, 2) = scale_s * Xi_s;

% =========================================================================
% SECTION 5: SOLVE CONJUGACY EQUATIONS
% =========================================================================

% Unstable manifold: solve Phi[mu,P](theta) = P(sigma) - T(P(mu*sigma))
[lambda_u, Pu] = unstableNewton_maps_1var_Billiard(coef, lambda_u, Pu, NFFT);

% Plot coefficient decay
coefAxis = 0:NFFT;
theMagnitudes = zeros(1, NFFT+1);
for n = 0:NFFT
    theMagnitudes(n+1) = log10(norm(Pu(:, n+1), inf));
end

figure('Color', 'w')
plot(coefAxis, theMagnitudes, 'b.')
title('Unstable Manifold - Coefficient Magnitudes (log10)')
xlabel('Taylor Degree')
ylabel('log_{10}(||coefficient||)')

% Stable manifold: solve Phi[mu,P](theta) = P(mu*theta) - T(P(theta))
[lambda_s, Ps] = stableNewton_maps_1var_Billiard(coef, lambda_s, Ps, NFFT);

% Plot coefficient decay
theMagnitudes = zeros(1, NFFT+1);
for n = 0:NFFT
    theMagnitudes(n+1) = log10(norm(Ps(:, n+1), inf));
end

figure('Color', 'w')
plot(coefAxis, theMagnitudes, 'b.')
title('Stable Manifold - Coefficient Magnitudes (log10)')
xlabel('Taylor Degree')
ylabel('log_{10}(||coefficient||)')

fprintf('Final eigenvalues: lambda_u = %.6f, lambda_s = %.6f\n', lambda_u, lambda_s);

% =========================================================================
% SECTION 6: EVALUATE MANIFOLD PARAMETERIZATIONS
% =========================================================================

numPoints = 1000;
thetas = linspace(-1, 1, numPoints);

Pu_image = zeros(2*per, numPoints);
Ps_image = zeros(2*per, numPoints);

for n = 1:numPoints
    theta = thetas(n);
    thisPointu = evaluate_taylor(Pu, theta);
    thisPoints = evaluate_taylor(Ps, theta);
    
    Pu_image(:, n) = real(thisPointu);
    Ps_image(:, n) = real(thisPoints);
end

% =========================================================================
% SECTION 7: COMPUTE GLOBAL MANIFOLDS
% =========================================================================

all_orbits_u1 = [];
all_orbits_s1 = [];

% Extract boundary points of parameterized manifolds
last_pointu = Pu_image(1:2, numPoints);
first_pointu = Pu_image(1:2, 1);
last_points = Ps_image(1:2, numPoints);
first_points = Ps_image(1:2, 1);

% Compute backward/forward orbit segments to establish fundamental domain
x01_j = BilliardMapINVOrbit(coef, last_pointu, per);
y01_j = BilliardMapINVOrbit(coef, first_pointu, per);
x11_j = BilliardMapOrbit(coef, f, last_points, per);
y11_j = BilliardMapOrbit(coef, f, first_points, per);

x01_all = x01_j(:, end);
y01_all = y01_j(:, end);
x11_all = x11_j(:, end);
y11_all = y11_j(:, end);

% Define fundamental domain ranges in theta coordinate
num_domain_points = 1000;
fund_u1 = linspace(x01_all(2), Pu_image(2, numPoints), num_domain_points);
fund_u2 = linspace(y01_all(2), Pu_image(2, 1), num_domain_points);
fund_s1 = linspace(x11_all(2), Ps_image(2, numPoints), num_domain_points);
fund_s2 = linspace(y11_all(2), Ps_image(2, 1), num_domain_points);

fund_u1_range = [min(fund_u1), max(fund_u1)];
fund_u2_range = [min(fund_u2), max(fund_u2)];
fund_s1_range = [min(fund_s1), max(fund_s1)];
fund_s2_range = [min(fund_s2), max(fund_s2)];

% Filter manifold points within fundamental domain
selected_indices_u = [];
Pu_image2 = zeros(2*per, numPoints);

for n = 1:numPoints
    theta = thetas(n);
    thisPointu2 = evaluate_taylor(Pu, theta);
    
    keep_this_n = false;
    for k = 1:per
        second_coord_u = real(thisPointu2(2*k));
        if (second_coord_u >= fund_u1_range(1) && second_coord_u <= fund_u1_range(2)) || ...
           (second_coord_u >= fund_u2_range(1) && second_coord_u <= fund_u2_range(2))
            keep_this_n = true;
            break;
        end
    end
    
    if keep_this_n
        selected_indices_u = [selected_indices_u, n];
        Pu_image2(:, n) = real(thisPointu2);
    end
end

selected_indices_s = [];
Ps_image2 = zeros(2*per, numPoints);

for n = 1:numPoints
    theta = thetas(n);
    thisPoints = evaluate_taylor(Ps, theta);
    
    keep_this_n = false;
    for k = 1:per
        second_coord_s = real(thisPoints(2*k));
        if (second_coord_s >= fund_s1_range(1) && second_coord_s <= fund_s1_range(2)) || ...
           (second_coord_s >= fund_s2_range(1) && second_coord_s <= fund_s2_range(2))
            keep_this_n = true;
            break;
        end
    end
    
    if keep_this_n
        selected_indices_s = [selected_indices_s, n];
        Ps_image2(:, n) = real(thisPoints);
    end
end

Pu_image2 = Pu_image2(:, selected_indices_u);
Ps_image2 = Ps_image2(:, selected_indices_s);

% Iterate fundamental domain forward/backward to generate global manifolds
for k = 1:length(selected_indices_u)
    for j = 1:per
        p1 = Pu_image2((2*j-1):(2*j), k);
        next1 = BilliardMapOrbit(coef, f, p1, num_orbit_steps);
        all_orbits_u1 = [all_orbits_u1, next1];
    end
end

for k = 1:length(selected_indices_s)
    for j = 1:per
        p2 = Ps_image2((2*j-1):(2*j), k);
        next2 = BilliardMapINVOrbit(coef, p2, num_orbit_steps);
        all_orbits_s1 = [all_orbits_s1, next2];
    end
end

% =========================================================================
% SECTION 8: VISUALIZATION AND OUTPUT
% =========================================================================

fig = figure('Color', 'w');
hold on
set(gca, 'FontSize', 20)

% Plot background orbits (green)
plot(mod(orbits(2, :), 1), orbits(1, :), 'g.')

% Plot periodic orbit (black)
plot(X(2:2:2*per), X(1:2:2*per), 'k.', 'MarkerSize', 20)

% Plot global manifolds (red = unstable, blue = stable)
plot(mod(all_orbits_u1(2, :), 1), all_orbits_u1(1, :), 'r.')
plot(mod(all_orbits_s1(2, :), 1), all_orbits_s1(1, :), 'b.')

% Plot parameterized manifold curves
for j = 1:per
    plot(mod(Pu_image(2*j, :), 1), Pu_image(2*j-1, :), 'r.')
    plot(mod(Ps_image(2*j, :), 1), Ps_image(2*j-1, :), 'b.')
end

xlabel('$\theta$', 'Interpreter', 'latex', 'FontSize', 20)
ylabel('$r$', 'Interpreter', 'latex', 'FontSize', 20)
set(gcf, 'Renderer', 'painters', 'PaperPositionMode', 'auto')
axis(bounds)

% Save figure and data
foldername = ['table_', num2str(tablenum)];
basefilename = ['manifold', 'tab', num2str(tablenum), ...
                'per', num2str(per), 'lev', num2str(level)];

filename = fullfile(foldername, [basefilename, '.png']);
print(fig, '-dpng', '-r300', filename)

filename_mat = fullfile(foldername, [basefilename, '.mat']);
save(filename_mat, 'orbits', 'X', 'all_orbits_u1', 'all_orbits_s1', ...
     'Pu_image', 'Ps_image', 'per');

fprintf('Figure saved: %s\n', filename);
fprintf('Data saved: %s\n', filename_mat);

end