% =========================================================================
% BILLIARD MAP PERIODIC ORBIT (STABLE AND UNSTABLE MANIFOLDS)
% =========================================================================
% 
% DESCRIPTION:
%   This script analyzes periodic orbits of billiard maps with perturbed
%   boundaries. It computes stable and unstable manifolds of periodic
%   points using numerical continuation methods.
%
% MAIN FEATURES:
%   - Supports 5 different table geometries (Tables 0-4, labeled A-E)
%   - Computes periodic orbits of various periods (2, 3, 4, 5, 10, 30)
%   - Visualizes stable/unstable manifolds around periodic points
%   - Uses Newton iteration for high-precision periodic point location
%
% REQUIRED DEPENDENCIES:
%   - BilliardIteration/ : Core iteration algorithms
%   - table_0/ through table_4/ : Table-specific data and functions
%   - BilliardMap.m : Core billiard map function
%   - MapManifold_Billiard.m : Manifold computation routine
%
% USAGE:
%   1. Uncomment ONE set of parameters below (table coefficients + orbit)
%   2. Run the script
%   3. Results will be visualized and saved to appropriate table directory
%
% =========================================================================

% Clean workspace
clearvars 
clear all
format long
close all

% Add required paths to MATLAB search path
addpath('BilliardIteration/', 'table_0/', 'table_1/', 'table_2/', ...
        'table_3/', 'table_4/')

% =========================================================================
% SECTION 1: TABLE GEOMETRY SELECTION
% =========================================================================
% Define boundary perturbation coefficients for the billiard table.
% The boundary is parameterized as:
%   x(s) = s + coefx(1)*sin(pi*s) + coefx(2)*sin(2*pi*s) + ...
%   y(t) = t + coefy(1)*sin(pi*t) + coefy(2)*sin(2*pi*t) + ...
%
% Available tables (uncomment ONE):
% -------------------------------------------------------------------------

% TABLE 0 (Label: A) - Moderate single-mode perturbation
% coefx = [1.1, 0.03]; 
% coefy = [1, 0.03];

% TABLE 1 (Label: D) - Stronger single-mode perturbation
% coefx = [2, 0.04]; 
% coefy = [1, 0.035];

% TABLE 2 (Label: C) - Two-mode perturbation
% coefx = [1.1, 0.08, 0.0002];  
% coefy = [1, 0.095, 0.0001];

% TABLE 3 (Label: B) - Moderate two-mode perturbation
% coefx = [1.1, 0.05, 0.00015];  
% coefy = [1, 0.035, 0.0001];

% TABLE 4 (Label: E) - Strong single-mode perturbation
coefx = [2, 0.05]; 
coefy = [1, 0.065];

% -------------------------------------------------------------------------
% Combine x and y coefficients into single array
coef = [coefx; coefy];

% =========================================================================
% SECTION 2: NUMERICAL PARAMETERS
% =========================================================================

numOrbits = 100;      % Number of intial orbit segments for phase space
numIterates = 400;    % Number of iterations per orbit segment
bounds = [-0.1, 1.1, -0.999, 0.999];  % Plot bounds: [xmin, xmax, ymin, ymax]
xplotinit = 0.23;     % Initial x-coordinate for plotting reference orbits

% =========================================================================
% SECTION 3: PERIODIC ORBIT SELECTION
% =========================================================================
% Uncomment ONE periodic orbit configuration below.
% Each configuration specifies:
%   V         : Initial guess for periodic point (x, y coordinates)
%   scale_u   : Scaling factor for unstable manifold computation
%   scale_s   : Scaling factor for stable manifold computation
%   numits    : Number of Newton iterations for refinement
%   tablenum  : Table index (0-4)
%   per       : Period of the orbit
%   level     : Orbit level (0 or 1 for symmetric pairs)
%
% IMPORTANT NOTES:
%   - Set n=1 in MapManifold_Billiard for Newton step during period 2
%   - Some configurations require adjusted NFFT values (noted below)
%   - "exact" indicates high-precision convergence achieved
% =========================================================================

% -------------------------------------------------------------------------
% PERIOD 2 ORBITS
% -------------------------------------------------------------------------
% Period-2 orbits exist for all tables and lie on the vertical symmetry axis

% TABLE 0 - Period 2 (EXACT)
% V = [0.0; 0.5];  
% scale_u = 0.48; 
% scale_s = 0.46; 
% numits = 6;
% tablenum = 0; 
% per = 2; 
% level = 0;

% TABLE 1 - Period 2 (Set NFFT = 10)
% V = [0.0; 0.5];  
% scale_u = 0.12; 
% scale_s = 0.1; 
% numits = 2;
% tablenum = 1; 
% per = 2; 
% level = 0;

% TABLE 2 - Period 2
% V = [0.0; 0.5];  
% scale_u = 0.32; 
% scale_s = 0.32; 
% numits = 6;
% tablenum = 2; 
% per = 2; 
% level = 0;

% TABLE 3 - Period 2
% V = [0.0; 0.5];  
% scale_u = 0.45; 
% scale_s = 0.43; 
% numits = 5;
% tablenum = 3; 
% per = 2; 
% level = 0;

% TABLE 4 - Period 2 (Set NFFT = 10) (EXACT)
% V = [0.0; 0.5];  
% scale_u = 0.12; 
% scale_s = 0.12; 
% numits = 2;
% tablenum = 4; 
% per = 2; 
% level = 0;

% -------------------------------------------------------------------------
% PERIOD 3 ORBITS
% -------------------------------------------------------------------------
% Period-3 orbits come in symmetric pairs (level 0 and level 1)

% TABLE 0 - Period 3, Level 0 (EXACT)
% V = [0.533341096440513; 0.853382728418630]; 
% scale_u = 0.25; 
% scale_s = 0.25; 
% numits = 10;
% tablenum = 0; 
% per = 3;
% level = 0;

% TABLE 0 - Period 3, Level 1 (EXACT)
% V = [-0.533341096440513; 0.853382728418630]; 
% scale_u = 0.25; 
% scale_s = 0.25; 
% numits = 10;
% tablenum = 0; 
% per = 3;
% level = 1;

% TABLE 2 - Period 3, Level 0 (EXACT)
% V = [0.445422324734639; 0.5]; 
% scale_u = 0.28; 
% scale_s = 0.28; 
% numits = 6;
% tablenum = 2; 
% per = 3;
% level = 0;

% TABLE 2 - Period 3, Level 1 (EXACT)
% V = [-0.445422324734638; 0.5]; 
% scale_u = 0.28; 
% scale_s = 0.28; 
% numits = 6;
% tablenum = 2; 
% per = 3;
% level = 1;

% TABLE 3 - Period 3, Level 0 (EXACT)
% V = [0.532766634761092; 0.334785202417012]; 
% scale_u = 0.27; 
% scale_s = 0.27; 
% numits = 7;
% tablenum = 3; 
% per = 3;
% level = 0;

% TABLE 3 - Period 3, Level 1 (EXACT)
% V = [-0.532766634761092; 0.334785202417012]; 
% scale_u = 0.27; 
% scale_s = 0.27; 
% numits = 7;
% tablenum = 3; 
% per = 3;
% level = 1;

% -------------------------------------------------------------------------
% PERIOD 4 ORBITS
% -------------------------------------------------------------------------

% TABLE 1 - Period 4, Level 0 (Set NFFT = 35)
% V = [0.5; 0];
% scale_u = 0.114; 
% scale_s = 0.114; 
% numits = 8; 
% tablenum = 1; 
% per = 4; 
% level = 0;

% TABLE 1 - Period 4, Level 1 (Set NFFT = 35)
% V = [-0.892289722890672; 0.221584381041448];
% scale_u = 0.114; 
% scale_s = 0.114; 
% numits = 8; 
% tablenum = 1; 
% per = 4; 
% level = 1;

% TABLE 2 - Period 4, Level 0 (EXACT)
% V = [0.707809782538882; 0];
% scale_u = 0.2; 
% scale_s = 0.2; 
% numits = 8; 
% tablenum = 2; 
% per = 4; 
% level = 0;

% TABLE 2 - Period 4, Level 1 (EXACT)
% V = [-0.733430851895140; 0.215220128706573];
% scale_u = 0.2; 
% scale_s = 0.2; 
% numits = 8; 
% tablenum = 2; 
% per = 4; 
% level = 1;

% TABLE 4 - Period 4, Level 0 (Set NFFT = 25)
% V = [0.89; 0.22];
% scale_u = 0.14; 
% scale_s = 0.115; 
% numits = 3; 
% tablenum = 4; 
% per = 4; 
% level = 0;

% TABLE 4 - Period 4, Level 1 (Set NFFT = 25)
% V = [-0.89; 0.22];
% scale_u = 0.14; 
% scale_s = 0.115; 
% numits = 3; 
% tablenum = 4; 
% per = 4; 
% level = 1;

% -------------------------------------------------------------------------
% PERIOD 5 ORBITS
% -------------------------------------------------------------------------

% TABLE 0 - Period 5, Level 0 (EXACT)
% V = [0.258508514723553; 0.430805087087808]; 
% scale_u = 0.29; 
% scale_s = 0.29; 
% numits = 11; 
% tablenum = 0; 
% per = 5; 
% level = 0;

% TABLE 0 - Period 5, Level 1 (EXACT)
% V = [-0.258508514723553; 0.430805087087808]; 
% scale_u = 0.29; 
% scale_s = 0.29; 
% numits = 11; 
% tablenum = 0; 
% per = 5; 
% level = 1;

% TABLE 1 - Period 5, Level 0 (Set NFFT = 30)
% V = [0.75; 0.44];
% scale_u = 0.125; 
% scale_s = 0.123; 
% numits = 10; 
% tablenum = 1; 
% per = 5; 
% level = 0;

% TABLE 1 - Period 5, Level 1 (Set NFFT = 30)
% V = [-0.75; 0.44];
% scale_u = 0.125; 
% scale_s = 0.123; 
% numits = 10; 
% tablenum = 1; 
% per = 5; 
% level = 1;

% TABLE 2 - Period 5, Level 0 (EXACT)
% V = [0.769488617362806; 0.499999999999998];
% scale_u = 0.17; 
% scale_s = 0.17; 
% numits = 25; 
% tablenum = 2; 
% per = 5; 
% level = 0;

% TABLE 2 - Period 5, Level 1 (EXACT)
% V = [-0.769488617362806; 0.499999999999998];
% scale_u = 0.17; 
% scale_s = 0.17; 
% numits = 25; 
% tablenum = 2; 
% per = 5; 
% level = 1;

% TABLE 3 - Period 5, Level 0 (EXACT)
% Note: Alternative initial guess found: [0.258508514723555; 0.430805087087806]
% V = [0.264053241259983; 0.422327043109421]; 
% scale_u = 0.21; 
% scale_s = 0.2; 
% numits = 10; 
% tablenum = 3; 
% per = 5; 
% level = 0;

% TABLE 3 - Period 5, Level 1 (EXACT)
% Note: Alternative initial guess found: [-0.258508514723555; 0.430805087087806]
% V = [-0.264053241259984; 0.422327043109420]; 
% scale_u = 0.21; 
% scale_s = 0.2; 
% numits = 10; 
% tablenum = 3; 
% per = 5; 
% level = 1;

% TABLE 4 - Period 5, Level 0 (Set NFFT = 25)
% V = [0.7; 0.46];
% scale_u = 0.125; 
% scale_s = 0.123; 
% numits = 7; 
% tablenum = 4; 
% per = 5; 
% level = 0;

% TABLE 4 - Period 5, Level 1 (Set NFFT = 25)
% V = [-0.7; 0.46];
% scale_u = 0.125; 
% scale_s = 0.123; 
% numits = 7; 
% tablenum = 4; 
% per = 5; 
% level = 1;

% -------------------------------------------------------------------------
% PERIOD 10 ORBITS
% -------------------------------------------------------------------------

% TABLE 0 - Period 10, Level 0 (EXACT)
% V = [-0.328258716973881; 0.700963992177044];
% scale_u = 0.39; 
% scale_s = 0.39; 
% numits = 10; 
% tablenum = 0; 
% per = 10; 
% level = 0;

% TABLE 1 - Period 10, Level 0 (EXACT)
% V = [0; 0.380151075279130];
% scale_u = 0.7;
% scale_s = 0.7; 
% numits = 20; 
% tablenum = 1; 
% per = 10; 
% level = 0;

% TABLE 2 - Period 10, Level 0 (EXACT)
% V = [0.165304542204223; 0.121962404960313];
% scale_u = 0.34; 
% scale_s = 0.34; 
% numits = 8; 
% tablenum = 2; 
% per = 10; 
% level = 0;

% TABLE 3 - Period 10, Level 0 (EXACT)
% V = [-0.174623726717599; 0.376735714816663];
% scale_u = 0.46; 
% scale_s = 0.48; 
% numits = 4; 
% tablenum = 3; 
% per = 10; 
% level = 0;

% TABLE 4 - Period 10, Level 0 (EXACT)
% Note: Alternative initial guess found: [0; 0.3]
V = [-0.410647076548679; 0.127496598101885]; 
scale_u = 0.65;
scale_s = 0.65; 
numits = 20; 
tablenum = 4; 
per = 10; 
level = 0;

% -------------------------------------------------------------------------
% PERIOD 30 ORBITS
% -------------------------------------------------------------------------

% TABLE 0 - Period 30, Level 0 (EXACT)
% V = [0.166688593308156; 0.510291508219422]; 
% scale_u = 0.15; 
% scale_s = 0.15; 
% numits = 20; 
% tablenum = 0; 
% per = 30;
% level = 0;

% TABLE 0 - Period 30, Level 1 (EXACT)
% V = [-0.166688593308155; 0.489708491780578]; 
% scale_u = 0.15; 
% scale_s = 0.15; 
% numits = 20; 
% tablenum = 0; 
% per = 30;
% level = 1;

% TABLE 3 - Period 30, Level 0 (EXACT)
% Note: Only 15 distinct points observed; verify if period is 30 or 15
% V = [0.487742098528984; 0.5]; 
% scale_u = 0.15; 
% scale_s = 0.15; 
% numits = 30; 
% tablenum = 3; 
% per = 30;
% level = 0;

% TABLE 3 - Period 30, Level 1 (EXACT)
% Note: Only 15 distinct points observed; verify if period is 30 or 15
% V = [-0.487742098528985; 0.5]; 
% scale_u = 0.15; 
% scale_s = 0.15; 
% numits = 30; 
% tablenum = 3; 
% per = 30;
% level = 1;

% =========================================================================
% SECTION 4: PERIODIC ORBIT COMPUTATION
% =========================================================================
% Build the full periodic orbit by iterating the map 'per' times starting
% from the initial guess V. Store all points in the orbit.

% Initialize array to store all points in the periodic orbit
X0 = zeros(2*per, 1);  % Each point has 2 coordinates (x, y)
X0(1:2) = V;           % Set first point to initial guess

% Iterate the map to generate the full periodic orbit
for j = 1:(per-1)
    twoj = 2*j;  % Index for j-th point
    
    % Apply billiard map to current point
    output1 = f(X0(twoj-1:twoj), coef);
    
    % Store next point in orbit
    X0(twoj+1:twoj+2) = output1(:, 1);
end

% =========================================================================
% SECTION 5: VERIFY PERIODICITY
% =========================================================================
% Check that the last iterate returns to the first point (periodic condition)

fV_return = f(X0(end-1:end), coef);  % Map final point
error0 = norm(X0(1:2) - fV_return(:, 1));  % Compute closure error

% Display initial error before Newton refinement
fprintf('Initial periodic orbit closure error: %.15e\n', error0);

% =========================================================================
% SECTION 6: MANIFOLD COMPUTATION AND VISUALIZATION
% =========================================================================
% Compute and plot stable/unstable manifolds of the periodic orbit

MapManifold_Billiard(coef, @f, @Df, numOrbits, numIterates, bounds, ...
    xplotinit, X0, per, scale_u, scale_s, numits, tablenum, level);

% =========================================================================
% HELPER FUNCTIONS
% =========================================================================

function output = f(V, coef, guess)
    % F - Billiard map forward iteration
    %
    % Evaluates the billiard map at point V with table coefficients coef.
    % Optional third argument 'guess' can be provided for optimization.
    %
    % INPUTS:
    %   V     : 2x1 vector [x; y] representing point on billiard table
    %   coef  : 2xN matrix of Fourier coefficients for boundary
    %   guess : (optional) initial guess for internal computations
    %
    % OUTPUTS:
    %   output : 2x3 matrix where first column is f(V)
    
    if ~exist('guess', 'var') 
        [output, ~] = BilliardMap(V, coef);
    else 
        [output, ~] = BilliardMap(V, coef, guess);
    end 
end

function Doutput = Df(V, coef, guess)
    % DF - Billiard map derivative (Jacobian)
    %
    % Computes the derivative of the billiard map at point V.
    % Returns 2x2 Jacobian matrix of partial derivatives.
    %
    % INPUTS:
    %   V     : 2x1 vector [x; y] representing point on billiard table
    %   coef  : 2xN matrix of Fourier coefficients for boundary
    %   guess : (optional) initial guess for internal computations
    %
    % OUTPUTS:
    %   Doutput : 2x2 Jacobian matrix [df_x/dx, df_x/dy; df_y/dx, df_y/dy]
    
    if ~exist('guess', 'var') 
        [~, Doutput] = BilliardMap(V, coef);
    else 
        [~, Doutput] = BilliardMap(V, coef, guess);
    end 
end

% =========================================================================
% END OF SCRIPT
% =========================================================================