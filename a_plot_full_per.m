% =========================================================================
% BILLIARD MAP MANIFOLD VISUALIZATION
% =========================================================================
%
% DESCRIPTION:
%   This script visualizes stable and unstable manifolds of periodic orbits
%   for billiard maps with perturbed boundaries. It loads precomputed
%   manifold data and overlays it on phase plane plots showing the billiard
%   dynamics in (theta, r) coordinates.
%
% MAIN FEATURES:
%   - Loads and combines multiple periodic orbit manifolds from .mat files
%   - Overlays manifolds on precomputed phase plane diagrams
%   - Color codes: Red = unstable manifolds, Blue = stable manifolds,
%                  Black = periodic orbit points
%   - Supports zooming to specific regions of interest
%   - Exports high-resolution figures (JPEG format)
%
% REQUIRED FILES:
%   - table_X/ directories (X = 0, 1, 2, 3, 4) containing:
%       * phaseplanex.fig : Base phase plane figure
%       * manifoldtab*.mat : Manifold data files with fields:
%           - X : Periodic orbit coordinates [x1; y1; x2; y2; ...]
%           - per : Period of the orbit
%           - all_orbits_u1 : Unstable manifold points
%           - all_orbits_s1 : Stable manifold points
%           - Pu_image : Iterated unstable manifold points
%           - Ps_image : Iterated stable manifold points
%
% OUTPUT:
%   - Saves visualization to: table_X/full_table_X.jpg
%
% COORDINATE SYSTEM:
%   - theta (horizontal axis): Angular coordinate, periodic with period 1
%   - r (vertical axis): Radial coordinate in [-1, 1]
%
% USAGE:
%   1. Set 'tablenum' variable to desired table (0, 1, 2, 3, or 4)
%   2. Optionally adjust axis limits for zooming
%   3. Optionally adjust output resolution (see DPI settings at end)
%   4. Run the script
%
% =========================================================================

% Clean workspace and set numerical display format
clearvars
format long
close all

% Add paths to table directories
addpath('table_0/', 'table_1/', 'table_2/', 'table_3/', 'table_4/')

% =========================================================================
% SECTION 1: CONFIGURATION
% =========================================================================

% Select which billiard table to visualize (0, 1, 2, 3, or 4)
tablenum = 4;  

% =========================================================================
% SECTION 2: LOAD BASE PHASE PLANE FIGURE
% =========================================================================
% Load the precomputed phase plane diagram for the selected table.
% This provides the background showing the overall dynamics.

% Standard phase plane (default)
% fig = openfig(sprintf('phaseplane%d.fig', tablenum));

% Alternative: Low alpha version (uncomment if needed)
fig = openfig(sprintf('phaseplaneLOWALPHA%d.fig', tablenum));

hold on;
set(gca, 'FontSize', 20);

% =========================================================================
% SECTION 3: LOAD AND PLOT MANIFOLD DATA
% =========================================================================
% Scan the table directory for all manifold data files and plot each one.
% Each file contains data for one periodic orbit and its manifolds.

% Construct path to table directory
folder_path = sprintf('table_%d/', tablenum);

% Find all manifold data files (pattern: manifoldtab*.mat)
mat_files = dir(fullfile(folder_path, 'manifoldtab*.mat'));

fprintf('Found %d manifold files in %s\n', length(mat_files), folder_path);

% Loop through each manifold file
for i = 1:length(mat_files)
    
    % Load manifold data
    data_path = fullfile(folder_path, mat_files(i).name);
    data = load(data_path);
    
    fprintf('Processing file %d/%d: %s (Period %d orbit)\n', ...
            i, length(mat_files), mat_files(i).name, data.per);
    
    % =====================================================================
    % Plot periodic orbit points (black dots)
    % =====================================================================
    % Extract x and y coordinates from interleaved array X
    % X format: [x1; y1; x2; y2; ...; x_per; y_per]
    x_coords = data.X(1:2:2*data.per);  % Odd indices: x-coordinates (r)
    y_coords = data.X(2:2:2*data.per);  % Even indices: y-coordinates (theta)
    
    plot(y_coords, x_coords, 'k.', 'MarkerSize', 20);
    
    % =====================================================================
    % Plot unstable manifold - initial segments (red dots)
    % =====================================================================
    % all_orbits_u1(1,:) contains r-coordinates
    % all_orbits_u1(2,:) contains theta-coordinates
    % Apply mod(theta, 1) to wrap theta to [0,1] due to periodicity
    plot(mod(data.all_orbits_u1(2, :), 1), data.all_orbits_u1(1, :), 'r.');
    
    % =====================================================================
    % Plot stable manifold - initial segments (blue dots)
    % =====================================================================
    % all_orbits_s1(1,:) contains r-coordinates
    % all_orbits_s1(2,:) contains theta-coordinates
    % Apply mod(theta, 1) to wrap theta to [0,1] due to periodicity
    plot(mod(data.all_orbits_s1(2, :), 1), data.all_orbits_s1(1, :), 'b.');
    
    % =====================================================================
    % Plot forward iterates of unstable/stable manifolds
    % =====================================================================
    % For each point in the periodic orbit, plot its associated manifold
    % iterates. This shows how manifolds evolve under forward iteration.
    for j = 1:data.per
        % Unstable manifold iterates at j-th periodic point (red)
        % Pu_image organized as [r1; theta1; r2; theta2; ...]
        plot(mod(data.Pu_image(2*j, :), 1), data.Pu_image(2*j-1, :), 'r.');
        
        % Stable manifold iterates at j-th periodic point (blue)
        % Ps_image organized as [r1; theta1; r2; theta2; ...]
        plot(mod(data.Ps_image(2*j, :), 1), data.Ps_image(2*j-1, :), 'b.');
    end
end

% =========================================================================
% SECTION 4: AXIS LABELS AND FORMATTING
% =========================================================================

% Set axis labels with LaTeX formatting
xlabel('$\theta$', 'Interpreter', 'latex', 'FontSize', 20)
ylabel('$r$', 'Interpreter', 'latex', 'FontSize', 20)

% =========================================================================
% SECTION 5: AXIS LIMITS (ZOOM CONTROL)
% =========================================================================
% Uncomment ONE of the following to set the viewing window:

% Full phase space view
  axis([0, 1, -1, 1])

% Zoomed view (default) - focuses on specific region of interest
% axis([0.45, 0.55, -0.1, 0.1]);

% =========================================================================
% SECTION 6: FIGURE SIZE AND RESOLUTION SETTINGS
% =========================================================================

% Set figure dimensions in inches (for print/export)
width_inches = 10;   % Width in inches
height_inches = 8;   % Height in inches

% Configure paper properties for printing/saving
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [width_inches height_inches]);
set(gcf, 'PaperPosition', [0 0 width_inches height_inches]);
set(gcf, 'PaperPositionMode', 'manual');

% Alternative renderer settings (uncomment if needed for compatibility)
% set(gcf, 'Renderer', 'painters');
% set(gcf, 'PaperPositionMode', 'auto');

% =========================================================================
% SECTION 7: SAVE OUTPUT FIGURE
% =========================================================================

% Construct output filename
output_filename = fullfile(folder_path, sprintf('full_table_%d.jpg', tablenum));

% Save figure as high-resolution JPEG
% Resolution guide:
%   -r300 : High quality (publication)  -> 10"×8" = 3000×2400 pixels
%   -r150 : Medium quality (default)    -> 10"×8" = 1500×1200 pixels
%   -r100 : Lower quality (web)         -> 10"×8" = 1000×800 pixels
print(fig, '-djpeg', '-r150', output_filename);

fprintf('\nVisualization saved to: %s\n', output_filename);
fprintf('Image dimensions: %d × %d pixels\n', ...
        round(width_inches*150), round(height_inches*150));

% =========================================================================
% ADDITIONAL NOTES
% =========================================================================
%
% RESOLUTION EXAMPLES:
%   10" × 8"  at 300 DPI (-r300) = 3000 × 2400 pixels (high quality)
%   6"  × 4"  at 300 DPI (-r300) = 1800 × 1200 pixels (medium size, high quality)
%   8"  × 6"  at 150 DPI (-r150) = 1200 × 900 pixels  (medium quality)
%   10" × 8"  at 150 DPI (-r150) = 1500 × 1200 pixels (default)
%
% COLOR CODING:
%   Red (r.)   : Unstable manifolds
%   Blue (b.)  : Stable manifolds
%   Black (k.) : Periodic orbit points
%
% COORDINATE INTERPRETATION:
%   theta : Measures position around boundary (mod 1, i.e., periodic)
%   r     : Measures "momentum" or velocity direction component
%
% TROUBLESHOOTING:
%   - If no manifolds appear: check that .mat files exist in table_X/
%   - If figure looks wrong: verify phaseplanex.fig file is correct version
%   - If colors don't match: ensure data structure fields are as expected
%   - If out of memory: reduce number of manifold points or DPI
%
% =========================================================================
% END OF SCRIPT
% =========================================================================