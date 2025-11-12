close all
clear all

% Configuration
numvals = 100;
numits = 2000;

% Table selection (uncomment one)
coefx = [1.1, 0.03];         coefy = [1, 0.03];           tablenum = 0;
% coefx = [2, 0.04];           coefy = [1, 0.035];          tablenum = 1;
% coefx = [1.1, 0.08, 0.0002]; coefy = [1, 0.095, 0.0001];  tablenum = 2;
% coefx = [1.1, 0.05, 0.00015]; coefy = [1, 0.035, 0.0001]; tablenum = 3;
% coefx = [2, 0.05];           coefy = [1, 0.065];          tablenum = 4;

coeff = [coefx; coefy];
rvec = linspace(-0.99, 0.99, numvals);
thetvec = [0, 0.3, 0.6, 0.9];

% Initialize arrays
rfinal = zeros(length(thetvec)*numvals, numits+1);
thfinal = zeros(length(thetvec)*numvals, numits+1);
rotfinal = zeros(length(thetvec)*numvals, numits+1);
digfinal = zeros(length(thetvec)*numvals, numits+1);

% Compute orbits
figure;
for i1 = 1:length(thetvec)
    theta = thetvec(i1);
    IJ = (i1-1)*numvals;
    
    for j = 1:numvals
        rit = zeros(numits+1, 1);
        thetait = zeros(numits+1, 1);
        rit(1) = rvec(j);
        thetait(1) = theta;
        
        for h = 1:numits
            output = RealF3([rit(h), thetait(h)], coeff);
            rit(h+1) = output(1, 1);
            thetait(h+1) = mod(output(2), 1);
        end
        
        thfinal(j+IJ, :) = thetait';
        rfinal(j+IJ, :) = rit';
        
        [rotn, diggn] = Extra_info2(thfinal(j+IJ, 1:end-1), numits);
        digfinal(j+IJ, :) = ((1-isnan(diggn))*diggn) * ones(1, numits+1);
        rotfinal(j+IJ, :) = min(mod(rotn, 1), 1) * ones(1, numits+1);
    end
end

% Prepare data for plotting
rplot = reshape(rfinal, [], 1);
thplot = reshape(thfinal, [], 1);
rotplot = reshape(rotfinal, [], 1);
digplot = reshape(digfinal, [], 1);
thplot2 = mod(thplot, 1);

% Plot with low alpha for transparency
trust = (digplot > 4.2);
dontrust = ~trust;

cmap = colormap(hsv);
cmap(1, :) = [0.9 0.9 0.9];

scatter(thplot2(trust), rplot(trust), 1, rotplot(trust), ...
    'MarkerFaceAlpha', 0.01, 'MarkerEdgeAlpha', 0.01); 
hold on
scatter(thplot2(dontrust), rplot(dontrust), 1, [0.9 0.9 0.9]);

caxis([0, 1])
colormap(cmap)
colorbar;
xlabel('\theta', 'FontSize', 20)
ylabel('r', 'FontSize', 20)

% Save figure
filename = ['phaseplaneLOWALPHA', num2str(tablenum)];
savefig(1, [filename, '.fig'])