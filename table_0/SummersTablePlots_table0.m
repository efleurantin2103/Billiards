coefx = [1.1,0.03]; coefy = [1,0.03]; %table 0 
 coeff = [coefx;coefy];
% 
% rvec = 0.3;
% thetvec = linspace(0,0.1,10);
% 
% numvals = 1;
% numits = 2; %%old info

%Summer's choice:
rvec = 0.3;
thetvec = (0.05);

numvals = 1;
numits = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%

plotthet = linspace(0,1,1000);
table_fn = Bill_Table(plotthet,coeff);
figure;
plot(table_fn(1,:),table_fn(2,:),'k'); hold on
axis("equal")
xlabel('\fontsize{20} X')
ylabel('\fontsize{20} Y')
rfinal = zeros(length(thetvec)*numvals,numits+1);
thfinal= zeros(length(thetvec)*numvals,numits+1);
rotfinal = zeros(length(thetvec)*numvals,numits+1);
digfinal = zeros(length(thetvec)*numvals,numits+1);


for i1 = 1:length(thetvec)
    theta = thetvec(i1);
    IJ = (i1-1)*numvals;
for j = 1:numvals
    rit = zeros(numits,1);
    thetait = zeros(numits,1);
    rit(1) = rvec(j); 
	thetait(1) = theta; 

	for k = 1:numits
		output = RealF3([rit(k);thetait(k)],coeff);
		rit(k+1) = output(1,1);
		thetait(k+1) = output(2,1);
    end
    coords = Bill_Table(thetait(1:end)',coeff);
    %plot(coords(1,:),coords(2,:),'b'); hold on
    plot(coords(1,:),coords(2,:)); hold on
   


    thfinal(j+IJ,:) = thetait(1:end)';
    rfinal(j+IJ,:) = rit(1:end)';

end


end


rplot = reshape(rfinal,length(thetvec)*(numits+1)*numvals,1);
thplot = reshape(thfinal,length(thetvec)*(numits+1)*numvals,1);
