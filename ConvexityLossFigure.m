
numep = 100;

epvect = linspace(0,10,numep); 

curvemin = zeros(numep,1); 
figure; 
hold on;
xlabel('\fontsize{20} x')
ylabel('\fontsize{20} y')
cmap = flip(winter(10),1);
set(gca(),'ColorOrder',cmap)

for k = 1:numep
	epsilon = epvect(k); 
	curvmin(k) = mincurvature(epsilon); 
end

fzero(@mincurvature2,3)

figure;
plot(epvect,curvmin,'linewidth',2)
xlabel('\fontsize{20} \epsilon')
ylabel('\fontsize{20} Minimum Signed Curvature')

figure(1);
print('-djpeg', '-r300', 'curvature1.jpg')

figure(2);
print('-djpeg', '-r300', 'curvature2.jpg')


function curvmin = mincurvature(epsilon)
	numpoints = 100; 
	t = linspace(0,1,numpoints);
	
	coefx = [1.1,epsilon*0.03];
	coefy = [1,epsilon*0.025];
	[bill,dbill,ddbill] =  ddBill_Table(t,coefx,coefy);	
	curvsign = (dbill(1,:).*ddbill(2,:)-dbill(2,:).*ddbill(1,:))./(dbill(1,:).^2+dbill(2,:).^2).^(3/2);
	curvmin = min(curvsign);
	plot(bill(1,:),bill(2,:),'linewidth',2)
end

function curvmin = mincurvature2(epsilon)
	numpoints = 100; 
	t = linspace(0,1,numpoints);
	
	coefx = [1.1,epsilon*0.03];
	coefy = [1,epsilon*0.03];
	[bill,dbill,ddbill] =  ddBill_Table(t,coefx,coefy);	
	curvsign = (dbill(1,:).*ddbill(2,:)-dbill(2,:).*ddbill(1,:))./(dbill(1,:).^2+dbill(2,:).^2).^(3/2);
	curvmin = min(curvsign);
end

function [z,dz,ddz] = ddBill_Table(t,coefx,coefy)

    tpi = 2*pi; 
    z   = zeros(2,length(t));
    dz  = zeros(2,length(t));
    ddz = zeros(2,length(t));

    for kk =1: length(coefx)        
        z(1,:)  = z(1,:) + coefx(kk)*cos(tpi*kk*t);
        z(2,:)  = z(2,:) + coefy(kk)*sin(tpi*kk*t);
        
        dz(1,:) = dz(1,:) + -1*tpi*kk*coefx(kk)*sin(tpi*kk*t);
        dz(2,:) = dz(2,:) +    tpi*kk*coefy(kk)*cos(tpi*kk*t);
		
		ddz(1,:) = ddz(1,:)+ -1*(tpi*kk)^2*coefx(kk)*cos(tpi*kk*t);
		ddz(2,:) = ddz(2,:)+ -1*(tpi*kk)^2*coefy(kk)*sin(tpi*kk*t);
    end
end