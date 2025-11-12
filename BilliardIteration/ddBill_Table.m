function [z,dz,ddz] = dBill_Table(t,coefx,coefy)

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