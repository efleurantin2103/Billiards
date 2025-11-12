function [z,dz] = Bill_Table(t,coeff)
%%INPUT 
%       t = parameter value corresponding to point on the table boundary.
%       coeff = 2xn matrix representing coefficients for x and y
%       coordiantes of table. (first row is x coordinates, second row is y
%       coordinates).
%%OUTPUT
%       z = cartesian coordinates for point on the table boundary. 
%       dz  = tangent line at point t on the table. 

    
    coefx = coeff(1,:);
    coefy = coeff(2,:);
    tpi = 2*pi; 
    z = zeros(2,length(t));
    dz = zeros(2,length(t));

    for kk =1: length(coefx)
        
        z(1,:) = z(1,:)+ coefx(kk)*cos(kk*tpi*t);
        z(2,:) = z(2,:) + coefy(kk)*sin(tpi*kk*t);
        dz(1,:) = dz(1,:) + -1*tpi*kk*coefx(kk)*sin(kk*tpi*t);
        dz(2,:) = dz(2,:)+ tpi*kk*coefy(kk)*cos(tpi*kk*t);

    end

    


end