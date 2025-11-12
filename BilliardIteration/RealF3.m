function [output,DF] = RealF3(input,coeff)
%% Input 
%        input = 2x1 vector [r1;th1] where r1 is the cosine of the angle made
%        between the balls trajectory and the tangent line at the point of
%        origin.
%        coeff = 2xn coefficient matrix with the first row representing the
%        coefficients of the x coordinate and the second row representing
%        the coefficients of the y coordinate for the billiard table. 
%% Output
%        output = 2x3 vector with first column [r2;th2] and second and
%        third column as gamma and tau valus fro r2v and v2r. 
%        DF = 2x2 Jacobian matrix representing d(r2,th2)/(d(r,th))
%
%% Example: [output,DF] = RealF3([0.1;0.3],[1.1,0.03;1,0.03])

th1 = input(2);
r1   = input(1);
output = zeros(2,3);
flag = 0; 

%% First Step: Find v given r.
[z1,dz1] = Bill_Table(th1,coeff);
[v1,angles1] = r2v3(r1,dz1/norm(dz1));
output(:,2) = angles1;
% Need Prior v for use in ellipse guess.
v = priorv(v1,dz1);

%% Second Step: Find th2 by finding intersection of v and billiard table.
% Use approximate ellipse as guess for perturbed ellipse case. 
[distguess,thguess] = nextstepell(th1,v,coeff(:,1));
guess = [distguess;thguess];   


L =@(s) z1 + s.*v1;
F =@(x) L(x(1)) - Bill_Table(x(2),coeff);
dF = @(x) [v1,-dBill_Table(x(2),coeff)];

% Use deflation to eliminate possibility of converging to point of origin.     
deflat =  @(x) (1/norm([x(1)-0;x(2)-th1])^2);
ddeflat = @(x) -2*[x(1),x(2)-th1]*deflat(x)^2;	
G = @(x) F(x)*deflat(x);
dG =@(x) F(x)*ddeflat(x)+deflat(x)*dF(x);

[newtval,~,kk] = Newtons(G,dG,guess);
    
    if kk > 999
        flag = 1;
        
    end

th2 = newtval(2);
%% Third Step: Find r2 from v and tangent of new point on table. 
[~,dz2] = Bill_Table(th2,coeff);	
[r2,angles2] = v2r3(v1,dz2); 

output(:,3) = angles2(1:2);
output(:,1) = [r2;th2];

%%%%%%%% DF Calculation %%%%%%%%%

 T = dBill_Table(th1,coeff);
 Tp = ddBill_Table(th1,coeff);
 U = dBill_Table(th2,coeff);
 Up = ddBill_Table(th2,coeff);
 sfinal = newtval(1);
% 
%% First Step    

 dv_dgamma = [-sin(angles1(1) + angles1(2)); cos(angles1(1) + angles1(2))];
 dv_dtau = [-sin(angles1(1) + angles1(2)); cos(angles1(1) + angles1(2))];
 dv_dgam_dtau = [dv_dgamma, dv_dtau];

 dG_drdth = [-1,0; 0, Tp(1)*sin(angles1(2))-Tp(2)*cos(angles1(2))];
 dG_dgam_dtau = [-sin(angles1(1)), 0; 0, T(1)*cos(angles1(2))+T(2)*sin(angles1(2))]; 
 dgamtau_dr_dth = -dG_dgam_dtau\dG_drdth;
 dv = dv_dgam_dtau*dgamtau_dr_dth;

%% Second Newtons method to find (s, theta2): Function H
 DHr = [sfinal*dv(:,1), T + sfinal*dv(:,2)];
 DHs = [v1, -1*U];
 ds = -1*DHs\DHr;
 dthetahat = ds(2,:);

%% Third Step 

 drhat_dgam_dtau = [-sin(angles2(3) - angles2(2)), sin(angles2(3) - angles2(2))];

 dK_dgam_dtau = [0 , U(1)*cos(angles2(2))+U(2)*sin(angles2(2));  v1(1)*cos(angles2(3))+v1(2)*sin(angles2(3)), 0];
 dK_drdth_v = [ 0, 0; dv(1,1)*sin(angles2(3))-dv(2,1)*cos(angles2(3)), dv(1,2)*sin(angles2(3))-dv(2,2)*cos(angles2(3))];

 dK_dth2 = [Up(1)*sin(angles2(2)) - Up(2)*cos(angles2(2)); 0];
 dK_drdth_U = dK_dth2 * dthetahat; 

 dKdrdth = dK_drdth_v + dK_drdth_U;

 dgamtau_dr_dth = -1*dK_dgam_dtau\dKdrdth;
 drhat = drhat_dgam_dtau*dgamtau_dr_dth;

 DF = [drhat; dthetahat]; 
end
%%% 		
function [tell,th2] = nextstepell(th1,v1,coeff)
	%Use the previous v so we can figure out the ellipse bounce instead of the actual bounce
	tpi = 2*pi;
	a = coeff(1);
    b = coeff(2);

	[~,vell] = findbouncevector(th1,v1,coeff);

	dell = vell./[a;b];
	tell = -2 * (dell(1)*cos(tpi*th1)+dell(2)*sin(tpi*th1))/(dell(1)^2+dell(2)^2);
	th2 = mod(atan2(sin(tpi*th1)+tell*dell(2),cos(tpi*th1)+tell*dell(1)),tpi)/tpi;
end
	
%%%% This may be unnecessary 
function [pts,v,r] = findbouncevector(th,v1,coeff)
    % Find normal vector to p2
    [pts,n1] = Bill_Table(th,coeff);
    n1 = n1/norm(n1);
    n = [-n1(2);n1(1)];
    v1 = v1/norm(v1);
    r = dot(v1,n1);
    v2 = v1-2*dot(n,v1)*n;
    v = v2/norm(v2);
end

    % Second derivative of table at theta
function ddz = ddBill_Table(t,coeff)
coefx = coeff(1,:);
coefy = coeff(2,:);
tpi = 2*pi; 
   ddz = zeros(2,1);

    for kk =1: length(coefx)
        ddz(1) = ddz(1) - (tpi*kk)^2*coefx(kk)*cos(kk*tpi*t);
        ddz(2) = ddz(2)- (tpi*kk)^2*coefy(kk)*sin(tpi*kk*t);
    end
end