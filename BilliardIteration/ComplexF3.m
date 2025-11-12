%Patrick Bishop, Jay Mireles James, Evelyn Sander 2024
function [output,DF] = ComplexF3(input,coeff,outputguess) 
%% Input: 
%        input = 2x1 vector [r1;th1] where r1 is the cosine of the angle made
%        between the balls trajectory and the tangent line at the point of
%        origin.
%        coeff = coefficient matrix for Bill_Table
%
%        outputguess = [r2guess, gamma1guess,gamma2guess; theta2guess, tau1guess, tau2guess
%			where 
%        [r2guess;th2guess] is the real approximation found by RealF. 
%        [gam1;tau1] is the guess for r2vcomplex.
%        [gam2; tau2] is the guess for v2rcomplex
%
%% OUTPUT
%       output = 2x3 vector with first column [r2;th2], second and third
%       column are currently the gamm and tau values for the first and
%       third step.
%       DF = 2x2 Jacobian matrix representing d(r2,th2)/(d(r,th))
%
%% EXAMPLE:  [output,DF] = ComplexF3([0.1;0.3],[1.1,0.03;1,0.03],[0.178453273522060,1.470628905633337, -1.330992626661454; 0.769948476908876, -2.801621532294791, 0.060389432633410])
%% EXAMPLE:  [output,DF] = ComplexF3([r1;theta1],[coefx;coefy],[r2guess, gamma1guess,gamma2guess; theta2guess, tau1guess, tau2guess])


%% Add tau and gamma into outputguess
r1 = input(1);
th1 = input(2);


r2guess = outputguess(1,1);
th2guess = outputguess(2,1);
output = zeros(2,3);
angleguess1 = outputguess(:,2);
angleguess2 = outputguess(:,3);


[pts,dz1] = Bill_Table(th1,coeff);
[guesspts] = Bill_Table(th2guess,coeff);

[v,angles1] = r2vCOMPLEX3(r1,dz1,angleguess1);
output(:,2) = angles1;

distguess =  norm(guesspts-pts); 
guess = [distguess;th2guess];    
    
L =@(s) pts + s*v;
F =@(x) L(x(1)) - Bill_Table(x(2),coeff);
dF = @(x) [v,-dBill_Table(x(2),coeff)];
       
[newtval,~,kk] = Newtons(F,dF,guess);

sfinal = newtval(1);

if kk > 999
    flag = 1
end
    
th2 = newtval(2);
[~,dz2] = Bill_Table(th2,coeff);	
[r2,angles2] = v2rComplex3(v,dz2,angleguess2);

output(:,3) = angles2(1:2);
output(:,1) = [r2;th2];


%%%%%%%% DF Calculation %%%%%%%%%

 T = dBill_Table(th1,coeff);
 Tp = ddBill_Table(th1,coeff);
 U = dBill_Table(th2,coeff);
 Up = ddBill_Table(th2,coeff);
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
 DHs = [v, -1*U];
 ds = -1*DHs\DHr;
 dthetahat = ds(2,:);

%% Third Step 

 drhat_dgam_dtau = [-sin(angles2(3) - angles2(2)), sin(angles2(3) - angles2(2))];

 dK_dgam_dtau = [0 , U(1)*cos(angles2(2))+U(2)*sin(angles2(2));  v(1)*cos(angles2(3))+v(2)*sin(angles2(3)), 0];
 dK_drdth_v = [ 0, 0; dv(1,1)*sin(angles2(3))-dv(2,1)*cos(angles2(3)), dv(1,2)*sin(angles2(3))-dv(2,2)*cos(angles2(3))];
 dK_dth2 = [Up(1)*sin(angles2(2)) - Up(2)*cos(angles2(2)); 0];
 dK_drdth_U = dK_dth2 * dthetahat; 

 dKdrdth = dK_drdth_v + dK_drdth_U;

 dgamtau_dr_dth = -1*dK_dgam_dtau\dKdrdth;
 drhat = drhat_dgam_dtau*dgamtau_dr_dth;

 DF = [drhat; dthetahat]; 
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

