function [v,angles] = r2vCOMPLEX3(r,dz,inputguess)
%%INPUT 
%       r = cos(gamma), cosine of angle made between the vector v and the tangent
%       line at the boundary .
%       dz = tangent line of boundary at initial point.
%       inputguess = 2x1 vector with [gamma_guess; tau_guess].
%%OUTPUT
%       v = vector in direction of the balls trajectory from initial point 
%       on boundary.  
%       angles = 2x1 vector with first values corresponding to the angle
%       gamma and the second value corresponding to the angle tau. 

gammaguess = inputguess(1);
tauguess = inputguess(2);

findgamma = @(gamma) cos(gamma)-r;
dfindgamma = @(gamma) -sin(gamma);

findtau = @(tau)  dz(2)*cos(tau)-dz(1)*sin(tau);
dfindtau = @(tau) -dz(2)*sin(tau)-dz(1)*cos(tau);


[gamma,stopv,~] = Newtons(findgamma,dfindgamma,gammaguess);
[tau,stopv,~] = Newtons(findtau,dfindtau,tauguess);

angles = [gamma;tau];
v = [cos(gamma+tau);sin(gamma+tau)];


end

