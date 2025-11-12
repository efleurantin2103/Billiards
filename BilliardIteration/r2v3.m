function [v,angles] = r2v3(r,dz)
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

gamma = acos(r);
tau = atan2(dz(2),dz(1));

angles = [gamma;tau];
v = [cos(gamma+tau);sin(gamma+tau)];

end




