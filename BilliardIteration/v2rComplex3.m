function [r,angles] = v2rComplex3(v,dz,inputguess);
%%INPUT 
%       v = trajectory of billiard ball.
%       dz = tangent line at the point on the boundary. 
%       inputguess = 2x1 vector with [gamma_guess; tau_guess].
%%OUTPUT 
%       r = cos(gamma_new), cosine of angle made between the vector v and the tangent
%       of the boundary at the new point of contact. 
%       angles = 2x1 vector with first values corresponding to the angle
%       gamma and the second value corresponding to the angle tau and the
%       third angle corresponding to rho. 

gammaguess = inputguess(1);
tauguess = inputguess(2);

findtau = @(tau) dz(2)*cos(tau) - dz(1)*sin(tau);
dfindtau = @(tau) -dz(2)*sin(tau) - dz(1)*cos(tau);
findgammaptau = @(gammaptau) v(2)*cos(gammaptau) - v(1)*sin(gammaptau);
dfindgammaptau = @(gammaptau) -v(2)*sin(gammaptau) - v(1)*cos(gammaptau);

[gammaptau,stopv,~] = Newtons(findgammaptau,dfindgammaptau,gammaguess+tauguess);
[tau,stopv,~] = Newtons(findtau,dfindtau,tauguess);

angles = [gammaptau-tau;tau;gammaptau];
r = cos(gammaptau-tau);


end

