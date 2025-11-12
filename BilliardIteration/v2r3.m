function [r,angles] = v2r3(v,dz)
%%INPUT 
%       v = trajectory of billiard ball. 
%       dz = tangent at point on billiard table
%%OUTPUT
%       r = cosine of angle between tangent and v
%       angles = 3x1 vector with first values corresponding to the angle
%       gamma and the second value corresponding to the angle tau and the
%       third angle corresponding to rho.
tau = atan2(dz(2),dz(1));
gamptau = atan2(v(2),v(1));

angles = [gamptau-tau;tau;gamptau];
r = cos(gamptau-tau);


