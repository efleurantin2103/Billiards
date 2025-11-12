function pv = priorv(v,dz)
% This function is used primarily for the nextstepell function only to help
% acquire a guess for finding the next point on the boundary. 

%%INPUT 
%       v = trajectory of ball
%       dz = tangent at point of origin on boundary. 
%%OUPUT 
%       pv = vector of ball leading up to point of origin (previous
%       trajectory before hitting first point). 

    n = [-dz(2);dz(1)];
    n = n/norm(n);
    v2 = v-2*dot(n,v)*n;
    pv = v2/norm(v2);
