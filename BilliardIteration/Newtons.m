%Multivariable Newton's Method
function [xfin,stop,k] = Newtons(F,dF,guess)  
%% INPUT 
%       F = function that you wish to find the zero for.
%       DF = derivative of F.
%       guess = initial guess for Newton's Method
%
%% OUPUT
%       xfin = zero for F. 
%       stop = the norm of the difference between the last and second to
%       last step. 
%       k = Iterate of Newton's Method at which the method is terminated. 

    x = guess;
    stop = 1;
    k = 1;

    while (stop > 10^(-14)) && (k<1000)
		diff = -dF(x)\F(x);
        x = x + diff; 
        stop = norm(diff,'inf');
        k = k+1;
     end
    
	extraits = 2;
    for j = 1:extraits
        diff = -1*dF(x)\F(x);
        x = x + diff; 
        stop = norm(diff,'inf');
        k = k+1;
     end
  
     xfin = x;
     
end