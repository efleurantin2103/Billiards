function [output,Df] = BilliardMap(V,coef,guess);

% function [output,Df] = BilliardMap(V,coef,[guess]);
% This function inputs 
%	V  = point
% coef = coefficients, 
% guess = optional guess  
% outputs 
% output = the billiard map point (2 vector)
% Df = derivative
% It uses the real version of billiards map when no guess
% and uses the complex implicit version when there is a guess.  

coefx = coef(1,:);
coefy = coef(2,:); 

if ~exist('guess','var')
	[output1,Df] = RealF3(V,[coefx;coefy]);
    output = output1(:,1);
else
	[output1,Df] = ComplexF3(V,[coefx;coefy],guess);
    output = output1(:,1);
end 





