function output = RealF3INV(input,coeff)
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

newinput = input;
newinput(1) = -newinput(1);

output = RealF3(newinput,coeff);

output(1,1) = -output(1,1);