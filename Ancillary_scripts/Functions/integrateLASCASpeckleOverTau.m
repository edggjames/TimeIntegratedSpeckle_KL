function [I_int] = integrateLASCASpeckleOverTau(M_1,M_2,T,t_c,H)
% function to calculate the time integrated expression of a LASCA speckle 
% pattern over a camera exposure period T

% M_1 and M_2 are two uncorrelated 2D complex valued matrices. 

% t_c is the decorrelation time of the sample

% This is equation 95 from Edward James' PhD thesis. 

if nargin == 4
    [a,b,c,d] = calcABCD(M_1,M_2);
elseif nargin == 5
    [a,b,c,d] = calcABCD(M_1,M_2,H);
end

% implement expression

W = T./t_c;
X = exp(-W);
Y = exp(W);
Z = exp(-2*W);

term_1 = (c.^2 + d.^2).*T + t_c*(a.^2 + b.^2 - c.^2 - d.^2);   
term_2 = (-a.^2 - b.^2 + c.^2 + d.^2).*X*t_c;
term_3 = 2*(a.*c + b.*d).*Z*t_c; 
term_4 = -1 + Y + Y .* sqrt(1-Y) .* atanh(sqrt(1-Y));
term_5 = sqrt(Z.*(-1+Y));
I_int = term_1+term_2-(term_3.*term_4./term_5);

I_int = I_int/T;
end

