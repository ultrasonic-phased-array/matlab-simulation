function y=fresnel_int(x)
% y = fresnel_int(x) computes the Fresnel integral defined as the integral
% from t = 0 to t = x of the function exp(i*pi*t^2/2). Uses the approximate
% expressions given by Abramowitz and Stegun, Handbook of Mathematical 
% Functions, Dover Publications, 1965, pp. 301-302.


%separate arguments into positive and negative values, change sign 
%of the negative values
xn =-x(x<0);      
xp=x(x >=0);

%compute cosine and sine integrals of the negative values, using the
%oddness property of the cosine and sign integrals

[cn,sn] =cs_int(xn);
cn= -cn;
sn = -sn;

%compute cosine and sine integrals of the positive values

[cp, sp]=cs_int(xp);

%combine cosine and sine integrals for positive and negative
%values and return the complex Fresnel integral
ct =[cn cp];
st =[sn sp];
y=ct+i*st;

%cs_int(xi) calculates approximations of the cosine and sine integrals
%for positive values of xi only(see Abramowitz and Stegun reference above) 
function [c, s]=cs_int(xi)
f =(1+0.926.*xi)./(2+1.792.*xi +3.104.*xi.^2);      % f function (see ref.)
g=1./(2+4.142.*xi+3.492.*xi.^2+6.67.*xi.^3);        % g function (see ref.)
c=0.5 +f.*sin(pi.*xi.^2./2) -g.*cos(pi.*xi.^2./2);  % cos integral approx.
s = 0.5 -f.*cos(pi.*xi.^2./2)-g.*sin(pi.*xi.^2./2); % sin integral approx.
