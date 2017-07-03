function p = fresnel_2D(b, f, c, x, z)
% p = fresnel_2D(b, f, c, x, z) calculates the normalized pressure 
% field at a point (x, z), (in mm), of a 1-D element of 
% length 2b(in mm), at a frequency, f,(in MHz)radiating
% into a fluid with wave speed, c, (in m/sec). This function uses the 
% fresnel_int function to calculate the Fresnel integral numerically.

% calculate wave number
kb =2000*pi*f*b/c;

% put (x, z) coordinates in normalized form 
xb=x/b;
zb=z/b;
% calculate term in Fresnel integral argument
arg = sqrt(kb./(pi*zb));

%calculate normalized pressure
p=sqrt(1/(2*i)).*exp(i*kb*zb).*(fresnel_int(arg.*(xb+1))...
    -fresnel_int(arg.*(xb -1)));