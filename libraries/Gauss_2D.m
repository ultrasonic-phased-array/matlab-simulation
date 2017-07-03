function p=Gauss_2D(b, f, c,  x, z)
% p = Gauss_2D(b,f,c,x,z) calculates the normalized pressure at a 
% point (x, z) (in mm)in a fluid whose wave speed is c (in m/sec)
% for a 1-D element of length 2b (in mm) radiating at a frequency, f, 
% (in MHz). The function uses a paraxial multi-Gaussian beam model and 
% 15 Gaussian coefficients developed by Wen and Breazeale that are
% contained in the MATLAB function gauss_c15.

% retrieve Wen and Breazeale coefficients
[A, B] = gauss_c15;

% calculate the wave number
kb = 2000*pi*f*b/c;

%normalize the (x,z) coordinates
xb = x/b;
zb = z/b;

%initialize the pressure to zero and then superimpose 15
%Gaussian beams to calculate the pressure wave field
p=0;
for nn = 1:15
    qb=zb-i*1000*pi*f*b./(B(nn)*c);
    qb0 = -i*1000*pi*f*b./(B(nn)*c);
    p=p+sqrt(qb0./qb).*A(nn).*exp(i*kb*xb.^2./(2*qb));
end