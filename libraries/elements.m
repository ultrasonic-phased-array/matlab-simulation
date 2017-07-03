function [ A, d, g, xc]=elements(f, c, dl, gd, N)
% [A,d,g,xc]=elements(f,c,dl,gd,N) calculates the 
% total length of an array,A (in mm),the element size,d=2b, 
% (in mm), the gap size, g, (in mm) and the location of the
% centroids of the array elements, xc, (in mm) for an array
% with N elements. The imputs are the frequency, f, (in MHz)
% the wave speed, c, (in m/sec), the element length divided 
% by the wavelength, dl, the gap size divided by the element
% length,gd, and the number of elements, N.

% dl is the element diameter, d, divided by the 
% wavelength,l, i.e. dl =d/l.
d=dl.*c./(1000*f);
%gd is the gap size, g, between elements as a fraction of the 
%element size, i.e. gd =g/d
g=gd.*d;
% A is the total aperture size of the array 
A = N*d + (N-1)*g;
% x= xc is the location of the centroid of each element
% where x = 0 is at the center of the array 
for nn = 1:N
    xc(nn) = (g+d)*((2*nn -1)/2 - N/2);
end
