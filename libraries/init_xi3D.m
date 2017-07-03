function [xi, P, Q] = init_xi3D(x,y,z)
% [xi, P, Q] = init_xi3D(x,y z) examines the sizes of the (x,y,z) variables
% (which specify points in the second medium, across a plane interface,
% to which a ray must travel from an array element or element segment)
% and returns a PxQ array of zero values to hold the distances xi 
% at which the ray intersects an interface, as well as the values (P,Q).
%  Eleven different combinations of sizes for (x,y,z) are
% allowed, which permits (x,y,z) to represent values in planes parallel
% to the x-,y-,or z-axes (three cases), or values along lines parallel to
%  the x-, y-,or z-axes (six cases since the line could be represented
% as row or column vectors),or values along an inclined line 
% in 3-D (two cases since the line could be represented as row or column 
% vectors). 

% get sizes of (x,y,z)
[nrx,ncx ]= size(x);
[nry,ncy]=size(y);
[nrz,ncz] = size(z);

% if x,z are equal size [nrx,ncx] matrices and y is single value, make 
% xi a [nrx, ncx] matrix
if nrx == nrz && ncx == ncz && nry ==1 && ncy ==1
    xi = zeros(nrx,ncx);
    P = nrx;
    Q = ncx;
% if x, y are equal size [nrx, ncx] matrices and z is a single value, make
% xi a [ nrx, ncx] matrix
elseif nrx == nry && ncx == ncy  && nrz ==1 && ncz ==1
    xi = zeros(nrx,ncx);
    P = nrx;
    Q = ncx;
% if y, z are equal size [nry,ncy] matrices and  x is a single value,  make 
% xi a [nry, ncy] matrix
elseif nry ==nrz && ncy == ncz  && nrx ==1 && ncx ==1
    xi=zeros(nry, ncy);
    P = nry;
    Q = ncy;
% if z is a [1,ncz] vector and x and y are single values, make
% xi a [1,ncz] vector
elseif nrz ==1 && ncz >1  && nrx ==1 && ncx ==1 && nry == 1  && ncy ==1
    xi =zeros(1, ncz);
    P = 1;
    Q = ncz;
% if z is a [nrz, 1] vector and x and y are single values, make
% xi a [nrz,1] vector
elseif ncz ==1 && nrz >1  && nrx ==1 && ncx == 1 && nry == 1 && ncy == 1
    xi =zeros(nrz,1);
    P= nrz;
    Q =1;
% if x is a [1,ncx] vector and y and z are single values, make
% xi a [1,ncx] vector
elseif nrx ==1 && ncx >1  && nry ==1 && ncy == 1 && nrz == 1 && ncz == 1
    xi =zeros(1,ncx);
    P= 1;
    Q = ncx;
% if x is a [nrx, 1] vector and  y and z are single values, make
% xi a [nrx, 1] vector
elseif ncx == 1 && nrx >1 && nry ==1 && ncy ==1 && nrz == 1 && ncz == 1
    xi = zeros(nrx, 1);
    P= nrx;
    Q =1;
% if y is a [1, ncy] vector and x and z are single values, make
% xi a [1, ncy] vector
elseif nry ==1 && ncy >1 && nrx ==1 && ncx == 1 && nrz == 1 && ncz == 1
    xi=zeros(1, ncy);
    P =1;
    Q = ncy;
 % if y is a [nry, 1] vector and x and z are single values, mke
 % xi a [nry,1] vector
elseif nry >1 && ncy ==1 && nrx ==1 && ncx ==1 && nrz ==1 && ncz ==1
    xi=zeros(nry, 1);
    P = nry;
    Q =1;
% if x, y, z are equal size [1, ncx] vectors, make
% xi a [ 1,ncx] vector
elseif nrx ==nry && ncx == ncy  && nrz == nrx &&  ncz == ncx && nrx ==1
    xi =zeros(1, ncx);
    P= 1;
    Q = ncx;
% if x, y, z are equal size [nrx,1] vectors, make
% xi a [nrx,1] vector
elseif nrx ==nry && ncx == ncy  && nrz == nrx &&  ncz == ncx && ncx ==1
    xi =zeros(nrx,1);
    P = nrx;
    Q = 1;
else error(' (x,y,z) combination given is not supported')
end
    
    