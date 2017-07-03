function [xi, P, Q] = init_xi(x,z)
% [xi,P,Q] =init_xi(x,z) examines the points(x,z), where x can be a row
% or column vector and z a scalar, or z a row or column vector and x a
% scalar, or both x and z can be equal sized scalars, vectors or matrices. 
% The dimensions (P,Q) of xi are chosen accordingly so that calls to 
% functions of (x, z, xi) can be made transparently and consistently 
% if one is evaluating that function along an axis, a line, or over a 2-D 
% array of points. An empty xi matrix of dimensions PxQ is returned, along
% with the dimensions P and Q.

% get sizes of x and z variables
[nrx, ncx] =size(x);
[nrz, ncz] = size(z);

% if x and z are equal sized  matrices,vectors,or scalars, xi is of the 
% same size
if nrx == nrz && ncx ==ncz
    xi=zeros(nrx, ncx);
    P=nrx;
    Q=ncx;
% if x is a column vector and z a scalar, xi is the same size column vector
elseif nrx > 1 && ncx ==1 && nrz ==1 && ncz ==1
    xi=zeros(nrx,1);
    P=nrx;
    Q=1;
% if z is a column vector and x a scalar, xi is the same size column vector
elseif nrz >1 && ncz == 1 && nrx ==1 && ncx ==1
    xi =zeros(nrz,1);
    P=nrz;
    Q=1;
% if x is a row vector and z a scalar, xi is the same size row vector
elseif ncx > 1 && nrx ==1 && nrz ==1 && ncz ==1
    xi=zeros(1, ncx);
    P=1;
    Q=ncx;
% if z is a row vector and x a scalar, xi is the same size row  vector
elseif ncz > 1 && nrz ==1 && nrx ==1 && ncx ==1
    xi=zeros(1,ncz);
    P=1;
    Q=ncz;
% other combinations are not supported 
else error('(x,z) must be (vector,scalar) pairs or equal matrices')
end
