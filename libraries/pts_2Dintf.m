function xi = pts_2Dintf( e, xn, angt, Dt0, c1,c2, x, z)
% xi = pts_2Dintf(e, xn, angt, Dt0, c1, c2, x, z) calculates the
% intersection of a ray from the center of a segment of an array element in
% one fluid to a point (x, z) (in mm) in a second fluid across a plane 
% interface, where e is the offset of the element from the center of the
% array and xn is the offset of the segment from the center of the element.
% (both in mm). The parameter angt is the angle (in degrees) that the array
% makes with respect to the x-axis (the interface) and Dt0 is the distance
% of the center of the array above the interface (in mm). The parameters
% c1, c2 are the wave speeds in the first and second medium, respectively, 
% (both in m/sec). This function uses the function init_xi(x,z) to examine
% the sizes of the (x,z) variables to decide on the corresponding number
% of rows and columns needed to calculate the locations xi (in mm) at
% which rays from the center of a segment to the points (x,z) intersect the
% interface. The function ferrari2 is then used with the appropriate input
% arguments to calculate the xi values (in mm). 


% calculate wave speed ratio
cr =c1/c2;

% based on sizes of (x, z), determine corresponding number of rows and
%columns (P,Q) needed for xi calculations and initialize xi as zeros.
[xi,P,Q]=init_xi(x,z);

% obtain sizes of (x,z) so appropriate arguments can be found in the calls
% to the function ferrari2 when making the xi calculations
[nrx, ncx] =size(x);
[nrz,ncz]=size(z);

% calculate xi locations using ferrari's method
for pp=1:P
    for qq=1:Q
    Dtn=Dt0+(e+xn)*sin(angt*pi/180);
    % if x is a point,and z is a row or column vector
    if nrx ==1 && ncx == 1
    Dxn= x -(e+xn)*cos(angt*pi/180);
    xi(pp,qq)=ferrari2(cr, z(pp,qq), Dtn,Dxn);
    % if z is a point, and x is a row or column vector
    elseif nrz ==1 && ncz ==1
     Dxn = x(pp,qq) -(e+xn)*cos(angt*pi/180);   
    xi(pp,qq)=ferrari2(cr, z, Dtn,Dxn);
    % if x and z are equal size PxQ matrices
    else
     Dxn = x(pp,qq) -(e+xn)*cos(angt*pi/180);   
    xi(pp,qq)=ferrari2(cr, z(pp,qq), Dtn,Dxn);        
    end
    
    end
end
