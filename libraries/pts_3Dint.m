function xi = pts_3Dint(ex, ey, xn, yn, angt, Dt0, c1, c2, x, y, z)
% xi = pts3Dint(ex, ey, xn,yn,angt,Dt0, c1,c2,x,y,z) calculates the
% distance, xi, (in mm) along the interface in the plane of incidence, 
% at which a ray from the center of an array element segment to a point 
% in the second medium intersects the interface. The parameters
% (ex, ey) are the element offsets (in mm)from the center of the 
% entire array to the center of the element in the x'- and y'-directions,
% respectively, and (xn,yn)are similarly the offsets as measured
% to the center of the element segment from the center of the
% element in the x'- and y'-directions (in mm). The parameter angt, is the
% angle of the array (in degrees) from the interface, and Dt0 is the
% distance (in mm) of the center of the array from the interface. (c1,c2)
% are the wave speeds in the first and second medium (in m/sec) and 
% (x,y,z) are the coordinates of end point of the ray in the second medium
% (all in mm). 

%calculate wave speed ratio
cr=c1/c2;
% determine size of array needed for xi calculations based on the sizes of
% the (x,y,z) variables) and also determine those sizes
[xi, P, Q ] = init_xi3D(x,y,z);  

[nrx,ncx] =size(x);
[nry, ncy] =size(y);
[nrz,ncz] =size(z);


% call ferrari2 function to compute xi with the arguments of that function
% determined by the sizes of the (x,y,z) variables.
De = Dt0 +(ex + xn)*sind(angt); 
for pp=1:P
    for qq = 1:Q
        
        % x and y are points, z is a row or column vector
        if nrx ==1 && ncx ==1 && nry ==1 && ncy ==1
            Db=sqrt((x-(ex +xn)*cosd(angt)).^2 +(y-(ey+yn)).^2);
            xi(pp,qq) =ferrari2(cr, z(pp,qq), De, Db);
        % y and z are points, x is a row or column vector
        elseif nry == 1 && ncy ==1 && nrz ==1 && ncz ==1
            Db=sqrt((x(pp,qq)-(ex +xn)*cosd(angt)).^2 +(y-(ey+yn)).^2) ;           
            xi(pp,qq) =ferrari2(cr, z, De, Db);
        % x and z are points, y is a row or column vector
        elseif nrx ==1 && ncx ==1 && nrz ==1 && ncz ==1
            Db=sqrt((x-(ex +xn)*cosd(angt)).^2 +(y(pp,qq)-(ey+yn)).^2);
            xi(pp,qq) =ferrari2(cr, z, De, Db);
        % y is a point, x and z are equal size PxQ matrices
        elseif nry ==1 && ncy ==1  && nrx == nrz && ncx == ncz
            Db=sqrt((x(pp,qq)-(ex +xn)*cosd(angt)).^2 +(y-(ey+yn)).^2);
            xi(pp,qq) = ferrari2(cr, z(pp,qq), De, Db);
            
          
        % z is a point, x and y are equal size PxQ matrices
        elseif nrz == 1 && ncz ==1 && nrx == nry && ncx == ncy
            Db=sqrt((x(pp,qq)-(ex +xn)*cosd(angt)).^2 +(y(pp,qq)-(ey+yn)).^2);
            xi(pp,qq) = ferrari2(cr, z, De, Db);
        % x is a point, y and z are equal size PxQ matrices
        elseif nrx ==1 && ncx ==1 && nry == nrz && ncy == ncz
            Db=sqrt((x-(ex +xn)*cosd(angt)).^2 +(y(pp,qq)-(ey+yn)).^2);
            xi(pp,qq) = ferrari2(cr, z(pp,qq), De, Db);
        % x, y, z are all equal size row or column vectors
        else
            Db=sqrt((x(pp,qq)-(ex +xn)*cosd(angt)).^2 +(y(pp,qq)-(ey+yn)).^2);
            xi(pp,qq) = ferrari2(cr, z(pp,qq),De, Db);
        end
    end
end