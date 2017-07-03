function xi= ferrari2(cr,DF,DT,DX)
% xi = ferrari2(cr, DF, DT, DX) solves for the intersection point, xi, on
% a plane interface along a Snells law ray path from a point located a
% distance DT (in mm)above the interface to a point located a 
% distance DF (in mm)below the interface.
%  Both DT and DF must be positive. DX (in mm)is the separation 
% distance between the points along the plane interface and can be positive
% or negative. cr = c1/c2 is the ratio of the wave speed in medium 1 to
% that of the wavespeed in medium 2.
% The intersection point,xi, is obtained by writing Snells law as a quartic
% equation in xi and solving the quartic with Ferrari's method. Of the
% four roots, two will be complex, one will be the wanted real solution 
% in the interval [0,DX] and one will be real but outside that interval. 
% reference: http://exampleproblems.com/wiki/index.php/Quartic_equation
% If the root returned by Ferrari's method lies inside the permissable
% interval, [0, DX], and is essentially real(set by a tolerance value 
% in line 76), the solution obtained by Ferrari's method is used.
% Otherwise, the MATLAB function fzero is used instead to find the 
% intersection point.

% if two media are identical, use explicit solution for the interface point
% along a straight ray
if abs(cr-1) < 10^(-6)
    xi = DX*DT/(DF+DT);
%otherwise, use Ferrari's method
else   
 cri=1/cr; % cri = c2/c1
 %define coefficients of quartic Ax^4 +Bx^3 +Cx^2 + Dx + E =0
 A = 1-cri^2;
 B = (2*(cri)^2*DX -2*DX)/DT;
 C = (DX^2 +DT^2 -(cri)^2*(DX^2 +DF^2))/(DT^2);
 D = -2*DX*DT^2/(DT^3);
 E= DX^2*DT^2/(DT^4);
%  begin Ferrari's solution
alpha = -3*B^2/(8*A^2) + C/A;
beta  = B^3/(8*A^3) - B*C/(2*A^2) + D/A;
gamma = -3*B^4/(256*A^4) + C*B^2/(16*A^3) - B*D/(4*A^2) + E/A;
% if beta =0 the quartic is a bi-quadratic whose solution is easier
if(beta == 0)
x(1) = -B/(4*A) + sqrt( (-alpha + sqrt(alpha^2-4*gamma))/2); 
x(2) = -B/(4*A) + sqrt( (-alpha - sqrt(alpha^2-4*gamma))/2);
x(3) = -B/(4*A) - sqrt( (-alpha + sqrt(alpha^2-4*gamma))/2);
x(4) = -B/(4*A) - sqrt( (-alpha - sqrt(alpha^2-4*gamma))/2);
% otherwise, proceed with Ferrari's method
else

P= -alpha^2/12 - gamma;
Q= -alpha^3/108 + alpha*gamma/3 - beta^2/8;
%

Rm=  Q/2 - sqrt(Q^2/4 + P^3/27);
%
U=Rm^(1/3);
%
if(U == 0)
    y=-5/6*alpha - U;
else
    y=-5/6*alpha - U + P/(3*U);
end
%
W=sqrt(alpha + 2*y );
%
x(1) = -B/(4*A) + 0.5*( + W + sqrt(-(3*alpha + 2*y + 2*beta/W )));
x(2) = -B/(4*A) + 0.5*( - W + sqrt(-(3*alpha + 2*y - 2*beta/W )));
x(3) = -B/(4*A) + 0.5*( + W - sqrt(-(3*alpha + 2*y + 2*beta/W )));
x(4) = -B/(4*A) + 0.5*( - W - sqrt(-(3*alpha + 2*y - 2*beta/W )));
end
% end of bi-quadratic solution or ferrari method with four roots

% find root that is real,lies in the interval [0, DX]
flag =0;
 for nn=1:4
     xr=real(x(nn));
     axi= DT*abs(imag(x(nn)));
     xt=xr*DT;
     tol = 10^(-6);
     % f is a function which should be zero if Snell's law
     % is satisfied and can also be used to check the 
     % accuracy of Ferrari's solution. Currently not used.
     % f =(DX-xt)*sqrt(xt^2+DT^2)-cri*xt*sqrt((DX-xt)^2+DF^2);

        if DX >=0 && (xt >=0 && xt<= DX)  && axi < tol
        xi = xr*DT;
        flag =1;
    
        elseif DX <0 && (xt <=0 && xt >= DX) && axi < tol
         
       xi = xr*DT;
       flag =1;
         
        end
 end      
        if flag == 0

        % if interface intersection value returned by Ferrari's 
        % method lies outside the permissable region or the 
        % tolerance on being real is not met, use fzero instead
  
           xi=fzero(@interface2,[0,DX], [], cr, DF, DT, DX);

        end

end
end