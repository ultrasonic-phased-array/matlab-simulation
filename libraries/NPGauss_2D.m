function p = NPGauss_2D(b,f,c,e,x,z)
% p = NPGauss_2D(b,f,c,e,x,z) calculates the normalized pressure
% of an element of length 2b (in mm), at a frequency, f, (in MHz),in a 
% fluid whose wave speed is c (in m/sec). The offset of the center of
%the element in the x-direction is e (in mm) and the pressure is
% calculated at a point (x,z) (in mm). The function uses a non-paraxial 
% expansion of a cylindrical wave and 10 Gaussians to model piston
% behavior of the element. 

% get the Gaussian coefficients of Wen and Breazeale
[A, B] = gauss_c10;

%define non-dimensional quantities
xb=x/b;
zb=z/b;
eb=e/b;
Rb=sqrt((xb-eb).^2 +zb.^2);
kb= 2000*pi*f*b/c;
Db = kb/2;
cosp=zb./Rb;


%calculate normalized pressure field from 10 Gaussians
p =0;
for nn= 1:10
    arg =(cosp.^2 +1i*B(nn).*Rb./Db);
    Dn = sqrt(arg);
    amp = A(nn).*exp(1i.*kb.*Rb)./Dn;
    p = p + amp.*exp(-1i.*kb.*(xb.^2)./(2.*Rb.*arg));
end