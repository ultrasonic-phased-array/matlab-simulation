function p = on_axis_foc2D(b, R, f, c, z)
% p = on_axis_foc2D(b,R, f,c,z) computes the on-axis normalized 
% pressure for a 1-D focused piston element of length 2b 
% and focal length R (in mm).
% The frequency is f (in MHz), b is the transducer half-length
% (in mm), c is the wave speed of the surrounding fluid 
% (in m/sec),and z is the on-axis distance (in mm). The
% paraxial approximation is used to write the pressure field in terms
% of a Fresnel integral. Note: the propagation term exp(ikz) is removed
% from the wave field calculation.

% ensure no division by zero at z =0
z = z +eps*(z == 0);

% define transducer wave number
kb = 2000*pi*f*b/c;

% define u and prevent division by zero
u =(1-z/R);
u = u + eps*( u == 0);

% argument of the Fresnel integral and denominator in on-axis pressure
% equation
x = sqrt((u.*kb.*b)./(pi.*z)).*( z <= R)+...
   sqrt((-u.*kb.*b)./(pi.*z)).*(z > R);
denom = sqrt(u).*(z <= R) + sqrt(-u).*( z > R);
Fr = fresnel_int(x).*( z <= R) + conj(fresnel_int(x)).*(z >R);

% compute normalized on-axis pressure (p/rho*c*v0) with
% the propagation phase term exp(ikz) removed. Use analytical
% values near the focus and the numerical Fresnel integral values
% away from the focus
p=(sqrt(2/i).*sqrt((b/R).*kb/pi)).*( abs(u) <= .005) + ...
   (sqrt(2/i).*Fr./denom).*(abs(u) > .005);




