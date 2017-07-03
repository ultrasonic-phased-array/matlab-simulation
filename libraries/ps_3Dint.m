function [vx,vy,vz] = ps_3Dint(lx,ly,f,mat,ex,ey,angt, Dt0,x,y,z, varargin )
% [vx,vy,vz] = ps_3Dint(lx,ly,f,mat,ex,ey,angt, Dt0, x,y,z,Ropt, Qopt)
% calculates the normalized velocity components (vx,vy,vz) of a rectangular 
% array element radiating waves through a planar fluid/solid interface. The 
% parameters (lx, ly) are the lengths of the element in the x'- and y'-
% directions, respectively (in mm), f is the frequency (in MHz), and mat is
% a vector mat = [d1, cp1, d2, cp2, cs2, type] where (d1, cp1) are the
% density (in gm/cm^3) and compressional wave speed (in m/sec) for the
% fluid and (d2, cp2, cs2) are similarly the density, P-wave speed, and
% S-wave speed for the solid, and type ='p' or 's' for a P-wave or
% S-wave, respectively, in the solid. The distances (ex, ey) are the 
% x'- and y'- coordinates of the centroid of the element relative to the
% center of the array (in mm). The parameters angt is the angle 
% (in degrees) the array makes with respect to the interface, and Dt0 
% is the distance of the center of the array above the interface (in mm).
% The parameters (x,y,z) specify the point(s) in the second medium at 
% which the fields are to be calculated (in mm), where x- and y- are 
% parallel to the interface and z is normal to the interface, pointing
% into the second medium. 
% Ropt and Qopt are optional arguments. Ropt specifies the number of
% segments to use in the x'-direction while Qopt specifies the number of
% segments in the y'-direction . If either Ropt or Qopt are not
% given as input arguments for a given direction then the function uses  
% one segmentper wavelength in that direction, based on the input 
% frequency, f, which must be a scalar when either Ropt or Qopt 
% are not given. 


%extract material densities, wave speeds, and the type of wave in the
%second medium from mat vector
d1 =mat(1);
cp1=mat(2);
d2 =mat(3);
cp2=mat(4);
cs2 =mat(5);
type =mat(6);

% wave speed in the first medium (a fluid) is for compressional waves
c1 =cp1;
% decide which wave speed to use in second medium for specified wave type
if strcmp(type, 'p')
    c2 =cp2;
elseif strcmp(type,'s')
    c2=cs2;
else error(' type must be ''p'' or ''s'' ')
end

%compute wave numbers for waves in first and second medium
k1=2000*pi*f/c1;
k2 =2000*pi*f/c2;

%if number of x-segments is specified then use
if nargin > 11
    R = varargin{1};
    
% else choose number of terms so each segment
% is a wave length  or less
else
    R=ceil(1000*f*lx/c1);
        if R < 1
             R=1;
        end
end

% if number of y-segments is specified then use
if nargin >12
    Q = varargin{2};
    
% else choose number of terms so that each segment
% is a wave length or less
else
    Q=ceil(1000*f*ly/c1);
        if Q < 1
        Q=1;
        end
end

% compute centroid locations of segments in x'- and y'-directions
% relative to the element centroid
xc=zeros(1,R);
yc=zeros(1,Q);
for rr=1:R
    xc(rr) = -lx/2 +(lx/R)*(rr-0.5);
end
for qq=1:Q
    yc(qq) = -ly/2 +(ly/Q)*(qq-0.5);
end

% calculate normalized velocity components as a sum over all the
% segments as an approximation of the Rayleigh-Sommerfeld
% integral
vx=0;
vy=0;
vz=0;

for rr = 1:R
    for qq = 1:Q
        % calculate distance xi along the interface for a ray from a 
        %segment to the specified point in the second medium
        Db = sqrt((x-(ex+xc(rr)).*cosd(angt)).^2 +(y-(ey+yc(qq))).^2);
        Ds = Dt0 + (ex +xc(rr)).*sind(angt);
        xi = pts_3Dint(ex,ey,xc(rr),yc(qq),angt,Dt0,c1,c2,x,y,z);
        
        % calculate incident and refracted angles along the ray,
        % including the special case when ray is at normal incidence
        if Db ==0
            ang1 =0;
        else
            ang1 = atand(xi./Ds);
        end
        
        if ang1 == 0
            ang2 =0;
        else
            ang2=atand((Db-xi)./z);
        end
        % calculate ray path lengths in each medium
        r1 =sqrt(Ds.^2 +xi.^2);
        r2=sqrt((Db-xi).^2 +z.^2);
        % calculate segment sizes in x'- and y'- directions
        dx=lx/R;
        dy =ly/Q;
        
        % calculate (x', y')components of unit vector along the ray in the
        % first medium
        if Db ==0
            uxt =-sind(angt);
            uyt = 0;
        else
        uxt=xi.*(x-(ex+xc(rr)).*cosd(angt)).*cosd(angt)./(Db.*r1) ...
            -Ds.*sind(angt)./r1;
         uyt = xi.*(y - (ey+yc(qq)))./(Db.*r1);
        end
        
        % calculate polarization components for P- and S-waves in the 
        % second medium, including special case of normal incidence 
        if Db == 0
            dpx =0;
            dpy=0;
            dpz=1;
            dsx =1;
            dsy =0;
            dsz=0;
        else
            dpx = (1-xi./Db).*(x-(ex+xc(rr)).*cosd(angt))./r2;
            dpy = (1 -xi./Db).*(y-(ey+yc(qq)))./r2;
            dpz=z./r2;
            dsx = sqrt(dpy.^2 +dpz.^2);
            dsy= -dpx.*dpy./dsx;
            dsz = -dpx.*dpz./dsx;
        end
        % choose polarization components to use based on wave type in the
        % second medium
        if strcmp(type, 'p' )
            px=dpx;
            py=dpy;
            pz =dpz;
        elseif strcmp(type, 's')
            px = dsx;
            py = dsy;
            pz =dsz;
        else error('wrong type')
        end
        % calculate transmission coefficients (based on velocity ratios)
        % for P- and S-waves and choose appropriate coefficient for the
        % specified wave type
        [tpp,tps]= T_fluid_solid(d1,cp1,d2,cp2,cs2, ang1);

        if strcmp(type,'p')
            T=tpp;
        elseif  strcmp(type, 's')
            T = tps;
        end
       % form up the directivity term
       argx = k1.*uxt.*dx/2;
       argx =argx +eps.*(argx == 0);
       argy = k1.*uyt.*dy/2;
       argy = argy + eps.*( argy == 0);
       dir = (sin(argx)./argx).*(sin(argy)./argy);
       % form up the denominator term
       D1 = r1 + r2.*(c2/c1).*(cosd(ang1)./cosd(ang2)).^2;
       D2 = r1 + r2.*(c2/c1);
       % put transmission coefficient, polarization, directivity, phase
       % term and denominator together to calculate velocity components. 
       vx = vx + T.*px.*dir.*exp(1i.*k1.*r1 +1i.*k2.*r2)./sqrt(D1.*D2);
       vy = vy + T.*py.*dir.*exp(1i.*k1.*r1 +1i.*k2.*r2)./sqrt(D1.*D2);  
       vz = vz + T.*pz.*dir.*exp(1i.*k1.*r1 +1i.*k2.*r2)./sqrt(D1.*D2);
    end
end
% include external factor for these components
vx = vx.*(-1i*k1*dx*dy)/(2*pi);  
vy = vy.*(-1i*k1*dx*dy)/(2*pi);  
vz = vz.*(-1i*k1*dx*dy)/(2*pi);  





