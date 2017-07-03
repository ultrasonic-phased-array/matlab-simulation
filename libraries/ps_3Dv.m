function p = ps_3Dv(lx,ly,f,c,ex,ey,x,y,z, varargin )
% p =ps_3Dv(lx, ly, f, c, ex, ey, x,y,z,Popt,Qopt) computes the normalized
% pressure, p, at a location (x,y,z) (in mm) in a fluid
% for a rectangular element of lengths (lx, ly)
% (in mm) along the x- and y-axes, respectively,at a frequency, f, (in MHz)
% ,and for a wave speed, c, (in m/sec) of the fluid. This
% function can used to describe an element in an array by
% specifying non-zero values for (ex,ey) (in mm), which are the offsets 
% of the center of the element along the x- and y-axes, respectively.
% The assumed harmonic time dependency is exp(-2i*pi*f*t)and
% the Rayleigh-Sommerfeld integral for a piston source is used
% as the beam model.
% Popt and Qopt are optional arguments. Popt specifies the number of
% segments to use in the x-direction while Qopt specifies the number of
% segments in the y-direction . If either Popt or Qopt are not
% given as input arguments for a given direction the function uses 
% one segment per wavelength in that direction, based on the input 
% frequency, f, which must be a scalar when either Popt or 
% Qopt are not given. 


%compute wave number
k=2000*pi*f/c;

%if number of x-segments is specified then use
if nargin > 9
    P = varargin{1};
    
% else choose number of terms so each segment
% length is at most a wave length
else
    P=ceil(1000*f*lx/c);
        if P < 1
             P=1;
        end
end

% if number of y-segments is specified then use
if nargin >10
    Q = varargin{2};
    
% else choose number of terms so that each segment
%is a wave length or less
else
    Q=ceil(1000*f*ly/c);
        if Q < 1
        Q=1;
        end
end

%compute centroid locations of segments in x- and y-directions
xc=zeros(1,P);
yc=zeros(1,Q);
for pp=1:P
    xc(pp) = -lx/2 +(lx/P)*(pp-0.5);
end
for qq=1:Q
    yc(qq) = -ly/2 +(ly/Q)*(qq-0.5);
end

% calculate normalized pressure as a sum over all the
% segments as an approximation of the Rayleigh-Sommerfeld
% integral
p=0;
for pp = 1:P
    for qq = 1:Q
        rpq=sqrt((x-xc(pp) -ex).^2 +(y-yc(qq)-ey).^2 +z.^2);
        ux= (x -xc(pp)-ex)./rpq;
        uy = (y-yc(qq)-ey)./rpq;
        ux =ux+eps*(ux == 0);
        uy =uy+eps*(uy == 0);     
        dirx = sin(k.*ux.*lx/(2*P))./(k.*ux.*lx/(2*P));
        diry =sin(k.*uy.*ly/(2*Q))./(k.*uy.*ly/(2*Q));
        p=p + dirx.*diry.*exp(1i*k.*rpq)./rpq;
    end
end
p = p.*(-1i*k*(lx/P)*(ly/Q))/(2*pi);  % include external factor






