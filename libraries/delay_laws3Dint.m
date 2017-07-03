function td = delay_laws3Dint(Mx,My,sx,sy,thetat, phi, ...
    theta2,DT0, DF, c1,c2, plt)
% td = delay_laws3Dint(Mx,My,sx,sy,thetat,phi, theta2, DT0,DF,c1,c2,plt)
% calculates the delay laws for steering and focusing a 2-D array 
% through a planar interface between two media in three dimensions. 
% (Mx, My)are the number of elements in the (x', y') directions, (sx, sy)
% are the pitches (in mm), and thetat is the the angle that array
% makes with the interface(in degrees).Steering and focusing to a point in
% the second medium is specified by giving the angles theta2 and 
% phi,(both in degrees). The height of the center of the array above
% the interface is DT0 (in mm).The wave speeds of the first and second 
% media are (c1, c2), respectively (in m/sec). The plt argument is a string
% ('y' or 'n')that specifies if a plot of the rays from the centroids
% of the elements to the point in the second medium is wanted ('y') or
% not ('n'). Plotting is not done if steering only (DF = inf) is specified.

% compute wave speed ratio
cr=c1/c2;

% compute element centroid locations
Mbx=(Mx-1)/2;
Mby=(My-1)/2;
mx=1:1:Mx;
ex=(mx-1-Mbx)*sx;
my=1:1:My;
ey=(my-1-Mby)*sy;

%initialize variables to be used
t=zeros(Mx,My);
Db=zeros(Mx,My);
De=zeros(1,Mx);
xi=zeros(Mx,My);

ang1 =asind(c1*sind(theta2)/c2); % ang1e in first medium (in degrees)

switch(DF)
    % steering only case, use linear steering law
    case inf
       
        ux= sind(ang1)*cosd(phi)*cosd(thetat) -cosd(ang1)*sind(thetat);
        uy =sind(ang1)*sind(phi);
        for m =1:Mx
            for n = 1:My
                t(m,n)= 1000*(ux*ex(m)+uy*ey(n))/c1;  %time in microsec
            end
        end
        td = abs(min(min(t))) +t; % make sure delay is positive
        
    % steering and focusing case
    otherwise
        % determine distances De, Db needed in arguments of ferrari2
        % function
    DQ=DT0*tand(ang1)+DF*tand(theta2);
    x=DQ*cosd(phi);
    y=DQ*sind(phi);
    for m=1:Mx
        for n = 1:My
            Db(m,n) = sqrt((x-ex(m)*cosd(thetat))^2 +(y-ey(n))^2);
        end 
    end
    De = DT0 +ex*sind(thetat);
    % use ferrari2 method to determine distance, xi, where a ray from an 
    % element to the point (x, y, DF)intesects the interface 
    % in the plane of incidence
    for m=1:Mx
        for n = 1:My
            xi(m,n) = ferrari2(cr,DF,De(m),Db(m,n));
        end
    end
    % use ray distances to calculate time advances (in microsec) 
    for m=1:Mx
        for n=1:My
            t(m,n) = 1000*sqrt(xi(m,n)^2 +De(m)^2)/c1 +...
                1000*sqrt(DF^2+(Db(m,n) -xi(m,n))^2)/c2;
        end
    end
    % turn time advances into delays and make all delays positive
    td =max(max(t)) -t;
    
    % plotting rays option
    if strcmp(plt, 'y')

      for m=1:Mx
          for n = 1:My
              xp(1,1) = ex(m)*cosd(thetat);
              zp(1,1)=DT0 +ex(m)*sind(thetat);
              yp(1,1) = ey(n);
              xp(2,1) = ex(m)*cosd(thetat) + xi(m,n)*(x-ex(m)*cosd(thetat))/Db(m,n);
              yp(2,1) = ey(n) + xi(m,n)*(y-ey(n))/Db(m,n);
              zp(2,1) =0;
              xp(3,1) = x;
              yp(3, 1) = y;
              zp(3,1) =-DF;
              plot3(xp,yp,zp)
              hold on
          end
      end
      hold off
    end
    %end plotting rays option
end
