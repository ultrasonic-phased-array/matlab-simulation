function td=delay_laws2D_int( M, s, angt,ang20,DT0,DF, c1, c2, plt)
% td = delay_laws2D_int(M,s,angt, an20, DT0, DF,c1,c2, plt) calculates 
% the delay laws for steering and focusing an array of 1-D elements
% through a planar interface between two media in two dimensions. The
% number of elements is M, the pitch is s (in mm), the angle that array
% makes with the interface is ang(in degrees). The height of the
% center of the array above the interface is DT0 (in mm). Steering and
% focusing to a point in the second medium are specified by giving the
% refracted angle, ang20, (in degrees) and the depth in the second
% medium, DF, (in mm). The wave speeds of the first and second media
% are (c1, c2), respectively (in m/sec). The plt argument is a string
% ('y' or 'n')that specifies if a plot of the rays from the centroids
% of the elements to the point in the second medium is wanted ('y') or
% not ('n')



cr = c1/c2; % wave speed ratio
Mb=(M-1)/2;
%compute location of element centroids, e
m=1:1:M;
e =(m-1-Mb)*s;
% computed parameters:
% ang10, incident angle of central ray, deg
% DX0, distance along interface from center of array to focal point, mm
% DT, heights of elements above interface, mm
% DX, distances along interface from elements to focal point, mm
        ang10 = asind((c1/c2).*sind(ang20));
        DX0 = DF.*tand(ang20) + DT0.*tand(ang10);
        DT = DT0 + e.*sind(angt);
        DX =DX0 - e.*cosd(angt);
switch (DF)
% steering only case, use linear law
    case inf
        if (ang10 -angt)>0
            td = 1000*(m-1)*s*sind(ang10-angt)/c1;
        else
            td = 1000*(M-m)*s*abs(sind(ang10-angt))/c1;
        end
% plotting rays option      
        if strcmp(plt,'y')
            for nn = 1:M
    xp2(1, nn) = e(nn)*cosd(angt);
    yp2(1, nn) = DT(nn);
    xp2(2, nn) = e(nn)*cosd(ang10-angt)/cosd(ang10) +DT0*tand(ang10);
    dm=e(nn)*cosd(ang10-angt)/cosd(ang10);
    if ang20 >0
        dM = e(M)*cosd(ang10-angt)/cosd(ang10);
    else
        dM =e(1)*cosd(ang10-angt)/cosd(ang10) ;
    end
    yp2(2, nn) = 0;
    xp2(3, nn) = xp2(2,nn) + (dM-dm)*sind(ang20)*sind(ang20);
    yp2(3, nn) = -(dM-dm)*sind(ang20)*cosd(ang20);
end
plot(xp2, yp2, 'b')
        end
% end plotting rays option

% steering and focusing case
    otherwise,
        
%solve for ray intersection locations on interface, xi,(in mm)
%and path lengths in medium 1 and medium 2, r1, r2 (mm)

        xi=zeros(1,M);
        r1=zeros(1,M);
        r2=zeros(1,M);
        for mm = 1:M
         xi(mm) = ferrari2(cr,DF,DT(mm),DX(mm));
        r1(mm) =sqrt(xi(mm)^2 +(DT0+e(mm)*sind(angt))^2);
        r2(mm) =sqrt( (xi(mm) +e(mm)*cosd(angt)-DX0)^2 +DF^2);
        end

% solve for time advances (in microsec), turn into delays, and
% make the delays ,td, positive
        t= 1000*r1/c1 +1000*r2/c2;
        td=max(t) -t;
% plotting rays option        
        if strcmp(plt, 'y')

            for nn = 1:M
    xp(1, nn) = e(nn)*cos(angt*pi/180);
    yp(1, nn) = DT(nn);
    xp(2, nn) = e(nn)*cos(angt*pi/180) +xi(nn);
    yp(2, nn) = 0;
    xp(3, nn) = DX0;
    yp(3, nn) = -DF;
            end
  plot(xp, yp, 'b')          
        end
%end plotting rays option   

end
end

