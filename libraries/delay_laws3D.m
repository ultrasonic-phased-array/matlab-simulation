function td= delay_laws3D(M, N, sx, sy, theta, phi, F, c)
% td = delay_laws3D(M,N,sx,sy,theta, phi, F, c) generates the time delays
% td (in microseconds) for a 2-D array of MxN elements in a single medium 
% with elements whose pitches are (sx,sy) in the x- and 
% y-directions, respectively (in mm). The steering direction is
% specified by the spherical coordinate angles (theta, phi) (both in
% degrees) and the focusing distance is specied by F (in mm). For steering
% only, F = inf. The wavespeed of medium is c (in m/sec).

% calculate locations of element centroids in  x- and y-directions
m=1:M;
n=1:N;
Mb =(M-1)/2;
Nb=(N-1)/2;
exm=(m-1-Mb)*sx;
eyn=(n-1-Nb)*sy;

%calculate delays (in microseconds)
switch(F)
% if steering only specified, use explicit steering law
    case(inf)
        for mm=1:M
            for nn=1:N
                dt(mm,nn)=1000*(exm(mm)*sind(theta)*cosd(phi) + ...
                    eyn(nn)*sind(theta)*sind(phi))/c;
            end
        end
        % make delays all positive
        td = abs(min(min(dt))) + dt;
% otherwise, if steering and focusing specified, use time delays to 
% the specified  point  
    otherwise,
        for mm=1:M
            for nn=1:N
            r(mm,nn) = sqrt((F*sind(theta)*cosd(phi) -exm(mm))^2 ...
            +(F*sind(theta)*sind(phi)-eyn(nn))^2 +F^2*(cosd(theta))^2);
            end
        end
        td = max(max(1000*r/c)) -1000*r/c;
end
end
            