function td=delay_laws2D(M, s, Phi, F, c)
% td = delay_laws2D(M,s,Phi,F,c) generates the time delay
% td (in microsec) for an array with M elements, pitch s  
% (in mm), where we want to steer the beam at the angle Phi 
% (in degrees) and focus it at the distance F (in mm) 
% in a single medium of wave speed c (in m/sec). For steering 
% at an angle Phi only the focal length, F, must be set equal 
% to inf.

Mb=(M-1)/2;
m=1:1:M  ;
em =s*((m-1)-Mb);  % location of centroids of elements

switch (F)
    % steering only case
    case inf
        if Phi > 0
            td=1000*s*sind(Phi)*(m-1)/c;
        else
            td=1000*s*sind(abs(Phi))*(M-m)/c;
        end
          
      %steering and focusing case 
    otherwise,
    r1=sqrt(F^2 +(Mb*s)^2 + 2*F*Mb*s*sind(Phi));
    rm = sqrt(F^2+em.^2 - 2*F*em*sind(Phi));
    rM=sqrt(F^2 +(Mb*s)^2 + 2*F*Mb*s*sind(abs(Phi)));
        if Phi > 0
            td=1000*(r1-rm)/c;
        else
            td=1000*(rM-rm)/c;
        end
end
  
