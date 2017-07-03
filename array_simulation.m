% This script solves for the normalized pressure wave field of a 2-D 
% array of rectangular elements radiating waves in a fluid using the
% MATLAB function ps_3Dv. Both time delay and apodization laws can 
% be specified for the array to steer it and focus it.

steps = 40; %steps for dynamic beam steering

% For video recording:
%writerObj = VideoWriter('out.avi');
%writerObj.FrameRate = 30; % How many frames per second.
%open(writerObj); 

for n=0:0.1:steps 
    %  ------------- give input parameters -------------------------
    lx = 5;   % element length in x-direction (mm)
    ly = 5;   % element length in y-direction (mm)
    gx=5;   % gap length in x-direction (mm)
    gy = 5; % gap length in y-direction (mm)
    f= 0.04;       % frequency (MHz)
    c = 340;   % wave speed (m/sec)
    L1 =19;      % number of elements in x-direction
    L2 =19;      % number of elements in y-direction
    theta =n;   % steering angle in theta direction (deg)
    phi =0;     %steering angle in phi direction (deg)
    Fl = 200;  % focal distance (mm)
    %weighting choices are 'rect','cos', 'Han', 'Ham', 'Blk', 'tri' 
    ampx_type ='rect';   % weighting coeffcients in x-direction
    ampy_type ='rect';   % weighting coefficients in y-direction

    % field points (x,y,z)to evaluate
    xs= linspace(-150,150, 200);
    zs= linspace(1, 300, 200);
    y=0;
    [x,z]=meshgrid(xs,zs);

    % ---------------- end input parameters ----------------------

    % calculate array pitches
    sx = lx+gx;
    sy = ly+gy;

    % compute centroid locations for the elements
    Nx = 1:L1;
    Ny = 1:L2;
    ex =(2*Nx -1-L1)*(sx/2);
    ey =(2*Ny -1 -L2)*(sy/2);

    % generate time delays, put in exponential 
    % and calculate amplitude weights
    td =delay_laws3D(L1,L2,sx,sy,theta,phi,Fl,c);
    delay = exp(1i.*2.*pi.*f.*td);

    Cx = discrete_windows(L1,ampx_type);
    Cy = discrete_windows(L2,ampy_type);


    % calculate normalized pressure
    p=0;
    for nn=1:L1
        for ll=1:L2
            p = p + Cx(nn)*Cy(ll)*delay(nn,ll)...
                *ps_3Dv(lx,ly,f,c,ex(nn),ey(ll),x,y,z);
        end
    end

    % ---------------- outputs --------------------------
    %plot results based on specification of (x,y,z) points
    subplot(2,1,1); imagesc(xs,zs,abs(p)); title("2D section view");  ylabel("Z in mm");

    %heatmap of individual phase delay
    phase = angle(delay);
    subplot(2,1,2); pcolor(phase); title("Phase of individual transceivers");

    colormap hot
    colorbar

    %Video recording settings:
    %set(gca, 'CameraViewAngle',cva);
    %set(gca, 'CameraUpVector',cuv);
    %set(gca, 'CameraTarget',ct);
    %set(gca, 'CameraPosition',cp);
    drawnow;

    %frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
    %writeVideo(writerObj, frame);
end
%close(writerObj); % Saves the movie.
        