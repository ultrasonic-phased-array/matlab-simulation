function amp=discrete_windows(M, type)
% amp = discrete_windows(M, type) returns the discrete apodization
% amplitudes for M elements of type 'cos' (cosine), 'Han' (Hanning)
% 'Ham' (Hamming), 'Blk' (Blackman), 'tri' (triangle),
% and 'rect' (a window with all ones, i.e. no apodization)

m=1:M;
switch type
    case 'cos'
        amp = sin(pi*(m-1)/(M-1));
    case 'Han'
        amp =(sin(pi*(m-1)/(M-1))).^2;
    case 'Ham'
        amp= 0.54 -0.46*cos(2*pi*(m-1)/(M-1));
    case 'Blk'
        amp=0.42 -0.5*cos(2*pi*(m-1)/(M-1)) + ...
            0.08*cos(4*pi*(m-1)/(M-1));
    case 'tri'
        amp =1 - abs(2*(m-1)/(M-1) -1);
    case 'rect'
        amp = ones(1,M);
    otherwise
        disp(' Wrong type. Choices are ''cos'', ''Han'', ''Ham'', ''Blk'', ''tri'', ''rect''  ')
       
end
