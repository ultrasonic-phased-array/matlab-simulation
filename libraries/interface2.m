function y =interface2(x, cr, df, dp, dpf)
% y = interface2(x, cr, df, dp, dpf) outputs the value of a function, y, 
% which is zero if the input argument,x,is the location along an interface
% where Snell's law is satisfied. The input parameter cr =c1/c2, where c1
% is the wave speed in medium one, and c2 is the wave speed in medium 2, 
% The other input parameters (df, dp, dpf) define a ray which goes from
% a point in medium one to the interface and then to a point in medium
% two, where df = DF is the depth of the point in medium two,
% dp = DT is the height of the point in medium one, and dpf = DX is the 
% separation distance between the points in medium one and two 
% (see Fig 5.4 in the text). The function y used here is c1 times the 
% function defined in Eq.(5.2.6) in the text.

% the function,y, 


y =x./sqrt(x.^2+dp^2)-cr*(dpf-x)./sqrt((dpf-x).^2 +df^2);
