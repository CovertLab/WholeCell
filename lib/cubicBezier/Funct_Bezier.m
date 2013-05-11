function [h, Qx, Qy] = Funct_Bezier(Px, Py, n, axesHandle)
% Px contains x-coordinates of control points [Px0,Px1,Px2,Px3]
% Py contains y-coordinates of control points [Py0,Py1,Py2,Py3]
% n is number of intervals

% Equation of Bezier Curve, utilizes Horner's rule for efficient computation.
% Q(t)=(-P0 + 3*(P1-P2) + P3)*t^3 + 3*(P0-2*P1+P2)*t^2 + 3*(P1-P0)*t + Px0

Px0 = Px(1);
Py0 = Py(1);
Px1 = Px(2);
Py1 = Py(2);
Px2 = Px(3);
Py2 = Py(3);
Px3 = Px(4);
Py3 = Py(4);

cx3 = -Px0 + 3 * (Px1 - Px2) + Px3;
cy3 = -Py0 + 3 * (Py1 - Py2) + Py3;
cx2 = 3 * (Px0 - 2 * Px1 + Px2); 
cy2 = 3 * (Py0 - 2 * Py1 + Py2);
cx1 = 3 * (Px1 - Px0);
cy1 = 3 * (Py1 - Py0);
cx0 = Px0;
cy0 = Py0;

dt = 1 / n;
Qx(1) = Px0; % Qx at t=0
Qy(1) = Py0; % Qy at t=0
for i = 1:n  
    t = i * dt;
    Qx(i + 1) = ((cx3 * t + cx2) * t + cx1) * t + cx0;
    Qy(i + 1) = ((cy3 * t + cy2) * t + cy1) * t + cy0;    
end

h = line(Qx, Qy, 'Parent', axesHandle);      %Qx(t) vs Qy(t)
% hold on
% plot(Px,Py,'ro','LineWidth',2); %control points
% plot(Px,Py,'g:','LineWidth',2); %control polygon
% 
% text(Px0+.1,Py0,'\bf{\it{P}_{o}}');
% text(Px1+.1,Py1,'\bf{\it{P}_{1}}');
% text(Px2+.1,Py2,'\bf{\it{P}_{2}}');
% text(Px3-.25,Py3,'\bf{\it{P}_{3}}');

% % % --------------------------------
% % % Author: Dr. Murtaza Khan
% % % Email : drkhanmurtaza@gmail.com
% % % --------------------------------


