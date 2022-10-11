function [xp,yp] = ellipse(x,y,ex,ey)
%x and y are the coordinates of the center of the ellipse
%a semi-major axis
%b semi-minor axis
%0.01 is the angle step
ang=0:0.01:2*pi; 
r = ey.*ex./sqrt((ex*cos(ang)).^2 + (ey*sin(ang)).^2);
xp = x + r.*sin(ang);
yp = y + r.*cos(ang);
end
