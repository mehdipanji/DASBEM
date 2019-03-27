function [distp] = DISTP(xp,yp,x1,y1,x2,y2,x3,y3)
% CALCULATE DISTANCE FROM COLLOCATION NODE TO INTEGRATION ELEMENT
d1 = (xp-x1)^2+(yp-y1)^2;
d2 = (xp-x2)^2+(yp-y2)^2;
d3 = (xp-x3)^2+(yp-y3)^2;
distp = 0.98*sqrt(min([d1,d2,d3]));