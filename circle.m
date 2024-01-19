function h = circle(x,y,r,c)
hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = fill(xunit, yunit,c);
end