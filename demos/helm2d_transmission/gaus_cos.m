function [f,fd,fdd] = gaus_cos(t,a,b,c,d)
f   = zeros(size(t));
fd  = zeros(size(t));
fdd = zeros(size(t));

f  = f  + a*exp(-t.^2/b^2).*cos(c*t+d);
fd = fd - a*exp(-t.^2/b^2).*(2*t.*cos(c*t+d)/b^2 + c*sin(c*t+d));
fdd= fdd + a*exp(-t.^2/b^2).*(4*c*t.*sin(c*t+d)/b^2 - ...
    (c^2+2/b^2-4*t.^2/b^2.*cos(c*t+d)));


end