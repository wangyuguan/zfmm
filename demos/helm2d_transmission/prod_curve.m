function [f,fd,fdd] = prod_curve(t,fun1,fun2,a,b,t0,t1)

    tsz = size(t);
    t = t(:);
    f  = zeros(2,numel(t));
    fd = zeros(2,numel(t));
    fdd= zeros(2,numel(t));

    f(1,:)  = t;
    fd(1,:) = 1;
    fdd(1,:)= 0;

    [f1,f1d,f1dd] = fun1(t);
    [f2,f2d,f2dd] = fun2(t);

    f(2,:)  = f1.*f2;
    fd(2,:) = f1.*f2d + f2.*f1d;
    fdd(2,:)= f1dd.*f2 + 2*f1d.*f2d + f1.*f2dd;
 

    if (nargin < 4)
        ifcmplx = false;
    else
        ifcmplx = true;
    end
    if (ifcmplx)
    phi = @(t,u,v,z) u*(t-z).*erfc(u*(t-z))*v - exp(-u^2*(t-z).^2)/sqrt(pi)*v;
    phid= @(t,u,v,z) u*erfc(u*(t-z))*v;
    phidd=@(t,u,v,z) -u*u*exp(-u^2*(t-z).^2)*2*v/sqrt(pi);
    f(1,:) = f(1,:).' + 1i*(phi(t,a,b,t0) - phi(t,-a,b,t1)); 
    fd(1,:)= fd(1,:).' + 1i*(phid(t,a,b,t0) - phid(t,-a,b,t1));
    fdd(1,:) =fdd(1,:).' + 1i*(phidd(t,a,b,t0) - phidd(t,-a,b,t1));
    end   

    f   = reshape(f,[2,tsz]);
    fd  = reshape(fd,[2,tsz]);
    fdd = reshape(fdd,[2,tsz]);
end
