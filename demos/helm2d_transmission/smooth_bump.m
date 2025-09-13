function [f,fd,fdd] = smooth_bump(t,a,b,t0,t1)
    b = b/abs(a);
    phi=  @(t,u,v,z) u*erfc(u*(t-z))*v;
    phid= @(t,u,v,z) -u^2*exp(-u^2*(t-z).^2)*2*v/sqrt(pi);
    phidd=@(t,u,v,z)  u^4*(t-z).*exp(-u^2*(t-z).^2)*4*v/sqrt(pi);
    f  = 2-((phi(t,a,b,t0) - phi(t,-a,b,t1))); 
    fd = -((phid(t,a,b,t0) - phid(t,-a,b,t1)));
    fdd= -((phidd(t,a,b,t0) - phidd(t,-a,b,t1)));
    %f  = 0*f/2;
    %fd = 0*fd/2;
    %fdd= 0*fdd/2;
    %f = ones(size(f));
    %fd= 0*zeros(size(f));
    %fdd=0*zeros(size(f));

    a*erfc(a*(t-t0))*b;
    
end

