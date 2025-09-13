
function [in,out] = is_interior(targets,delta,xlims,eps,type)

if strcmpi(type, 'flat')
    fun1 = @(t) gaus_cos(t,0,2,1,0);
elseif strcmpi(type, 'bump')
    fun1 = @(t) gaus_cos(t,1,4,4,0);
elseif strcmpi(type, 'gotham')
    as = 2*[1,-0.5,0.3];
    bs = [2*pi,pi/1.3,sqrt(2)*pi];
    cs = [0.8,1.5,3.2];
    fun1 = @(t) cos_wig(t,as,bs,cs);
end
a0 = 2;
b0 = 1;
t0 = xlims(1) + delta;
t1 = xlims(2) - delta;
toff = sqrt(-log(eps))/a0;
tt0 = t0 + toff;
tt1 = t1 - toff;
fun2 = @(t) smooth_bump(t,a0,b0,tt0,tt1);
fr = @(t) prod_curve(t, fun1, fun2);
% tr = real(chnkrtotal.r(1,:));
% fr_tr = fr(tr);
% figure(2)
% plot(tr, fr_tr(2,:))

fr_xtarg =  fr(targets(1,:));
in = targets(2,:)<fr_xtarg(2,:);
out = ~in;

end