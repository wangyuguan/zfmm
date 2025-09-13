clear
close('all')
run('../../../chunkie/startup.m')
addpath('../../zfmm2d/matlab_mac/')

zk = 2*pi;
ri = 1.3;
zks = zk*[1; ri];

% zks(1) is for the point source below the interface, so it defines the
% solution above the interface


% zks(2) is for the point source above the interface, so it defines the
% solution below the interface

xlims = [-10, 10];
slope = .25;
a_i = 3;
delta = 1;

%% determine end point based on accuracy
toff = 1.5*sqrt(-log(eps))/a_i;
xc = xlims(2) + toff;
im_max = -log((1e-12))/real(zk);
re_max = xc + max(toff, im_max/slope);
xend = re_max - xlims(2);


cstruct = [];
cstruct.type = 'bump';
cstruct.disc = 'adap';
% maxchunklen = 2*pi/abs(zk);
cstruct.maxchunklen =  .25;

[chnkrl, chnkrm, chnkrr, f] = get_boundary_curves(cstruct, xlims, ...
      delta, a_i, xend, slope);
chnkrs = [chnkrl, chnkrm, chnkrr];
chnkrtotal = merge(chnkrs);


figure 
plot(chnkrtotal)
title('surface')



%% point sources
src_out = [];
src_out.r = [0.1;-4]; 


src_in = [];
src_in.r = [-3;6];


%% target points
rmin = min(chnkrtotal); rmax = max(chnkrtotal);
xl = rmax(1)-rmin(1);
yl = rmax(2)-rmin(2);
nx = 200;
xtarg = linspace(xlims(1),xlims(2),nx); 
ytarg = linspace(xlims(1),xlims(2),nx);
[xxtarg,yytarg] = meshgrid(xtarg,ytarg);
targets = zeros(2,length(xxtarg(:)));
targets(1,:) = xxtarg(:); 
targets(2,:) = yytarg(:);


[in,out] = is_interior(targets,delta,xlims,eps,cstruct.type);

targets_in = [];
targets_in.r = targets(:,in);
targets_out = [];
targets_out.r = targets(:,out);



%% compute ground truth 
Sk_in = kernel('helm', 's', zks(1));
Sk_out = kernel('helm', 's', zks(2));

Spk_in = kernel('helm', 'sp', zks(1));
Spk_out = kernel('helm', 'sp', zks(2));

% zztarg = zeros(size(xxtarg));
% ztargin = Sk_in.eval(src_in, targets_in);
% ztargout = Sk_out.eval(src_out, targets_out);
% zztarg(in) = ztargin;
% zztarg(out) = ztargout;
% 
% maxin = max(abs(zztarg(:)));
% maxsc = max(abs(zztarg(:)));
% maxtot = max(abs(zztarg(:)));
% maxu = max(max(maxin,maxsc),maxtot);

% figure 
% h=pcolor(xxtarg,yytarg,imag(zztarg));
% set(h,'EdgeColor','none')
% clim([-maxu,maxu])
% colormap(redblue);
% colorbar
% hold on
% plot(chnkrtotal,'k','LineWidth',3)
% axis equal tight
% hold on 
% scatter(src_out.r(1),src_out.r(2),'filled')
% hold on 
% scatter(src_in.r(1),src_in.r(2),'filled')
% set(gca, "box","off","Xtick",[],"Ytick",[]);
% title('ground truth')


%% Compute rhs and solution
rhs1 = Sk_in.eval(src_in, chnkrtotal)-Sk_out.eval(src_out, chnkrtotal);
rhs1 = rhs1(:);
rhs2 = Spk_in.eval(src_in, chnkrtotal)-Spk_out.eval(src_out, chnkrtotal);
rhs2 = rhs2(:);
rhs = zeros(2*chnkrtotal.npt,1);
rhs(1:2:end) = rhs1;
rhs(2:2:end) = rhs2;


%% setup kernels
skdiff = kernel('helmdiff', 's', zks);
dkdiff = kernel('helmdiff', 'd', zks);
skpdiff = kernel('helmdiff', 'sp', zks);
dkpdiff = kernel('helmdiff', 'dp', zks);

fkern = kernel([dkdiff, (-1)*skdiff; dkpdiff, (-1)*skpdiff]);

opts.corrections = true;
sysmatloc = chunkermat(chnkrtotal, fkern, opts);
Ax = @(x) solver_apply_A_fmm(x,zks,chnkrtotal,sysmatloc);
tic
sol = gmres(Ax,rhs,[],1e-12,200);
% load('sol_analytical_test.mat','sol')
t = toc 


%% post processing 
dkin = kernel('helm', 'd', zks(1));
skin = kernel('helm', 's', zks(1));
fkernin = kernel([dkin, (-1)*skin]);

dkout = kernel('helm', 'd', zks(2));
skout = kernel('helm', 's', zks(2));
fkernout = kernel([dkout, (-1)*skout]);



nx = 1000;
xtarg = linspace(-10,10,nx); 
ytarg = linspace(-2.5,2.5,nx);
[xxtarg,yytarg] = meshgrid(xtarg,ytarg);
targets = zeros(2,length(xxtarg(:)));
targets(1,:) = xxtarg(:); 
targets(2,:) = yytarg(:);
[in,out] = is_interior(targets,delta,xlims,eps,cstruct.type);
zzsoln = zeros(size(xxtarg));
zsolnin =  post_process_fmm(zks(1),targets(:,in),chnkrtotal,sol);
zsolnout = post_process_fmm(zks(2),targets(:,out),chnkrtotal,sol);
zzsoln(in) = zsolnin;
zzsoln(out) = zsolnout;


targets_in = [];
targets_in.r = targets(:,in);
targets_out = [];
targets_out.r = targets(:,out);
zztarg = zeros(size(xxtarg));
ztargin = Sk_in.eval(src_in, targets_in);
ztargout = Sk_out.eval(src_out, targets_out);
zztarg(in) = ztargin;
zztarg(out) = ztargout;


%% plot solutions
% figure 
% h=pcolor(xxtarg,yytarg,imag(zzsoln));
% set(h,'EdgeColor','none')
% clim([-maxu,maxu])
% colormap(redblue);
% colorbar
% hold on
% plot(chnkrtotal,'k','LineWidth',3)
% axis equal tight
% hold on 
% scatter(src_out.r(1),src_out.r(2),'filled')
% hold on 
% scatter(src_in.r(1),src_in.r(2),'filled')
% set(gca, "box","off","Xtick",[],"Ytick",[]);
% title('soln')


%% plot relative errors 

err = zzsoln-zztarg;
rel_err = abs(err)./abs(zztarg);
log10err = log10(abs(rel_err));
idx = find(abs(real(chnkrtotal.r(1,:)))<10);

figure('Units','inches','Position',[1 1 10 5])
h=pcolor(xxtarg,yytarg,log10err);
set(h,'EdgeColor','none')
colormap(redblue);
colorbar
cb = colorbar;
xlabel('$x_1$','FontSize',24,Interpreter='latex')
ylabel('$x_2$','FontSize',24,Interpreter='latex')
% cb.Ticks = 0 : -1 : -10;
% set(gca, "box","off","Xtick",[],"Ytick",[]);
hold on
plot(real(chnkrtotal.r(1,idx)),real(chnkrtotal.r(2,idx)),'k','LineWidth',3)
ax = gca; 
ax.FontSize = 14;
ax.XLabel.FontSize = 24;
ax.YLabel.FontSize = 24;
pos = ax.Position;
ax.Position = [pos(1)+0.05 pos(2) pos(3)*0.9 pos(4)];



%% plot the magnitude of the dipole & charge
npt = size(chnkrtotal.r,2)*size(chnkrtotal.r,3);
figure('Units','inches','Position',[1 1 10 5])
semilogy(real(chnkrtotal.r(1,:)),(abs(sol(1:2:end))),'LineWidth',3,'Color','red','LineStyle','-.');
hold on 
semilogy(real(chnkrtotal.r(1,:)),(abs(sol(2:2:end))),'LineWidth',3,'Color','blue','LineStyle','-');
grid on 
xlabel('$x_1$','FontSize',24,Interpreter='latex')
ax = gca; 
ax.FontSize = 14;
ax.XLabel.FontSize = 24;
ax.YLabel.FontSize = 24;
legend({'dipole','charge'},FontSize=20,Interpreter='latex')


