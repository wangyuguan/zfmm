clc 
clear all 


addpath ../../matlab_mac/ 


rng(1)
n = 10000;
slope = 0.2; 
zk = complex(2*pi);
zsrc = rand(2,n);
zsrc = complex(zsrc+1j*slope*zsrc);
ztarg = zsrc;
isep = 1;

charge = complex(rand(n,1));
dipstr = complex(rand(n,1));
dipvec = complex(rand(2,n));

ifcharge = 1;
ifdipole = 1;
ifpgh = 2;

tic 
[pot_ex, grad_ex] = zh2devaldirect(zk,n,ztarg,n,zsrc,charge,dipstr,dipvec,ifcharge,ifdipole,ifpgh); 
toc 


eps = 1e-12;
tic 
[pot, grad] = zhfmm2d(eps,zk,n,zsrc,ifcharge,charge,ifdipole,dipstr,dipvec,n,ztarg,ifpgh,isep);
toc

err_pot = norm(pot(:)-pot_ex(:))/norm(pot_ex(:))
err_grad = norm(grad(:)-grad_ex(:))/norm(grad_ex(:))