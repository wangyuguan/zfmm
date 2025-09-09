clc 
clear all 


addpath ../../matlab_linux/ 



n = 20000;
slope = 0.3; 
zsrc = rand(3,n);
zsrc = complex(zsrc+1j*slope*zsrc);
ztarg = zsrc;
isep = 1;

charge = complex(rand(n,1));
dipstr = complex(rand(n,1));
dipvec = complex(rand(3,n));

ifcharge = 1;
ifdipole = 1;
ifpgh = 2;

tic 
[pot_ex, grad_ex] = zl3devaldirect(n,ztarg,n,zsrc,charge,dipstr,dipvec,ifcharge,ifdipole,ifpgh); 
toc 


eps = 1e-12;
tic 
[pot, grad] = zlfmm3d(eps,n,zsrc,ifcharge,charge,ifdipole,dipstr,dipvec,n,ztarg,isep,ifpgh);
toc

err_pot = norm(pot(:)-pot_ex(:))/norm(pot_ex(:))
err_grad = norm(grad(:)-grad_ex(:))/norm(grad_ex(:))