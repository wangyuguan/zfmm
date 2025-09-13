function Ax = solver_apply_A_fmm(x,zks,chnkrtotal,sysmatloc)

ns = chnkrtotal.npt;
nt = chnkrtotal.npt;


xd = x(1:2:end);
xs = x(2:2:end);


zk1 = complex(zks(1));
zk2 = complex(zks(2));


dipstr = complex(chnkrtotal.wts(:).*xd(:));
charge = complex(chnkrtotal.wts(:).*xs(:));
zsrc = chnkrtotal.r;
zsrc = complex(reshape(zsrc,[2,ns]));
ztarg = chnkrtotal.r;
ztarg = complex(reshape(ztarg,[2,ns]));
dipvec = chnkrtotal.n;
dipvec = complex(reshape(dipvec,[2,ns]));


ifcharge = 1;
ifdipole = 1; 
ifpgh = 2;
eps = 1e-12;
isep = 1;
[pot1,grad1] = zhfmm2d(eps,zk1,ns,zsrc,ifcharge,charge,ifdipole,dipstr,dipvec,nt,ztarg,ifpgh,isep);
pot1 = -pot1;
grad1 = -grad1;

[pot2,grad2] = zhfmm2d(eps,zk2,ns,zsrc,ifcharge,charge,ifdipole,dipstr,dipvec,nt,ztarg,ifpgh,isep);
pot2 = -pot2;
grad2 = -grad2;


Ax = zeros(size(x));
Ax(1:2:end) = pot1-pot2;
Ax(2:2:end) = sum(dipvec.*(grad1-grad2),1);
Ax = Ax + x + sysmatloc*x;
end