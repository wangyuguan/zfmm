

function pot = post_process_fmm(zk,targets,chnkrtotal,soln)
ns = chnkrtotal.npt;
nt = size(targets,2);
ztarg = complex(targets);
zsrc = complex(reshape(chnkrtotal.r,[2,ns]));
zk = complex(zk);
w = chnkrtotal.wts;
dipvec = chnkrtotal.n;
dipvec = complex(reshape(dipvec,[2,ns]));
dipstr = soln(1:2:end);
dipstr = complex(dipstr(:).*w(:));
charge = soln(2:2:end);
charge = complex(charge(:).*w(:));
isep = 1;
ifdipole = 1; 
ifcharge = 1;
ifpgh = 1;
eps = 1e-13;
[pot,~] = zhfmm2d(eps,zk,ns,zsrc,ifcharge,charge,ifdipole,dipstr,dipvec,nt,ztarg,ifpgh,isep);
pot = -pot;
end

