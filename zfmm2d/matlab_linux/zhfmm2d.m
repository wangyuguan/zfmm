function [pot,grad] = zhfmm2d(eps,zk,ns,zsrc,ifcharge,charge,ifdipole,dipstr,dipvec,nt,ztarg,ifpgh,ifprint,isep)

pot = complex(zeros(nt,1));
grad = complex(zeros(2,nt));

mex_id_ = 'zhfmm2d(i double[x], i dcomplex[x], i int[x], i dcomplex[xx], i int[x], i dcomplex[x], i int[x], i dcomplex[x], i dcomplex[xx], io dcomplex[x], io dcomplex[xx], i int[x], i dcomplex[xx], i int[x], i int[x], i int[x])';
[pot, grad] = zhfmm2dmex(mex_id_, eps, zk, ns, zsrc, ifcharge, charge, ifdipole, dipstr, dipvec, pot, grad, nt, ztarg, ifpgh, ifprint, isep, 1, 1, 1, 2, ns, 1, ns, 1, ns, 2, ns, nt, 2, nt, 1, 2, nt, 1, 1, 1);

end