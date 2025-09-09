function [pot,grad] = zlfmm2d(eps,ns,zsrc,ifcharge,charge,ifdipole,dipstr,dipvec,nt,ztarg,ifpgh,isep)

pot = complex(zeros(nt,1));
grad = complex(zeros(2,nt));
ifprint = 0;

mex_id_ = 'zlfmm2d(i double[x], i int[x], i dcomplex[xx], i int[x], i dcomplex[x], i int[x], i dcomplex[x], i dcomplex[xx], io dcomplex[x], io dcomplex[xx], i int[x], i dcomplex[xx], i int[x], i int[x], i int[x])';
[pot, grad] = zlfmm2dmex(mex_id_, eps, ns, zsrc, ifcharge, charge, ifdipole, dipstr, dipvec, pot, grad, nt, ztarg, ifpgh, isep, ifprint, 1, 1, 2, ns, 1, ns, 1, ns, 2, ns, nt, 2, nt, 1, 2, nt, 1, 1, 1);

end