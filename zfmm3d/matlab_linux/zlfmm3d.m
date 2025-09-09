function [pot,grad] = zlfmm3d(eps,ns,zsrc,ifcharge,charge,ifdipole,dipstr,dipvec,nt,ztarg,isep,ifpgh)

pot = complex(zeros(nt,1));
grad = complex(zeros(3,nt));
ifprint = 0;

mex_id_ = 'zlfmm3d(i double[x], i int[x], i dcomplex[xx], i int[x], i dcomplex[x], i int[x], i dcomplex[x], i dcomplex[xx], io dcomplex[x], io dcomplex[xx], i int[x], i dcomplex[xx], i int[x], i int[x], i int[x])';
[pot, grad] = zlfmm3dmex(mex_id_, eps, ns, zsrc, ifcharge, charge, ifdipole, dipstr, dipvec, pot, grad, nt, ztarg, ifpgh, isep, ifprint, 1, 1, 3, ns, 1, ns, 1, ns, 3, ns, nt, 3, nt, 1, 3, nt, 1, 1, 1);

end