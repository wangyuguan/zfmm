function [pot,grad] = zl2devaldirect(nt,ztarg,ns,zsrc,charge,dipstr,dipvec,ifcharge,ifdipole,ifpgh)

pot = complex(zeros(nt,1));
grad = complex(zeros(2,nt));

mex_id_ = 'zl2devaldirect(i int[x], i dcomplex[xx], i int[x], i dcomplex[xx], i dcomplex[x], i dcomplex[x], i dcomplex[xx], i int[x], io dcomplex[x], io dcomplex[xx], i int[x], i int[x])';
[pot, grad] = zl2devaldirectmex(mex_id_, nt, ztarg, ns, zsrc, charge, dipstr, dipvec, ifcharge, pot, grad, ifdipole, ifpgh, 1, 2, nt, 1, 2, ns, ns, ns, 2, ns, 1, nt, 2, nt, 1, 1);

end