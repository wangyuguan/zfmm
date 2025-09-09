function [pot,grad] = zl3devaldirect(nt,ztarg,ns,zsrc,charge,dipstr,dipvec,ifcharge,ifdipole,ifpgh)

pot = complex(zeros(nt,1));
grad = complex(zeros(3,nt));

mex_id_ = 'zl3devaldirect(i int[x], i dcomplex[xx], i int[x], i dcomplex[xx], i dcomplex[x], i dcomplex[x], i dcomplex[xx], i int[x], io dcomplex[x], io dcomplex[xx], i int[x], i int[x])';
[pot, grad] = zl3devaldirectmex(mex_id_, nt, ztarg, ns, zsrc, charge, dipstr, dipvec, ifcharge, pot, grad, ifdipole, ifpgh, 1, 3, nt, 1, 3, ns, ns, ns, 3, ns, 1, nt, 3, nt, 1, 1);

end