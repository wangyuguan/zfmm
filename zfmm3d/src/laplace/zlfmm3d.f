C
C   subroutines for 3d complex helmholtz fast multipole method  
C
C       For i=1,...,nt,  evaluate:
C
C                              ns
C       u(xi) = (1/(4*pi)) * \sum charge(j)*1/R(xi,yj) 
C                              j=1
C          - (1/(4*pi)) * dipstr(j)*<dipvec(:,j), \nabla_y R(xi,yj)>
C                 
C               where R(x,y)=sqrt((x1-y1)^2+(x2-y2)^2+(z2-z2)^2)
C
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine zlfmm3d(eps,ns,zsrc,ifcharge,charge,ifdipole,dipstr,
     1 dipvec,pot,grad,nt,ztarg,ifpgh,isep,ifprint)
  
c----------------------------------------------
      implicit real *8 (a-h,o-z)
      complex *16 zsrc(3,ns),ztarg(3,nt),dipvec(3,ns)
      complex *16 charge(ns),dipstr(ns),pot(nt),grad(3,nt)
      dimension src(3,ns),targ(3,nt)
      dimension srcsort(3,ns)
      dimension targsort(3,nt)
      dimension zc(3)
      complex *16 zsrcsort(3,ns),chargesort(ns),dipstrsort(ns)
      complex *16 ztargsort(3,nt),dipvecsort(3,ns)
      integer *8 ltree,iptr(8)
      integer, allocatable :: itree(:)
      real *8, allocatable :: centers(:,:),boxsize(:)
      complex *16, allocatable :: zcsrc(:,:),zctarg(:,:)
      integer, allocatable :: isrc(:),isrcse(:,:)
      integer, allocatable :: itarg(:),itargse(:,:)
      real *8, allocatable :: rscales(:)
      integer, allocatable :: nterms(:)
      complex *16 ima 

      pot = 0 
      grad = 0 
      ima = (0.0d0,1.0d0)
      zsrcsort = 0
      chargesort = 0 
      dipstrsort = 0 
      dipvecsort = 0
      ztargsort = 0

      if ((ifcharge.lt.0).or.(ifcharge.gt.1)) then 
        ifcharge = 1  
      endif 

      if ((ifdipole.lt.0).or.(ifdipole.gt.1)) then 
        ifdipole = 0  
      endif

      if ((ifpgh.lt.0).or.(ifpgh.gt.2)) then 
        ifpgh = 1
      endif 

C
CCC   get real locations 
C
      do i = 1,ns
        src(1:3,i) = real(zsrc(1:3,i))
      enddo   

      do i = 1,nt
        targ(1:3,i) = real(ztarg(1:3,i))
      enddo 
C
CC      call the tree memory management
CC      code to determine number of boxes,
CC      number of levels and length of tree
C

      ifpghtarg = 0

      call lndiv(eps,ns,nt,ifcharge,ifdipole,ifpgh,
     1   ifpghtarg,ndiv,idivflag)
      if (ifprint.eq.1) then 
          call prinf('ndiv=*',ndiv,1)
          call prinf('idivflag=*',idivflag,1)
      endif 
      nlmax = 51
      nlevels = 0
      nboxes = 0
      ltree = 0
      nlmin = 0
      iper = 0
      ifunif = 0
      mnbors = (2*isep+1)**3
      call pts_tree_mem(src, ns, targ, nt, idivflag, ndiv, nlmin,            
     1     nlmax, ifunif, iper, isep, nlevels, nboxes, ltree)
      
      if (ifprint.eq.1) then 
          call prinf('nboxes=*',nboxes,1)
          call prinf('nlevels=*',nlevels,1)
          call prinf('ltree=*',ltree,1)
          call prinf('ndiv=*',ndiv,1)
      endif 

      allocate(itree(ltree))
      allocate(centers(3,nboxes))
      allocate(boxsize(0:nlevels))

      iptr = 0 
      itree = 0 
      centers = 0 
      boxsize = 0 

c
c       call the tree code
c
      call pts_tree_build(src,ns,targ,nt,idivflag,ndiv,nlmin,
     1  nlmax,ifunif,iper,isep,nlevels,nboxes,ltree,itree,iptr,
     2  centers,boxsize)

C c        * iptr(1) - laddr
C c        * iptr(2) - ilevel
C c        * iptr(3) - iparent
C c        * iptr(4) - nchild
C c        * iptr(5) - ichild
C c        * iptr(6) - ncoll
C c        * iptr(7) - coll
C c        * iptr(8) - ltree

      allocate(itarg(nt))
      allocate(itargse(2,nboxes))
      allocate(isrc(ns))
      allocate(isrcse(2,nboxes))
      call pts_tree_sort(nt,targ,itree,ltree,nboxes,nlevels,
     1   iptr,centers,itarg,itargse)
      call pts_tree_sort(ns,src,itree,ltree,nboxes,nlevels,
     1   iptr,centers,isrc,isrcse)

      
C
CC     reorder sources, charges, and targets
C
      call dreorderf(3,ns,src,srcsort,isrc)
      call zreorderf(3,ns,zsrc,zsrcsort,isrc)
      if (ifcharge.eq.1) call zreorderf(1,ns,charge,chargesort,isrc)
      call dreorderf(3,nt,targ,targsort,itarg)
      call zreorderf(3,nt,ztarg,ztargsort,itarg)
      if (ifdipole.eq.1) call zreorderf(1,ns,dipstr,dipstrsort,isrc)
      if (ifdipole.eq.1) call zreorderf(3,ns,dipvec,dipvecsort,isrc)

C
CC      compute the centers of the complex boxes 
C
      allocate(zcsrc(3,nboxes),zctarg(3,nboxes))
      zcsrc = 0 
      zctarg = 0
      do ibox=1,nboxes 
C
CC          compute complex scource boxcenters 
C
        if (isrcse(1,ibox).le.isrcse(2,ibox)) then 
          npts = isrcse(2,ibox)-isrcse(1,ibox)+1 
          zc = 0 
          do i=isrcse(1,ibox),isrcse(2,ibox)
            zc = zc+aimag(zsrcsort(1:3,i))
          enddo 
          zc = zc/npts 
          zcsrc(1:3,ibox) = centers(1:3,ibox)+ima*zc
        endif        


C            call prin2('zcsrc=*',zcsrc(1:3,ibox),6)
C
CC          compute complex target boxcenters 
C

        if (itargse(1,ibox).le.itargse(2,ibox)) then 
          npts = itargse(2,ibox)-itargse(1,ibox)+1 
          zc = 0 
          do i=itargse(1,ibox),itargse(2,ibox)
            zc = zc+aimag(ztargsort(1:3,i))
          enddo 
          zc = zc/npts
          zctarg(1:3,ibox) = centers(1:3,ibox)+ima*zc
        endif 

      enddo 
C
CC      rescaling parameters
C

      b0 = boxsize(0)
      if (ifprint.eq.1) call prin2('b0=*',b0,1)
      b0inv = 1.0d0/b0
      b0inv2 = b0inv**2
      call drescale(3*ns,srcsort,b0inv)
      call drescale(3*nt,targsort,b0inv)
      call zrescale(3*ns,zsrcsort,b0inv)
      call zrescale(3*nt,ztargsort,b0inv)
      if (ifcharge.eq.1) call zrescale(ns,chargesort,b0inv)
      if (ifdipole.eq.1) call zrescale(3*ns,dipvecsort,b0inv2)
      call drescale(3*nboxes,centers,b0inv)
      call zrescale(3*nboxes,zcsrc,b0inv)
      call zrescale(3*nboxes,zctarg,b0inv)
      call drescale(nlevels+1,boxsize,b0inv)

      if (ifprint.eq.1) call prin2('boxsize=*',boxsize,nlevels+1)


c
cc      compute scaling factor for multipole/local expansions
c       and lengths of multipole and local expansions
c
      allocate(rscales(0:nlevels))
      rscales = 0 
      do ilev = 0,nlevels
        rscales(ilev) = boxsize(ilev)
        if (rscales(ilev).gt.1) then 
          rscales(ilev) = 1
        endif 
      enddo

      if (ifprint.eq.1) then 
        call prin2('rscales=*',rscales,nlevels+1)
      endif 

      allocate(nterms(0:nlevels))
      nterms = 0 
      nmax = 0
      do i=0,nlevels
        call zl3dterms(eps, nterms(i))
        if(nterms(i).gt.nmax) nmax = nterms(i)
      enddo

      if (ifprint.eq.1) then 
        call prinf('nmax=*',nmax,1)
        call prinf('nterms=*',nterms,nlevels+1)
      endif 

      call zlfmm3dmain(ns,zsrcsort,chargesort,dipstrsort,dipvecsort,
     1 nt,ztargsort,isrc,itarg,isrcse,itargse,itree,ltree,iptr,
     2 nlevels,nboxes,boxsize,rscales,nterms,nmax,centers,zcsrc,
     3 zctarg,itree(iptr(1)),itree(iptr(2)),itree(iptr(3)),
     4 itree(iptr(4)),itree(iptr(5)),itree(iptr(6)),
     5 pot,grad,ifcharge,ifdipole,ifpgh,ifprint,isep)

      if (ifpgh.ge.2) call zrescale(3*nt,grad,b0inv)


      end 

C
C
C
C
C
      subroutine zlfmm3dmain(ns,zsrcsort,chargesort,dipstrsort,
     1 dipvecsort,nt,ztargsort,isrc,itarg,isrcse,itargse,itree,ltree,
     2 iptr,nlevels,nboxes,boxsize,rscales,nterms,nmax,centers,zcsrc,
     3 zctarg,laddr,ilevel,iparent,nchild,ichild,ncol,pot,grad,
     4 ifcharge,ifdipole,ifpgh,ifprint,isep)
      implicit real *8 (a-h,o-z)
      complex *16 zsrcsort(3,ns),chargesort(ns)
      complex *16 dipstrsort(ns),dipvecsort(3,ns)
      dimension isrc(ns),itarg(nt)
      dimension isrcse(2,nboxes), itargse(2,nboxes)
      complex *16 ztargsort(3,nt),pot(nt),grad(3,nt)
      integer *8 ltree, iptr(8),iaddr(2,nboxes),lmptot
      dimension ilevel(nboxes),nchild(nboxes),ncol(nboxes)
      dimension iparent(nboxes),ichild(8,nboxes)
      dimension itree(ltree)
      dimension boxsize(0:nlevels), rscales(0:nlevels)
      dimension nterms(0:nlevels)
      dimension centers(3,nboxes)
      complex *16 zcsrc(3,nboxes),zctarg(3,nboxes)
      integer  laddr(2,0:nlevels)
      integer, allocatable :: list1(:,:),list2(:,:),list3(:,:)
      integer, allocatable :: list4(:,:)
      integer, allocatable :: nlist1s(:),nlist2s(:),nlist3s(:)
      integer, allocatable :: nlist4s(:)
      complex *16, allocatable :: rmlexp(:)
      complex *16 potsort(nt),gradsort(3,nt)
c      complex *16 mp(0:nmax,-nmax:nmax,nboxes)
c      complex *16 loc(0:nmax,-nmax:nmax,nboxes)
      real *8, allocatable :: amat(:,:)

C
CC       initialize all variables 
C
      pot = 0 
      potsort = 0 
      mp = 0 
      loc = 0 
      grad = 0 
      gradsort = 0 

      if (ifprint.eq.1) then 
        do ilev = 0, nlevels
          call prinf('laddr=*',laddr(1,ilev),2)
        enddo
      endif 

C
CC       compute list info
C
      mnlist1 = 0
      mnlist2 = 0
      mnlist3 = 0
      mnlist4 = 0
      iper = 0 
      mnbors = (2*isep+1)**3

      call computemnlists(nlevels,nboxes,itree(iptr(1)),boxsize,
     1  centers,itree(iptr(3)),itree(iptr(4)),
     1  itree(iptr(5)),isep,itree(iptr(6)),mnbors,
     1  itree(iptr(7)),iper,mnlist1,mnlist2,mnlist3,mnlist4)

      if (ifprint.eq.1) then 
        call prinf('mnlist1=*',mnlist1,1)
        call prinf('mnlist2=*',mnlist2,1)
        call prinf('mnlist3=*',mnlist3,1)
        call prinf('mnlist4=*',mnlist4,1)
      endif 
        
      allocate(list1(mnlist1,nboxes),nlist1s(nboxes))
      allocate(list2(mnlist2,nboxes),nlist2s(nboxes))
      allocate(list3(mnlist3,nboxes),nlist3s(nboxes))
      allocate(list4(mnlist4,nboxes),nlist4s(nboxes))

      call computelists(nlevels,nboxes,itree(iptr(1)),boxsize,
     1  centers,itree(iptr(3)),itree(iptr(4)),
     1  itree(iptr(5)),isep,itree(iptr(6)),mnbors,
     1  itree(iptr(7)),iper,nlist1s,mnlist1,list1,nlist2s,
     1  mnlist2,list2,nlist3s,mnlist3,list3,
     1  nlist4s,mnlist4,list4)

      

      nd = 1
      call mpalloc(nd,itree(iptr(1)),iaddr,nlevels,lmptot,nboxes,
     1    nterms)

c      call prinf('iaddr=*',iaddr(1,nboxes),1)
c      call prinf('iaddr=*',iaddr(2,nboxes),1)
c      call prinf('lmptot=*',lmptot,1)
c      stop 

      allocate(rmlexp(lmptot))
      rmlexp = 0
      

      if (nlevels.gt.0) nmax=nterms(0)*2
      allocate(amat(0:nmax,-nmax:nmax))
      amat = 0 
      call cpu_time(t1)
      call get_amat(nmax,amat)
      call cpu_time(t2)
      if (ifprint.eq.1) then 
        call prin2('computing binomial coefficients spends * secs', 
     1 t2-t1, 1)
      endif 
c
cc       ... step 1, locate all charges, assign them to boxes, and
cc       form multipole expansions
c 
      call cpu_time(t1)
C$        t1=omp_get_wtime()
      do ilev=1,nlevels
        nterm = nterms(ilev)
        rscale = rscales(ilev)
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(istart,iend,npts)
C$OMP$SCHEDULE(DYNAMIC)
        do ibox=laddr(1,ilev),laddr(2,ilev) 
          istart = isrcse(1,ibox)
          iend = isrcse(2,ibox)
          npts = iend-istart+1
          if (nchild(ibox).eq.0 .and. npts.gt.0) then  
            call zl3dformmp(zsrcsort(1,istart),npts,chargesort(istart),
     1 dipstrsort(istart),dipvecsort(1,istart),zcsrc(1,ibox),rscale,
     2 nterm,rmlexp(iaddr(1,ibox)),ifcharge,ifdipole)   
          endif 
        enddo 
C$OMP END PARALLEL DO
      enddo
      call cpu_time(t2)
C$        t2=omp_get_wtime()
      if (ifprint.eq.1) call prin2('step 1 spends * secs', t2-t1, 1)

c
cc       ... step 2, for every parant box, form multipole expansions
cc       from their children 
c

      call cpu_time(t1)
      do ilev=nlevels-1,1,-1 
        nterm = nterms(ilev)
        rscale = rscales(ilev)
        width = boxsize(ilev)
        nterm1 = nterms(ilev+1)
        rscale1 = rscales(ilev+1)
        width1 = boxsize(ilev+1)
C$        t1=omp_get_wtime()
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,j,jbox)
C$OMP$SCHEDULE(DYNAMIC)
        do ibox=laddr(1,ilev),laddr(2,ilev)
          if (nchild(ibox).gt.0 .and.
     1 isrcse(2,ibox).ge.isrcse(1,ibox)) then
            do j = 1,nchild(ibox)
              jbox = ichild(j,ibox)
              if (isrcse(2,jbox).ge.isrcse(1,jbox)) then 
                call zl3dmpmp(zcsrc(1,jbox),nterm1,
     1 rmlexp(iaddr(1,jbox)),rscale1,zcsrc(1,ibox),nterm,
     2 rmlexp(iaddr(1,ibox)),rscale,nmax,amat)
              endif 
            enddo 
          endif  
        enddo 
C$OMP END PARALLEL DO
      enddo 
      call cpu_time(t2)
C$        t2=omp_get_wtime()
      if (ifprint.eq.1) call prin2('step 2 spends * secs', t2-t1, 1)

C
CC       ... step 3, for all childless boxes, directly compute interactions 
CC           with all boxes in list1 
C
        
      call cpu_time(t1)
C$        t1=omp_get_wtime()
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,istart,iend,nptstarg,ilist)
C$OMP$PRIVATE(jbox,jstart,jend,nptssrc)
C$OMP$SCHEDULE(DYNAMIC)
      do ibox=1,nboxes
        istart = itargse(1,ibox)
        iend = itargse(2,ibox)
        nptstarg = iend-istart+1
        if (nptstarg.gt.0 .and.
     1     nlist1s(ibox).gt.0 .and. nchild(ibox).eq.0) then 
          do ilist=1,nlist1s(ibox)
            jbox = list1(ilist,ibox)
            jstart = isrcse(1,jbox)
            jend = isrcse(2,jbox)
            nptssrc = jend-jstart+1
            if (nptssrc.gt.0) then 
              call zl3devaldirect(nptstarg,ztargsort(1,istart),
     1 nptssrc,zsrcsort(1,jstart),chargesort(jstart),
     2 dipstrsort(jstart),dipvecsort(1,jstart),ifcharge,
     3 potsort(istart),gradsort(1,istart),ifdipole,ifpgh)
            endif 
          enddo 
        endif 
      enddo 
C$OMP END PARALLEL DO
      call cpu_time(t2)
C$        t2=omp_get_wtime()
      if (ifprint.eq.1) call prin2('step 3 spends * secs', t2-t1, 1)
C
CC       ... step 4, for each box, convert the multipole expansions of all 
CC             boxes in list2 into the its local expansion  
C

      call cpu_time(t1)
C$        t1=omp_get_wtime()
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,rscale,nterm,width,ilist,jbox)
C$OMP$SCHEDULE(DYNAMIC)
      do ibox=1,nboxes
        rscale = rscales(ilevel(ibox))
        nterm = nterms(ilevel(ibox))
        width = boxsize(ilevel(ibox))
        if ((itargse(2,ibox).ge.itargse(1,ibox)) .and.
     1 (nlist2s(ibox).gt.0)) then 
          do ilist=1,nlist2s(ibox) 
            jbox = list2(ilist,ibox)
            if (isrcse(2,jbox).ge.isrcse(1,jbox)) then   
              call zl3dmploc(zcsrc(1,jbox),nterm,
     1 rmlexp(iaddr(1,jbox)),rscale,zctarg(1,ibox),nterm,
     2 rmlexp(iaddr(2,ibox)),rscale,nmax,amat)   
            endif 
          enddo 
        endif 
      enddo 
C$OMP END PARALLEL DO
      call cpu_time(t2)
C$        t2=omp_get_wtime()
      if (ifprint.eq.1) call prin2('step 4 spends * secs', t2-t1, 1)

c
cc       ... step 5, for each childless box, evaluate the multipole 
cc             expansion from boxes in list3 
c
      
      call cpu_time(t1)
C$        t1=omp_get_wtime()
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox, istart, iend, npts)
C$OMP$PRIVATE(ilist, jbox)
C$OMP$PRIVATE(rscale, nterm)
C$OMP$SCHEDULE(DYNAMIC)
      do ibox=1,nboxes
        istart = itargse(1,ibox)
        iend = itargse(2,ibox)
        npts = iend-istart+1
        if (npts.gt.0 .and. nchild(ibox).eq.0 .and.
     1 nlist3s(ibox).gt.0) then 
          do ilist = 1,nlist3s(ibox)
            jbox = list3(ilist,ibox)
            if (isrcse(2,jbox).ge.isrcse(1,jbox)) then 
              rscale = rscales(ilevel(jbox))
              nterm = nterms(ilevel(jbox))
              call zl3devalmp(npts,ztargsort(1,istart),nterm,
     1 zcsrc(1,jbox),rscale,rmlexp(iaddr(1,jbox)),
     2 potsort(istart),gradsort(1,istart),ifpgh)
            endif 
          enddo 
        endif 
      enddo 
C$OMP END PARALLEL DO
      call cpu_time(t2)
C$        t2=omp_get_wtime()
      if (ifprint.eq.1) call prin2('step 5 spends * secs', t2-t1, 1)

C
CC       ... step 6, for each box, form its local expansion from boxes 
CC             in list4 
C

      call cpu_time(t1)
C$        t1=omp_get_wtime()
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,nterm,rscale)
C$OMP$PRIVATE(ilist,jbox,jstart,jend,npts)
C$OMP$SCHEDULE(DYNAMIC)
      do ibox=1,nboxes
        nterm = nterms(ilevel(ibox))
        rscale = rscales(ilevel(ibox))
        if (itargse(2,ibox).ge.itargse(1,ibox) .and.
     1 nlist4s(ibox).gt.0) then
          do ilist=1,nlist4s(ibox)
            jbox = list4(ilist,ibox)
            jstart = isrcse(1,jbox)
            jend = isrcse(2,jbox)
            npts = jend-jstart+1
            if (npts.gt.0) then 
              call zl3dformloc(zsrcsort(1,jstart),npts,
     1 chargesort(jstart),dipstrsort(jstart),dipvecsort(1,jstart),
     2 zctarg(1,ibox),rscale,nterm,rmlexp(iaddr(2,ibox)),
     3 ifcharge,ifdipole)
            endif 
          enddo 
        endif 
      enddo 
C$OMP END PARALLEL DO
      call cpu_time(t2)
C$        t2=omp_get_wtime()
      if (ifprint.eq.1) call prin2('step 6 spends * secs', t2-t1, 1)

C
CC       ... step 7, for each parent box, shift its local expansion 
CC             to its children 
C

      call cpu_time(t1)
C$        t1=omp_get_wtime()
      do ilev=1,nlevels-1
        nterm = nterms(ilev)
        rscale = rscales(ilev)
        width = boxsize(ilev)
        nterm1 = nterms(ilev+1)
        rscale1 = rscales(ilev+1)
        width1 = boxsize(ilev+1)
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,j,jbox)
C$OMP$SCHEDULE(DYNAMIC)
        do ibox=laddr(1,ilev),laddr(2,ilev)
          if (itargse(2,ibox).ge.itargse(1,ibox) .and.
     1 (nlist2s(ibox).gt.0 .or.  nlist4s(ibox).gt.0)) then
            do j=1,nchild(ibox)
              jbox = ichild(j,ibox)
              if (itargse(2,jbox).ge.itargse(1,jbox)) then 
                call zl3dlocloc(zctarg(1,ibox),nterm,
     1 rmlexp(iaddr(2,ibox)),rscale,zctarg(1,jbox),nterm1,
     2 rmlexp(iaddr(2,jbox)),rscale1,nmax,amat)
              endif 
            enddo 
          endif 
        enddo 
C$OMP END PARALLEL DO
      enddo 
      call cpu_time(t2)
C$        t2=omp_get_wtime()
      if (ifprint.eq.1) call prin2('step 7 spends * secs', t2-t1, 1)

C
CC       ... step 8, for each childless box, add its local expansion 
CC          to the potential 
C
      call cpu_time(t1)
C$        t1=omp_get_wtime()
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,istart,iend,npts,rscale,nterm)
C$OMP$SCHEDULE(DYNAMIC)
      do ibox=1,nboxes 
        istart = itargse(1,ibox)
        iend = itargse(2,ibox)
        npts = iend-istart+1
        if (nchild(ibox).eq.0 .and. npts.gt.0) then
          rscale = rscales(ilevel(ibox))
          nterm = nterms(ilevel(ibox))
          call zl3devalloc(npts,ztargsort(1,istart),nterm,
     1 zctarg(1,ibox),rscale,rmlexp(iaddr(2,ibox)),
     2 potsort(istart),gradsort(1,istart),ifpgh)
        endif 
      enddo 
C$OMP END PARALLEL DO
      call cpu_time(t2)
C$        t2=omp_get_wtime()

      if (ifprint.eq.1) call prin2('step 8 spends * secs', t2-t1, 1)
      call zreorderi(1,nt,potsort,pot,itarg) 
      if (ifpgh.ge.2) call zreorderi(3,nt,gradsort,grad,itarg) 
      
      end 
C
C
C
C
c 
      subroutine zl3devaldirect(nt,ztarg,ns,zsrc,charge,dipstr,dipvec,
     1 ifcharge,pot,grad,ifdipole,ifpgh)
      implicit real *8 (a-h,o-z)
      complex *16 ztarg(3,nt),pot(nt),grad(3,nt)
      complex *16 zsrc(3,ns),charge(ns),dipstr(ns)
      complex *16 dipvec(3,ns)
      complex *16 rz3,y1mx1,y2mx2,y3mx3,dot1,dot2,dot3,dot 
      complex *16 rz,rz_pi4_inv,rz3_pi4_inv,rz5_pi4_inv
c      complex *16 x1my1,x2my2,x3my3

      done = 1 
      pi = atan(done)*4
      pi4 = pi*4 
      do i=1,nt 
        do j=1,ns 
c          call zdist(3,ztargsort(1,i),zsrcsort(1,j),rz)
          y1mx1 = zsrc(1,j)-ztarg(1,i)
          y2mx2 = zsrc(2,j)-ztarg(2,i)
          y3mx3 = zsrc(3,j)-ztarg(3,i)
          rz = sqrt(y1mx1**2+y2mx2**2+y3mx3**2)
          rz_pi4_inv = done/(pi4*rz)
          if ((ifdipole.eq.1).or.(ifpgh.ge.2)) then 
            rz3_pi4_inv = done/(pi4*rz**3)
          endif 
          if (abs(rz).gt.1E-16) then 

            if (ifcharge.eq.1) then 
              pot(i) = pot(i)+charge(j)*rz_pi4_inv
              if (ifpgh.ge.2) then
                grad(1,i) = grad(1,i)+charge(j)*y1mx1*rz3_pi4_inv
                grad(2,i) = grad(2,i)+charge(j)*y2mx2*rz3_pi4_inv
                grad(3,i) = grad(3,i)+charge(j)*y3mx3*rz3_pi4_inv
              endif 
            endif 

            if (ifdipole.eq.1) then
              dot1 = -y1mx1*dipvec(1,j)
              dot2 = -y2mx2*dipvec(2,j)
              dot3 = -y3mx3*dipvec(3,j)
              dot = dot1+dot2+dot3 
              pot(i) = pot(i)-dipstr(j)*dot*rz3_pi4_inv 
              if (ifpgh.ge.2) then 
                dot = dot*3/(pi4*rz**5)
                grad(1,i) = grad(1,i)-dipstr(j)*
     1                       (dipvec(1,j)*rz3_pi4_inv+dot*y1mx1)
                grad(2,i) = grad(2,i)-dipstr(j)*
     1                       (dipvec(2,j)*rz3_pi4_inv+dot*y2mx2)
                grad(3,i) = grad(3,i)-dipstr(j)*
     1                       (dipvec(3,j)*rz3_pi4_inv+dot*y3mx3)
              endif        
            endif 

          endif 
        enddo
      enddo 
      end 
C        
C
C
C
C       
      subroutine zl3dformmp(zsrc,npts,charge,dipstr,dipvec,zcsrc,rscale,
     1 nterm,mp,ifcharge,ifdipole)
      implicit real *8 (a-h,o-z)
      complex *16 charge(npts),dipstr(npts)
      complex *16 zsrc(3,npts),dipvec(3,npts)
      complex *16 zcsrc(3),mp(0:nterm,-nterm:nterm)
      complex *16 z1,z2,z3,rz,rzxy,ctheta,stheta,cphi,sphi,eiphi,zr
      complex *16 ima
      complex *16 e_zr(3),e_th(3),e_ph(3),dot_zr,dot_th,dot_ph
      complex *16 dzr,dth,dph,dot_th_zr,dot_ph_zr
      complex *16, allocatable :: yl(:,:),dyl(:,:)
      complex *16, allocatable :: rzall(:),eiphiall(:)
      real *8, allocatable :: rat1(:,:),rat2(:,:)

      allocate(rat1(0:nterm,0:nterm),rat2(0:nterm,0:nterm))
      allocate(yl(0:nterm,0:nterm),dyl(0:nterm,0:nterm))
      allocate(rzall(0:nterm),eiphiall(-nterm:nterm))

      call ylgndrini(nterm, rat1, rat2)

      ima = (0.0d0,1.0d0)
      done = 1
      pi4 = atan(done)*4*4 
      pi4_sqrt_inv = done/sqrt(pi4)
      rzall(0) = 1 
      eiphiall(0) = 1 
      

      do i=1,npts 
        call zcart2sph(zsrc(1:3,i),zcsrc,zr,ctheta,stheta,cphi,sphi,
     1 eiphi)
        rz = zr/rscale


        do n=1,nterm 
          rzall(n) = rzall(n-1)*rz 
          eiphiall(n) = eiphiall(n-1)*eiphi 
          eiphiall(-n) = 1/eiphiall(n)
        enddo 


        if (ifdipole.eq.0) then 
          call zylgndrf(nterm,ctheta,yl,rat1,rat2)
        else 
          call zylgndr2sf(nterm,ctheta,yl,dyl,rat1,rat2)
c          call zylgndr2f(nterm,ctheta,yl,dyl,rat1,rat2)
        endif 

        do n=0,nterm
          do m=0,n 
            yl(n,m) = yl(n,m)/(2*n+done)*pi4_sqrt_inv*((-1)**m) 
            if (ifdipole.eq.1) then 
              dyl(n,m) = dyl(n,m)/(2*n+done)*pi4_sqrt_inv*((-1)**m)
            endif 
          enddo 
        enddo  

        do n=0,nterm 
          mp(n,0) = mp(n,0)+charge(i)*rzall(n)*yl(n,0)
          do m=1,n 
            if (ifdipole.eq.0) then 
              mp(n,m) = mp(n,m)+charge(i)*rzall(n)*
     1                 yl(n,abs(m))*eiphiall(-m)
              mp(n,-m) = mp(n,-m)+charge(i)*rzall(n)*
     1                 yl(n,abs(m))*eiphiall(m)
            else 
              mp(n,m) = mp(n,m)+charge(i)*rzall(n)*stheta*
     1                 yl(n,abs(m))*eiphiall(-m)
              mp(n,-m) = mp(n,-m)+charge(i)*rzall(n)*stheta*
     1                 yl(n,abs(m))*eiphiall(m)
            endif 
          enddo 
        enddo 

c        do n=0,nterm 
c          do m=-n,n
c
c            if (ifcharge.eq.1) then 
c              mp(n,m) = mp(n,m)+charge(i)*rzall(n)*yl(n,abs(m))*
c     1                             eiphiall(-m)
c            endif 
c
c            if (ifdipole.eq.1) then 
c
c              if (n.eq.0) then 
c                dzr = 0 
c              else 
c                dzr = n*rzall(n-1)/rscale*yl(n,abs(m))*eiphiall(-m)
c              endif  
c
c              dth = -rzall(n)*dyl(n,abs(m))*stheta*eiphiall(-m)
c              dph = -ima*m*rzall(n)*yl(n,abs(m))*eiphiall(-m)
c
c              mp(n,m) = mp(n,m)-dipstr(i)*
c     1                         (dot_zr*dzr+dot_th*dth+dot_ph*dph)
c            endif 
c            
c          enddo 
c        enddo 


        if (ifdipole.eq.1) then 

          e_zr(1) = stheta*cphi
          e_zr(2) = stheta*sphi
          e_zr(3) = ctheta

          e_th(1) = ctheta*cphi
          e_th(2) = ctheta*sphi
          e_th(3) = -stheta

          e_ph(1) = -sphi
          e_ph(2) = cphi
          e_ph(3) = 0

          dot_zr = 0 
          dot_th = 0 
          dot_ph = 0 
          do idim=1,3 
            dot_zr = dot_zr+dipvec(idim,i)*e_zr(idim)
            dot_th = dot_th+dipvec(idim,i)*e_th(idim)
            dot_ph = dot_ph+dipvec(idim,i)*e_ph(idim)
          enddo
c          call prin2('dot_zr=*',dot_zr,2)
c          call prin2('dot_th=*',dot_th,2)
c          call prin2('dot_ph=*',dot_ph,2)

c          dot_th = dot_th/zr
c          dot_ph = dot_ph/zr/stheta
          dot_th_zr = dot_th/zr 
          dot_ph_zr = dot_ph/zr

          do n=0,nterm 

            if (n.eq.0) then 
              dzr = 0 
            else 
              dzr = n*rzall(n-1)/rscale*yl(n,0)
            endif 
            dth = -rzall(n)*dyl(n,0)*stheta
            mp(n,0) = mp(n,0)-dipstr(i)*(dot_zr*dzr+dot_th_zr*dth)

            do m=1,n

              dzr = n*rzall(n-1)/rscale*yl(n,abs(m))
     1                    *stheta*eiphiall(-m)
              dth = -rzall(n)*dyl(n,abs(m))*eiphiall(-m)
              dph = -ima*m*rzall(n)*yl(n,abs(m))*eiphiall(-m)
              mp(n,m) = mp(n,m)-dipstr(i)*
     1              (dot_zr*dzr+dot_th_zr*dth+dot_ph_zr*dph)

              dzr = n*rzall(n-1)/rscale*yl(n,abs(m))
     1                     *stheta*eiphiall(m)
              dth = -rzall(n)*dyl(n,abs(m))*eiphiall(m)
              dph = ima*m*rzall(n)*yl(n,abs(m))*eiphiall(m)
              mp(n,-m) = mp(n,-m)-dipstr(i)*
     1             (dot_zr*dzr+dot_th_zr*dth+dot_ph_zr*dph)
            enddo 
          enddo 
        endif 
      enddo 

      end 
C
C
C
C
C
      subroutine zl3dformloc(zsrc,npts,charge,dipstr,dipvec,zctarg,
     1 rscale,nterm,loc,ifcharge,ifdipole)
      implicit real *8 (a-h,o-z)
      complex *16 charge(npts),dipstr(npts)
      complex *16 zsrc(3,npts),dipvec(3,npts)
      complex *16 zctarg(3),loc(0:nterm,-nterm:nterm)
      complex *16 z1,z2,z3,rz,rzxy,ctheta,stheta,cphi,sphi,eiphi,zr
      complex *16 ima
      complex *16 e_zr(3),e_th(3),e_ph(3),dot_zr,dot_th,dot_ph
      complex *16 dzr,dth,dph,dot_th_zr,dot_ph_zr
      complex *16, allocatable :: yl(:,:),dyl(:,:)
      complex *16, allocatable :: rzall(:),eiphiall(:)
      real *8, allocatable :: rat1(:,:),rat2(:,:)

      allocate(rat1(0:nterm,0:nterm),rat2(0:nterm,0:nterm))
      allocate(yl(0:nterm,0:nterm),dyl(0:nterm,0:nterm))
      allocate(rzall(-nterm-1:nterm+1),eiphiall(-nterm:nterm))

      call ylgndrini(nterm, rat1, rat2)

      ima = (0.0d0,1.0d0)
      done = 1
      pi4 = atan(done)*4*4 
      pi4_sqrt_inv = done/sqrt(pi4)
      rzall(0) = 1 
      eiphiall(0) = 1 

      do i=1,npts 
        call zcart2sph(zsrc(1:3,i),zctarg,zr,ctheta,stheta,cphi,sphi,
     1 eiphi)
        rz = zr/rscale

        do n=1,nterm 
          rzall(n) = rzall(n-1)*rz 
          eiphiall(n) = eiphiall(n-1)*eiphi 
          eiphiall(-n) = 1/eiphiall(n)
          rzall(-n) = 1/rzall(n)
        enddo 
        rzall(nterm+1) = rzall(nterm)*rz
        rzall(-nterm-1) = 1/rzall(nterm+1)

        if (ifdipole.eq.0) then 
          call zylgndrf(nterm,ctheta,yl,rat1,rat2)
        else 
          call zylgndr2sf(nterm,ctheta,yl,dyl,rat1,rat2)
        endif 

        do n=0,nterm 
          do m=0,n 
            yl(n,m) = yl(n,m)*pi4_sqrt_inv*((-1)**m)/(2*n+done)
            if (ifdipole.eq.1) then 
              dyl(n,m) = dyl(n,m)*pi4_sqrt_inv*((-1)**m)/(2*n+done)
            endif 
          enddo 
        enddo 

c        do n=0,nterm 
c          do m=-n,n
c
c            if (ifcharge.eq.1) then 
c              loc(n,m) = loc(n,m)+charge(i)*yl(n,abs(m))
c     1             *eiphiall(-m)*rzall(-n-1)
c            endif 
c
c            if (ifdipole.eq.1) then 
c              dzr = -(n+1)*yl(n,abs(m))*eiphiall(-m)*rzall(-n-1)/zr 
c              dth = -dyl(n,abs(m))*stheta*eiphiall(-m)*rzall(-n-1)
c              dph = -m*ima*yl(n,abs(m))*eiphiall(-m)*rzall(-n-1)
c
c              loc(n,m) = loc(n,m)-dipstr(i)*
c     1                         (dot_zr*dzr+dot_th*dth+dot_ph*dph)
c            endif 
c
c          enddo 
c        enddo 

        do n=0,nterm 
          loc(n,0) = loc(n,0)+charge(i)*rzall(-n-1)*yl(n,0)
          do m=1,n 
            if (ifdipole.eq.0) then 
              loc(n,m) = loc(n,m)+charge(i)*rzall(-n-1)*
     1                 yl(n,abs(m))*eiphiall(-m)
              loc(n,-m) = loc(n,-m)+charge(i)*rzall(-n-1)*
     1                 yl(n,abs(m))*eiphiall(m)
            else 
              loc(n,m) = loc(n,m)+charge(i)*rzall(-n-1)*stheta*
     1                 yl(n,abs(m))*eiphiall(-m)
              loc(n,-m) = loc(n,-m)+charge(i)*rzall(-n-1)*stheta*
     1                 yl(n,abs(m))*eiphiall(m)
            endif 
          enddo 
        enddo 

        if (ifdipole.eq.1) then 

          e_zr(1) = stheta*cphi
          e_zr(2) = stheta*sphi
          e_zr(3) = ctheta

          e_th(1) = ctheta*cphi
          e_th(2) = ctheta*sphi
          e_th(3) = -stheta

          e_ph(1) = -sphi
          e_ph(2) = cphi
          e_ph(3) = 0

          dot_zr = 0 
          dot_th = 0 
          dot_ph = 0 
          do idim=1,3 
            dot_zr = dot_zr+dipvec(idim,i)*e_zr(idim)
            dot_th = dot_th+dipvec(idim,i)*e_th(idim)
            dot_ph = dot_ph+dipvec(idim,i)*e_ph(idim)
          enddo

          dot_th_zr = dot_th/zr
          dot_ph_zr = dot_ph/zr

          do n=0,nterm 

            dzr = -(n+1)*yl(n,0)*rzall(-n-1)/zr 
            dth = -dyl(n,0)*stheta*rzall(-n-1)
            loc(n,0) = loc(n,0)-dipstr(i)*
     1          (dot_zr*dzr+dot_th_zr*dth)

            do m=1,n 
              dzr = -(n+1)*yl(n,abs(m))*stheta*
     1                 eiphiall(-m)*rzall(-n-1)/zr 
              dth = -dyl(n,abs(m))*eiphiall(-m)*rzall(-n-1)
              dph = -m*ima*yl(n,abs(m))*eiphiall(-m)*rzall(-n-1)
              loc(n,m) = loc(n,m)-dipstr(i)*
     1               (dot_zr*dzr+dot_th_zr*dth+dot_ph_zr*dph)


              dzr = -(n+1)*yl(n,abs(m))*stheta*
     1                 eiphiall(m)*rzall(-n-1)/zr 
              dth = -dyl(n,abs(m))*eiphiall(m)*rzall(-n-1)
              dph = m*ima*yl(n,abs(m))*eiphiall(m)*rzall(-n-1)
              loc(n,-m) = loc(n,-m)-dipstr(i)*
     1               (dot_zr*dzr+dot_th_zr*dth+dot_ph_zr*dph)
            enddo 

          enddo 

        endif 

      enddo 
      end 
C
C
C
C
C
      subroutine zl3devalmp(nt,ztarg,nterm,zc,rscale,mp,pot,grad,ifpgh)
      implicit real *8 (a-h,o-z)
      complex *16 ztarg(3,nt),pot(nt),grad(3,nt)
      complex *16 zc(3),mp(0:nterm,-nterm:nterm)
      complex *16 z1,z2,z3,rz,rzxy,ctheta,stheta,cphi,sphi,eiphi,zr
      complex *16 ima
      complex *16 e_zr(3),e_th(3),e_ph(3),e_th_zr(3),e_ph_zr(3)
      complex *16 dzr,dth,dph
      complex *16, allocatable :: yl(:,:),dyl(:,:)
      complex *16, allocatable :: rzall(:),eiphiall(:)
      real *8, allocatable :: rat1(:,:),rat2(:,:)

      allocate(yl(0:nterm,0:nterm),dyl(0:nterm,0:nterm))
      allocate(rzall(0:nterm+1),eiphiall(-nterm:nterm))
      allocate(rat1(0:nterm,0:nterm),rat2(0:nterm,0:nterm))

      done = 1
      ima = (0.0d0,1.0d0)
      pi4 = atan(done)*4*4 
      pi4_sqrt_inv = done/sqrt(pi4)
      rzall(0) = 1 
      eiphiall(0) = 1 

      call ylgndrini(nterm,rat1,rat2)

      do i=1,nt 
        call zcart2sph(ztarg(1:3,i),zc,zr,ctheta,stheta,cphi,sphi,
     1 eiphi)
        rz = zr/rscale

        do n=1,nterm 
          rzall(n) = rzall(n-1)*rz 
          eiphiall(n) = eiphiall(n-1)*eiphi 
          eiphiall(-n) = 1/eiphiall(n)
        enddo 
        rzall(nterm+1) = rzall(nterm)*rz

        if (ifpgh.eq.1) then 
          call zylgndrf(nterm,ctheta,yl,rat1,rat2)
        else 
          call zylgndr2sf(nterm,ctheta,yl,dyl,rat1,rat2)
          e_zr(1) = stheta*cphi
          e_zr(2) = stheta*sphi
          e_zr(3) = ctheta

          e_th(1) = ctheta*cphi
          e_th(2) = ctheta*sphi
          e_th(3) = -stheta
          e_th_zr = e_th/zr 

          e_ph(1) = -sphi
          e_ph(2) = cphi
          e_ph(3) = 0
          e_ph_zr = e_ph/zr 
        endif 

        do n=0,nterm 
          do m=0,n 
            yl(n,m) = yl(n,m)*pi4_sqrt_inv*((-1)**m)/rscale
            if (ifpgh.ge.2) then 
              dyl(n,m) = dyl(n,m)*pi4_sqrt_inv*((-1)**m)/rscale
            endif 
          enddo 
        enddo 

c        do n=0,nterm 
c          do m=-n,n 
c            pot(i) = pot(i)+mp(n,m)*yl(n,abs(m))*eiphiall(m)
c     1                                       /rzall(n+1)
c            if (ifpgh.ge.2) then 
c              dzr = -(n+1)*yl(n,abs(m))*eiphiall(m)/rzall(n+1)/zr
c              dth = -dyl(n,abs(m))*stheta*eiphiall(m)/rzall(n+1)
c              dph = m*ima*yl(n,abs(m))*eiphiall(m)/rzall(n+1)
c
c              grad(1:3,i) = grad(1:3,i)+mp(n,m)*
c     1         (dzr*e_zr+dth*e_th/zr+dph*e_ph/zr/stheta)
c            endif 
c          enddo 
c        enddo 
c      enddo 

        do n=0,nterm 
          pot(i) = pot(i)+mp(n,0)*yl(n,0)/rzall(n+1)
          do m=1,n 
            if (ifpgh.eq.1) then
              pot(i) = pot(i)+yl(n,abs(m))*
     1  (mp(n,m)*eiphiall(m)+mp(n,-m)*eiphiall(-m))/rzall(n+1)
            else 
              pot(i) = pot(i)+yl(n,abs(m))*stheta*
     1  (mp(n,m)*eiphiall(m)+mp(n,-m)*eiphiall(-m))/rzall(n+1)
            endif 
          enddo 
        enddo 

        if (ifpgh.ge.2) then 
          do n=0,nterm 
            grad(1:3,i) = grad(1:3,i)+mp(n,0)*
     1         (-(n+1)*yl(n,0)/rzall(n+1)/zr*e_zr
     1             -dyl(n,0)*stheta/rzall(n+1)*e_th_zr)
            do m=1,n 
              dzr = -(n+1)*yl(n,abs(m))*stheta*eiphiall(m)
     1             /rzall(n+1)/zr
              dth = -dyl(n,abs(m))*eiphiall(m)/rzall(n+1)
              dph = m*ima*yl(n,abs(m))*eiphiall(m)/rzall(n+1)

              grad(1:3,i) = grad(1:3,i)+mp(n,m)*
     1         (dzr*e_zr+dth*e_th_zr+dph*e_ph_zr)


              dzr = -(n+1)*yl(n,abs(m))*stheta*eiphiall(-m)
     1              /rzall(n+1)/zr
              dth = -dyl(n,abs(m))*eiphiall(-m)/rzall(n+1)
              dph = -m*ima*yl(n,abs(m))*eiphiall(-m)/rzall(n+1)

              grad(1:3,i) = grad(1:3,i)+mp(n,-m)*
     1         (dzr*e_zr+dth*e_th_zr+dph*e_ph_zr)
            enddo 
          enddo 
        endif 
      enddo 
      end 
C
C
C
C
C
      subroutine zl3devalloc(nt,ztarg,nterm,zc,rscale,loc,pot,grad,
     1 ifpgh)
      implicit real *8 (a-h,o-z)
      complex *16 ztarg(3,nt),pot(nt),grad(3,nt)
      complex *16 zc(3),loc(0:nterm,-nterm:nterm)
      complex *16 z1,z2,z3,rz,rzxy,ctheta,stheta,cphi,sphi,eiphi,zr
      complex *16 ima
      complex *16 e_zr(3),e_th(3),e_ph(3),e_th_zr(3),e_ph_zr(3)
      complex *16 dzr,dth,dph
      complex *16, allocatable :: yl(:,:),dyl(:,:)
      complex *16, allocatable :: rzall(:),eiphiall(:)
      real *8, allocatable :: rat1(:,:),rat2(:,:)

      allocate(yl(0:nterm,0:nterm),dyl(0:nterm,0:nterm))
      allocate(rzall(0:nterm),eiphiall(-nterm:nterm))
      allocate(rat1(0:nterm,0:nterm),rat2(0:nterm,0:nterm))

      done = 1
      ima = (0.0d0,1.0d0)
      pi4 = atan(done)*4*4 
      pi4_sqrt_inv = done/sqrt(pi4)
      rzall(0) = 1 
      eiphiall(0) = 1 

      call ylgndrini(nterm,rat1,rat2)

      do i=1,nt 
        call zcart2sph(ztarg(1:3,i),zc,zr,ctheta,stheta,cphi,sphi,
     1 eiphi)
        rz = zr/rscale

        do n=1,nterm 
          rzall(n) = rzall(n-1)*rz 
          eiphiall(n) = eiphiall(n-1)*eiphi 
          eiphiall(-n) = 1/eiphiall(n)
        enddo 

        if (ifpgh.eq.1) then 
          call zylgndrf(nterm,ctheta,yl,rat1,rat2)
        else 
          call zylgndr2sf(nterm,ctheta,yl,dyl,rat1,rat2)
          e_zr(1) = stheta*cphi
          e_zr(2) = stheta*sphi
          e_zr(3) = ctheta

          e_th(1) = ctheta*cphi
          e_th(2) = ctheta*sphi
          e_th(3) = -stheta
          e_th_zr = e_th/zr 

          e_ph(1) = -sphi
          e_ph(2) = cphi
          e_ph(3) = 0
          e_ph_zr = e_ph/zr 
        endif 

        do n=0,nterm 
          do m=0,n 
            yl(n,m) = yl(n,m)*pi4_sqrt_inv*((-1)**m)/rscale
            if (ifpgh.ge.2) then 
              dyl(n,m) = dyl(n,m)*pi4_sqrt_inv*((-1)**m)/rscale
            endif 
          enddo 
        enddo 

c        do n=0,nterm 
c          do m=-n,n 
c            pot(i) = pot(i)+loc(n,m)*rzall(n)*yl(n,abs(m))*eiphiall(m)
c
c            if (ifpgh.ge.2) then 
c              if (n.eq.0) then 
c                dzr = 0 
c              else  
c                dzr = n*rzall(n-1)/rscale*yl(n,abs(m))*eiphiall(m)
c              endif 
c              dth = -rzall(n)*dyl(n,abs(m))*stheta*eiphiall(m)
c              dph = m*ima*rzall(n)*yl(n,abs(m))*eiphiall(m)
c
c              grad(1:3,i) = grad(1:3,i)+loc(n,m)*
c     1         (dzr*e_zr+dth*e_th/zr+dph*e_ph/zr/stheta)
c            endif      
c          enddo 
c        enddo 

        do n=0,nterm 
          pot(i) = pot(i)+loc(n,0)*rzall(n)*yl(n,0)
          do m=1,n 
            if (ifpgh.eq.1) then 
              pot(i) = pot(i)+yl(n,abs(m))*
     1  (loc(n,m)*eiphiall(m)+loc(n,-m)*eiphiall(-m))*rzall(n)
            else 
              pot(i) = pot(i)+yl(n,abs(m))*stheta*
     1  (loc(n,m)*eiphiall(m)+loc(n,-m)*eiphiall(-m))*rzall(n)
            endif 
          enddo 
        enddo 
      

        if (ifpgh.ge.2) then 
          do n=0,nterm 
            if (n.eq.0) then 
              dzr = 0 
            else 
              dzr = n*rzall(n-1)/rscale*yl(n,0)
            endif  
            dth = -rzall(n)*dyl(n,0)*stheta
            grad(1:3,i) = grad(1:3,i)+loc(n,0)*
     1         (dzr*e_zr+dth*e_th_zr)


            do m=1,n

              if (n.eq.0) then 
                dzr = 0 
              else 
                dzr = n*rzall(n-1)/rscale*yl(n,abs(m))
     1                        *stheta*eiphiall(m)
              endif 
              dth = -rzall(n)*dyl(n,abs(m))*eiphiall(m)
              dph = m*ima*rzall(n)*yl(n,abs(m))*eiphiall(m)
              grad(1:3,i) = grad(1:3,i)+loc(n,m)*
     1         (dzr*e_zr+dth*e_th_zr+dph*e_ph_zr)




              if (n.eq.0) then 
                dzr = 0 
              else 
                dzr = n*rzall(n-1)/rscale*yl(n,abs(m))
     1                        *stheta*eiphiall(-m)
              endif 
              dth = -rzall(n)*dyl(n,abs(m))*eiphiall(-m)
              dph = -m*ima*rzall(n)*yl(n,abs(m))*eiphiall(-m)
              grad(1:3,i) = grad(1:3,i)+loc(n,-m)*
     1         (dzr*e_zr+dth*e_th_zr+dph*e_ph_zr)
            enddo 
          enddo 
        endif 
      enddo 
      end 
C
C
C
C
C
      subroutine zl3dterms(eps, nterms)
c***********************************************************************
c
c     Determine number of terms in mpole expansions.
c
c     The method is based on examining the decay of \rho^n / r^{n+1}
c     for \rho a worst case source and r a worst case target.
c
c     INPUT:
c
c     eps   tolerance
c
c     OUTPUT:
c
c     nterms  required expansion order
c-----------------------------------------------------------------------
c
      implicit real *8 (a-h,o-z)
C      integer *8 j,nterms
C      real *8 xtemp1,eps
      real *8 jfun
c
      z1 = 1.5d0
      z2 = sqrt(3d0)/2.d0
c
      jfun = z2
      hfun = 1.0d0/(z1**2)
      nterms = 1
      do j = 2, 1000
        hfun = hfun/z1
        jfun = jfun*z2
        xtemp1 = jfun*hfun
        if(xtemp1 .lt. eps)then
          nterms = j
          return
        endif
      enddo
      return
      end

