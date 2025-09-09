      subroutine zhfmm3d(eps,zk,ns,zsrc,ifcharge,charge,ifdipole,
     1 dipstr,dipvec,pot,grad,nt,ztarg,ifpgh,isep,ifprint)
c----------------------------------------------
      implicit real *8 (a-h,o-z)
      complex *16 zk,zkfmm,ima 
      complex *16 zsrc(3,ns),ztarg(3,nt)
      complex *16 charge(ns),pot(nt),grad(3,nt)
      complex *16 dipstr(ns),dipvec(3,ns)
      dimension src(3,ns),targ(3,nt)
      dimension srcsort(3,ns)
      dimension targsort(3,nt)
      dimension zc(3)
      complex *16 zsrcsort(3,ns),chargesort(ns),dipstrsort(ns)
      complex *16 dipvecsort(3,ns)
      complex *16 ztargsort(3,nt)
      integer *8 ltree,iptr(8)
      integer, allocatable :: itree(:)
      real *8, allocatable :: centers(:,:),boxsize(:)
      complex *16, allocatable :: zcsrc(:,:),zctarg(:,:)
      integer, allocatable :: isrc(:),isrcse(:,:)
      integer, allocatable :: itarg(:),itargse(:,:)
      real *8, allocatable :: rscales(:)
      integer, allocatable :: nterms(:)

      done = 1 
      pi = atan(done)*4 
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

      ifpghtarg = 1
      idivflag = 2

      call hndiv(eps,ns,nt,ifcharge,ifdipole,ifpgh,
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
      call pts_tree_sort(nt,targ,itree,ltree,nboxes,nlevels,iptr,
     1   centers,itarg,itargse)
      call pts_tree_sort(ns,src,itree,ltree,nboxes,nlevels,iptr,
     1   centers,isrc,isrcse)



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
      b0inv = 1.0d0/b0
      b0inv2 = b0inv**2
      if (ifprint.eq.1) then 
        call prin2('b0=*',b0,1)
        call prin2('b0inv=*',b0inv,1)
        call prin2('b0inv2=*',b0inv2,1)
      endif 

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
      zkfmm = zk*b0
c
cc      compute scaling factor for multipole/local expansions
c       and lengths of multipole and local expansions
c
      allocate(rscales(0:nlevels))
      rscales = 0 
      do ilev = 0,nlevels
        rscales(ilev) = boxsize(ilev)*abs(zkfmm)
        if(rscales(ilev).gt.1) then 
          rscales(ilev) = 1
        endif 
      enddo

      if (ifprint.eq.1) call prin2('rscales=*',rscales,nlevels+1)
c      print *, rscales(0:nlevels)


      allocate(nterms(0:nlevels))
      nterms = 0 
      nmax = 0
      do i=0,nlevels
        call h3dterms(boxsize(i),zkfmm,eps,nterms(i),isep)
        nterms(i) = nterms(i)+n0
        if(nterms(i).gt.nmax) nmax = nterms(i)
        if(i.gt.0) then 
          if(nterms(i).gt.nterms(i-1)) nterms(i)=nterms(i-1)
        endif 
      enddo



      if (ifprint.eq.1) then 
        call prinf('nmax=*',nmax,1)
        call prinf('nterms=*',nterms,nlevels+1)
      endif 

      


      call zhfmm3dmain(zkfmm,ns,zsrcsort,chargesort,dipstrsort,
     1 dipvecsort,nt,ztargsort,isrc,itarg,isrcse,itargse,itree,ltree,
     2 iptr,nlevels,nboxes,boxsize,rscales,nterms,nmax,centers,zcsrc,
     3 zctarg,itree(iptr(1)),itree(iptr(2)),itree(iptr(3)),
     4 itree(iptr(4)),itree(iptr(5)),itree(iptr(6)),itree(iptr(7)),
     5 pot,grad,ifcharge,ifdipole,ifpgh,ifprint,isep)

      if (ifpgh.ge.2) call zrescale(3*nt,grad,b0inv)

      end 
C
C
C
C
C
      subroutine zhfmm3dmain(zk,ns,zsrcsort,chargesort,dipstrsort,
     1 dipvecsort,nt,ztargsort,isrc,itarg,isrcse,itargse,itree,ltree,
     2 iptr,nlevels,nboxes,boxsize,rscales,nterms,nmax,centers,zcsrc,
     3 zctarg,laddr,ilevel,iparent,nchild,ichild,ncol,coll,pot,grad,
     4 ifcharge,ifdipole,ifpgh,ifprint,isep)
c   Helmholtz FMM in C^3: evaluate all pairwise particle
      implicit real *8 (a-h,o-z)
      complex *16 zk,ima,rz
      complex *16 zsrcsort(3,ns),chargesort(ns),dipstrsort(ns)
      complex *16 dipvecsort(3,ns)
      dimension isrc(ns),itarg(nt)
      dimension isrcse(2,nboxes), itargse(2,nboxes)
      complex *16 ztargsort(3,nt),pot(nt),potsort(nt)
      complex *16 gradsort(3,nt),grad(3,nt)
      integer *8 ltree,iptr(8),iaddr(2,nboxes),lmptot,isep
      dimension ilevel(nboxes),nchild(nboxes),ncol(nboxes)
      dimension iparent(nboxes),ichild(8,nboxes)
c      integer coll(mnbors,nboxes)
      dimension itree(ltree)
      dimension boxsize(0:nlevels), rscales(0:nlevels)
      dimension nterms(0:nlevels)
      dimension centers(3,nboxes)
      complex *16 zcsrc(3,nboxes),zctarg(3,nboxes)
      integer laddr(2,0:nlevels)
      integer, allocatable :: list1(:,:),list2(:,:),list3(:,:)
      integer, allocatable :: list4(:,:)
      integer, allocatable :: nlist1s(:),nlist2s(:),nlist3s(:)
      integer, allocatable :: nlist4s(:)
      complex *16, allocatable :: rmlexp(:)

      ima = (0.0d0,1.0d0)
      done = 1 
      pi = 4*atan(done)
C
CC       initialize all variables 
C
      pot = 0 
      potsort = 0 
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
    

C     allocate space for multipole and local expansion
      nd = 1 
      call mpalloc(nd,itree(iptr(1)),iaddr,nlevels,lmptot,nboxes,
     1    nterms)
      allocate(rmlexp(lmptot))
      rmlexp = 0

      
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
C$OMP$PRIVATE(ibox,istart,iend,npts)
C$OMP$SCHEDULE(DYNAMIC)
        do ibox=laddr(1,ilev),laddr(2,ilev) 
          istart = isrcse(1,ibox)
          iend = isrcse(2,ibox)
          npts = iend-istart+1 
          if (nchild(ibox).eq.0 .and. npts.gt.0) then
            call zh3dformmp(zk,rscale,zsrcsort(1,istart),npts,
     1 chargesort(istart),dipstrsort(istart),dipvecsort(1,istart),
     2 zcsrc(1,ibox),nterm,rmlexp(iaddr(1,ibox)),ifcharge,ifdipole)
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
                call zh3dmpmp(zcsrc(1,jbox),
     1 zcsrc(1,ibox),nterm1,nterm,rscale1,rscale,width1,width,zk,
     2 rmlexp(iaddr(1,jbox)),rmlexp(iaddr(1,ibox)))
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
     1  nlist1s(ibox).gt.0. .and. nchild(ibox).eq.0) then 
          do ilist=1,nlist1s(ibox)
            jbox = list1(ilist,ibox)
            jstart = isrcse(1,jbox)
            jend = isrcse(2,jbox)
            nptssrc = jend-jstart+1
            if (nptssrc.gt.0) then 
              call zh3devaldirect(zk,nptstarg,ztargsort(1,istart),
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
        if (itargse(2,ibox).ge.itargse(1,ibox) .and.
     1 nlist2s(ibox).gt.0) then 
          do ilist=1,nlist2s(ibox) 
            jbox = list2(ilist,ibox)
            if (isrcse(2,jbox).ge.isrcse(1,jbox)) then
              call zh3dmploc(zcsrc(1,jbox),zctarg(1,ibox),nterm,nterm,
     1 rscale,rscale,width,width,zk,rmlexp(iaddr(1,jbox)),
     2 rmlexp(iaddr(2,ibox)))
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
C$OMP$PRIVATE(ibox,istart,iend,npts)
C$OMP$PRIVATE(ilist,jbox,jstart,jend)
C$OMP$PRIVATE(rscale, nterm)
C$OMP$SCHEDULE(DYNAMIC)
      do ibox=1,nboxes
        istart = itargse(1,ibox)
        iend = itargse(2,ibox)
        npts = iend-istart+1
        if (npts.gt.0 .and.
     1          nchild(ibox).eq.0 .and. nlist3s(ibox).gt.0) then 
          do ilist = 1,nlist3s(ibox)
            jbox = list3(ilist,ibox)
            if (isrcse(2,jbox).ge.isrcse(1,jbox)) then 
              rscale = rscales(ilevel(jbox))
              nterm = nterms(ilevel(jbox))
              call zh3devalmp(zk,npts,ztargsort(1,istart),rscale,
     1 nterm,zcsrc(1,jbox),rmlexp(iaddr(1,jbox)),potsort(istart),
     2 gradsort(1,istart),ifpgh)
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
C$OMP$PRIVATE(ibox,nterm,rscale,ilist,jbox,jstart,jend,npts)
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
              call zh3dformloc(zk,rscale,zsrcsort(1,jstart),npts,
     1 chargesort(jstart),dipstrsort(jstart),dipvecsort(1,jstart),
     2 zctarg(1,ibox),nterm,rmlexp(iaddr(2,ibox)),ifcharge,ifdipole)
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
     1    (nlist2s(ibox).gt.0).or.(nlist4s(ibox).gt.0)) then
            do j=1,nchild(ibox)
              jbox = ichild(j,ibox)
              if (itargse(2,jbox).ge.itargse(1,jbox)) then 
                call zh3dlocloc(zctarg(1,ibox),
     1 zctarg(1,jbox),nterm,nterm1,rscale,rscale1,width,width1,zk,
     2 rmlexp(iaddr(2,ibox)),rmlexp(iaddr(2,jbox)))
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
          call zh3devalloc(zk,npts,ztargsort(1,istart),rscale,
     1 nterm,zctarg(1,ibox),rmlexp(iaddr(2,ibox)),potsort(istart),
     2 gradsort(1,istart),ifpgh)
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
C
      subroutine zh3devaldirect(zk,nt,ztarg,ns,zsrc,charge,dipstr,
     1 dipvec,ifcharge,pot,grad,ifdipole,ifpgh)
      implicit real *8 (a-h,o-z)
      complex *16 zk,rz,h0,ima,z 
      complex *16 ztarg(3,nt),pot(nt),grad(3,nt)
      complex *16 zsrc(3,ns),charge(ns),dipstr(ns)
      complex *16 dipvec(3,ns)
      complex *16 x1my1,x2my2,x3my3 
      complex *16 dot1,dot2,dot3,dot 
      complex *16 tmp1,tmp2,tmp3,tmp4
      complex *16 rz3_pi4_inv

      ima = (0.0d0,1.0d0)
      done = 1 
      pi = atan(done)*4
      pi4  = pi*4 
        
      do i=1,nt 
        do j=1,ns 
          x1my1 =  ztarg(1,i)-zsrc(1,j)
          x2my2 =  ztarg(2,i)-zsrc(2,j)
          x3my3 =  ztarg(3,i)-zsrc(3,j)
          rz = sqrt(x1my1**2+x2my2**2+x3my3**2)
          z = ima*zk*rz
          if ((abs(rz).gt.1E-16)) then         
            if (ifcharge.eq.1) then 
              h0 = exp(z)/(pi4*rz)  
              pot(i) = pot(i)+charge(j)*h0 
              if (ifpgh.ge.2) then 
                rz3_pi4_inv = 1/(pi4*rz**3)
                h0 = exp(z)*(done-z)*rz3_pi4_inv
                grad(1,i) = grad(1,i)+charge(j)*(-x1my1)*h0
                grad(2,i) = grad(2,i)+charge(j)*(-x2my2)*h0
                grad(3,i) = grad(3,i)+charge(j)*(-x3my3)*h0
              endif 
            endif 

            if (ifdipole.eq.1) then 
              dot1 = x1my1*dipvec(1,j)
              dot2 = x2my2*dipvec(2,j)
              dot3 = x3my3*dipvec(3,j)

              dot = dot1+dot2+dot3
              tmp3 = exp(z)
              tmp4 = (done-z)
              rz3_pi4_inv = 1/(pi4*rz**3)
              pot(i) = pot(i)-dipstr(j)*tmp3*dot*tmp4*rz3_pi4_inv

              if (ifpgh.ge.2) then 
                tmp1 = dot*tmp3*rz3_pi4_inv;
                tmp1 = tmp1*(zk**2-3*tmp4/rz**2);
                tmp2 = tmp4*tmp3*rz3_pi4_inv;

                grad(1,i) = grad(1,i)
     1                -dipstr(j)*(dipvec(1,j)*tmp2+x1my1*tmp1)
                grad(2,i) = grad(2,i)
     1                -dipstr(j)*(dipvec(2,j)*tmp2+x2my2*tmp1)
                grad(3,i) = grad(3,i)
     1                -dipstr(j)*(dipvec(3,j)*tmp2+x3my3*tmp1)

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
      subroutine zh3dformmp(zk,rscale,zsrc,npts,charge,dipstr,dipvec,
     1 zcsrc,nterm,mp,ifcharge,ifdipole)
      implicit real *8 (a-h,o-z)
      complex *16 zk,ima,zrk,pi4_inv_i_zk,z1
      complex *16 zsrc(3,npts),charge(npts)
      complex *16 dipstr(npts),dipvec(3,npts)
      complex *16 zcsrc(3),mp(0:nterm,-nterm:nterm)
      complex *16 ctheta,stheta,cphi,sphi,eiphi,zr
      complex *16 e_zr(3),e_th(3),e_ph(3),dot_zr,dot_th,dot_ph
      complex *16 dot_th_zr,dot_ph_zr 
      complex *16, allocatable :: eiphis(:)
      complex *16, allocatable :: fjs(:),fjder(:)
      complex *16, allocatable :: yl(:,:), dyl(:,:)
      real *8, allocatable :: rat1(:,:),rat2(:,:)
      integer, allocatable :: m1(:)

      allocate(eiphis(-nterm:nterm))
      allocate(fjs(0:nterm+100),fjder(0:nterm+100))
      allocate(yl(0:nterm,0:nterm),dyl(0:nterm,0:nterm))
      allocate(rat1(0:nterm,0:nterm),rat2(0:nterm,0:nterm))
      allocate(m1(0:nterm))

      m1(0) = 1 
      do n=1,nterm
        m1(n) = m1(n-1)*(-1) 
      enddo 

      call ylgndrini(nterm,rat1,rat2)
      done = 1
      ima = (0.0d0,1.0d0)
      ifder = 0 
      if (ifdipole.eq.1) ifder = 1 
      pi = atan(done)*4
      pi4 = pi*4 
      pi4_inv = 1/pi4 
      pi4_inv_i_zk = pi4_inv*ima*zk 

      do i=1,npts 
        
        call zcart2sph(zsrc(1:3,i),zcsrc,zr,ctheta,stheta,cphi,sphi,
     1 eiphi)
        eiphis(0) = 1 
        do m=1,nterm 
          eiphis(m) = eiphis(m-1)*eiphi
        enddo 

        do m=-nterm,-1 
          eiphis(m) = 1/eiphis(-m)
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
        endif 

        zrk = zr*zk 
        call besseljs3d(nterm,zrk,rscale,fjs,ifder,fjder)

        if (ifdipole.eq.0) then  
          call zylgndr(nterm,ctheta,yl)
          do n=0,nterm
            do m=0,n 
              yl(n,m) = yl(n,m)*pi4_inv_i_zk*m1(m)
            enddo 
          enddo 
        else 
          call zylgndr2sf(nterm,ctheta,yl,dyl,rat1,rat2)
          do n=0,nterm
            do m=0,n 
              yl(n,m) = yl(n,m)*pi4_inv_i_zk*m1(m)
              dyl(n,m) = dyl(n,m)*pi4_inv_i_zk*m1(m)
            enddo 
          enddo 
        endif 

        if (ifcharge.eq.1) then 
          if (ifdipole.eq.0) then 
            do n=0,nterm 
              mp(n,0) = mp(n,0)+charge(i)*fjs(n)*yl(n,0)
              do m=1,n 
                z1 = charge(i)*fjs(n)*yl(n,m)
                mp(n,m) = mp(n,m)+z1*eiphis(-m)
                mp(n,-m) = mp(n,-m)+z1*eiphis(m)
              enddo 
            enddo 
          else 
            do n=0,nterm 
              mp(n,0) = mp(n,0)+charge(i)*fjs(n)*yl(n,0)
              do m=1,n
                z1 =  charge(i)*fjs(n)*yl(n,m)*stheta
                mp(n,m) = mp(n,m)+z1*eiphis(-m)
                mp(n,-m) = mp(n,-m)+z1*eiphis(m)
              enddo 
            enddo 
          endif 
        endif 

        if (ifdipole.eq.1) then 
          do n=0,nterm 
            mp(n,0) = mp(n,0)-
     1    dipstr(i)*zk*fjder(n)*yl(n,0)*dot_zr
            mp(n,0) = mp(n,0)+
     1    dipstr(i)*fjs(n)*dyl(n,0)*stheta*dot_th_zr 
            do m=1,n 
              z1 = dipstr(i)*zk*fjder(n)*yl(n,m)*stheta*dot_zr
              mp(n,m) = mp(n,m)-z1*eiphis(-m)
              mp(n,-m) = mp(n,-m)-z1*eiphis(m)
              z1 = dipstr(i)*fjs(n)*dyl(n,m)*dot_th_zr 
              mp(n,m) = mp(n,m)+z1*eiphis(-m)
              mp(n,-m) = mp(n,-m)+z1*eiphis(m)
              z1 = dipstr(i)*m*ima*fjs(n)*yl(n,m)*dot_ph_zr 
              mp(n,m) = mp(n,m)+z1*eiphis(-m)
              mp(n,-m) = mp(n,-m)-z1*eiphis(m)
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
      subroutine zh3dformloc(zk,rscale,zsrc,npts,charge,dipstr,dipvec,
     1 zctarg,nterm,loc,ifcharge,ifdipole)
      implicit real *8 (a-h,o-z)
      complex *16 zk,ima,zrk,z1,pi4_inv_i_zk
      complex *16 zsrc(3,npts),charge(npts)
      complex *16 dipstr(npts),dipvec(3,npts)
      complex *16 zctarg(3),loc(0:nterm,-nterm:nterm)
      complex *16 ctheta,stheta,cphi,sphi,eiphi,zr
      complex *16 e_zr(3),e_th(3),e_ph(3),dot_zr,dot_th,dot_ph
      complex *16 dot_th_zr,dot_ph_zr 
      complex *16, allocatable :: eiphis(:)
      complex *16, allocatable :: fhs(:),fhder(:)
      complex *16, allocatable :: yl(:,:),dyl(:,:)
      real *8, allocatable :: rat1(:,:),rat2(:,:)
      integer, allocatable :: m1(:)

      allocate(eiphis(-nterm:nterm))
      allocate(fhs(0:nterm+100),fhder(0:nterm+100))
      allocate(yl(0:nterm,0:nterm),dyl(0:nterm,0:nterm))
      allocate(rat1(0:nterm,0:nterm),rat2(0:nterm,0:nterm))
      allocate(m1(0:nterm))

      m1(0) = 1 
      do n=1,nterm
        m1(n) = m1(n-1)*(-1) 
      enddo 

      call ylgndrini(nterm,rat1,rat2)
      done = 1 
      ima = (0.0d0,1.0d0)
      ifder = 0 
      if (ifdipole.eq.1) ifder = 1 
      pi = atan(done)*4
      pi4 = pi*4 
      pi4_inv = 1/pi4 
      pi4_inv_i_zk = pi4_inv*ima*zk

      do i=1,npts 
        call zcart2sph(zsrc(1:3,i),zctarg,zr,ctheta,stheta,cphi,sphi,
     1 eiphi)

        eiphis(0) = 1 
        do m=1,nterm 
          eiphis(m) = eiphis(m-1)*eiphi
        enddo 

        do m=-nterm,-1 
          eiphis(m) = 1/eiphis(-m)
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
        endif 

        zrk = zk*zr 
        call h3dall(nterm,zrk,rscale,fhs,ifder,fhder)

        if (ifdipole.eq.0) then  
          call zylgndrf(nterm,ctheta,yl,rat1,rat2)
          do n=0,nterm
            do m=0,n 
              yl(n,m) = yl(n,m)*pi4_inv_i_zk*m1(m)
            enddo 
          enddo 
        else 
          call zylgndr2sf(nterm,ctheta,yl,dyl,rat1,rat2)
          do n=0,nterm
            do m=0,n 
              yl(n,m) = yl(n,m)*pi4_inv_i_zk*m1(m)
              dyl(n,m) = dyl(n,m)*pi4_inv_i_zk*m1(m)
            enddo 
          enddo 
        endif 

        if (ifcharge.eq.1) then 
          if (ifdipole.eq.0) then 
            do n=0,nterm 
              loc(n,0) = loc(n,0)+charge(i)*fhs(n)*yl(n,0)
              do m=1,n
                z1 =  charge(i)*fhs(n)*yl(n,m)
                loc(n,m) = loc(n,m)+z1*eiphis(-m)
                loc(n,-m) = loc(n,-m)+z1*eiphis(m)
              enddo 
c              do m=-n,n 
c                loc(n,m) = loc(n,m)+
c     1              charge(i)*fhs(n)*yl(n,abs(m))*eiphis(-m)
c              enddo 
            enddo 
          else 
            do n=0,nterm
              loc(n,0) = loc(n,0)+charge(i)*fhs(n)*yl(n,0)
              do m=1,n 
                z1 = charge(i)*fhs(n)*yl(n,m)*stheta
                loc(n,m) = loc(n,m)+z1*eiphis(-m)
                loc(n,-m) = loc(n,-m)+z1*eiphis(m)
              enddo  
c              do m=-n,n
c                if (m.eq.0) then  
c                  loc(n,0) = loc(n,0)+charge(i)*fhs(n)*yl(n,0)
c                else 
c                  loc(n,m) = loc(n,m)+
c     1       charge(i)*fhs(n)*yl(n,abs(m))*stheta*eiphis(-m)
c                endif 
c              enddo 
            enddo 
          endif 
        endif 


        if (ifdipole.eq.1) then 
          do n=0,nterm 
c            do m=-n,n 
c              if (m.eq.0) then 
c                loc(n,0) = loc(n,0)-
c     1    dipstr(i)*zk*fhder(n)*yl(n,0)*dot_zr
c                loc(n,0) = loc(n,0)+
c     1    dipstr(i)*fhs(n)*dyl(n,0)*stheta*dot_th/zr 
c              else 
c                loc(n,m) = loc(n,m)-
c     1    dipstr(i)*zk*fhder(n)*yl(n,abs(m))*stheta*eiphis(-m)*dot_zr
c                loc(n,m) = loc(n,m)+
c     1    dipstr(i)*fhs(n)*dyl(n,abs(m))*eiphis(-m)*dot_th/zr 
c                loc(n,m) = loc(n,m)+
c     1    dipstr(i)*m*ima*fhs(n)*yl(n,abs(m))*eiphis(-m)*dot_ph/zr 
c              endif 
c            enddo 
          loc(n,0) = loc(n,0)-
     1       dipstr(i)*zk*fhder(n)*yl(n,0)*dot_zr
          loc(n,0) = loc(n,0)+
     1       dipstr(i)*fhs(n)*dyl(n,0)*stheta*dot_th/zr 

            do m=1,n  
              z1 = dipstr(i)*zk*fhder(n)*yl(n,m)*stheta*dot_zr
              loc(n,m) = loc(n,m)-z1*eiphis(-m)
              loc(n,-m) = loc(n,-m)-z1*eiphis(m)
              z1 = dipstr(i)*fhs(n)*dyl(n,m)*dot_th/zr 
              loc(n,m) = loc(n,m)+z1*eiphis(-m)
              loc(n,-m) = loc(n,-m)+z1*eiphis(m)
              z1 = dipstr(i)*m*ima*fhs(n)*yl(n,m)*dot_ph/zr 
              loc(n,m) = loc(n,m)+z1*eiphis(-m)
              loc(n,-m) = loc(n,-m)-z1*eiphis(m)
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
      subroutine zh3devalmp(zk,nt,ztarg,rscale,nterm,zc,mp,pot,grad,
     1 ifpgh)
      implicit real *8 (a-h,o-z)
      complex *16 zk,ima,zrk,z1,z2 
      complex *16 ztarg(3,nt),grad(3,nt),pot(nt)
      complex *16 zc(3),mp(0:nterm,-nterm:nterm)
      complex *16 ctheta,stheta,cphi,sphi,eiphi,zr
      complex *16 e_zr(3),e_th(3),e_ph(3),e_th_zr(3),e_ph_zr(3)

      complex *16, allocatable :: eiphis(:)
      complex *16, allocatable :: fhs(:),fhder(:)
      complex *16, allocatable :: yl(:,:),dyl(:,:)
      real *8, allocatable :: rat1(:,:),rat2(:,:)
      integer, allocatable :: m1(:)

      allocate(eiphis(-nterm:nterm))
      allocate(fhs(0:nterm+100),fhder(0:nterm+100))
      allocate(yl(0:nterm,0:nterm),dyl(0:nterm,0:nterm))
      allocate(rat1(0:nterm,0:nterm),rat2(0:nterm,0:nterm))
      allocate(m1(0:nterm))


      call ylgndrini(nterm,rat1,rat2)
      ima = (0.0d0,1.0d0)
      done = 1 
      ifder = 0 
      if (ifpgh.ge.2) ifder = 1 
      m1(0) = 1 
      do i=1,nterm 
        m1(i)=m1(i-1)*(-1)
      enddo 

      do i=1,nt 
        call zcart2sph(ztarg(1:3,i),zc,zr,ctheta,stheta,cphi,sphi,
     1 eiphi)

        eiphis(0) = 1 
        do m=1,nterm 
          eiphis(m) = eiphis(m-1)*eiphi
        enddo 

        do m=-nterm,-1 
          eiphis(m) = 1/eiphis(-m)
        enddo 

        if (ifpgh.eq.2) then
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

        zrk = zr*zk 
        call h3dall(nterm,zrk,rscale,fhs,ifder,fhder)

        if (ifpgh.eq.1) then 
          call zylgndrf(nterm,ctheta,yl,rat1,rat2)
          do n=0,nterm
            do m=0,n 
              yl(n,m) = yl(n,m)*m1(m)
            enddo 
          enddo 
        else 
          call zylgndr2sf(nterm,ctheta,yl,dyl,rat1,rat2)
          do n=0,nterm
            do m=0,n 
              yl(n,m) = yl(n,m)*m1(m)
              dyl(n,m) = dyl(n,m)*m1(m)
            enddo 
          enddo 
        endif 

        do n=0,nterm 
          pot(i) = pot(i)+mp(n,0)*fhs(n)*yl(n,0)
          do m=1,n 
            if (ifpgh.eq.1) then
              z1 = fhs(n)*yl(n,m)
              pot(i) = pot(i)+
     1            z1*(mp(n,m)*eiphis(m)+mp(n,-m)*eiphis(-m))
            else 
              z1 = fhs(n)*yl(n,m)*stheta
              pot(i) = pot(i)+
     1            z1*(mp(n,m)*eiphis(m)+mp(n,-m)*eiphis(-m))
            endif 
          enddo 
c          do m=-n,n 
c            if (m.eq.0) then 
c              pot(i) = pot(i)+mp(n,0)*fhs(n)*yl(n,0)
c            else 
c              if (ifpgh.eq.1) then
c                pot(i) = pot(i)+mp(n,m)*fhs(n)*yl(n,abs(m))
c     1                           *eiphis(m)
c              else 
c              pot(i) = pot(i)+mp(n,m)*fhs(n)*yl(n,abs(m))
c     1                           *stheta*eiphis(m)
c              endif 
c            endif 
c          enddo
        enddo 

        

        if (ifpgh.ge.2) then 
          do n=0,nterm 
c            do m=-n,n 
c              if (m.eq.0) then 
c                grad(1:3,i) = grad(1:3,i)
c     1                   +mp(n,0)*zk*fhder(n)*yl(n,0)*e_zr
c                grad(1:3,i) = grad(1:3,i)
c     1              -mp(n,0)*fhs(n)*dyl(n,0)*stheta*e_th/zr
c              else 
c                grad(1:3,i) = grad(1:3,i)+
c     1    mp(n,m)*zk*fhder(n)*yl(n,abs(m))*stheta*eiphis(m)*e_zr
c                grad(1:3,i) = grad(1:3,i)-
c     1    mp(n,m)*fhs(n)*dyl(n,abs(m))*eiphis(m)*e_th/zr
c                grad(1:3,i) = grad(1:3,i)+
c     1    m*ima*mp(n,m)*fhs(n)*yl(n,abs(m))*eiphis(m)*e_ph/zr
c              endif 
c            enddo 
            grad(1:3,i) = grad(1:3,i)+mp(n,0)*
     1     (zk*fhder(n)*yl(n,0)*e_zr-fhs(n)*dyl(n,0)*stheta*e_th_zr)
            do m=1,n
              z1 = mp(n,m)*eiphis(m)
              z2 = mp(n,-m)*eiphis(-m)
              grad(1:3,i) = grad(1:3,i)+
     1     zk*fhder(n)*yl(n,m)*stheta*e_zr*(z1+z2)
              grad(1:3,i) = grad(1:3,i)-
     1     fhs(n)*dyl(n,m)*e_th_zr*(z1+z2)
              grad(1:3,i) = grad(1:3,i)+
     1     m*ima*fhs(n)*yl(n,m)*e_ph_zr*(z1-z2)
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
      subroutine zh3devalloc(zk,nt,ztarg,rscale,nterm,zc,loc,pot,grad,
     1 ifpgh)
      implicit real *8 (a-h,o-z)
      complex *16 zk,ima,zrk,z1,z2
      complex *16 ztarg(3,nt),pot(nt)
      complex *16 zc(3),loc(0:nterm,-nterm:nterm)
      complex *16 grad(3,nt)
      complex *16 ctheta,stheta,cphi,sphi,eiphi,zr
      complex *16 e_zr(3),e_th(3),e_ph(3),e_th_zr(3),e_ph_zr(3)

      complex *16, allocatable :: eiphis(:)
      complex *16, allocatable :: yl(:,:),dyl(:,:)
      complex *16, allocatable :: fjs(:),fjder(:)
      real *8, allocatable :: rat1(:,:),rat2(:,:)
      integer, allocatable :: m1(:)

      allocate(eiphis(-nterm:nterm))
      allocate(yl(0:nterm,0:nterm),dyl(0:nterm,0:nterm))
      allocate(fjs(0:nterm+100),fjder(0:nterm+100))
      allocate(rat1(0:nterm,0:nterm),rat2(0:nterm,0:nterm))
      allocate(m1(0:nterm))

      call ylgndrini(nterm,rat1,rat2)
      ima = (0.0d0,1.0d0)
      done = 1 
      ifder = 0 
      if (ifpgh.ge.2) ifder = 1 

      m1(0) = 1 
      do i=1,nterm 
        m1(i) = m1(i-1)*(-1)
      enddo 
      

      do i=1,nt 
        call zcart2sph(ztarg(1:3,i),zc,zr,ctheta,stheta,cphi,sphi,
     1 eiphi)
        eiphis(0) = 1 
        do m=1,nterm 
          eiphis(m) = eiphis(m-1)*eiphi
        enddo 

        do m=-nterm,-1 
          eiphis(m) = 1/eiphis(-m)
        enddo 

        if (ifpgh.eq.2) then
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

        zrk = zr*zk 
        call besseljs3d(nterm,zrk,rscale,fjs,ifder,fjder)

        if (ifpgh.eq.1) then 
          call zylgndrf(nterm,ctheta,yl,rat1,rat2)
          do n=0,nterm
            do m=0,n 
              yl(n,m) = yl(n,m)*m1(m)
            enddo 
          enddo 
        else 
          call zylgndr2sf(nterm,ctheta,yl,dyl,rat1,rat2)
          do n=0,nterm
            do m=0,n 
              yl(n,m) = yl(n,m)*m1(m)
              dyl(n,m) = dyl(n,m)*m1(m)
            enddo 
          enddo 
        endif 

        do n=0,nterm 
          pot(i) = pot(i)+loc(n,0)*fjs(n)*yl(n,0)
          do m=1,n 
            if (ifpgh.eq.1) then
              z1 = fjs(n)*yl(n,m)
              pot(i) = pot(i)+
     1            z1*(loc(n,m)*eiphis(m)+loc(n,-m)*eiphis(-m))
            else 
              z1 = fjs(n)*yl(n,m)*stheta
              pot(i) = pot(i)+
     1            z1*(loc(n,m)*eiphis(m)+loc(n,-m)*eiphis(-m))
            endif  
          enddo 
c          do m=-n,n 
c            if (m.eq.0) then 
c              pot(i) = pot(i)+loc(n,0)*fjs(n)*yl(n,0)
c            else 
c              if (ifpgh.eq.1) then 
c                pot(i) = pot(i)+loc(n,m)*fjs(n)*yl(n,abs(m))
c     1                           *eiphis(m)
c              else 
c                pot(i) = pot(i)+loc(n,m)*fjs(n)*yl(n,abs(m))
c     1                           *stheta*eiphis(m)
c              endif 
c            endif 
c          enddo
        enddo 

        if (ifpgh.ge.2) then 
          do n=0,nterm 
c            do m=-n,n 
c              if (m.eq.0) then 
c                grad(1:3,i) = grad(1:3,i)
c     1                   +loc(n,0)*zk*fjder(n)*yl(n,0)*e_zr
c                grad(1:3,i) = grad(1:3,i)
c     1              -loc(n,0)*fjs(n)*dyl(n,0)*stheta*e_th_zr
c              else 
c                grad(1:3,i) = grad(1:3,i)+
c     1 loc(n,m)*zk*fjder(n)*yl(n,abs(m))*stheta*eiphis(m)*e_zr
c                grad(1:3,i) = grad(1:3,i)-
c     1 loc(n,m)*fjs(n)*dyl(n,abs(m))*eiphis(m)*e_th_zr
c                grad(1:3,i) = grad(1:3,i)+
c     1    m*ima*loc(n,m)*fjs(n)*yl(n,abs(m))*eiphis(m)*e_ph_zr
c              endif 
c            enddo 

            grad(1:3,i) = grad(1:3,i)+loc(n,0)*
     1     (zk*fjder(n)*yl(n,0)*e_zr-fjs(n)*dyl(n,0)*stheta*e_th_zr)
            do m=1,n
              z1 = loc(n,m)*eiphis(m)
              z2 = loc(n,-m)*eiphis(-m)
              grad(1:3,i) = grad(1:3,i)+
     1     zk*fjder(n)*yl(n,m)*stheta*e_zr*(z1+z2)
              grad(1:3,i) = grad(1:3,i)-
     1     fjs(n)*dyl(n,m)*e_th_zr*(z1+z2)
              grad(1:3,i) = grad(1:3,i)+
     1     m*ima*fjs(n)*yl(n,m)*e_ph_zr*(z1-z2)
            enddo
          enddo 
        endif 
      enddo 

      end 
