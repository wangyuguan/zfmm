C
C   subroutines for 2d complex laplace fast multipole method  
C
C       For i=1,...,nt,  evaluate:
C
C                ns
C       u(xi) = \sum charge(j)*log(R(xi,yj))/2/pi
C                j=1
C               - dipstr(j)*<dipvec(:,j), \nabla_y log(R(xi,yj))>/2/pi
C                 
C               where R(x,y)=sqrt((x1-y1)^2+(x2-y2)^2)
C
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine zlfmm2d(eps,ns,zsrc,ifcharge,charge,ifdipole,dipstr,
     1 dipvec,pot,grad,nt,ztarg,ifpgh,isep,ifprint)
c----------------------------------------------
c   INPUT PARAMETERS:
c   eps          : FMM precision requested (real *8)
c   ns           : number of sources (integer)
c   zsrc(2,ns)   : source locations (complex *16)
c   ifcharge     : flag for include single layer (1 or 0)
c   charge(ns)   : charge strengths (complex *16)
c   ifdipole     : flag for double layer (1 or 0)
c   dipstr(ns)   : dipole strengths (complex *16)
c   dipvec(2,ns) : dipole vectors (complex *16)
c   nt           : number of targets (integer)
c   ztarg(2,nt)  : complex target locations (complex *16)
c   ifpgh        : ifpghtarg = 1, only potential is computed at targets
c                 ifpghtarg = 2, potential and gradient are computed at targets
c   ifprint      : if print properties (1 or 0)

c   OUTPUT PARAMETERS
c   pot(nt)      : potential at the target locations (complex *16)
c   grad(2,nt)   : gradient of potential at the target locations (complex *16)
c----------------------------------------------
      implicit real *8 (a-h,o-z)
      complex *16 ima 
      complex *16 zsrc(2,ns),ztarg(2,nt)
      complex *16 charge(ns),pot(nt),grad(2,nt)
      complex *16 dipstr(ns),dipvec(2,ns)
      complex *16 dipstrsort(ns),dipvecsort(2,ns)
      dimension src(2,ns),targ(2,nt)
      dimension srcsort(2,ns),targsort(2,nt)
      dimension zc(2)
      complex *16 zsrcsort(2,ns),chargesort(ns),ztargsort(2,nt)
      dimension iptr(8)
      integer ltree
      integer, allocatable :: itree(:)
      real *8, allocatable :: centers(:,:),boxsize(:),rscales(:)
      complex *16, allocatable :: zcsrc(:,:),zctarg(:,:)
      integer, allocatable :: isrc(:),isrcse(:,:),itarg(:),itargse(:,:)
      integer, allocatable :: nterms(:)

      done = 1 
      pi = atan(done)*4 
      ima = (0.0d0,1.0d0)

      if ((ifcharge.lt.0).or.(ifcharge.gt.1)) then 
        ifcharge = 1  
      endif 

      if ((ifdipole.lt.0).or.(ifdipole.gt.1)) then 
        ifdipole = 0 
      endif

      if ((ifpgh.lt.0).or.(ifpgh.gt.2)) then 
        ifpgh = 1
      endif 

c   get real locations 
      do i = 1,ns
        src(1:2,i) = real(zsrc(1:2,i))
      enddo   

      do i = 1,nt
        targ(1:2,i) = real(ztarg(1:2,i))
      enddo   

c
cc      call the tree memory management
c       code to determine number of boxes,
c       number of levels and length of tree
c

      ifpghtarg = 0 
      
        
      call lndiv2d(eps,ns,nt,ifcharge,ifdipole,ifpgh,ifpghtarg,ndiv,
     1 idivflag)

c
         
      if (ifprint.eq.1) then 
        call prinf('ndiv=*',ndiv,1)
        call prinf('idivflag=*',idivflag,1)
      endif 


      ltree = 0
      nlmin = 0
      nlmax = 51
      ifunif = 0
      iper = 0

      call pts_tree_mem(src,ns,targ,nt,idivflag,ndiv,nlmin,nlmax,
     1 ifunif,iper,isep,nlevels,nboxes,ltree)
      

      if (ifprint.eq.1) then 
        call prinf('nboxes=*',nboxes,1)
        call prinf('nlevels=*',nlevels,1)
      endif 

      allocate(itree(ltree))
      allocate(centers(2,nboxes))
      allocate(boxsize(0:nlevels))

c
c       call the tree code
c

      call pts_tree_build(src,ns,targ,nt,idivflag,ndiv,nlmin,nlmax,
     1 ifunif,iper,isep,nlevels,nboxes,ltree,itree,iptr,centers,
     2 boxsize)

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
     1 centers,itarg,itargse)
      call pts_tree_sort(ns,src,itree,ltree,nboxes,nlevels,iptr,
     1 centers,isrc,isrcse)


c      do ibox=1,nboxes 
c        print *, isrcse(1,ibox), isrcse(2,ibox)
c      enddo 


c
cc     reorder sources charges, and targets
c

      call dreorderf(2,ns,src,srcsort,isrc)
      call zreorderf(2,ns,zsrc,zsrcsort,isrc)
      if (ifcharge.eq.1) then 
        call zreorderf(1,ns,charge,chargesort,isrc)
      endif 
      call dreorderf(2,nt,targ,targsort,itarg)
      call zreorderf(2,nt,ztarg,ztargsort,itarg)
      if (ifdipole.eq.1) then 
        call zreorderf(1,ns,dipstr,dipstrsort,isrc)
        call zreorderf(2,ns,dipvec,dipvecsort,isrc)
      endif 

c
c      compute the centers of the complex boxes 
c

      allocate(zcsrc(2,nboxes),zctarg(2,nboxes))
      zcsrc = 0 
      zctarg = 0
      do ibox=1,nboxes 
c
c          compute complex scource boxcenters 
c
        if (isrcse(1,ibox).le.isrcse(2,ibox)) then 
          npts = isrcse(2,ibox)-isrcse(1,ibox)+1 
          zc = 0 
          do i=isrcse(1,ibox),isrcse(2,ibox)
            zc = zc+aimag(zsrcsort(1:2,i))
          enddo 
          zc = zc/npts 
          zcsrc(1:2,ibox) = centers(1:2,ibox)+ima*zc 
        endif             
c
c          compute complex target boxcenters 
c

        if (itargse(1,ibox).le.itargse(2,ibox)) then 
          npts = itargse(2,ibox)-itargse(1,ibox)+1 
          zc = 0 
          do i=itargse(1,ibox),itargse(2,ibox)
              zc = zc+aimag(ztargsort(1:2,i))
          enddo 
          zc = zc/npts
          zctarg(1:2,ibox) = centers(1:2,ibox)+ima*zc
C                 print *, zctarg(1:2,ibox), centers(1:2,ibox)
        endif 
      enddo 


c
cc      compute scaling factor for multipole/local expansions
c       and lengths of multipole and local expansions
c

      allocate(rscales(0:nlevels),nterms(0:nlevels))

      nmax = 0
      ier = 0
      do i=0,nlevels
        rscales(i) = boxsize(i)
        call l2dterms(eps,nterms(i),ier)
        if(nterms(i).gt.nmax) nmax = nterms(i)
      enddo 
      if (ifprint.eq.1) then 
        call prinf('nmax=*',nmax,1)
        call prin2('rscales=*',rscales,nlevels+1)
        call prinf('nterms=*',nterms,nlevels+1)
      endif 



      call zlfmm2dmain(ns,zsrcsort,chargesort,dipstrsort,dipvecsort,
     1 nt,ztargsort,isrc,itarg,isrcse,itargse,itree,ltree,iptr,nlevels,
     2 nboxes,boxsize,rscales,nterms,nmax,centers,zcsrc,zctarg,    
     3 itree(iptr(1)),itree(iptr(2)),itree(iptr(3)),itree(iptr(4)),
     4 itree(iptr(5)),itree(iptr(6)),isep,pot,grad,ifcharge,
     5 ifdipole,ifpgh,ifprint)
     
      pot = pot/2/pi 
      grad = grad/2/pi 
      end 
c
c
c
c
c
      subroutine zlfmm2dmain(ns,zsrcsort,chargesort,dipstrsort,
     1 dipvecsort,nt,ztargsort,isrc,itarg,isrcse,itargse,itree,ltree,
     2 iptr,nlevels,nboxes,boxsize,rscales,nterms,nmax,centers,zcsrc,
     3 zctarg,laddr,ilevel,iparent,nchild,ichild,ncol,isep,pot,grad,
     4 ifcharge,ifdipole,ifpgh,ifprint)
      implicit real *8 (a-h,o-z)
      complex *16 zsrcsort(2,ns), chargesort(ns)
      complex *16 dipstrsort(ns),dipvecsort(2,ns)
      dimension isrc(ns),itarg(nt),isrcse(2,nboxes),itargse(2,nboxes)
      complex *16 ztargsort(2,nt),pot(nt),grad(2,nt)
      integer ltree 
      dimension ilevel(nboxes),nchild(nboxes),ncol(nboxes)
      dimension iparent(nboxes),ichild(4,nboxes)
      dimension itree(ltree), iptr(8) 
      dimension boxsize(0:nlevels),rscales(0:nlevels)
      dimension nterms(0:nlevels),centers(2,nboxes)
      complex *16 zcsrc(2,nboxes),zctarg(2,nboxes)
      integer laddr(2,0:nlevels)
      integer, allocatable :: list1(:,:),list2(:,:),list3(:,:)
      integer, allocatable :: list4(:,:)
      integer, allocatable :: nlist1s(:),nlist2s(:),nlist3s(:)
      integer, allocatable :: nlist4s(:)
      complex *16 potsort(nt),gradsort(2,nt)
      complex *16 mp(-nmax:nmax,nboxes),loc(-nmax:nmax,nboxes)
      dimension carray(0:2*nmax,0:2*nmax)

c
c        initialize all variables 
c
      pot = 0 
      potsort = 0 
      grad = 0 
      gradsort = 0 
      mp = 0 
      loc = 0 

      ldc = 2*nmax 
      call init_carray(carray,ldc)

      if (ifprint.eq.1) then 
        call prinf('ifprint=*',ifprint,1)
        call prinf('nmax=*',nmax,1)
      endif 

c
c        compute list info
c
      iper = 0 
      mnbors = (2*isep+1)**2
      call computemnlists(nlevels,nboxes,itree(iptr(1)),boxsize,
     1  centers,itree(iptr(3)),itree(iptr(4)),itree(iptr(5)),isep,
     1  itree(iptr(6)),mnbors,itree(iptr(7)),iper,mnlist1,mnlist2,
     1  mnlist3,mnlist4)

      if (ifprint.eq.1) then 
        call prinf('mnlist1=*',mnlist1,1)
        call prinf('mnlist2=*',mnlist2,1)
        call prinf('mnlist3=*',mnlist3,1)
        call prinf('mnlist4=*',mnlist4,1)
      endif 

      allocate(nlist1s(nboxes),list1(mnlist1,nboxes))
      allocate(nlist2s(nboxes),list2(mnlist2,nboxes))
      allocate(nlist3s(nboxes),list3(mnlist3,nboxes))
      allocate(nlist4s(nboxes),list4(mnlist4,nboxes))
        
      call computelists(nlevels,nboxes,itree(iptr(1)),boxsize,centers,
     1 itree(iptr(3)),itree(iptr(4)),itree(iptr(5)),isep,
     1 itree(iptr(6)),mnbors,itree(iptr(7)),iper,nlist1s,mnlist1,list1,
     1 nlist2s,mnlist2,list2,nlist3s,mnlist3,list3,nlist4s,mnlist4,
     1 list4)

C
CC      list1 of ibox: only for childless box, include itself, childless colleagues 
CC             and colleagues' children adjacent to ibox 
CC      list2 of ibox:: parent's colleagues' children that are well separated 
CC      list3 of ibox:: only for childless box, colleagues' children non-adjacent to ibox
CC
CC      list4 of ibox:: inverse of list 3 


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
          if (nchild(ibox).eq.0.and.npts.gt.0) then
            call zl2dformmp(rscale,zsrcsort(1,istart),npts,
     1 chargesort(istart),dipstrsort(istart),dipvecsort(1,istart),
     2 zcsrc(1,ibox),nterm,mp(-nterm,ibox),ifcharge,ifdipole)
          endif 
        enddo 
C$OMP END PARALLEL DO
      enddo 
      call cpu_time(t2)
      if (ifprint.eq.1) call prin2('from mp spends * seconds',t2-t1,1)


c
cc       ... step 2, for every parant box, form multipole expansions
cc       from their children 
c

      call cpu_time(t1)
C$        t1=omp_get_wtime()
      do ilev=nlevels-1,1,-1 
        nterm = nterms(ilev)
        rscale = rscales(ilev)
        nterm1 = nterms(ilev+1)
        rscale1 = rscales(ilev+1)
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,j,jbox)
C$OMP$SCHEDULE(DYNAMIC)
        do ibox=laddr(1,ilev),laddr(2,ilev)
          if (nchild(ibox).gt.0 .and. 
     1  isrcse(2,ibox).ge.isrcse(1,ibox)) then
            do j = 1,nchild(ibox)
              jbox = ichild(j,ibox)
              if (isrcse(2,jbox).ge.isrcse(1,jbox)) then 
                call zl2dmpmp(zcsrc(1,ibox),zcsrc(1,jbox),rscale,
     1 rscale1,nterm,nterm1,mp(-nterm,ibox),mp(-nterm1,jbox),carray,
     1 ldc)   
              endif 
            enddo 
          endif  
        enddo 
C$OMP END PARALLEL DO
      enddo 
      call cpu_time(t2)
C$        t2=omp_get_wtime()
      if (ifprint.eq.1) call prin2('m2m spends * seconds,',t2-t1,1)

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
        if (nptstarg.gt.0 .and. nlist1s(ibox).gt.0. .and.
     1  nchild(ibox).eq.0) then 

          do ilist=1,nlist1s(ibox)
            jbox = list1(ilist,ibox)
            jstart = isrcse(1,jbox)
            jend = isrcse(2,jbox)
            nptssrc = jend-jstart+1

            if (nptssrc.gt.0) then 
              call zl2devaldirect(nptstarg,ztargsort(1,istart),
     1              nptssrc,zsrcsort(1,jstart),chargesort(jstart),
     2            dipstrsort(jstart),dipvecsort(1,jstart),ifcharge,
     3            potsort(istart),gradsort(1,istart),ifdipole,ifpgh)
            endif 
          enddo 
        endif 
      enddo 
C$OMP END PARALLEL DO
      call cpu_time(t2)
C$        t2=omp_get_wtime()
      if (ifprint.eq.1) call prin2('direct eval spends * seconds',
     1 t2-t1, 1)


C
CC       ... step 4, for each box, convert the multipole expansions of all 
CC             boxes in list2 into the its local expansion  
C

      call cpu_time(t1)
C$        t1=omp_get_wtime()
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,rscale,nterm)
C$OMP$PRIVATE(ilist,jbox)
C$OMP$SCHEDULE(DYNAMIC)
      do ibox=1,nboxes
        rscale = rscales(ilevel(ibox))
        nterm = nterms(ilevel(ibox))
        if (itargse(2,ibox).ge.itargse(1,ibox)
     1 .and. nlist2s(ibox).gt.0) then 
          do ilist=1,nlist2s(ibox) 
            jbox = list2(ilist,ibox)
            if (isrcse(2,jbox).ge.isrcse(1,jbox)) then        
              call zl2dmploc(zctarg(1,ibox),zcsrc(1,jbox),rscale,
     1 rscale,nterm,nterm,loc(-nterm,ibox),mp(-nterm,jbox),carray,ldc)                     
            endif 
          enddo 
        endif 
      enddo 
C$OMP END PARALLEL DO
      call cpu_time(t2)
C$        t2=omp_get_wtime()
      if (ifprint.eq.1) call prin2('m2l spends * seconds', t2-t1, 1)


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
        if ((npts.gt.0).and.(nchild(ibox).eq.0).and.
     1 (nlist3s(ibox).gt.0)) then 
          do ilist = 1,nlist3s(ibox)
            jbox = list3(ilist,ibox)
            if (isrcse(2,jbox).ge.isrcse(1,jbox)) then 
              rscale = rscales(ilevel(jbox))
              nterm = nterms(ilevel(jbox))
              call zl2devalmp(npts,ztargsort(1,istart),rscale,
     1 nterm,zcsrc(1,jbox),mp(-nterm,jbox),potsort(istart),
     2 gradsort(1,istart),ifpgh)
            endif 
          enddo 
        endif 
      enddo 
C$OMP END PARALLEL DO
      call cpu_time(t2)
C$        t2=omp_get_wtime()
      if (ifprint.eq.1) call prin2('eval mp spends * seconds',t2-t1,1)

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
        if ((itargse(2,ibox).ge.itargse(1,ibox))
     1 .and.(nlist4s(ibox).gt.0)) then
          do ilist=1,nlist4s(ibox)
            jbox = list4(ilist,ibox)
            jstart = isrcse(1,jbox)
            jend = isrcse(2,jbox)
            npts = jend-jstart+1
            if (npts.gt.0) then 
              call zl2dformloc(rscale,zsrcsort(1,jstart),npts,
     1 chargesort(jstart),dipstrsort(jstart),dipvecsort(1,jstart),
     2 zctarg(1,ibox),nterm,loc(-nterm,ibox),ifcharge,ifdipole)
            endif 
          enddo 
        endif 
      enddo 
C$OMP END PARALLEL DO
      call cpu_time(t2)
C$        t2=omp_get_wtime() 
      if (ifprint.eq.1) call prin2('form loc spends * seconds',t2-t1,1)

C
CC       ... step 7, for each parent box, shift its local expansion 
CC             to its children 
C

      call cpu_time(t1)
C$        t1=omp_get_wtime()
      do ilev=1,nlevels-1
        nterm = nterms(ilev)
        rscale = rscales(ilev)
        nterm1 = nterms(ilev+1)
        rscale1 = rscales(ilev+1)
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox)
C$OMP$PRIVATE(j,jbox)
C$OMP$SCHEDULE(DYNAMIC)
        do ibox=laddr(1,ilev),laddr(2,ilev)
          if (itargse(2,ibox).ge.itargse(1,ibox) .and.
     1   (nlist2s(ibox).gt.0 .or. nlist4s(ibox).gt.0)) then
            do j=1,nchild(ibox)
              jbox = ichild(j,ibox)
              if (itargse(2,jbox).ge.itargse(1,jbox)) then 
                call zl2dlocloc(zctarg(1,jbox),zctarg(1,ibox),rscale1,
     1 rscale,nterm1,nterm,loc(-nterm1,jbox),loc(-nterm,ibox),carray,
     1 ldc)
              endif 
            enddo 
          endif 
        enddo 
C$OMP END PARALLEL DO
      enddo 
      call cpu_time(t2)
C$        t2=omp_get_wtime()
      if (ifprint.eq.1) call prin2('l2l * seconds',t2-t1,1)

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
          call zl2devalloc(npts,ztargsort(1,istart),rscale,nterm,
     1 zctarg(1,ibox),loc(-nterm,ibox),potsort(istart),
     1 gradsort(1,istart),ifpgh)
        endif 
      enddo 
C$OMP END PARALLEL DO
      call cpu_time(t2)
C$        t2=omp_get_wtime()
      if (ifprint.eq.1) call prin2('eval loc * seconds',t2-t1,1)


C
CC      resort the potential into the original order 
C
      call zreorderi(1,nt,potsort,pot,itarg) 
      if (ifpgh.gt.1) call zreorderi(2,nt,gradsort,grad,itarg) 

      end 
c
c
c
c
c
      subroutine zl2devaldirect(nt,ztargsort,ns,zsrcsort,chargesort,
     1 dipstrsort,dipvecsort,ifcharge,potsort,gradsort,ifdipole,ifpgh)
      implicit real *8 (a-h,o-z)
      complex *16 rz,rz_square
      complex *16 tmp,tmp1,tmp2
      complex *16 dipstrsort(ns),dipvecsort(2,ns)
      complex *16 ztargsort(2,nt),potsort(nt),gradsort(2,nt)
      complex *16 zsrcsort(2,ns),chargesort(ns)
      complex *16 y1mx1,y2mx2,x1my1,x2my2 

      do i=1,nt 
        do j=1,ns 
          y1mx1 = zsrcsort(1,j)-ztargsort(1,i)
          y2mx2 = zsrcsort(2,j)-ztargsort(2,i)
          x1my1 = -y1mx1 
          x2my2 = -y2mx2 
          rz_square = y1mx1**2+y2mx2**2
          rz = sqrt(rz_square)


          if (abs(rz).gt.1E-16) then 

            if (ifcharge.eq.1) then 
              potsort(i) = potsort(i)+chargesort(j)*log(rz)
            endif 

            if (ifdipole.eq.1) then 
              potsort(i) = potsort(i)-dipstrsort(j)*
     1   (dipvecsort(1,j)*y1mx1+dipvecsort(2,j)*y2mx2)/rz_square
            endif 

            if(ifpgh.gt.1)then 

              if (ifcharge.eq.1) then 
                gradsort(1,i) = gradsort(1,i)+chargesort(j)*
     1       x1my1/rz_square
                gradsort(2,i) = gradsort(2,i)+chargesort(j)*
     1       x2my2/rz_square
              endif 

              if (ifdipole.eq.1) then 
                gradsort(1,i) = gradsort(1,i)-dipstrsort(j)*
     1              (dipvecsort(1,j)*(x1my1**2-x2my2**2)+
     1            dipvecsort(2,j)*2*x1my1*x2my2)/rz_square**2 
                gradsort(2,i) = gradsort(2,i)-dipstrsort(j)*
     1              (dipvecsort(1,j)*2*x1my1*x2my2+
     1      dipvecsort(2,j)*(x2my2**2-x1my1**2))/rz_square**2 
              endif 

            endif 

          endif 
        enddo
      enddo  
      end
c
c
c
c
c
      subroutine zl2dformmp(rscale,zsrcsort,npts,chargesort,
     1 dipstrsort,dipvecsort,zcsrc,nterm,mp,ifcharge,ifdipole)
      implicit real *8 (a-h,o-z)
      complex *16 zsrcsort(2,npts),chargesort(npts)
      complex *16 dipstrsort(npts),dipvecsort(2,npts)
      complex *16 zcsrc(2),mp(-nterm:nterm)
      complex *16 ez,rz,z,ima,dr,dth 
      complex *16 vr(2),ve(2),pr,pe
      complex *16 zall(0:nterm),ezall(0:nterm)
      

      ima = (0.0d0,1.0d0)
      done = 1 
      zall(0) = 1 
      ezall(0) = 1 

      do i=1,npts 

        call zcart2polar(zsrcsort(1,i),zcsrc,rz,ez)
        z = rz/rscale 

        do m=1,nterm 
          zall(m) = zall(m-1)*z
          ezall(m) = ezall(m-1)*ez 
        enddo 


        if (ifdipole.eq.1) then 
          vr = (zsrcsort(1:2,i)-zcsrc)/rz
          ve(1) = -vr(2)
          ve(2) = vr(1)

          pr = vr(1)*dipvecsort(1,i)+vr(2)*dipvecsort(2,i)
          pe = ve(1)*dipvecsort(1,i)+ve(2)*dipvecsort(2,i)
        endif 

        if (ifcharge.eq.1) then 
          mp(0) = mp(0)+chargesort(i)
        endif 

        do m=1,nterm 

          if (ifcharge.eq.1) then 
            mp(m) = mp(m)-zall(m)/ezall(m)*chargesort(i)/m/2
            mp(-m) = mp(-m)-zall(m)*ezall(m)*chargesort(i)/m/2
          endif 

          if (ifdipole.eq.1) then 
            dr = -zall(m-1)/ezall(m)/2/rscale 
            dth = ima*zall(m)/ezall(m)/2

            mp(m) = mp(m)-dipstrsort(i)*(dr*pr+dth*pe/rz)

            dr = -zall(m-1)*ezall(m)/2/rscale
            dth = -ima*zall(m)*ezall(m)/2

            mp(-m) = mp(-m)-dipstrsort(i)*(dr*pr+dth*pe/rz)
          endif 

        enddo 
        

      enddo 

      end 
c
c
c
c
c
      subroutine zl2dformloc(rscale,zsrcsort,npts,chargesort,
     1   dipstrsort,dipvecsort,zcsrc,nterm,loc,ifcharge,ifdipole)
      implicit real *8 (a-h,o-z)
      complex *16 zsrcsort(2,npts),chargesort(npts)
      complex *16 dipstrsort(npts),dipvecsort(2,npts)
      complex *16 zcsrc(2),loc(-nterm:nterm)
      complex *16 ez,rz,z,ima,dr,dth 
      complex *16 vr(2),ve(2),pr,pe
      complex *16 zall(0:nterm),ezall(0:nterm)


      ima = (0.0d0,1.0d0)
      done = 1 
      zall(0) = 1 
      ezall(0) = 1 

      do i=1,npts 

        call zcart2polar(zsrcsort(1,i),zcsrc,rz,ez)
        z = rz/rscale 

        do m=1,nterm 
          zall(m) = zall(m-1)*z
          ezall(m) = ezall(m-1)*ez 
        enddo 

        if (ifdipole.eq.1) then 
          vr = (zsrcsort(1:2,i)-zcsrc)/rz
          ve(1) = -vr(2)
          ve(2) = vr(1)

          pr = vr(1)*dipvecsort(1,i)+vr(2)*dipvecsort(2,i)
          pe = ve(1)*dipvecsort(1,i)+ve(2)*dipvecsort(2,i)

          loc(0) = loc(0)-dipstrsort(i)*pr/rz

        endif 

        if (ifcharge.eq.1) then 
          loc(0) = loc(0)+chargesort(i)*log(rz)
        endif 

        
        do m=1,nterm 

          if (ifcharge.eq.1) then 
            loc(m) = loc(m)-done/zall(m)/ezall(m)*chargesort(i)/m/2
            loc(-m) = loc(-m)-done/zall(m)*ezall(m)*chargesort(i)/m/2
          endif 

          if (ifdipole.eq.1) then 
            dr = done/zall(m)/rz/ezall(m)/2
            dth = ima/zall(m)/ezall(m)/2

            loc(m) = loc(m)-dipstrsort(i)*(dr*pr+dth*pe/rz)


            dr = done/zall(m)/rz*ezall(m)/2
            dth = -ima/zall(m)*ezall(m)/2

            loc(-m) = loc(-m)-dipstrsort(i)*(dr*pr+dth*pe/rz)
          endif 

        enddo 
        

      enddo 

      end 
c
c
c
c
c
      subroutine zl2devalmp(nt,ztargsort,rscale,nterm,zc,mp,potsort,
     1 gradsort,ifpgh)
      implicit real *8 (a-h,o-z)
      complex *16 zk,zr,ez,z,ima,dr,dth 
      complex *16 ztargsort(2,nt),potsort(nt),gradsort(2,nt)
      complex *16 zc(2),mp(-nterm:nterm)
      complex *16 vr(2),ve(2),vtmp(2)
      complex *16 zall(0:nterm),ezall(0:nterm)

      ima = (0.0d0,1.0d0)
      zall(0) = 1 
      ezall(0) = 1 
      
      do i=1,nt 
        call zcart2polar(ztargsort(1,i),zc,zr,ez)
        z = zr/rscale 


        do m=1,nterm 
          zall(m) = zall(m-1)*z
          ezall(m) = ezall(m-1)*ez 
        enddo 


        potsort(i) = potsort(i) + mp(0)*log(zr)

        if (ifpgh.gt.1) then 
          vr = (ztargsort(1:2,i)-zc)/zr 
          ve(1) = -vr(2) 
          ve(2) = vr(1)
          dr = 1/zr
          gradsort(1:2,i) = gradsort(1:2,i)+mp(0)*dr*vr
        endif 

        do m=1,nterm 
          potsort(i) = potsort(i)+mp(m)/zall(m)*ezall(m)
          potsort(i) = potsort(i)+mp(-m)/zall(m)/ezall(m)

          if(ifpgh.gt.1)then 
            dr = -m/zall(m)*ezall(m)/zr 
            dth = ima*m/zall(m)*ezall(m)
            vtmp = dr*vr+dth*ve/zr
            gradsort(1:2,i) = gradsort(1:2,i)+mp(m)*vtmp

            dr = -m/zall(m)/ezall(m)/zr 
            dth = -ima*m/zall(m)/ezall(m)
            vtmp = dr*vr+dth*ve/zr
            gradsort(1:2,i) = gradsort(1:2,i)+mp(-m)*vtmp
          endif 
        enddo 

      enddo 

      end 
c
c
c
c
c
      subroutine zl2devalloc(nt,ztargsort,rscale,nterm,zc,loc,potsort,
     1     gradsort,ifpgh)
      implicit real *8 (a-h,o-z)
      complex *16 zk,zr,ez,z,ima,dr,dth 
      complex *16 zall(0:nterm),ezall(0:nterm)
      complex *16 ztargsort(2,nt),potsort(nt),gradsort(2,nt)
      complex *16 zc(2),loc(-nterm:nterm)
      complex *16 vr(2),ve(2),vtmp(2)

      ima = (0.0d0,1.0d0)
      zall(0) = 1 
      ezall(0) = 1 

      do i=1,nt 
        call zcart2polar(ztargsort(1,i),zc,zr,ez)
        z = zr/rscale 

        do m=1,nterm 
          zall(m) = zall(m-1)*z
          ezall(m) = ezall(m-1)*ez 
        enddo 


        if (ifpgh.gt.1) then 
          vr = (ztargsort(1:2,i)-zc)/zr 
          ve(1) = -vr(2) 
          ve(2) = vr(1)
        endif 

        potsort(i) = potsort(i) + loc(0)


        do m=1,nterm 
          potsort(i) = potsort(i)+loc(m)*zall(m)*ezall(m)
          potsort(i) = potsort(i)+loc(-m)*zall(m)/ezall(m)


          if(ifpgh.gt.1)then 
            dr = m*zall(m-1)*ezall(m)/rscale
            dth = ima*m*zall(m)*ezall(m)
            vtmp = dr*vr+dth*ve/zr
            gradsort(1:2,i) = gradsort(1:2,i)+loc(m)*vtmp

            dr = m*zall(m-1)/ezall(m)/rscale
            dth = -ima*m*zall(m)/ezall(m)
            vtmp = dr*vr+dth*ve/zr
            gradsort(1:2,i) = gradsort(1:2,i)+loc(-m)*vtmp
          endif 
        enddo 

      enddo 

      end 