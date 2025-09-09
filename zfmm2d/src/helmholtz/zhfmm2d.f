C
C   subroutines for 2d complex helmholtz fast multipole method  
C
C       For i=1,...,nt,  evaluate:
C
C                ns
C       u(xi) = \sum charge(j)*1j*H0(zk*R(xi,yj))/4 
C                j=1
C               - dipstr(j)*1j*<dipvec(:,j), \nabla_y H0(zk*R(xi,yj))>/4
C                 
C               where R(x,y)=sqrt((x1-y1)^2+(x2-y2)^2)
C
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine zhfmm2d(eps,zk,ns,zsrc,ifcharge,charge,ifdipole,
     1 dipstr,dipvec,pot,grad,nt,ztarg,ifpgh,isep,ifprint)
c----------------------------------------------
c   INPUT PARAMETERS:
c   eps          : real FMM precision requested 
c   zk           : Helmholtz parameter 
c   ns           : number of sources 
c   zsrc(2,ns)   : source locations 
c   ifcharge     : flag for charges 
c   charge(ns)   : charge strengths 
c   ifdipole     : flag for dipoles
c   dipstr(ns)   : dipole strengths 
c   dipvec(2,ns) : dipole vectors
c   nt           : number of targets 
c   ztarg(2,nt)  : complex target locations 
c   ifpgh        : ifpghtarg = 1, only potential is computed at targets
c                  ifpghtarg = 2, potential and gradient are computed at targets
c                  ifpghtarg = 3 (not available yet), potential, gradient and hessian are computed at targets
c   ifprint      : if print properties 

c   OUTPUT PARAMETERS
c   pot(nt)      : potential at the target locations 
c   grad(2,nt)   : gradient of potential at the target locations 
c----------------------------------------------
      implicit real *8 (a-h,o-z)
      complex *16 zk,ima 
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
      call hndiv2d(eps,ns,nt,ifcharge,ifdipole,ifpgh,ifpghtarg,ndiv,
     1 idivflag)
         
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
        call prinf('ltree=*',ltree,1)
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
        rscales(i) = min(abs(zk*boxsize(i)/(2.0d0*pi)),1.0d0)
        call h2dterms(boxsize(i),zk,eps,nterms(i),ier,isep)
        if(nterms(i).gt.nmax) nmax = nterms(i)
      enddo 
      if (ifprint.eq.1) then 
        call prinf('nmax=*',nmax,1)
        call prin2('rscales=*',rscales,nlevels+1)
        call prinf('nterms=*',nterms,nlevels+1)
      endif 

      call zhfmm2dmain(zk,ns,zsrcsort,chargesort,dipstrsort,dipvecsort,
     1 nt,ztargsort,isrc,itarg,isrcse,itargse,itree,ltree,iptr,nlevels,
     2 nboxes,boxsize,rscales,nterms,nmax,centers,zcsrc,zctarg,    
     3 itree(iptr(1)),itree(iptr(2)),itree(iptr(3)),itree(iptr(4)),
     4 itree(iptr(5)),itree(iptr(6)),isep,pot,grad,ifcharge,
     5 ifdipole,ifpgh,ifprint)

      pot = pot*ima/4 
      grad = grad*ima/4

      end 
C
C
C
C
C
      subroutine zhfmm2dmain(zk,ns,zsrcsort,chargesort,dipstrsort,
     1 dipvecsort,nt,ztargsort,isrc,itarg,isrcse,itargse,itree,ltree,
     2 iptr,nlevels,nboxes,boxsize,rscales,nterms,nmax,centers,zcsrc,
     3 zctarg,laddr,ilevel,iparent,nchild,ichild,ncol,isep,pot,grad,
     4 ifcharge,ifdipole,ifpgh,ifprint)
c   Helmholtz FMM in C^2: evaluate all pairwise particle
      implicit real *8 (a-h,o-z)
      complex *16 zk 
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
      integer laddr(2,0:nlevels),iaddr(4,nboxes)
      integer, allocatable :: list1(:,:),list2(:,:),list3(:,:)
      integer, allocatable :: list4(:,:)
      integer, allocatable :: nlist1s(:),nlist2s(:),nlist3s(:)
      integer, allocatable :: nlist4s(:)
      complex *16 potsort(nt),gradsort(2,nt)
      integer lmptot
      complex *16, allocatable :: rmlexp(:)

c
c        initialize all variables 
c
      pot = 0 
      potsort = 0 
      grad = 0 
      gradsort = 0 

      if (ifprint.eq.1) then 
        call prinf('ifprint=*',ifprint,1)
        call prin2('zk=*',zk,2)
        call prinf('nmax=*',nmax,1)

        do ilev=0,nlevels
          call prinf('laddr=*',laddr(1:2,ilev),2)
        enddo
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



C     allocate space for multipole and local expansions 
      nd = 1
      call h2dmpalloc(nd,itree(iptr(1)),iaddr,nlevels,lmptot,nterms)
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
            call zh2dformmp(zk,rscale,zsrcsort(1,istart),npts,
     1 chargesort(istart),dipstrsort(istart),dipvecsort(1,istart),
     2 zcsrc(1,ibox),nterm,rmlexp(iaddr(1,ibox)),ifcharge,ifdipole)
          endif 
        enddo 
C$OMP END PARALLEL DO
      enddo 
      call cpu_time(t2)
C$        t2=omp_get_wtime()
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
          if (nchild(ibox).gt.0
     1 .and. isrcse(2,ibox).ge.isrcse(1,ibox)) then
            do j = 1,nchild(ibox)
              jbox = ichild(j,ibox)
              if (isrcse(2,jbox).ge.isrcse(1,jbox)) then 
                call zh2dmpmp(zk,zcsrc(1,ibox),zcsrc(1,jbox),rscale,
     1 rscale1,nterm,nterm1,rmlexp(iaddr(1,ibox)),
     2 rmlexp(iaddr(1,jbox)))
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
        if ((nptstarg.gt.0).and.(nlist1s(ibox).gt.0.).and.(
     1  nchild(ibox).eq.0)) then 

          do ilist=1,nlist1s(ibox)
            jbox = list1(ilist,ibox)
            jstart = isrcse(1,jbox)
            jend = isrcse(2,jbox)
            nptssrc = jend-jstart+1

            if (nptssrc.gt.0) then 
              call zh2devaldirect(zk,nptstarg,ztargsort(1,istart),
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
        if ((itargse(2,ibox).ge.itargse(1,ibox))
     1 .and.(nlist2s(ibox).gt.0)) then 
          do ilist=1,nlist2s(ibox) 
            jbox = list2(ilist,ibox)
            if (isrcse(2,jbox).ge.isrcse(1,jbox)) then        
              call zh2dmploc(zk,zctarg(1,ibox),zcsrc(1,jbox),rscale,
     1 nterm,rmlexp(iaddr(2,ibox)),rmlexp(iaddr(1,jbox)))                 
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
              call zh2devalmp(zk,npts,ztargsort(1,istart),rscale,
     1 nterm,zcsrc(1,jbox),rmlexp(iaddr(1,jbox)),potsort(istart),
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
            npts  = jend-jstart+1
            if (npts.gt.0) then 
              call zh2dformloc(zk,rscale,zsrcsort(1,jstart),npts,
     1 chargesort(jstart),dipstrsort(jstart),dipvecsort(1,jstart),
     2 zctarg(1,ibox),nterm,rmlexp(iaddr(2,ibox)),ifcharge,ifdipole)
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
                call zh2dlocloc(zk,zctarg(1,jbox),zctarg(1,ibox),
     1 rscale1,rscale,nterm1,nterm,rmlexp(iaddr(2,jbox)),
     2 rmlexp(iaddr(2,ibox)))
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
          call zh2devalloc(zk,npts,ztargsort(1,istart),rscale,nterm,
     1 zctarg(1,ibox),rmlexp(iaddr(2,ibox)),potsort(istart),
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
C
C
C
C
C
      subroutine zh2dformmp(zk,rscale,zsrcsort,npts,chargesort,
     1 dipstrsort,dipvecsort,zcsrc,nterm,mp,ifcharge,ifdipole)
      implicit real *8 (a-h,o-z)
      complex *16 zk,zsrcsort(2,npts),chargesort(npts)
      complex *16 dipstrsort(npts),dipvecsort(2,npts)
      complex *16 zcsrc(2),mp(-nterm:nterm)
      complex *16 ez,rz,z,ima,dr,dth 
      complex *16 fjs(0:(nterm+100)),fjder(0:(nterm+100))
      complex *16 vr(2),ve(2),pr,pe
      complex *16 ezs(-nterm:nterm)

      ima = (0.0d0,1.0d0)
      ifder = 0 
      if (ifdipole.eq.1) then 
        ifder = 1 
      endif 

      ezs(0) = 1 
      do i=1,npts 

        call zcart2polar(zsrcsort(1,i),zcsrc,rz,ez)
        z = rz*zk

        do m=1,nterm 
          ezs(m) = ezs(m-1)*ez
        enddo 

        do m=-nterm,-1 
          ezs(m) = 1/ezs(-m)
        enddo 

        call jbessel2d(nterm,z,rscale,fjs,ifder,fjder)
c
        if (ifdipole.eq.1) then 
          vr = (zsrcsort(1:2,i)-zcsrc)/rz
          ve(1) = -vr(2)
          ve(2) = vr(1)


          pr = vr(1)*dipvecsort(1,i)+vr(2)*dipvecsort(2,i)
          pe = ve(1)*dipvecsort(1,i)+ve(2)*dipvecsort(2,i)
        endif 

        do m=0,nterm 

          if (ifcharge.eq.1) then 
            mp(m) = mp(m)+fjs(m)*ezs(-m)*chargesort(i)
          endif 
c
          if (ifdipole.eq.1) then 
            dr = ezs(-m)*fjder(m)*zk
            dth = fjs(m)*(ima*(-m))*ezs(-m)

            mp(m) = mp(m)-dipstrsort(i)*(dr*pr+dth*pe/rz)
          endif 

        enddo 
        
        do m=-nterm,-1

          if (ifcharge.eq.1) then 
            mp(m) = mp(m)+((-1)**m)*fjs(-m)*ezs(-m)*chargesort(i)  
          endif   
c
c
          if (ifdipole.eq.1) then 
            dr = ezs(-m)*fjder(-m)*zk
            dth = fjs(-m)*(ima*(-m))*ezs(-m)
            
            mp(m) = mp(m)-dipstrsort(i)*((-1)**m)*(dr*pr+dth*pe/rz)
          endif 

        enddo 

      enddo 

      end 
C
C
C
C
C
      subroutine zh2dformloc(zk,rscale,zsrcsort,npts,chargesort,
     1 dipstrsort,dipvecsort,zctarg,nterm,loc,ifcharge,ifdipole)
      implicit real *8 (a-h,o-z)
      complex *16 zk,zsrcsort(2,npts),chargesort(npts)
      complex *16 dipstrsort(npts),dipvecsort(2,npts)
      complex *16 zctarg(2),loc(-nterm:nterm)
      complex *16 ez,rz,z,ima,dr,dth 
      complex *16 fhs(0:nterm+100),fhder(0:nterm+100)
      complex *16 vr(2),ve(2),pr,pe
      complex *16 ezs(-nterm:nterm)

      ima = (0.0d0,1.0d0)
      ifder = 0 
      if (ifdipole.eq.1) then 
        ifder = 1 
      endif 

      ezs(0) = 1

      do i=1,npts 

        call zcart2polar(zsrcsort(1,i),zctarg,rz,ez)
        z = rz*zk

        do m = 1,nterm 
          ezs(m) = ezs(m-1)*ez 
        enddo 

        do m = -nterm,-1 
          ezs(m) = 1/ezs(-m)
        enddo  

        call h2dall(nterm,z,rscale,fhs,ifder,fhder)

c
        if (ifdipole.eq.1) then 
          vr = (zsrcsort(1:2,i)-zctarg)/rz
          ve(1) = -vr(2)
          ve(2) = vr(1)
          pr = vr(1)*dipvecsort(1,i)+vr(2)*dipvecsort(2,i)
          pe = ve(1)*dipvecsort(1,i)+ve(2)*dipvecsort(2,i)
        endif 



        do m=0,nterm 
          if (ifcharge.eq.1) then 
            loc(m) = loc(m)+fhs(m)*(ez**(-m))*chargesort(i)
          endif 

          if (ifdipole.eq.1) then 
            dr = (ezs(-m))*fhder(m)*zk
            dth = fhs(m)*(ima*(-m))*(ezs(-m))
            loc(m) = loc(m)-dipstrsort(i)*(dr*pr+dth*pe/rz)
          endif 
        enddo 
        
        do m=-nterm,-1

          if (ifcharge.eq.1) then 
            loc(m) = loc(m)+((-1)**m)*fhs(-m)*(ezs(-m))*chargesort(i)
          endif 

          if (ifdipole.eq.1) then 
            dr = (ezs(-m))*fhder(-m)*zk
            dth = fhs(-m)*(ima*(-m))*(ezs(-m))
            
            loc(m) = loc(m)-dipstrsort(i)*((-1)**m)*(dr*pr+dth*pe/rz)
          endif 

        enddo 

      enddo 

      end 

C
C
C
C
C
      subroutine zh2devaldirect(zk,nt,ztargsort,ns,zsrcsort,chargesort,
     1 dipstrsort,dipvecsort,ifcharge,potsort,gradsort,ifdipole,ifpgh)
      implicit real *8 (a-h,o-z)
      complex *16 zk,rz,h0,h1,h2,dh1,z
      complex *16 tmp,tmp1,tmp2,gradtmp 
      complex *16 dipstrsort(ns),dipvecsort(2,ns)
      complex *16 ztargsort(2,nt),potsort(nt),gradsort(2,nt)
      complex *16 zsrcsort(2,ns),chargesort(ns)
      complex *16 y1mx1,y2mx2,x1my1,x2my2 
      complex *16 hvec(0:10)

      do i=1,nt 
        do j=1,ns 
c          call zdist(2,ztargsort(1,i),zsrcsort(1,j),rz)
          y1mx1 = zsrcsort(1,j)-ztargsort(1,i)
          y2mx2 = zsrcsort(2,j)-ztargsort(2,i)
          x1my1 = -y1mx1 
          x2my2 = -y2mx2 
          rz = sqrt(y1mx1**2+y2mx2**2)
          z = rz*zk 

          call hank101(z,h0,h1)

          if (abs(z).gt.1E-16 .and. aimag(z).gt.-32) then 

            if (ifcharge.eq.1) then 
              potsort(i) = potsort(i)+chargesort(j)*h0
            endif 

            if (ifdipole.eq.1) then 
              potsort(i) = potsort(i)+dipstrsort(j)*dipvecsort(1,j)*h1
     1           *zk*y1mx1/rz 
              potsort(i) = potsort(i)+dipstrsort(j)*dipvecsort(2,j)*h1
     1           *zk*y2mx2/rz 
            endif 

            if(ifpgh.gt.1)then 

              if (ifcharge.eq.1) then 
                gradsort(1:2,i) = gradsort(1:2,i)-chargesort(j)*zk*h1*
     1    (ztargsort(1:2,i)-zsrcsort(1:2,j))/rz 
              endif 

              if (ifdipole.eq.1) then 
                h2 = (2/z)*h1-h0
                dh1 = (h0-h2)/2

                tmp = dipvecsort(1,j)*y1mx1+dipvecsort(2,j)*y2mx2
                tmp1 = zk*dh1*x1my1/(rz**2)-h1*x1my1/(rz**3)

                tmp2 = zk*dh1*x2my2/(rz**2)-h1*x2my2/(rz**3)

                gradsort(1,i) = gradsort(1,i)+dipstrsort(j)*zk*tmp1*tmp
     1  -dipstrsort(j)*(zk*h1/rz)*dipvecsort(1,j)
                
                gradsort(2,i) = gradsort(2,i)+dipstrsort(j)*zk*tmp2*tmp
     1  -dipstrsort(j)*(zk*h1/rz)*dipvecsort(2,j)
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
      subroutine zh2devalmp(zk,nt,ztargsort,rscale,nterm,zc,mp,potsort,
     1     gradsort,ifpgh)
      implicit real *8 (a-h,o-z)
      complex *16 zk,zr,ei,z,ima,dr,dth 
      complex *16 ztargsort(2,nt),potsort(nt),gradsort(2,nt)
      complex *16 zc(2),mp(-nterm:nterm)
      complex *16 fhs(0:nterm+100),fhder(0:nterm+100)
      complex *16 vr(2),ve(2)
      complex *16 eis(-nterm:nterm)
      

      ifder = 0 
      if (ifpgh.gt.1) then 
        ifder = 1 
      endif 

      ima = (0.0d0,1.0d0)
      eis(0) = 1
      do i=1,nt 
        call zcart2polar(ztargsort(1,i),zc,zr,ei)

        do m=1,nterm 
          eis(m) = eis(m-1)*ei 
        enddo 

        do m = -nterm,-1  
          eis(m) = 1/eis(-m)
        enddo 

        z = zr*zk
        call h2dall(nterm,z,rscale,fhs,ifder,fhder)

        do m=0,nterm 
          potsort(i) = potsort(i)+mp(m)*fhs(m)*eis(m)

          if(ifpgh.gt.1)then 
            vr = (ztargsort(1:2,i)-zc)/zr 
            ve(1) = -vr(2) 
            ve(2) = vr(1)
            dr = eis(m)*fhder(m)*zk
            dth = fhs(m)*(ima*m)*eis(m)
            gradsort(1:2,i) = gradsort(1:2,i)+mp(m)*(dr*vr+dth*ve/zr)
          endif 
        enddo 

        do m=-nterm,-1 
          potsort(i) = potsort(i)+mp(m)*((-1)**m)*fhs(-m)*eis(m)

          if(ifpgh.gt.1)then 
            dr = eis(m)*fhder(-m)*zk
            dth = fhs(-m)*(ima*m)*eis(m)
            gradsort(1:2,i) = gradsort(1:2,i)+
     1           ((-1)**m)*mp(m)*(dr*vr+dth*ve/zr)
          endif 
        enddo 
      enddo 

      end 
C
C
C
C
C
      subroutine zh2devalloc(zk,nt,ztargsort,rscale,nterm,zc,loc,
     1  potsort,gradsort,ifpgh)
      implicit real *8 (a-h,o-z)
      complex *16 zk,zr,ei,z,ima,dr,dth 
      complex *16 ztargsort(2,nt),potsort(nt),gradsort(2,nt)
      complex *16 zc(2),loc(-nterm:nterm)
      complex *16 fjs(0:(nterm+100)),fjder(0:(nterm+100))
      complex *16 vr(2),ve(2)
      complex *16 eis(-nterm:nterm)

      ima = (0.0d0,1.0d0)
      ifder = 0 
      if (ifpgh.gt.1) then 
        ifder = 1 
      endif 
    
      z0 = 0
      eis(0) = 1 
      do i=1,nt 
        call zcart2polar(ztargsort(1,i),zc,zr,ei)

        do m=1,nterm 
          eis(m) = eis(m-1)*ei 
        enddo 

        do m = -nterm,-1  
          eis(m) = 1/eis(-m)
        enddo 

        z = zr*zk
        call jbessel2d(nterm,z,rscale,fjs,ifder,fjder)

        do m=0,nterm 
          potsort(i) = potsort(i)+loc(m)*fjs(m)*(ei**m)

          if (ifpgh.gt.1) then 
            vr = (ztargsort(1:2,i)-zc)/zr 
            ve(1) = -vr(2)
            ve(2) = vr(1)
            dr = eis(m)*fjder(m)*zk
            dth = fjs(m)*(ima*m)*eis(m)
            gradsort(1:2,i) = gradsort(1:2,i)+loc(m)*(dr*vr+dth*ve/zr)
          endif 
        enddo 

        do m=-nterm,-1 
          potsort(i) = potsort(i)+loc(m)*((-1)**m)*fjs(-m)*eis(m)

          if(ifpgh.gt.1)then 
            dr = eis(m)*fjder(-m)*zk
            dth = fjs(-m)*(ima*m)*eis(m)
            gradsort(1:2,i) = gradsort(1:2,i)+
     1           ((-1)**m)*loc(m)*(dr*vr+dth*ve/zr)
          endif 
        enddo 
      enddo 

      end 
c
c
c
c
c
c
      subroutine h2dmpalloc(nd,laddr,iaddr,nlevels,lmptot,
     1                          nterms)
c     This subroutine determines the size of the array
c     to be allocated for the multipole expansions
c     iaddr(1,i) points to the starting location of the multipole
c     expansion of box i 
c     iaddr(2,i) points to the local
c     expansion of box i
c     iaddr(3,i) points to the outgoing diag form of box i
c     iaddr(4,i) points to the incoming diag form of box i
c  
c     Input arguments
c     nd          in: integer
c                 number of expansions
c
c     laddr       in: Integer(2,0:nlevels)
c                 indexing array providing access to boxes at each
c                 level
c
c     nlevels     in: Integer
c                 Total numner of levels
c     
c     nterms      in: Integer(0:nlevels)
c                 Number of terms requried in expansions at each
c                 level
c
c------------------------------------------------------------------
c     Output arguments
c     iaddr       out: Integer(4,nboxes)
c                 Points the multipole and local expansions in box i
c 
c     lmptot      out: Integer
c                 Total length of expansions array required
c------------------------------------------------------------------

      implicit none
      integer nlevels,nterms(0:nlevels),nd,nsig,nt1,nt2,next235
      integer iaddr(4,*), lmptot, laddr(2,0:nlevels)
      integer ibox,i,iptr,istart,nn,itmp
      real *8 ddn,dn
c
      istart = 1
      do i = 0,nlevels
         nn = (2*nterms(i)+1)*2*nd
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,itmp)
         do ibox = laddr(1,i),laddr(2,i)
c     Allocate memory for the multipole expansion         
c
            itmp = ibox - laddr(1,i)
            iaddr(1,ibox) = istart + itmp*nn
         enddo
C$OMP END PARALLEL DO         
         istart = istart + (laddr(2,i)-laddr(1,i)+1)*nn
      enddo
c
      do i=0,nlevels
         nn = (2*nterms(i)+1)*2*nd
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,itmp)
         do ibox = laddr(1,i),laddr(2,i)
c     Allocate memory for the local expansion
c
            itmp = ibox - laddr(1,i)
            iaddr(2,ibox) = istart + itmp*nn 
         enddo
C$OMP END PARALLEL DO         
         istart = istart + (laddr(2,i)-laddr(1,i)+1)*nn
      enddo
c
      do i=0,nlevels
         dn = 2*(nterms(i)+nterms(i)) + 1
         nn = 2*nd*next235(dn)
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,itmp)
         do ibox = laddr(1,i),laddr(2,i)
c     Allocate memory for the local expansion
c
            itmp = ibox - laddr(1,i)
            iaddr(3,ibox) = istart + itmp*nn 
         enddo
C$OMP END PARALLEL DO         
         istart = istart + (laddr(2,i)-laddr(1,i)+1)*nn
      enddo
c
      do i=0,nlevels
         dn = 2*(nterms(i)+nterms(i)) + 1
         nn = 2*nd*next235(dn)
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,itmp)
         do ibox = laddr(1,i),laddr(2,i)
c     Allocate memory for the local expansion
c
            itmp = ibox - laddr(1,i)
            iaddr(4,ibox) = istart + itmp*nn 
         enddo
C$OMP END PARALLEL DO         
         istart = istart + (laddr(2,i)-laddr(1,i)+1)*nn
      enddo
c
      lmptot = istart

      return
      end