      implicit real *8 (a-h,o-z)
      complex *16 ima 
      integer unit,iostat

      real *8, allocatable :: charge_real(:),charge_imag(:)
      real *8, allocatable :: dipstr_real(:),dipstr_imag(:)
      real *8, allocatable :: dipvec_real(:,:),dipvec_imag(:,:)
      real *8, allocatable :: zsrc_real(:,:),zsrc_imag(:,:)
      real *8, allocatable :: ztarg_real(:,:),ztarg_imag(:,:)

      complex *16, allocatable :: charge(:)
      complex *16, allocatable :: dipstr(:)
      complex *16, allocatable :: dipvec(:,:)
      complex *16, allocatable :: zsrc(:,:)
      complex *16, allocatable :: ztarg(:,:) 

      complex *16, allocatable :: pot(:),pot0(:)
      complex *16, allocatable :: grad(:,:),grad0(:,:)

      call prini(6,13)

      ima = (0.0d0,1.0d0)
c     source and target location will be real if this is off 
      ifcplx = 1
      isep = 1
      done = 1 
      pi = atan(done)*4
      unit = 10


      open(unit,file='test_data_2d.bin',form='unformatted',
     1 access='stream',status='old')
     
      read(unit) npt 
      call prinf('nt=*',npt,1)

      ns = npt 
      nt = npt 

      allocate(dipstr_real(ns),dipstr_imag(ns))
      allocate(dipvec_real(2,ns),dipvec_imag(2,ns))
      allocate(zsrc_real(2,ns),zsrc_imag(2,ns))
      allocate(ztarg_real(2,nt),ztarg_imag(2,nt))
      allocate(charge_real(ns),charge_imag(ns))


      allocate(dipstr(ns))
      allocate(dipvec(2,ns))
      allocate(zsrc(2,ns))
      allocate(ztarg(2,nt))
      allocate(charge(ns))
      allocate(pot(nt),pot0(nt))
      allocate(grad(2,nt),grad0(2,nt))


      read(unit) zsrc_real  

      read(unit) zsrc_imag

      read(unit) ztarg_real

      read(unit) ztarg_imag

      read(unit) charge_real

      read(unit) charge_imag

      read(unit) dipvec_real

      read(unit) dipvec_imag

      read(unit) dipstr_real

      read(unit) dipstr_imag



      close(unit)


      do i=1,ns  
        dipstr(i) = dipstr_real(i)+ima*dipstr_imag(i)
        charge(i) = charge_real(i)+ima*charge_imag(i)
        do j=1,2
          dipvec(j,i) = dipvec_real(j,i)+ima*dipvec_imag(j,i)
          zsrc(j,i) = zsrc_real(j,i)+ifcplx*ima*zsrc_imag(j,i)
        enddo 
      enddo 

      do i=1,nt 
        do j=1,2
          ztarg(j,i) = zsrc(j,i)
        enddo 
      enddo 


      eps = 1E-12
      ifcharge = 1
      ifdipole = 0
      ifpgh = 1
      ifprint = 1
      pot = 0 
      grad = 0
      call cpu_time(t1)
C$        t1=omp_get_wtime()
      call zlfmm2d(eps,ns,zsrc,ifcharge,charge,ifdipole,dipstr,
     1 dipvec,pot,grad,nt,ztarg,ifpgh,isep,ifprint)
      call cpu_time(t2)
C$        t2=omp_get_wtime()
      call prin2('pot=*',pot,24)
      call prin2('fmm spends *',t2-t1,1)


      if (1.eq.1) then 
        pot0 = 0  
        grad0 = 0
        call cpu_time(t1)
C$        t1=omp_get_wtime()
        call zl2devaldirect(nt,ztarg,ns,zsrc,charge,dipstr,dipvec,
     1 ifcharge,pot0,grad0,ifdipole,ifpgh)
        call cpu_time(t2)
C$        t2=omp_get_wtime()
        call prin2('pot0=*',pot0,24)
        call prin2('direct method spends *',t2-t1,1)

        charge_norm = 0 
        pot0_norm = 0
        pot_err = 0 
        do i=1,nt
          pot0_norm = pot0_norm+abs(pot0(i))**2
          pot_err = pot_err+abs(pot0(i)-pot(i))**2
          charge_norm = charge_norm+abs(charge(i))**2
        enddo 
        pot0_norm = sqrt(pot0_norm)
        charge_norm = sqrt(charge_norm)
        pot_err = sqrt(pot_err)
        rel_err = pot_err/(pot0_norm+charge_norm)
        call prin2('pot0_norm=*',pot0_norm,1)
        call prin2('charge_norm=*',charge_norm,1)
        call prin2('pot_err=*',pot_err,1)
        call prin2('rel err of pot=*',rel_err,1)
      endif 


      end 