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

      isep = 1
      ima = (0.0d0,1.0d0)
      done = 1 

      
      unit=10
      open(unit,file='test_data.bin',form='unformatted',
     1 access='stream',status='old')
     
      read(unit) npt 
      ns = npt 
      nt = npt 
      allocate(dipstr_real(ns),dipstr_imag(ns))
      allocate(dipvec_real(3,ns),dipvec_imag(3,ns))
      allocate(zsrc_real(3,ns),zsrc_imag(3,ns))
      allocate(ztarg_real(3,nt),ztarg_imag(3,nt))
      allocate(charge_real(ns),charge_imag(ns))


      allocate(dipstr(ns))
      allocate(dipvec(3,ns))
      allocate(zsrc(3,ns))
      allocate(ztarg(3,nt))
      allocate(charge(ns))
      allocate(pot(nt),pot0(nt))
      allocate(grad(3,nt),grad0(3,nt))

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
        do j=1,3 
          dipvec(j,i) = dipvec_real(j,i)+ima*dipvec_imag(j,i)
          zsrc(j,i) = zsrc_real(j,i)+ima*zsrc_imag(j,i)
        enddo 
      enddo 


      do i=1,nt 
        do j=1,3 
          ztarg(j,i) = ztarg_real(j,i)+ima*ztarg_imag(j,i)
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
      call zlfmm3d(eps,ns,zsrc,ifcharge,charge,ifdipole,dipstr,
     1 dipvec,pot,grad,nt,ztarg,ifpgh,isep,ifprint)
      call cpu_time(t2)
C$        t2=omp_get_wtime()
      call prin2('pot=*',pot,24)
      call prin2('fmm spends *',t2-t1,1)



      pot0 = 0 
      grad0 = 0
      call cpu_time(t1)
C$        t1=omp_get_wtime()
      call zl3devaldirect(nt,ztarg,ns,zsrc,charge,dipstr,dipvec,
     1 ifcharge,pot0,grad0,ifdipole,ifpgh)
      call cpu_time(t2)
C$        t2=omp_get_wtime()
      call prin2('pot0=*',pot0,24)
      call prin2('direct method spends *',t2-t1,1)

 
c      call getl2err(1,npt,pot0,pot,err)
c      call prin2('rel err of pot=*',err,1)

      err = 0
      err_i = 0  
      err_max = 0 
      i_max = 0
      pot0_norm = 0
      charge_norm = 0 
      do i=1,nt
        err_i = abs(pot0(i)-pot(i))
        err = err+err_i**2
        if (err_i.gt.err_max) then 
          i_max = i 
          err_max = err_i 
        endif 
        pot0_norm = pot0_norm+abs(pot0(i))**2
        charge_norm=charge_norm+abs(charge(i))**2
      enddo
      err = sqrt(err) 
      charge_norm = sqrt(charge_norm)
      call prin2('err=*',err,1)
      pot0_norm = sqrt(pot0_norm)
      call prin2('pot0_norm=*',pot0_norm,1)
      call prin2('rel err=*',err/(pot0_norm+charge_norm),1)
      call prin2('charge norm=*',charge_norm,1)

      end 