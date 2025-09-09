      subroutine zl3dmpmp(x0y0z0,nterm0,mp0,rscale0,xyz,nterm,mp,
     1  rscale,nmax,amat)
      implicit real *8 (a-h,o-z)
      complex *16 x0y0z0(3),xyz(3)
      complex *16 mp0(0:nterm0,-nterm0:nterm0)
      complex *16 mp(0:nterm,-nterm:nterm)
      complex *16 rz,ctheta,eiphi
      dimension amat(0:nmax,-nmax:nmax)
      complex *16 xc1(3),xc2(3),zshift
      complex *16 cbeta,sbeta
      complex *16 calpha,salpha 

      complex *16, allocatable :: mp0z(:,:)
      complex *16, allocatable :: mp0zy(:,:)
      complex *16, allocatable :: mpzy(:,:)
      complex *16, allocatable :: mpz(:,:)
      complex *16, allocatable :: mpn(:,:)

      allocate(mp0z(0:nterm0,-nterm0:nterm0))
      allocate(mp0zy(0:nterm0,-nterm0:nterm0))
      allocate(mpzy(0:nterm,-nterm:nterm))
      allocate(mpz(0:nterm,-nterm:nterm))
      allocate(mpn(0:nterm,-nterm:nterm))


c      call lmp_ini(nterm0,mp0z)
c      call lmp_ini(nterm0,mp0zy)
c      call lmp_ini(nterm,mpzy)
c      call lmp_ini(nterm,mpz)
c      call lmp_ini(nterm,mpn)
      mp0z = 0
      mp0zy = 0
      mpzy = 0 
      mpz = 0
      mpn = 0

      e2 = abs(xyz(2)-x0y0z0(2))
      if (e2.le.1E-16) then 
        xc1 = xyz-x0y0z0 
        cbeta = 1 
        sbeta = 0
        do n=0,nterm0 
          do m=-n,n 
            mp0z(n,m) = mp0(n,m) 
          enddo 
        enddo 
      else 
        call z3dzrotini(xyz,x0y0z0,xc1,cbeta,sbeta)
        call z3dzrot(nterm0,mp0,mp0z,cbeta,sbeta)
      endif 

      e1 = abs(xc1(1))
      if (e1.le.1E-16) then
        xc2 = xc1 
        calpha = 1 
        salpha = 0 
        do n=0,nterm0
          do m=-n,n
            mp0zy(n,m) = mp0z(n,m) 
          enddo 
        enddo 
      else 
        call z3dyrotini(xc1,xc2,calpha,salpha)
        call z3dyrot(nterm0,mp0z,mp0zy,calpha,salpha)
      endif  

      zshift = -xc2(3)
      call zl3dmpmp_zshift(zshift,nterm0,mp0zy,rscale0,
     1     nterm,mpzy,rscale,nmax,amat)

      e1 = abs(xc1(1))
      if (e1.le.1E-16) then
        do n=0,nterm 
          do m=-n,n
              mpz(n,m) = mpzy(n,m)
          enddo 
        enddo    
      else 
        salpha = -salpha 
        call z3dyrot(nterm,mpzy,mpz,calpha,salpha)
      endif  


      if (e2.le.1E-16) then 
        do n=0,nterm 
          do m=-n,n
            mpn(n,m) = mpz(n,m) 
          enddo 
        enddo 
      else 
        sbeta = -sbeta 
        call z3dzrot(nterm,mpz,mpn,cbeta,sbeta)
      endif 

      do n=0,nterm 
        do m=-n,n 
          mp(n,m) = mp(n,m)+mpn(n,m)
        enddo 
      enddo 

      end
C
C
C
C
C
      subroutine zl3dmploc(x0y0z0,nterm0,mp0,rscale0,xyz,nterm,loc,
     1 rscale,nmax,amat)
      implicit real *8 (a-h,o-z)
      complex *16 x0y0z0(3),xyz(3)
      complex *16 mp0(0:nterm0,-nterm0:nterm0)
      complex *16 loc(0:nterm,-nterm:nterm)
      complex *16 rz,ctheta,eiphi
      dimension amat(0:nmax,-nmax:nmax)
      complex *16 xc1(3),xc2(3),zshift
      complex *16 cbeta,sbeta
      complex *16 calpha,salpha

      complex *16, allocatable :: mp0z(:,:)
      complex *16, allocatable :: mp0zy(:,:)
      complex *16, allocatable :: loczy(:,:)
      complex *16, allocatable :: locz(:,:)
      complex *16, allocatable :: locn(:,:) 

      allocate(mp0z(0:nterm0,-nterm0:nterm0))
      allocate(mp0zy(0:nterm0,-nterm0:nterm0))
      allocate(loczy(0:nterm,-nterm:nterm))
      allocate(locz(0:nterm,-nterm:nterm))
      allocate(locn(0:nterm,-nterm:nterm))

      mp0z = 0
      mp0zy = 0 
      loczy = 0 
      locz = 0 
      locn = 0 

c      call cpu_time(t1)
      e2 = abs(xyz(2)-x0y0z0(2))
      if (e2.le.1E-16) then 
        xc1 = xyz-x0y0z0 
        cbeta = 1 
        sbeta = 0
        do n=0,nterm0 
          do m=-n,n 
            mp0z(n,m) = mp0(n,m) 
          enddo 
        enddo 
      else 
        call z3dzrotini(xyz,x0y0z0,xc1,cbeta,sbeta)
c        call lmp_ini(nterm0,mp0z)
        call z3dzrot(nterm0,mp0,mp0z,cbeta,sbeta)
      endif 
c      call cpu_time(t2)
c      call prin2('z rotate spends *',t2-t1,1)

c      call cpu_time(t1)
      e1 = abs(xc1(1))
      if (e1.le.1E-16) then
        xc2 = xc1 
        calpha = 1 
        salpha = 0 
        do n=0,nterm0
          do m=-n,n
            mp0zy(n,m) = mp0z(n,m) 
          enddo 
        enddo 
      else 
        call z3dyrotini(xc1,xc2,calpha,salpha)
c        mp0zy = 0
        call z3dyrot(nterm0,mp0z,mp0zy,calpha,salpha)
      endif  
c      call cpu_time(t2)
c      call prin2('y rotate spends *',t2-t1,1)


c      call cpu_time(t1)
      zshift = -xc2(3)
      call lmp_ini(nterm,loczy)
      call zl3dmploc_zshift(zshift,nterm0,mp0zy,rscale0,
     1     nterm,loczy,rscale,nmax,amat)
c      call cpu_time(t2)
c      call prin2('z shift spends *',t2-t1,1)

      e1 = abs(xc1(1))
      if (e1.le.1E-16) then
        do n=0,nterm 
          do m=-n,n
            locz(n,m) = loczy(n,m)
          enddo 
        enddo    
      else 
        salpha = -salpha 
c        locz = 0
        call z3dyrot(nterm,loczy,locz,calpha,salpha)
      endif  


      if (e2.le.1E-16) then 
        do n=0,nterm 
          do m=-n,n
            locn(n,m) = locz(n,m) 
          enddo 
        enddo 
      else 
        sbeta = -sbeta 
c        locn = 0
        call z3dzrot(nterm,locz,locn,cbeta,sbeta)
      endif 

      do n=0,nterm 
        do m=-n,n 
          loc(n,m) = loc(n,m)+locn(n,m)
        enddo 
      enddo 

      end
C
C
C
C
C
      subroutine zl3dlocloc(x0y0z0,nterm0,loc0,rscale0,xyz,nterm,loc,
     1 rscale,nmax,amat)
      implicit real *8 (a-h,o-z)
      complex *16 x0y0z0(3),xyz(3)
      complex *16 loc0(0:nterm0,-nterm0:nterm0)
      complex *16 loc(0:nterm,-nterm:nterm)
      complex *16 rz,ctheta,eiphi
      dimension amat(0:nmax,-nmax:nmax)
      complex *16 xc1(3),xc2(3),zshift
      complex *16 cbeta,sbeta
      complex *16 calpha,salpha 

      complex *16, allocatable :: loc0z(:,:)
      complex *16, allocatable :: loc0zy(:,:)
      complex *16, allocatable :: loczy(:,:)
      complex *16, allocatable :: locz(:,:)
      complex *16, allocatable :: locn(:,:)

      allocate(loc0z(0:nterm0,-nterm0:nterm0))
      allocate(loc0zy(0:nterm0,-nterm0:nterm0))
      allocate(loczy(0:nterm,-nterm:nterm))
      allocate(locz(0:nterm,-nterm:nterm))
      allocate(locn(0:nterm,-nterm:nterm))

c      call lmp_ini(nterm0,loc0z)
c      call lmp_ini(nterm0,loc0zy)
c      call lmp_ini(nterm,loczy)
c      call lmp_ini(nterm,locz)
c      call lmp_ini(nterm,locn)
      loc0z = 0 
      loc0zy = 0 
      loczy = 0 
      locz = 0 
      locn = 0 

      e2 = abs(xyz(2)-x0y0z0(2))
      if (e2.le.1E-16) then 
        xc1 = xyz-x0y0z0 
        cbeta = 1 
        sbeta = 0
        do n=0,nterm0 
          do m=-n,n 
            loc0z(n,m) = loc0(n,m) 
          enddo 
        enddo 
      else 
        call z3dzrotini(xyz,x0y0z0,xc1,cbeta,sbeta)
        call z3dzrot(nterm0,loc0,loc0z,cbeta,sbeta)
      endif 

      e1 = abs(xc1(1))
      if (e1.le.1E-16) then
        xc2 = xc1 
        calpha = 1 
        salpha = 0 
        do n=0,nterm0
          do m=-n,n
            loc0zy(n,m) = loc0z(n,m) 
          enddo 
        enddo 
      else 
        call z3dyrotini(xc1,xc2,calpha,salpha)
        call z3dyrot(nterm0,loc0z,loc0zy,calpha,salpha)
      endif  

      zshift = -xc2(3)
      call zl3dlocloc_zshift(zshift,nterm0,loc0zy,rscale0,
     1     nterm,loczy,rscale,nmax,amat)

      e1 = abs(xc1(1))
      if (e1.le.1E-16) then
        do n=0,nterm 
          do m=-n,n
            locz(n,m) = loczy(n,m)
          enddo 
        enddo    
      else 
        salpha = -salpha 
        call z3dyrot(nterm,loczy,locz,calpha,salpha)
      endif  


      if (e2.le.1E-16) then 
        do n=0,nterm 
          do m=-n,n
            locn(n,m) = locz(n,m) 
          enddo 
        enddo 
      else 
        sbeta = -sbeta 
        call z3dzrot(nterm,locz,locn,cbeta,sbeta)
      endif 

      do n=0,nterm 
        do m=-n,n 
          loc(n,m) = loc(n,m)+locn(n,m)
        enddo 
      enddo 

      end
c C
c C
c C
c C
C
      subroutine zl3dmpmp_zshift(zshift,nterm0,mp0,rscale0,nterm,mp,
     1 rscale,nmax,amat)
      implicit real *8 (a-h,o-z)
      complex *16 zshift
      complex *16 mp0(0:nterm0,-nterm0:nterm0)
      complex *16 mp(0:nterm,-nterm:nterm)
      complex *16 rz,ctheta
      complex *16 ima,tmp
      dimension amat(0:nmax,-nmax:nmax)

      real *8, allocatable :: rscale0s(:),rscales(:)
      complex *16, allocatable :: rzall(:)
      real *8, allocatable :: yl(:)

      allocate(rscale0s(0:nterm0),rscales(0:nterm))
      allocate(rzall(0:nterm),yl(0:nterm))

      done = 1 
      ima = (0.0d0,1.0d0)
      pi = atan(done)*4
      pi4 = pi*4 
      pi4_sqrt = dsqrt(pi4) 

      rscale0s(0) = 1 
      rscales(0) = 1 
      do n=1,nterm0 
        rscale0s(n) = rscale0s(n-1)*rscale0 
      enddo 

      do n=1,nterm
        rscales(n) = rscales(n-1)*rscale 
      enddo 

      if (abs(sqrt(zshift**2)-zshift)>1E-16) then 
        rz = -zshift 
        do n=0,nterm 
          yl(n) = dsqrt(2*n+1.0d0) * (-1)**n 
        enddo 
      else 
        rz = zshift
        do n=0,nterm 
          yl(n) = dsqrt(2*n+1.0d0) 
        enddo 
      endif 


      rz = rz/rscale0
      rzall(0) = 1 
      do n=1,nterm 
        rzall(n) = rzall(n-1)*rz
      enddo 

      do j=0,nterm 
        tmpj = pi4_sqrt*(rscale0s(j)/rscales(j))
        do k=-j,j 
          do n=0,j 
            if ((abs(k).le.(j-n)).and.((j-n).le.nterm0))then 
              tmp = mp0(j-n,k)*amat(n,0)*amat(j-n,k)/amat(j,k)
              tmp = tmp*rzall(n)*yl(n)/(2*n+done)
              mp(j,k) = mp(j,k)+tmp* tmpj
            endif 
          enddo 
        enddo 
      enddo 

      end

C
C
C
C
C
      subroutine zl3dmploc_zshift(zshift,nterm0,mp0,rscale0,nterm,loc,
     1  rscale,nmax,amat)
      implicit real *8 (a-h,o-z)
      complex *16 zshift
      complex *16 mp0(0:nterm0,-nterm0:nterm0)
      complex *16 loc(0:nterm,-nterm:nterm)
      complex *16 rz,eiphi,ctheta
      complex *16 ima,tmp,tmpj
      dimension amat(0:nmax,-nmax:nmax)

      real *8, allocatable :: rscale0s(:),rscales(:),yl(:)
      integer, allocatable :: m1_power(:)
      complex *16, allocatable :: rzall(:)

      allocate(m1_power(0:nterm+nterm0))
      allocate(rscale0s(0:nterm+1),rscales(0:nterm+1))
      allocate(rzall(0:nterm+nterm0),yl(0:nterm+nterm0))

      ima = (0.0d0,1.0d0)
      done  = 1
      pi = atan(done)*4
      pi4 = pi*4 
      pi4_sqrt = sqrt(pi4)

      m1_power(0) =1 
      do m=1,nterm+nterm0 
        m1_power(m)=(-1)*m1_power(m-1)
      enddo 

      if (abs(sqrt(zshift**2)-zshift)>1E-16) then 
        rz = -zshift 
        do n=0,nterm + nterm0
          yl(n) = dsqrt(2*n+1.0d0) * (-1)**n 
        enddo 
      else 
        rz = zshift
        do n=0,nterm + nterm0
          yl(n) = dsqrt(2*n+1.0d0) 
        enddo 
      endif 

      rz = rz/rscale0
      rzall(0) = 1 
      do n=1,nterm+nterm0
        rzall(n) = rzall(n-1)*rz
      enddo 


      do j=0,nterm 
        tmpj = pi4_sqrt/rz/(2*j+done)
        do k=-j,j 
          do n=0,nterm0
            if ((abs(k).le.n))then
              tmp = mp0(n,k)*amat(n,k)*amat(j,k)/amat(j+n,0)
              tmp = tmp*m1_power(abs(n+k))/rzall(n+j)
              loc(j,k) = loc(j,k)+tmp*yl(n+j)*tmpj
            endif 
          enddo 
        enddo 
      enddo 

C
      end 
c C
c C
c C
c C
C
      subroutine zl3dlocloc_zshift(zshift,nterm0,loc0,rscale0,nterm,
     1 loc,rscale,nmax,amat)
      implicit real *8 (a-h,o-z)
      complex *16 zshift
      complex *16 loc0(0:nterm0,-nterm0:nterm0)
      complex *16 loc(0:nterm,-nterm:nterm)
      complex *16 rz,ctheta,eiphi
      complex *16 ima,tmp
      dimension amat(0:nmax,-nmax:nmax)
      
      integer, allocatable :: m1_power(:)
      real *8, allocatable :: rscale0s(:),rscales(:)
      complex *16, allocatable :: rzall(:)
      real *8 , allocatable :: yl(:)

      allocate(m1_power(0:nterm+nterm0))
      allocate(rscale0s(0:nterm0+1),rscales(0:nterm0+1))
      allocate(rzall(0:nterm0+nterm),yl(0:nterm0))

      ima = (0.0d0,1.0d0)
      done = 1 
      pi = atan(done)*4 
      pi4 = pi*4 
      pi4_sqrt = sqrt(pi4)
      eiphi = 1 

      m1_power(0) =1 
      do m=1,nterm+nterm0 
        m1_power(m)=(-1)*m1_power(m-1)
      enddo 

      rscale0s(0) = 1 
      rscales(0) = 1 
      do n=1,nterm0+1
        rscale0s(n) = rscale0s(n-1)*rscale0 
      enddo 

      do n=1,nterm0+1
        rscales(n) = rscales(n-1)*rscale 
      enddo 

      if (abs(sqrt(zshift**2)-zshift)>1E-16) then 
        rz = -zshift 
        do n=0,nterm0 
          yl(n) = dsqrt(2*n+1.0d0) * (-1)**n 
        enddo 
      else 
        rz = zshift
        do n=0,nterm0 
          yl(n) = dsqrt(2*n+1.0d0) 
        enddo 
      endif 


      rz = rz/rscale
      rzall(0) = 1 
      do m=1,nterm+nterm0 
        rzall(m) = rzall(m-1)*rz
      enddo 


      do j=0,nterm
        tmpj =  pi4_sqrt/(2*j+1.0d0)
        do k=-j,j 
          do n=j,nterm0
              tmp = loc0(n,k)*(amat(n-j,0)*amat(j,k)/amat(n,k))
              tmp = tmp*m1_power(abs(n+j))*rzall(n-j)*yl(n-j)                                                   
              tmp = tmp*(2*n+1.0d0)/(2*(n-j)+1.0d0)
              tmp = tmp*rscales(n+1)/rscale0s(n+1)
              loc(j,k) = loc(j,k)+tmp*tmpj 
          enddo 
        enddo 
      enddo 
      end 

C
C
C
C
      subroutine get_amat(nmax,amat)
      implicit real *8 (a-h,o-z)
      dimension amat(0:nmax,-nmax:nmax)
      
      real *8, allocatable :: facts(:)
      
      allocate(facts(0:2*nmax))
      done = 1
      pi = atan(done)*4
      pi4 = pi*4
      pi4_inv = 1/pi4  
      facts(0) = 1

      do i=1,2*nmax
        facts(i) =  dsqrt(i+0.0d0)*facts(i-1)
      enddo 

      amat = 0 
      do n=0,nmax 
        const = ((-1)**n)*dsqrt((2*n+1)*pi4_inv)
        do m=-n,n 
          amat(n,m) = const/facts(n-m)/facts(n+m)
        enddo 
      enddo 
      end 