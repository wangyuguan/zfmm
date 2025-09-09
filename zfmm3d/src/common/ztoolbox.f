C-----------------------------------------------------------------------
      subroutine zreorderf(ndim,n,arr,arrsort,iarr) 
c
cc       this subroutine sorts the array arr and stores
c        it in arrsort using the sorting order defined by
c        iarr
c
c        arrsort(j,i) = arr(j,iarr(i)), j =1,2,\ldots ndim
c                                       i=1,2,\ldots n
c
      implicit real *8 (a-h,o-z)
      complex *16 arr(ndim,n),arrsort(ndim,n)
      dimension iarr(n)

      do i=1,n
        do idim=1,ndim
          arrsort(idim,i) = arr(idim,iarr(i))
        enddo
      enddo

      end 
c
c
c
c
c
      subroutine zrescale(n,a,r)
      implicit none
      integer i,n
      real *8 r
      complex *16 a(n)

C$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,n
        a(i) = a(i)*r
      enddo
C$OMP END PARALLEL DO

      


      return
      end
C
C
C
C
C
C-----------------------------------------------------------------------
      subroutine zreorderi(ndim,n,arr,arrsort,iarr) 
c
cc       this subroutine sorts the array arr and stores
c        it in arrsort using the inverse of the
c        sorting order defined by
c        iarr.
c
c        Note that this subroutine is the inverse of 
c        dreorderf
c
c        arrsort(j,iarr(i)) = arr(j,i), j =1,2,\ldots ndim
c                                       i=1,2,\ldots n
c
      implicit real *8 (a-h,o-z)
      complex *16 arr(ndim,n),arrsort(ndim,n)
      dimension iarr(n)

      do i=1,n
        do idim=1,ndim
          arrsort(idim,iarr(i)) = arr(idim,i)
        enddo
      enddo

      end 
c
c
c
c
c
C-----------------------------------------------------------------------
      subroutine zdist(ndim,x,y,rz)
      implicit real *8 (a-h,o-z)
      complex *16 x(ndim),y(ndim),rz 
      rz = 0 
      do i=1,ndim 
        rz = rz+(x(i)-y(i))**2
      enddo 
      rz = sqrt(rz)
      end
C
C
C
C
C
C-----------------------------------------------------------------------
      subroutine zcart2polar(xy,x0y0,rz,ez)
      implicit real *8 (a-h,o-z)
      complex *16 xy(2),x0y0(2),z1,z2 
      complex *16 rz,ez,ima 
      ima = (0.0d0,1.0d0)
      z1 = xy(1)-x0y0(1)
      z2 = xy(2)-x0y0(2)
      rz = sqrt(z1**2+z2**2)
      ez = z1/rz+ima*(z2/rz) 
      end
c
c
c
c
c
      subroutine zsqrt(z2,z)
      implicit real *8 (a-h,o-z)
      complex *16 z2,z,ima 
      
      pi = atan(1.0d0)*4 
      ima = (0.0d0,1.0d0)
      z2_real = real(z2)
      z2_imag = aimag(z2)
      r = dsqrt(abs(z2))
      phase = atan2(z2_imag,z2_real)
      if (phase.eq.-pi) phase = pi 
      phase = phase*0.5d0 
      z = r*cos(phase)+ima*r*sin(phase)

      end 
c
c
c
c
c
c
C-----------------------------------------------------------------------
      subroutine zcart2sph(xyz,x0y0z0,rz,ctheta,stheta,cphi,sphi,
     1 eiphi)
      implicit real *8 (a-h,o-z)
      complex *16 xyz(3),x0y0z0(3),z1,z2,z3 
      complex *16 rz,ima,ctheta,stheta,cphi,sphi,eiphi,zrxy,zrxy2
      complex *16 rz2,stheta2
      ima = (0.0d0,1.0d0)
      z1 = xyz(1)-x0y0z0(1)
      z2 = xyz(2)-x0y0z0(2)
      z3 = xyz(3)-x0y0z0(3)
      rz2 = z1**2+z2**2+z3**2
      call zsqrt(rz2,rz)
      ctheta = z3/rz 
      stheta2 = z1**2+z2**2
      call zsqrt(stheta2,stheta)
      stheta = stheta/rz  
      zrxy2 = z1**2+z2**2
      call zsqrt(zrxy2,zrxy)
      if (abs(zrxy).le.1E-16) then
        cphi = 1 
        sphi = 0 
        eiphi = 1 
      else 
        cphi = z1/zrxy 
        sphi = z2/zrxy 
        eiphi = cphi+ima*sphi
      endif 

      end
c
c
c
c
c
C-----------------------------------------------------------------------
      subroutine getl2err(ndim,n,u0,u,err)
      implicit real *8 (a-h,o-z)
      complex *16 u0(ndim,n),u(ndim,n)
      err = 0
      u0_norm  = 0 
      do j=1,n
          do i=1,ndim
              err = err+abs(u0(i,j)-u(i,j))**2 
              u0_norm = u0_norm + abs(u0(i,j))**2 
          enddo  
      enddo 
      err = sqrt(err)/sqrt(u0_norm)
      end 
C
C
C
C
C
      subroutine fftshift(n,y,yshift)
      implicit real *8 (a-h,o-z)
      complex *16 y(n)
      complex *16 yshift(n)

      m = (n+1)/2
      yshift = 0
      do i=1,n-m
        yshift(i) = y(m+i)
      enddo 

      do i=1,m
        yshift(n-m+i) = y(i) 
      enddo 

      end 
c
c
c
c
c
      subroutine ifftshift(n,y,yshift)
      implicit real *8 (a-h,o-z)
      complex *16 y(n)
      complex *16 yshift(n)

      m = (n+1)/2
      yshift = 0
      do i=1,m
        yshift(i) = y(n-m+i) 
      enddo 

      do i=1,m
        yshift(m+i) = y(i) 
      enddo 

      end 
c
c
c
c
c
c
      subroutine phi_fun(n,x,phi)
      implicit real *8 (a-h,o-z)
      dimension x(n)
      dimension phi(n)
      done = 1 
      pi = atan(done)*4
      pi_sqrt_inv = done/sqrt(pi)
      do i=1,n 
        phi(i) = x(i)*erfc(x(i))-exp(-x(i)**2)/pi_sqrt_inv
      enddo 
      end 
c
c
c
c
c
      subroutine path_cmpl(n,x,z,a,b,ell)
      implicit real *8 (a-h,o-z)
      dimension x(n),z(n)
      dimension x1(n),x2(n),phi1(n),phi2(n)


      phi1 = 0 
      phi2 = 0 
      x1 = 0
      x2 = 0 

      do i=1,n
        x1(i) = b*(x(i)+ell) 
        x2(i) = -b*(x(i)-ell) 
      enddo 

      call phi_fun(n,x1,phi1)
      call phi_fun(n,x2,phi2)

      do i=1,n
        z(i) = (phi1(i)-phi2(i))*a
      enddo 

      end 
c
c
c
c
c
      subroutine myzylgndr2sf(nmax, ctheta, stheta, y, d, rat1, rat2)
      implicit real *8 (a-h,o-z)
      complex *16 ctheta, stheta
      complex *16 z, y(0:nmax,0:nmax), d(0:nmax,0:nmax), u,u2
      dimension rat1(0:nmax,0:nmax), rat2(0:nmax,0:nmax)
c      u=-sqrt(1-z*z)
      z = ctheta
      u2 = 1-z**2
      if (abs(u2).ge.1E-12) then 
        call zsqrt(u2,u)
        u = -u 
      else 
        u = -stheta
      endif 
c      u = sqrt(u2)
      
      y(0,0)=1
      d(0,0)=0
c
c       ... first, evaluate standard Legendre polynomials
c
      m=0
      if (m.lt.nmax)  y(m+1,m)=z*y(m,m)*rat1(m+1,m)
      if (m.lt.nmax)  d(m+1,m)=(z*d(m,m)+y(m,m))*rat1(m+1,m)
      do n=m+2, nmax
        y(n,m)=rat1(n,m)*z*y(n-1,m)-rat2(n,m)*y(n-2,m)
        d(n,m)=rat1(n,m)*(z*d(n-1,m)+y(n-1,m))-rat2(n,m)*d(n-2,m)
      enddo
c
c       ... then, evaluate scaled associated Legendre functions
c
      do m=1, nmax
c
        if (m.eq.1)  y(m,m)=y(m-1,m-1)*(-1)*rat1(m,m)
        if (m.gt.1)  y(m,m)=y(m-1,m-1)*u*rat1(m,m)
        if (m.gt.0)  d(m,m)=y(m,m)*(-m)*z
c
        if (m.lt.nmax)  y(m+1,m)=z*y(m,m)*rat1(m+1,m)
        if (m.lt.nmax)  
     $      d(m+1,m)=(z*d(m,m)+(1-z**2)*y(m,m))*rat1(m+1,m)
        do n=m+2, nmax
            y(n,m)=rat1(n,m)*z*y(n-1,m)-rat2(n,m)*y(n-2,m)
            d(n,m)=rat1(n,m)*(z*d(n-1,m)+(1-z**2)*y(n-1,m))-
     $         rat2(n,m)*d(n-2,m)
         enddo
      enddo
c
      do n=0, nmax
	      do m=0, n
          y(n,m)=y(n,m)*sqrt(2*n+1.0d0)
          d(n,m)=d(n,m)*sqrt(2*n+1.0d0)
         enddo
      enddo
      return
      end




