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
C
C
C
C
C
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
C
C
C
C
C
c C-----------------------------------------------------------------------
c         subroutine eval_mp(z,zc,kwave,rscale,nterm,mp,u)
c         implicit real *8 (a-h,o-z)
c         complex *16 z(2),zc(2)
c         complex *16 kwave 
c         complex *16 mp(-nterm:nterm)
c         complex *16 fhs(0:nterm),fhder(0:1)
c         complex *16 zr,ei 
c         complex *16 u

c C         u = 0 
c         ifder0 = 0
c         call zcart2polar(z,zc,zr,ei)
c         call h2dall(nterm,zr*kwave,rscale,fhs(0:nterm),
c      1              ifder0,fhder)
        
c         do m=0,nterm 
c             u = u+mp(m)*fhs(m)*(ei**m)
c         enddo 

c         do m=-nterm,-1 
c             u = u+mp(m)*((-1)**m)*fhs(-m)*(ei**m)
c         enddo 
c         end
c C
c C
c C
c C
c C
c C-----------------------------------------------------------------------
c         subroutine eval_loc(z,zc,kwave,rscale,nterm,loc,u)
c         implicit real *8 (a-h,o-z)
c         complex *16 z(2),zc(2)
c         complex *16 kwave 
c         complex *16 loc(-nterm:nterm)
c         complex *16 fjs(0:nterm),fjder(0:1)
c         complex *16 zr,ei
c         complex *16 u
        
c C         u = 0 
c         ifder0 = 0
c         call zcart2polar(z,zc,zr,ei)
c         call jbessel2d(nterm,zr*kwave,rscale,fjs(0:nterm),
c      1              ifder0,fjder)
c C         call prin2('fjs=*',fjs(0),nterm+1)
        
c         do m=0,nterm 
c             u = u+loc(m)*fjs(m)*(ei**m)
c         enddo 

c         do m=-nterm,-1 
c             u = u+loc(m)*((-1)**m)*fjs(-m)*(ei**m)
c         enddo 
c         end
C-----------------------------------------------------------------------
C
C
C
C
C
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
c      call prin2('norm=*',u0_norm,1)
      end 
C
C
C
C
C








