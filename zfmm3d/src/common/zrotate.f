      subroutine z3dzrot(nterm,lmp,lmpn,cbeta,sbeta)
      implicit real *8 (a-h,o-z)
      complex *16 eibeta,cbeta,sbeta
      complex *16 lmp(0:nterm,-nterm:nterm)
      complex *16 lmpn(0:nterm,-nterm:nterm)
      complex *16 ima 
      complex *16, allocatable :: eibetas(:)

      allocate(eibetas(-nterm:nterm))
      ima = (0.0d0,1.0d0)
      eibeta = cbeta+ima*sbeta 

      eibetas(0) = 1 
      do m=1,nterm 
        eibetas(m) = eibetas(m-1)*eibeta 
      enddo 

      do m=-nterm,-1 
        eibetas(m) = 1/eibetas(-m)
      enddo 

      do n=0,nterm
        do m=-n,n
          lmpn(n,m) = lmp(n,m)*eibetas(m)
        enddo 
      enddo 
      end 
c
c
c
c
      subroutine z3dzrotini(centern,center,xc,cbeta,sbeta)
      implicit real *8 (a-h,o-z)
      complex *16 centern(3),center(3),xc(3)
      complex *16 cbeta,sbeta 
      complex *16 rxy,rxy2 

      xc(1) = centern(1)-center(1)
      xc(2) = centern(2)-center(2)
      xc(3) = centern(3)-center(3)

      rxy2 = xc(1)**2+xc(2)**2
      call zsqrt(rxy2,rxy)
      cbeta = -xc(1)/rxy 
      sbeta = -xc(2)/rxy 

      xc(1) = -rxy 
      xc(2) = 0 
      end
C
C
C
C
C
      subroutine z3dyrotini(xc,xcn,calpha,salpha)
      implicit real *8 (a-h,o-z)
      complex *16 xc(3)
      complex *16 xcn(3)
      complex *16 calpha,salpha 
      complex *16 rxz,rxz2 
      rxz2 = xc(1)**2+xc(3)**2
      call zsqrt(rxz2,rxz)    
c      rxz = sqrt(xc(1)**2+xc(3)**2)  
      xcn(1) = 0 
      xcn(2) = 0 
      xcn(3) = rxz 
      calpha = xc(3)/rxz 
      salpha = xc(1)/rxz 

c      if (1.eq.0) then 
c        if (abs(salpha).eq.0) then 
c          if (real(calpha).gt.0) then 
c            calpha = 1 
c          else 
c            calpha = -1 
c          endif 
c        endif 
c
c        if (abs(calpha).eq.0) then 
c          if (real(salpha).gt.0) then 
c            salpha = 1 
c          else 
c            salpha = -1 
c          endif 
c        endif 
c      endif 

      end
c
c
c
c
c
      subroutine z3dyrot(nterm,lmp,lmpn,calpha,salpha)
      implicit real *8 (a-h,o-z)
      complex *16 lmp(0:nterm,-nterm:nterm)
      complex *16 lmpn(0:nterm,-nterm:nterm)
      complex *16 calpha,salpha 
      complex *16 x1,x2,x3,zrxy,zrxy2
      complex *16 ctheta,stheta,stheta2
c      complex *16 ctheta_quad,stheta_quad
      complex *16 cphi,sphi,eiphi
      complex *16 ynm,dynm
      complex *16 df_dth,df_dph,dth_dx,dth_dz,dph_dx,df_dx,df_dz 
      complex *16 ima,z1,z2,z3,z4,z5,z6  
      complex *16 ztmp1,ztmp2,ztmp3,ztmp4

      complex *16, allocatable :: ffti(:), ffto(:), ffti_shift(:)
      real *8, allocatable :: wsave(:)
      complex *16, allocatable :: yls(:,:), dyls(:,:), eiphis(:)
      complex *16, allocatable :: fdat(:,:), dfdat(:,:)
      real *8, allocatable :: rat1(:,:), rat2(:,:)  
      integer, allocatable :: m1_power(:)


      allocate(ffti(2*nterm+1+100),ffti_shift(2*nterm+1+100))
      allocate(ffto(2*nterm+1+100)) 
c      allocate(ffti(-nterm:nterm),ffto(-nterm:nterm))
      allocate(wsave(8*(3*nterm+1)+100) )
      allocate(yls(0:nterm,0:nterm) )
      allocate(dyls(0:nterm,0:nterm) )
      allocate(eiphis(-nterm:nterm) )
      allocate(fdat(0:nterm,-nterm:nterm) )
      allocate(dfdat(0:nterm,-nterm:nterm) )
      allocate(rat1(0:nterm,0:nterm),rat2(0:nterm,0:nterm) ) 
      allocate(m1_power(0:nterm))

      ima = (0.0d0,1.0d0)
      done = 1 
      fdat = 0
      dfdat = 0 
      call lmp_ini(nterm,lmpn)
      pi = atan(done)*4
      pi4 = 4*pi 
      pi4_sqrt = dsqrt(pi4)
      pi4_sqrt_inv = 1.0d0/dsqrt(pi4)
      ffti = 0 
      ffto = 0 
      ffti_shift = 0

      m1_power(0) = 1
      do m=1,nterm
        m1_power(m) = m1_power(m-1)*(-1)
      enddo 
      rat1 = 0 
      rat2 = 0
      call ylgndrini(nterm, rat1, rat2)

c      call cpu_time(t1)
      do j=-nterm,nterm 
        ph = j*2*pi/(2*nterm+1)
        x1 = cos(ph)*calpha 
        x2 = sin(ph)

        zrxy2 = x1**2+x2**2
        call zsqrt(zrxy2,zrxy)
        ctheta = -cos(ph)*salpha 
        call zsqrt(1-(cos(ph)**2)*(1-calpha**2),stheta)


        if (ph.eq.0) then 
          if (abs(stheta+calpha).le.abs(stheta-calpha)) then 
            stheta = -calpha 
          else 
            stheta = calpha 
          endif 
        endif 
        
        if (abs(zrxy).ge.1E-16) then 
          cphi = x1/zrxy 
          sphi = x2/zrxy 
          eiphi = cphi+ima*sphi
        else 
          cphi = 1 
          sphi = 0 
          eiphi = 1
        endif 

        eiphis(0) = 1 
        do m=1,nterm
          eiphis(m) = eiphis(m-1)*eiphi
        enddo 

        do m=-nterm,-1 
          eiphis(m) = 1/eiphis(-m)
        enddo  

        call myzylgndr2sf(nterm,ctheta,stheta,yls,dyls,rat1,rat2)

        z1 = pi4_sqrt_inv*stheta
        z2 = (stheta*calpha-ctheta*cphi*salpha)*pi4_sqrt_inv
        z3 = sphi*salpha*pi4_sqrt_inv*ima

C
c       evaluate at the original coordinates 
c
c
        do n=0,nterm 
          fdat(n,j) = fdat(n,j)+lmp(n,0)*yls(n,0)*pi4_sqrt_inv
          dfdat(n,j) = dfdat(n,j)-lmp(n,0)*dyls(n,0)*z2*stheta

          do m=1,n 

            z4 = lmp(n,m)*eiphis(m)
            z5 = lmp(n,-m)*eiphis(-m)
            z6 = z4+z5 

            fdat(n,j) = fdat(n,j)+z1*z6*m1_power(m)*yls(n,m)
            dfdat(n,j) = dfdat(n,j)+m1_power(m)*(-dyls(n,m)*z2*z6
     1                      +yls(n,m)*z3*m*z4-yls(n,m)*z3*m*z5)
          enddo 
        enddo  
      enddo
c      call cpu_time(t2)
c      call prin2('eval spends *', t2-t1, 1)

      ctheta = 0 
      eiphi = 1
      factor = pi4_sqrt/(2*nterm+1)
      call zylgndr2sf(nterm,ctheta,yls,dyls,rat1,rat2)
      call zffti(2*nterm+1,wsave)
      do n=0,nterm 
        do j=-nterm,nterm 
          ffti(j+nterm+1) = fdat(n,j)+dfdat(n,j)
        enddo 
        
        call ifftshift(2*nterm+1,ffti,ffti_shift)
        call zfftf(2*nterm+1,ffti_shift,wsave)
        call fftshift(2*nterm+1,ffti_shift,ffto)

        lmpn(n,0)=ffto(nterm+1)/(yls(n,0)-dyls(n,0))*factor
        do m=1,n
          z1 = m1_power(m)*(yls(n,m)-dyls(n,m))/factor
          lmpn(n,m) = ffto(m+nterm+1)/z1 
          lmpn(n,-m) = ffto(-m+nterm+1)/z1 
        enddo 
      enddo 

      end
c
c
c
c
c
      subroutine lmp_ini(nterm,lmp)
      implicit real *8 (a-h,o-z)
      complex *16 lmp(0:nterm,-nterm:nterm)

      do n=0,nterm 
        do m=-n,n 
          lmp(n,m) = 0
        enddo 
      enddo 
      
      end 
c
c
c
c
c c
c         subroutine z3dyrot_slow(nterm,lmp,lmpn,calpha,salpha)
c         implicit real *8 (a-h,o-z)
c         complex *16 lmp(0:nterm,-nterm:nterm)
c         complex *16 lmpn(0:nterm,-nterm:nterm)
c         complex *16 calpha,salpha 
c         complex *16 x1,x2,x3,zrxy 
c         complex *16 ctheta,stheta
c         complex *16 cphi,sphi,eiphi
c         complex *16 ynm,dynm
c         complex *16 df_dth,df_dph 
c         complex *16 dth_dx,dth_dz 
c         complex *16 dph_dx
c         complex *16 df_dx,df_dz 
c         complex *16 ima 
c         complex *16 ffti(2*nterm+1+10),ffti_shift(2*nterm+1+10)
c         complex *16 ffto(2*nterm+1+10)
c         dimension wsave(4*(3*nterm+1)+1000)
c         complex *16 yls(0:nterm,0:nterm)
c         complex *16 dyls(0:nterm,0:nterm)
c         complex *16 eiphis(-nterm:nterm)
c         complex *16 fdat(0:nterm,-nterm:nterm)
c         complex *16 dfdat(0:nterm,-nterm:nterm)
c         dimension rat1(0:nterm,0:nterm),rat2(0:nterm,0:nterm)  


c         ima = (0.0d0,1.0d0)
c         fdat = 0
c         dfdat = 0 
c         lmpn = 0
c         fdat = 0 
c         dfdat = 0 
c         pi = atan(1.0d0)*4
c         pi4 = 4*pi 
c         pi4_sqrt = sqrt(pi4)
c         pi4_sqrt_inv = 1.0d0/sqrt(pi4)

c         call ylgndrini(nterm, rat1, rat2)
c         call cpu_time(t1)
c         do n=0,nterm 
c           do j=-n,n 
c             ph = j*2*pi/(2*n+1)
c             x1 = cos(ph)*calpha 
c             x2 = sin(ph)

c             zrxy = sqrt(x1**2+x2**2)
c             ctheta = -cos(ph)*salpha 
c             stheta = sqrt(1-ctheta**2)

c             if (abs(zrxy).ge.1E-16) then 
c               cphi = x1/zrxy 
c               sphi = x2/zrxy 
c               eiphi = cphi+ima*sphi
c             else 
c               cphi = 1 
c               sphi = 0 
c               eiphi = 1
c             endif 

c             eiphis(0) = 1 
c             do m=1,n
c               eiphis(m) = eiphis(m-1)*eiphi
c             enddo 

c             do m=-n,-1 
c               eiphis(m) = 1/eiphis(-m)
c             enddo  

c             call zylgndr2sf(n,ctheta,yls(0:n,0:n),dyls(0:n,0:n),
c      1           rat1(0:n,0:n),rat2(0:n,0:n))
     
c             fdat(n,j) = fdat(n,j)+lmp(n,0)*yls(n,0)*pi4_sqrt_inv
c             do m=1,n 
c               ynm =  ((-1)**m)*yls(n,abs(m))*pi4_sqrt_inv*stheta
c               fdat(n,j) = fdat(n,j)+lmp(n,m)*ynm*eiphis(m)+
c      1                     lmp(n,-m)*ynm*eiphis(-m)
c             enddo 


c             do m=-n,n
c               ynm = ((-1)**m)*yls(n,abs(m))*eiphis(m)*pi4_sqrt_inv
c               dynm = -((-1)**m)*dyls(n,abs(m))*eiphis(m)*pi4_sqrt_inv
c               if (m.eq.0) then 
c                 df_dth = dynm*stheta                                  
c                 df_dph = 0
c               else 
c                 df_dth = dynm
c                 df_dph = ynm*ima*m
c               endif 

c c              df_dx = df_dth*ctheta*cphi-df_dph*sphi
c c              df_dz = -df_dth*stheta
              
c c              df_dth = -(df_dx*salpha+df_dz*calpha)

c c              df_dth = df_dth*stheta*calpha - 
c c     1            (df_dth*ctheta*cphi-df_dph*sphi)*salpha

c               df_dth = df_dth*(stheta*calpha-ctheta*cphi*salpha)+
c      1           df_dph*sphi*salpha
              
c               dfdat(n,j) = dfdat(n,j)+lmp(n,m)*df_dth
c             enddo 
c           enddo 
c         enddo 
c         call cpu_time(t2)
c         call prin2('eval spends *', t2-t1, 1)

c         ctheta = 0 
c         eiphi = 1



c         call zylgndr2sf(nterm,ctheta,yls(0:nterm,0:nterm),
c      1          dyls(0:nterm,0:nterm),rat1,rat2)

c         do n=0,nterm 
c           nn = 2*n+1
c           call zffti(nn,wsave)

c           do m=-n,n 
c             ffti(m+n+1) = fdat(n,m)+dfdat(n,m)
c           enddo 
          
c           call ifftshift(nn,ffti,ffti_shift)
c           call zfftf(nn,ffti_shift,wsave)
c           call fftshift(nn,ffti_shift,ffto)

c           do m=-n,n 
c             ffto(m+n+1) = ffto(m+n+1)/(2*n+1)
c             lmpn(n,m) = ffto(m+n+1)/(((-1)**m)*yls(n,abs(m))
c      1          -((-1)**m)*dyls(n,abs(m)))*pi4_sqrt
c           enddo 
c         enddo 

c         end
