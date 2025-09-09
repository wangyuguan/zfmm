cccc    shift from mp1 tp loc2 
      subroutine zh3dmploc(zc1,zc2,nterm1,nterm2,rscale1,
     1   rscale2,width1,width2,kwave,mp1,loc2)
      implicit real *8 (a-h,o-z)
      complex *16 zc1(3),zc2(3)
      complex *16 kwave 
      complex *16 mp1(0:nterm1,-nterm1:nterm1)
      complex *16 loc2(0:nterm2,-nterm2:nterm2)
      complex *16 xc1(3),xc2(3)
      complex *16 zshift
      complex *16 cbeta,sbeta
      complex *16 calpha,salpha 


      complex *16, allocatable :: mp1z(:,:)
      complex *16, allocatable :: mp1zy(:,:)
      complex *16, allocatable :: loc2zy(:,:)
      complex *16, allocatable :: loc2z(:,:)
      complex *16, allocatable :: loc2n(:,:)


      allocate(mp1z(0:nterm1,-nterm1:nterm1))
      allocate(mp1zy(0:nterm1,-nterm1:nterm1))
      allocate(loc2zy(0:nterm2,-nterm2:nterm2))
      allocate(loc2z(0:nterm2,-nterm2:nterm2))
      allocate(loc2n(0:nterm2,-nterm2:nterm2))

      e2 = abs(zc2(2)-zc1(2))
      if (abs(e2).le.1E-16) then 
        xc1(1:3) = zc2(1:3)-zc1(1:3)
        cbeta = 1 
        sbeta = 0 
        do n=0,nterm1 
          do m=-n,n 
            mp1z(n,m) = mp1(n,m)
          enddo 
        enddo
      else 
        call z3dzrotini(zc2,zc1,xc1,cbeta,sbeta)
        call z3dzrot(nterm1,mp1,mp1z,cbeta,sbeta)
      endif 

      e1 = abs(xc1(1))
      if (abs(e1).le.1E-16) then  
        xc2(1:3) = xc1(1:3)
        calpha = 1 
        salpha = 0  
        do n=0,nterm1 
          do m=-n,n 
            mp1zy(n,m) = mp1z(n,m)
          enddo 
        enddo
      else 
        call z3dyrotini(xc1,xc2,calpha,salpha)
        call z3dyrot(nterm1,mp1z,mp1zy,calpha,salpha)
      endif

      zshift = xc2(3)
      nquad = nterm2*2.5
      r = width1/2*sqrt(3.0d0)

      loc2zy = 0
      call zh3dzshiftmploc(r,nquad,zshift,nterm1,nterm2,
     1  rscale1,rscale2,kwave,mp1zy,loc2zy)

      if (abs(e1).le.1E-16) then 
        do n=0,nterm2
          do m=-n,n 
            loc2z(n,m) = loc2zy(n,m)
          enddo 
        enddo
      else 
        salpha = -salpha
        call z3dyrot(nterm2,loc2zy,loc2z,calpha,salpha)
      endif 

      if (abs(e2).le.1E-16) then 
        do n=0,nterm2
          do m=-n,n 
            loc2n(n,m) = loc2z(n,m)
          enddo 
        enddo
      else 
        sbeta = -sbeta 
        call z3dzrot(nterm2,loc2z,loc2n,cbeta,sbeta)
      endif 
      do n=0,nterm2 
        do m=-n,n 
          loc2(n,m) = loc2(n,m)+loc2n(n,m)
        enddo 
      enddo 
      end
C
C
C
C
C

      subroutine zh3dmpmp(zc1,zc2,nterm1,nterm2,rscale1,
     1   rscale2,width1,width2,kwave,mp1,mp2)
cccc    shift from mp1 tp mp2 
      implicit real *8 (a-h,o-z)
      complex *16 zc1(3),zc2(3)
      complex *16 kwave 
      complex *16 mp1(0:nterm1,-nterm1:nterm1)
      complex *16 mp2(0:nterm2,-nterm2:nterm2)
      complex *16 xc1(3),xc2(3)
      complex *16 zshift
      complex *16 cbeta,sbeta
      complex *16 calpha,salpha 

      complex *16, allocatable :: mp1z(:,:)
      complex *16, allocatable :: mp1zy(:,:)
      complex *16, allocatable :: mp2zy(:,:)
      complex *16, allocatable :: mp2z(:,:)
      complex *16, allocatable :: mp2n(:,:)


      allocate(mp1z(0:nterm1,-nterm1:nterm1))
      allocate(mp1zy(0:nterm1,-nterm1:nterm1))
      allocate(mp2zy(0:nterm2,-nterm2:nterm2))
      allocate(mp2z(0:nterm2,-nterm2:nterm2))
      allocate(mp2n(0:nterm2,-nterm2:nterm2))

cccc    shift from zc1 to zc2 

      e2 = abs(zc2(2)-zc1(2))
      if (abs(e2).le.1E-16) then 
        xc1(1:3) = zc2(1:3)-zc1(1:3)
        cbeta = 1 
        sbeta = 0 
        do n=0,nterm1 
          do m=-n,n 
              mp1z(n,m) = mp1(n,m)
          enddo 
        enddo
      else 
        call z3dzrotini(zc2,zc1,xc1,cbeta,sbeta)
        call z3dzrot(nterm1,mp1,mp1z,cbeta,sbeta)
      endif 

      e1 = abs(xc1(1))
      if (abs(e1).le.1E-16) then 
        xc2(1:3) = xc1(1:3)
        calpha = 1 
        salpha = 0  
        do n=0,nterm1 
          do m=-n,n 
            mp1zy(n,m) = mp1z(n,m)
          enddo 
        enddo
      else  
        call z3dyrotini(xc1,xc2,calpha,salpha)
        call z3dyrot(nterm1,mp1z,mp1zy,calpha,salpha)
      endif 

      zshift = xc2(3)
      nquad = nterm2*2.5
      r = width1*sqrt(3.0d0)

      call zh3dzshiftmpmp(r,nquad,zshift,nterm1,nterm2,rscale1,
     1  rscale2,kwave,mp1zy,mp2zy)
      if (abs(e1).le.1E-16) then 
        do n=0,nterm2
          do m=-n,n 
            mp2z(n,m) = mp2zy(n,m)
          enddo 
        enddo
      else 
        salpha = -salpha
        call z3dyrot(nterm2,mp2zy,mp2z,calpha,salpha)
      endif 

      if (abs(e2).le.1E-16) then 
        do n=0,nterm2
          do m=-n,n 
            mp2n(n,m) = mp2z(n,m)
          enddo 
        enddo
      else 
        sbeta = -sbeta 
        call z3dzrot(nterm2,mp2z,mp2n,cbeta,sbeta)
      endif 

      do n=0,nterm2 
        do m=-n,n 
          mp2(n,m) = mp2(n,m)+mp2n(n,m)
        enddo 
      enddo 

      end
C
C
C
C
C
      subroutine zh3dlocloc(zc1,zc2,nterm1,nterm2,rscale1,
     1   rscale2,width1,width2,kwave,loc1,loc2)
cccc    shift from loc1 to loc2 
      implicit real *8 (a-h,o-z)
      complex *16 zc1(3),zc2(3)
      complex *16 kwave 
      complex *16 loc1(0:nterm1,-nterm1:nterm1)
      complex *16 loc2(0:nterm2,-nterm2:nterm2)
      complex *16 xc1(3),xc2(3)
      complex *16 zshift
      complex *16 cbeta,sbeta
      complex *16 calpha,salpha 

      complex *16, allocatable :: loc1z(:,:)
      complex *16, allocatable :: loc1zy(:,:)
      complex *16, allocatable :: loc2zy(:,:)
      complex *16, allocatable :: loc2z(:,:)
      complex *16, allocatable :: loc2n(:,:)

      allocate(loc1z(0:nterm1,-nterm1:nterm1))
      allocate(loc1zy(0:nterm1,-nterm1:nterm1))
      allocate(loc2zy(0:nterm2,-nterm2:nterm2))
      allocate(loc2z(0:nterm2,-nterm2:nterm2))
      allocate(loc2n(0:nterm2,-nterm2:nterm2))

      loc1z = 0 
      loc1zy = 0 
      loc2zy = 0 
      loc2z = 0 
      loc2n = 0

      e2 = abs(zc2(2)-zc1(2))
      if (abs(e2).le.1E-16) then 
        xc1(1:3) = zc2(1:3)-zc1(1:3)
        cbeta = 1 
        sbeta = 0 
        do n=0,nterm1 
          do m=-n,n 
            loc1z(n,m) = loc1(n,m)
          enddo 
        enddo
      else 
        call z3dzrotini(zc2,zc1,xc1,cbeta,sbeta)
        call z3dzrot(nterm1,loc1,loc1z,cbeta,sbeta)
      endif 


      e1 = abs(xc1(1))
      if (abs(e1).le.1E-16) then  
        xc2(1:3) = xc1(1:3)
        calpha = 1 
        salpha = 0  
        do n=0,nterm1 
          do m=-n,n 
            loc1zy(n,m) = loc1z(n,m)
          enddo 
        enddo
      else 
        call z3dyrotini(xc1,xc2,calpha,salpha)
        call z3dyrot(nterm1,loc1z,loc1zy,calpha,salpha)
      endif 


      zshift = xc2(3)
      nquad = nterm2*2.5
      r = width2*sqrt(3.0d0)

      call zh3dzshiftlocloc(r,nquad,zshift,nterm1,nterm2,rscale1,
     1  rscale2,kwave,loc1zy,loc2zy)


      if (abs(e1).le.1E-16) then 
        do n=0,nterm2
          do m=-n,n 
            loc2z(n,m) = loc2zy(n,m)
          enddo 
        enddo
      else 
        salpha = -salpha
        call z3dyrot(nterm2,loc2zy,loc2z,calpha,salpha)
      endif 

      if (abs(e2).le.1E-16) then 
        do n=0,nterm2
          do m=-n,n 
            loc2n(n,m) = loc2z(n,m)
          enddo 
        enddo
      else 
        sbeta = -sbeta 
        call z3dzrot(nterm2,loc2z,loc2n,cbeta,sbeta)
      endif 

      do n=0,nterm2 
        do m=-n,n 
          loc2(n,m) = loc2(n,m)+loc2n(n,m)
        enddo 
      enddo 

      end 
C
C
C
C
C
cccc    shift from mp1 tp loc2 along z-axis
      subroutine zh3dzshiftmploc(r,nquad,zshift,nterm1,nterm2,
     1      rscale1,rscale2,kwave,mp1,loc2)
      implicit real *8 (a-h,o-z)
      complex *16 zshift 
      complex *16 mp1(0:nterm1,-nterm1:nterm1)
      complex *16 loc2(0:nterm2,-nterm2:nterm2)
      complex *16 kwave 
      dimension u(2),v(2)
      
      complex *16 xi,zi,ri,zk
      complex *16 cthi,sthi,eiphi 
      complex *16 ri2,sthi2
      complex *16 df_r,df_th,df_x,df_y,df 
      complex *16 z1,z2,z3,cthi_ri,sthi_ri
      complex *16, allocatable:: ynm(:,:)
      complex *16, allocatable:: dynm(:,:)
      complex *16, allocatable:: fhs(:),fhder(:)
      complex *16, allocatable:: fjs(:),fjder(:),fjsum(:)
      complex *16, allocatable:: yls(:,:)
      complex *16, allocatable:: dyls(:,:)
      complex *16, allocatable :: fdat(:,:),dfdat(:,:)
      real *8, allocatable :: x(:),whts(:)
      real *8, allocatable :: rat1(:,:),rat2(:,:)
      integer, allocatable :: minus1(:)

      nterm = max(nterm1,nterm2)

      allocate(fhs(0:(nterm+100)))
      allocate(fhder(0:(nterm+100)))
      allocate(fjs(0:(nterm+100)))
      allocate(fjsum(0:(nterm+100)))
      allocate(fjder(0:(nterm+100)))
      allocate(yls(0:nterm,0:nterm))
      allocate(dyls(0:nterm,0:nterm))
      allocate(fdat(1:nquad,-nterm1:nterm1))
      allocate(dfdat(1:nquad,-nterm1:nterm1))
      allocate(x(nquad),whts(nquad))
      allocate(minus1(0:nterm))
      allocate(rat1(0:nterm1,0:nterm1))
      allocate(rat2(0:nterm1,0:nterm1))

      call ylgndrini(nterm1,rat1,rat2)


      minus1(0) = 1 
      do i=1,nterm 
        minus1(i) = minus1(i-1)*(-1)
      enddo 

      itype = 1
      ifder = 1 
      pi = atan(1.0d0)*4
      pi4 = pi*4
      pi4_sqrt_inv = 1/dsqrt(pi4)
      pi2_inv = 1/(2*pi)
      call legeexps(itype,nquad,x,u,v,whts)

      fdat = 0 
      dfdat = 0
      do i=1,nquad 
        yi = sqrt(1-x(i)**2)
        xi = r*yi 
        zi = zshift+r*x(i)
        ri2 = (r*yi)**2+zi**2
        call zsqrt(ri2,ri)
        cthi = zi/ri 
        sthi2 = 1-cthi**2
        call zsqrt(sthi2,sthi)
        cthi_ri = cthi/ri
        sthi_ri = sthi/ri 

        zk = ri*kwave 
        call zylgndr2sf(nterm1,cthi,yls(0:nterm1,0:nterm1),
     1     dyls(0:nterm1,0:nterm1),rat1,rat2)
        do n=0,nterm1 
          do m=0,n 
            yls(n,m) = yls(n,m)*minus1(m)*pi4_sqrt_inv
            dyls(n,m) = dyls(n,m)*minus1(m)*pi4_sqrt_inv
          enddo 
        enddo 

        call h3dall(nterm1,zk,rscale1,fhs,ifder,fhder)

        do n=0,nterm1
          fhder(n) = kwave*fhder(n) 
        enddo 

        do n=0,nterm1 
c          do m=-n,n
c            if (m.ne.0) then  
c              fdat(i,m) = fdat(i,m)+mp1(n,m)*fhs(n)*minus1(abs(m))
c     1                           *yls(n,abs(m))*sthi
c            else 
c              fdat(i,m) = fdat(i,m)+mp1(n,m)*fhs(n)*minus1(abs(m))
c     1                         *yls(n,abs(m))
c            endif 
c          enddo
          fdat(i,0) = fdat(i,0)+mp1(n,0)*fhs(n)*yls(n,0)
          df_r = fhder(n)*yls(n,0)
          df_th = -fhs(n)*dyls(n,0)*sthi
          df_x = df_r*sthi+df_th*cthi_ri
          df_y = df_r*cthi-df_th*sthi_ri 
          df = df_x*yi+df_y*x(i)
          dfdat(i,0) = dfdat(i,0)+mp1(n,0)*df

          do m=1,n 
            z1 = fhs(n)*yls(n,m)*sthi
            fdat(i,m) = fdat(i,m)+mp1(n,m)*z1 
            fdat(i,-m) = fdat(i,-m)+mp1(n,-m)*z1 

            df_r = fhder(n)*yls(n,m)*sthi
            df_th = -fhs(n)*dyls(n,m)
            df_x = df_r*sthi+df_th*cthi_ri 
            df_y = df_r*cthi-df_th*sthi_ri 
            df = df_x*yi+df_y*x(i)
            dfdat(i,m) = dfdat(i,m)+mp1(n,m)*df
            dfdat(i,-m) = dfdat(i,-m)+mp1(n,-m)*df
          enddo 
        enddo 
            
c        do n=0,nterm1
c          do m=-n,n 
c            if (m.ne.0) then 
c              df_r = fhder(n)*minus1(abs(m))
c     1                             *yls(n,abs(m))*sthi
c            else 
c              df_r = fhder(n)*minus1(abs(m))*yls(n,abs(m))
c            endif 
c
c            if (m.eq.0) then 
c              df_th = -fhs(n)*minus1(abs(m))
c     1                               *dyls(n,abs(m))*sthi
c            else 
c              df_th = -fhs(n)*minus1(abs(m))*dyls(n,abs(m))
c            endif 
c     
c            df_x = df_r*sthi+df_th*cthi/ri 
c            df_y = df_r*cthi-df_th*sthi/ri 
c
c            df = df_x*sqrt(1.0d0-x(i)**2)+df_y*x(i)
c            dfdat(i,m) = dfdat(i,m)+mp1(n,m)*df
c          enddo
c        enddo 
      enddo 

      zk = kwave*r 
      call besseljs3d(nterm2,zk,rscale2,fjs,ifder,fjder)
      do n=0,nterm2
        fjder(n) = kwave*fjder(n)
        fjsum(n) = (fjs(n)**2+fjder(n)**2)*pi2_inv 
      enddo 

      loc2 = 0
      do i=1,nquad 
        cthi = x(i)
        call zylgndrf(nterm2,cthi,yls(0:nterm2,0:nterm2),rat1,rat2)
        do n=0,nterm2 
          do m=0,n 
            yls(n,m) = yls(n,m)*minus1(m)*pi4_sqrt_inv
          enddo 
        enddo 

        do n=0,nterm2 
c          do m=-n,n 
c            loc2(n,m) = loc2(n,m)+whts(i)*yls(n,abs(m))
c     1  *(fjs(n)*fdat(i,m)+fjder(n)*dfdat(i,m))/fjsum(n)
c          enddo 
          loc2(n,0) = loc2(n,0)+whts(i)*yls(n,0)
     1       *(fjs(n)*fdat(i,0)+fjder(n)*dfdat(i,0))/fjsum(n)
          do m=1,n 
            z1 = yls(n,m)*fjs(n)/fjsum(n)
            z2 = yls(n,m)*fjder(n)/fjsum(n)
            loc2(n,m) = loc2(n,m)+whts(i)*(z1*fdat(i,m)+z2*dfdat(i,m))
            loc2(n,-m) = loc2(n,-m)+whts(i)*
     1                   (z1*fdat(i,-m)+z2*dfdat(i,-m))
          enddo 
        enddo 
      enddo 


      end
C
C
C
C
C
cccc    shift from mp1 tp mp2 along z-axis
      subroutine zh3dzshiftmpmp(r,nquad,zshift,nterm1,nterm2,
     1  rscale1,rscale2,kwave,mp1,mp2)
      implicit real *8 (a-h,o-z)
      complex *16 zshift 
      complex *16 mp1(0:nterm1,-nterm1:nterm1)
      complex *16 mp2(0:nterm2,-nterm2:nterm2)
      complex *16 kwave 
      dimension u(2),v(2)
      complex *16 xi,yi,zi,ri,zk 
      complex *16 cthi,sthi,eiphi 
      complex *16 df_r,df_th,df_x,df_y,df 
      complex *16 ri2,sthi2
      complex *16 z1,z2,z3,cthi_ri,sthi_ri

      complex *16, allocatable:: fhs(:),fhder(:),fhsum(:)
      complex *16, allocatable:: fdat(:,:)
      complex *16, allocatable:: dfdat(:,:)
      complex *16, allocatable:: yls(:,:)
      complex *16, allocatable:: dyls(:,:)
      real *8, allocatable :: x(:),whts(:)
      real *8, allocatable :: rat1(:,:),rat2(:,:)
      integer, allocatable :: minus1(:)

      nterm = max(nterm1,nterm2)

      allocate(fhs(0:nterm+100))
      allocate(fhder(0:nterm+100))
      allocate(fhsum(0:(nterm+100)))
      allocate(fdat(1:nquad,-nterm:nterm))
      allocate(dfdat(1:nquad,-nterm:nterm))
      allocate(yls(0:nterm,0:nterm))
      allocate(dyls(0:nterm,0:nterm))
      allocate(x(nquad),whts(nquad))
      allocate(minus1(0:nterm))
      allocate(rat1(0:nterm,0:nterm))
      allocate(rat2(0:nterm,0:nterm))


      call ylgndrini(nterm,rat1,rat2)
      minus1(0) = 1 
      do i=1,nterm 
        minus1(i) = minus1(i-1)*(-1)
      enddo 

      itype = 1
      ifder = 1 
      pi = atan(1.0d0)*4
      pi4 = pi*4
      pi4_sqrt_inv = 1/dsqrt(pi4)
      pi2_inv = 1/(2*pi)
      call legeexps(itype,nquad,x,u,v,whts)

      fdat = 0 
      dfdat = 0
      do i=1,nquad 
        yi = sqrt(1-x(i)**2)
        zi = zshift+r*x(i)
        ri2 = (r*yi)**2+zi**2
        call zsqrt(ri2,ri)
        cthi = zi/ri 
        sthi2 = 1-cthi**2
        call zsqrt(sthi2,sthi)
        cthi_ri = cthi/ri
        sthi_ri = sthi/ri 

        eiphi = 1 
        zk = ri*kwave 
        yls = 0 
        dyls = 0
        call zylgndr2sf(nterm,cthi,yls(0:nterm,0:nterm),
     1     dyls(0:nterm,0:nterm),rat1,rat2)
        do n=0,nterm1 
          do m=0,n 
            yls(n,m) = yls(n,m)*minus1(m)*pi4_sqrt_inv
            dyls(n,m) = dyls(n,m)*minus1(m)*pi4_sqrt_inv
          enddo 
        enddo 
        call h3dall(nterm,zk,rscale1,fhs,ifder,fhder)

        do n=0,nterm
          fhder(n) = kwave*fhder(n) 
        enddo 

        do n=0,nterm1
          fdat(i,0) = fdat(i,0)+mp1(n,0)*fhs(n)*yls(n,0)
          df_r = fhder(n)*yls(n,0)
          df_th = -fhs(n)*dyls(n,0)*sthi
          df_x = df_r*sthi+df_th*cthi_ri
          df_y = df_r*cthi-df_th*sthi_ri 
          df = df_x*yi+df_y*x(i)
          dfdat(i,0) = dfdat(i,0)+mp1(n,0)*df

          do m=1,n 
            z1 = fhs(n)*yls(n,m)*sthi
            fdat(i,m) = fdat(i,m)+mp1(n,m)*z1 
            fdat(i,-m) = fdat(i,-m)+mp1(n,-m)*z1 

            df_r = fhder(n)*yls(n,m)*sthi
            df_th = -fhs(n)*dyls(n,m)
            df_x = df_r*sthi+df_th*cthi_ri 
            df_y = df_r*cthi-df_th*sthi_ri 
            df = df_x*yi+df_y*x(i)
            dfdat(i,m) = dfdat(i,m)+mp1(n,m)*df
            dfdat(i,-m) = dfdat(i,-m)+mp1(n,-m)*df
          enddo 
        enddo
      enddo 

      zk = kwave*r 
      call h3dall(nterm2,zk,rscale2,fhs,ifder,fhder)
      do n=0,nterm2
        fhder(n) = kwave*fhder(n)
        fhsum(n) = (fhs(n)**2+fhder(n)**2)*pi2_inv 
      enddo 

      mp2 = 0
      do i=1,nquad 
        cthi = x(i)
        eiphi = 1 
        yls = 0 
        call zylgndrf(nterm2,cthi,yls(0:nterm2,0:nterm2),
     1   rat1(0:nterm2,0:nterm2),rat2(0:nterm2,0:nterm2))
        do n=0,nterm2 
          do m=0,n 
            yls(n,m) = yls(n,m)*minus1(m)*pi4_sqrt_inv
          enddo 
        enddo 

        do n=0,nterm2
c          do m=-n,n 
c            mp2(n,m) = mp2(n,m)+2*pi*whts(i)*
c     1          ((-1)**m)*yls(n,abs(m))
c     1          *(fhs(n)*fdat(i,m)+fhder(n)*dfdat(i,m))
c     1          /(fhs(n)**2+fhder(n)**2)
c          enddo 
          mp2(n,0) = mp2(n,0)+whts(i)*yls(n,0)
     1       *(fhs(n)*fdat(i,0)+fhder(n)*dfdat(i,0))/fhsum(n)
          do m=1,n 
            z1 = yls(n,m)*fhs(n)/fhsum(n)
            z2 = yls(n,m)*fhder(n)/fhsum(n)
            mp2(n,m) = mp2(n,m)+whts(i)*(z1*fdat(i,m)+z2*dfdat(i,m))
            mp2(n,-m) = mp2(n,-m)+whts(i)*
     1                   (z1*fdat(i,-m)+z2*dfdat(i,-m))
          enddo
        enddo 
      enddo 

      end 
C
C
C
C
C
cccc    shift from loc1 tp loc2 
      subroutine zh3dzshiftlocloc(r,nquad,zshift,nterm1,nterm2,
     1  rscale1,rscale2,kwave,loc1,loc2)
      implicit real *8 (a-h,o-z)
      complex *16 zshift 
      complex *16 loc1(0:nterm1,-nterm1:nterm1)
      complex *16 loc2(0:nterm2,-nterm2:nterm2)
      complex *16 kwave 
      dimension u(2),v(2)
      complex *16 xi,zi,ri,zk,yi 
      complex *16 cthi,sthi,eiphi 
      complex *16 ri2,sthi2
      complex *16 df_r,df_th,df_x,df_y,df 
      complex *16 z1,z2,z3,cthi_ri,sthi_ri
      real *8, allocatable :: x(:),whts(:)
      complex *16, allocatable:: fjs(:),fjder(:),fjsum(:)
      complex *16, allocatable:: fdat(:,:),dfdat(:,:)
      complex *16, allocatable:: yls(:,:)
      complex *16, allocatable:: dyls(:,:)
      real *8, allocatable :: rat1(:,:),rat2(:,:)
      integer, allocatable :: minus1(:)


      nterm = max(nterm1,nterm2)
      allocate(fjs(0:nterm+100))
      allocate(fjder(0:nterm+100))
      allocate(fjsum(0:(nterm+100)))
      allocate(fdat(1:nquad,-nterm:nterm))
      allocate(dfdat(1:nquad,-nterm:nterm))
      allocate(x(nquad),whts(nquad))
      allocate(yls(0:nterm,0:nterm))
      allocate(dyls(0:nterm,0:nterm))
      allocate(minus1(0:nterm))
      allocate(rat1(0:nterm,0:nterm))
      allocate(rat2(0:nterm,0:nterm))

      call ylgndrini(nterm,rat1,rat2)

      minus1(0) = 1 
      do i=1,nterm 
        minus1(i) = minus1(i-1)*(-1)
      enddo 


      itype = 1
      ifder = 1 
      pi = atan(1.0d0)*4
      pi4 = pi*4
      pi4_sqrt_inv = 1/sqrt(pi4)
      pi2_inv = 1/(2*pi)
      call legeexps(itype,nquad,x,u,v,whts)
      fdat = 0 
      dfdat = 0
      do i=1,nquad 
        yi = sqrt(1-x(i)**2)
        xi = r*sqrt(1-x(i)**2)
        zi = zshift+r*x(i)
        ri2 = xi**2+zi**2
        call zsqrt(ri2,ri)
        cthi = zi/ri 
        sthi2 = 1-cthi**2
        call zsqrt(sthi2,sthi)
        cthi_ri = cthi/ri 
        sthi_ri = sthi/ri 

        zk = ri*kwave 
        call zylgndr2sf(nterm1,cthi,yls(0:nterm1,0:nterm1),
     1    dyls(0:nterm1,0:nterm1),rat1(0:nterm1,0:nterm1),
     2    rat2(0:nterm1,0:nterm1))
        do n=0,nterm1 
          do m=0,n 
            yls(n,m) = yls(n,m)*pi4_sqrt_inv*minus1(m)
            dyls(n,m) = dyls(n,m)*pi4_sqrt_inv*minus1(m)
          enddo 
        enddo 

        call besseljs3d(nterm1,zk,rscale1,fjs,ifder,fjder)
        do n=0,nterm1
            fjder(n) = kwave*fjder(n) 
        enddo 

        do n=0,nterm1 
          fdat(i,0) = fdat(i,0)+loc1(n,0)*fjs(n)*yls(n,0)
          df_r = fjder(n)*yls(n,0)
          df_th = -fjs(n)*dyls(n,0)*sthi
          df_x = df_r*sthi+df_th*cthi_ri
          df_y = df_r*cthi-df_th*sthi_ri 
          df = df_x*yi+df_y*x(i)
          dfdat(i,0) = dfdat(i,0)+loc1(n,0)*df

          do m=1,n 
            z1 = fjs(n)*yls(n,m)*sthi
            fdat(i,m) = fdat(i,m)+loc1(n,m)*z1 
            fdat(i,-m) = fdat(i,-m)+loc1(n,-m)*z1 

            df_r = fjder(n)*yls(n,m)*sthi
            df_th = -fjs(n)*dyls(n,m)
            df_x = df_r*sthi+df_th*cthi_ri 
            df_y = df_r*cthi-df_th*sthi_ri 
            df = df_x*yi+df_y*x(i)
            dfdat(i,m) = dfdat(i,m)+loc1(n,m)*df
            dfdat(i,-m) = dfdat(i,-m)+loc1(n,-m)*df
          enddo 
        enddo 
      enddo 

      zk = kwave*r 
      call besseljs3d(nterm2,zk,rscale2,fjs,ifder,fjder)
      do n=0,nterm2
        fjder(n) = kwave*fjder(n)
        fjsum(n) = (fjs(n)**2+fjder(n)**2)*pi2_inv 
      enddo 

      loc2 = 0
      do i=1,nquad 
        cthi = x(i)
        call zylgndrf(nterm2,cthi,yls(0:nterm2,0:nterm2),
     1 rat1(0:nterm2,0:nterm2),rat2(0:nterm2,0:nterm2))
        do n=0,nterm2 
          do m=0,n 
            yls(n,m) = minus1(m)*yls(n,m)*pi4_sqrt_inv
          enddo 
        enddo 

        do n=0,nterm2 
          loc2(n,0) = loc2(n,0)+whts(i)*yls(n,0)
     1       *(fjs(n)*fdat(i,0)+fjder(n)*dfdat(i,0))/fjsum(n)
          do m=1,n 
            z1 = yls(n,m)*fjs(n)/fjsum(n)
            z2 = yls(n,m)*fjder(n)/fjsum(n)


            loc2(n,m) = loc2(n,m)+whts(i)*(z1*fdat(i,m)+z2*dfdat(i,m))
            loc2(n,-m) = loc2(n,-m)+whts(i)*(z1*fdat(i,-m)
     1               +z2*dfdat(i,-m))
          enddo 
        enddo 
      enddo 
      end