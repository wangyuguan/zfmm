        subroutine cdist(x,y,zr)
        implicit real *8 (a-h,o-z)
        complex *16 x(3),y(3)
        complex *16 zr 
        zr = sqrt((x(1)-y(1))**2+(x(2)-y(2))**2+(x(3)-y(3))**2)
        end
c
c
c
c
c
        subroutine cinput(z,zc,zr,ctheta,eiphi)
        implicit real *8 (a-h,o-z)
        complex *16 z(3),zc(3)
        complex *16 zr,ctheta,eiphi 
        complex *16 rxy
        complex *16 ima 
        complex *16 x1,x2,x3 
        ima = (0.0d0,1.0d0)
        x1 = z(1)-zc(1)
        x2 = z(2)-zc(2)
        x3 = z(3)-zc(3)
        zr = sqrt(x1**2+x2**2+x3**2)
        rxy = sqrt(x1**2+x2**2)
        ctheta = x3/zr 
        eiphi = x1/rxy+ima*x2/rxy
        end 
c
c
c
c
c
        subroutine zspharm(nmax,ctheta,eiphi,y)
        implicit real *8 (a-h,o-z)
        complex *16 ctheta
        complex *16 eiphi
        complex *16 y(0:nmax,-nmax:nmax)
        complex *16 yl(0:nmax,0:nmax)

        pi = atan(1.0d0)*4.0d0
        call zylgndr(nmax,ctheta,yl)

        y = 0
        do n=0,nmax 
            do m=-n,n
                y(n,m) = ((-1)**m)*yl(n,abs(m))*(eiphi**m)
                y(n,m) = y(n,m)/sqrt(4.0d0*pi)
            enddo 
        enddo 

        end
C
C
C
C
C
        subroutine unzspharm(nmax,ctheta,eiphi,y)
        implicit real *8 (a-h,o-z)
        complex *16 ctheta
        complex *16 eiphi
        complex *16 y(0:nmax,-nmax:nmax)
        complex *16 yl(0:nmax,0:nmax)

        pi = atan(1.0d0)*4.0d0
        call zylgndr(nmax,ctheta,yl)

        y = 0
        done = 1 
        do n=0,nmax 
            do m=-n,n
                y(n,m) = yl(n,abs(m))*(eiphi**m)/sqrt(done*(2*n+1))
            enddo 
        enddo 

        end
C
C
C
C
C
        subroutine zspharmf(nmax,ctheta,eiphi,y,rat1,rat2)
        implicit real *8 (a-h,o-z)
        complex *16 ctheta
        complex *16 eiphi
        complex *16 y(0:nmax,-nmax:nmax)
        complex *16 yl(0:nmax,0:nmax)
        dimension rat1(0:nmax,0:nmax),rat2(0:nmax,0:nmax)

        pi = atan(1.0d0)*4.0d0
        call zylgndrf(nmax,ctheta,yl,rat1,rat2)

        y = 0
        do n=0,nmax 
            do m=-n,n
                y(n,m) = ((-1)**m)*yl(n,abs(m))*(eiphi**m)
                y(n,m) = y(n,m)/sqrt(4.0d0*pi)
            enddo 
        enddo 

        end
c
c
c
c
c
        subroutine zspharm2(nmax,ctheta,eiphi,y,dy)
        implicit real *8 (a-h,o-z)
        complex *16 ctheta,eiphi
        complex *16 y(0:nmax,-nmax:nmax)
        complex *16 dy(0:nmax,-nmax:nmax)
        complex *16 yl(0:nmax,0:nmax)
        complex *16 dyl(0:nmax,0:nmax)
        complex *16 stheta 

        stheta = sqrt(1-ctheta**2)
        pi = atan(1.0d0)*4.0d0

        call zylgndr2s(nmax,ctheta,yl,dyl)

        do n=0,nmax 
            do m=-n,n 
                y(n,m) = ((-1)**m)*yl(n,abs(m))*(eiphi**m)
                y(n,m) = y(n,m)/sqrt(4.0d0*pi)
                if (m.ne.0) then 
                    y(n,m) = y(n,m)*stheta;
                endif 

                dy(n,m) = ((-1)**m)*dyl(n,abs(m))*(eiphi**m)
                dy(n,m) = -dy(n,m)/sqrt(4.0d0*pi)

                if (m.eq.0)then 
                    dy(n,m) = dy(n,m)*stheta
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
        subroutine zspharm2f(nmax,ctheta,eiphi,y,dy,rat1,rat2)
        implicit real *8 (a-h,o-z)
        complex *16 ctheta,eiphi
        complex *16 y(0:nmax,-nmax:nmax)
        complex *16 dy(0:nmax,-nmax:nmax)
        complex *16 yl(0:nmax,0:nmax)
        complex *16 dyl(0:nmax,0:nmax)
        complex *16 stheta 
        dimension rat1(0:nmax,0:nmax),rat2(0:nmax,0:nmax)

        stheta = sqrt(1-ctheta**2)
        pi = atan(1.0d0)*4.0d0

        call zylgndr2sf(nmax,ctheta,yl,dyl,rat1,rat2)

        do n=0,nmax 
            do m=-n,n 
                y(n,m) = ((-1)**m)*yl(n,abs(m))*(eiphi**m)
                y(n,m) = y(n,m)/sqrt(4.0d0*pi)
                if (m.ne.0) then 
                    y(n,m) = y(n,m)*stheta;
                endif 

                dy(n,m) = ((-1)**m)*dyl(n,abs(m))*(eiphi**m)
                dy(n,m) = -dy(n,m)/sqrt(4.0d0*pi)

                if (m.eq.0)then 
                    dy(n,m) = dy(n,m)*stheta
                endif 
            enddo
        enddo 
        end 
c
c
c
c
c
        subroutine zspharm2s(nmax,ctheta,eiphi,y,dy)
        implicit real *8 (a-h,o-z)
        complex *16 ctheta,eiphi
        complex *16 y(0:nmax,-nmax:nmax)
        complex *16 dy(0:nmax,-nmax:nmax)
        complex *16 yl(0:nmax,0:nmax)
        complex *16 dyl(0:nmax,0:nmax)
        complex *16 stheta 

        stheta = sqrt(1-ctheta**2)
        pi = atan(1.0d0)*4.0d0

        yl = 0 
        dyl = 0
        call zylgndr2s(nmax,ctheta,yl,dyl)

        do n=0,nmax 
            do m=-n,n 
                y(n,m) = ((-1)**m)*yl(n,abs(m))*(eiphi**m)
                y(n,m) = y(n,m)/sqrt(4.0d0*pi)

                dy(n,m) = ((-1)**m)*dyl(n,abs(m))*(eiphi**m)
                dy(n,m) = -dy(n,m)/sqrt(4.0d0*pi)
            enddo
        enddo 
        end 
C
C
C
C
C
        subroutine zspharm2sf(nmax,ctheta,eiphi,y,dy,rat1,rat2)
        implicit real *8 (a-h,o-z)
        complex *16 ctheta,eiphi
        complex *16 y(0:nmax,-nmax:nmax)
        complex *16 dy(0:nmax,-nmax:nmax)
        complex *16 yl(0:nmax,0:nmax)
        complex *16 dyl(0:nmax,0:nmax)
        complex *16 stheta 
        dimension rat1(0:nmax,0:nmax),rat2(0:nmax,0:nmax)

        stheta = sqrt(1-ctheta**2)
        pi = atan(1.0d0)*4.0d0

        yl = 0 
        dyl = 0
        call zylgndr2sf(nmax,ctheta,yl,dyl,rat1,rat2)

        do n=0,nmax 
            do m=-n,n 
                y(n,m) = ((-1)**m)*yl(n,abs(m))*(eiphi**m)
                y(n,m) = y(n,m)/sqrt(4.0d0*pi)

                dy(n,m) = ((-1)**m)*dyl(n,abs(m))*(eiphi**m)
                dy(n,m) = -dy(n,m)/sqrt(4.0d0*pi)
            enddo
        enddo 
        end 