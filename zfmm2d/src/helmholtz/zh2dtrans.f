C
C   subroutines for translating the multipole and local 
C   expansions for 2d helmholtz using Graf's addition formulas 
C
C
      subroutine zh2dlocloc(zk,zc1,zc2,rscale1,rscale2,nterm1,nterm2,
     1 loc1,loc2)  
c----------------------------------------------        
c     translate a local expansion loc2 at zc2 to a local expansion at zc1
c     and add to loc1      
c
c----------------------------------------------
c   INPUT PARAMETERS:
c   zk           : Helmholtz parameter 
c   zc1(2)       : expansion center of target local expansion 
c   zc2(2)       : expansion center of original local expansion 
c   rscale1      : scaling parameter of target local expansion 
c   rscale2      : scaling parameter of original local expansion 
c   nterm1       : truncation limit of target local expansion 
c   nterm2       : truncation limit of original local expansion  
c   loc2(-nterm2:nterm2)   : coefficients of original local expansion 
c
c
c   OUTPUT PARAMETERS
c   loc1(-nterm1:nterm1)   : coefficients of target local expansion 
c----------------------------------------------
      implicit real *8 (a-h,o-z)
      complex *16 zk,rz,ez,ztmp,z
      complex *16 zc1(2),zc2(2)
      complex *16 loc1(-nterm1:nterm1),loc2(-nterm2:nterm2)
      complex *16 fjs(0:(nterm1+nterm2+100))
      complex *16 fjder(0:(nterm1+nterm2+100))
      complex *16 ezs(-(nterm1+nterm2):(nterm1+nterm2))
      dimension rscales1(0:(nterm1+nterm2))
      dimension rscales2(0:(nterm1+nterm2))

      nterm = nterm1+nterm2

      call zcart2polar(zc1,zc2,rz,ez)
      ifder0 = 0 
      

      ezs(0) = 1 
      do m=1,nterm 
        ezs(m) = ezs(m-1)*ez 
      enddo

      do m=-nterm,-1 
        ezs(m) = 1/ezs(-m) 
      enddo

      z = rz*zk
      call jbessel2d(nterm,z,rscale1,fjs,ifder0,fjder)

      rscales1(0) = 1 
      do m=1,nterm1+nterm2 
        rscales1(m) = rscales1(m-1)*rscale1
      enddo 
      rscales2(0) = 1 
      do m=1,nterm1+nterm2
        rscales2(m) = rscales2(m-1)*rscale2
      enddo 

      do m=-nterm1,nterm1 
        do n=-nterm2,nterm2 
c          rfrac = rscale1**(abs(m)+abs(m-n))/rscale2**(abs(n))
          rfrac = rscales1(abs(m))*rscales1(abs(m-n))/rscales2(abs(n))
          if ((n-m).ge.0)then 
            ztmp = fjs(n-m)*rfrac
          else 
            ztmp = ((-1)**(n-m))*fjs(m-n)*rfrac
          endif 
          ztmp = loc2(n)*ztmp*(ezs(n-m))
          loc1(m) = loc1(m)+ztmp
        enddo 
      enddo 

      end
C
C
C
C
C
      subroutine zh2dmploc(zk,zc1,zc2,rscale,nterm,loc,mp)
c----------------------------------------------        
c     convert a multipole expansion mp at zc2 to a local expansion at zc1
c     and add to loc 
c----------------------------------------------
c   INPUT PARAMETERS:
c   zk           : Helmholtz parameter 
c   zc1(2)       : expansion center of target local expansion 
c   zc2(2)       : expansion center of original multipole expansion 
c   rscale      : scaling parameter of multipole and local expansion
c   nterm        : truncation limit of multipole and local expansions 
c   mp(-nterm:nterm)   : coefficients of original multipole expansion 
c
c
c   OUTPUT PARAMETERS
c   loc(-nterm:nterm)   : coefficients of target local expansion 
c----------------------------------------------
      implicit real *8 (a-h,o-z)
      complex *16 zk,rz,ez,ztmp
      complex *16 zc1(2),zc2(2)
      complex *16 mp(-nterm:nterm),loc(-nterm:nterm)
      complex *16 fhs(0:2*nterm+100),fhder(0:2*nterm+100)
      complex *16 ezs(-(2*nterm):(2*nterm))
      complex *16 rscales(0:2*nterm) 

      call zcart2polar(zc1,zc2,rz,ez)
      ezs(0) = 1 
      do m=1,2*nterm 
        ezs(m) = ezs(m-1)*ez 
      enddo

      do m=-2*nterm,-1 
        ezs(m) = 1/ezs(-m) 
      enddo

      ifder0 = 0 
      call h2dall(2*nterm,rz*zk,rscale,fhs,ifder0,fhder)

      rscales(0) = 1 
      do m=1,2*nterm 
        rscales(m) = rscales(m-1)*rscale 
      enddo 

      do m=-nterm,nterm
        do n=-nterm,nterm
          if ((n-m).ge.0)then
            ztmp = fhs(n-m)
          else  
            ztmp = fhs(m-n)*((-1)**(m-n))
          endif 
          ztmp = ztmp*mp(n)*(ezs(n-m))

c          rscale_pow = rscale**(abs(m)+abs(n)-abs(m-n))
          ztmp = ztmp*rscales(abs(m))*rscales(abs(n))/rscales(abs(m-n))
          loc(m) = loc(m)+ztmp
        enddo
      enddo 

      end
c
c
c
c
c c
c       subroutine zh2dmploc1(zk,zc1,zc2,rscale1,rscale2,nterm1,nterm2,
c      1 loc,mp)
c c----------------------------------------------        
c c     convert a multipole expansion mp at zc2 to a local expansion at zc1
c c     and add to loc 
c c----------------------------------------------
c c   INPUT PARAMETERS:
c c   zk           : Helmholtz parameter 
c c   zc1(2)       : expansion center of target local expansion 
c c   zc2(2)       : expansion center of original multipole expansion 
c c   rscale      : scaling parameter of multipole and local expansion
c c   nterm        : truncation limit of multipole and local expansions 
c c   mp(-nterm:nterm)   : coefficients of original multipole expansion 
c c
c c
c c   OUTPUT PARAMETERS
c c   loc(-nterm:nterm)   : coefficients of target local expansion 
c c----------------------------------------------
c       implicit real *8 (a-h,o-z)
c       complex *16 zk,rz,ez,ztmp
c       complex *16 zc1(2),zc2(2)



c       allocate() mp(-nterm:nterm),loc(-nterm:nterm)
c       allocate() fhs(0:2*nterm+100),fhder(0:2*nterm+100)
c       allocate() ezs(-(2*nterm):(2*nterm))
c       allocate(rscales1(0:nterm1)) 
c       allocate(rscales2(0:nterm2)) 


c       nterm = max(nterm1,nterm2)

c       call zcart2polar(zc1,zc2,rz,ez)
c       ezs(0) = 1 
c       do m=1,2*nterm 
c         ezs(m) = ezs(m-1)*ez 
c       enddo

c       do m=-2*nterm,-1 
c         ezs(m) = 1/ezs(-m) 
c       enddo

c       ifder0 = 0 
c       call h2dall(2*nterm,rz*zk,rscale,fhs(0:2*nterm),ifder0,fhder)

c       rscales1(0) = 1 
c       do m=1,nterm1
c         rscales1(m) = rscales1(m-1)*rscale1 
c       enddo 

c       rscales2(0) = 1 
c       do m=1,nterm2
c         rscales2(m) = rscales2(m-1)*rscale2
c       enddo 

c       do m=-nterm,nterm
c         do n=-nterm,nterm
c           if ((n-m).ge.0)then
c             ztmp = fhs(n-m)
c           else  
c             ztmp = fhs(m-n)*((-1)**(m-n))
c           endif 
c           ztmp = ztmp*mp(n)*(ezs(n-m))

c c          rscale_pow = rscale**(abs(m)+abs(n)-abs(m-n))
c           ztmp = ztmp*rscales(abs(m)+abs(n)-abs(m-n))
c           loc(m) = loc(m)+ztmp
c         enddo
c       enddo 

c       end
c C
C
C
C
C
      subroutine zh2dmpmp(zk,zc1,zc2,rscale1,rscale2,nterm1,nterm2,mp1,
     1 mp2)
c----------------------------------------------        
c     translate a multipole expansion mp2 at zc2 to a multipole expansion at zc1
c     and add to mp1 
c----------------------------------------------
c   INPUT PARAMETERS:
c   zk           : Helmholtz parameter 
c   mp1(2)       : expansion center of target multipole expansion 
c   mp2(2)       : expansion center of original multipole expansion 
c   rscale1      : scaling parameter of target multipole expansion 
c   rscale2      : scaling parameter of original multipole expansion 
c   nterm1       : truncation limit of target multipole expansion 
c   nterm2       : truncation limit of original multipole expansion  
c   mp2(-nterm2:nterm2)   : coefficients of original multipole expansion 
c
c
c   OUTPUT PARAMETERS
c   mp1(-nterm1:nterm1)   : coefficients of target multipole expansion 
c----------------------------------------------
      implicit real *8 (a-h,o-z)
      complex *16 zc1(2),zc2(2),zk
      complex *16 mp1(-nterm1:nterm1),mp2(-nterm2:nterm2)
      complex *16 ima
      complex *16 fjs(0:(nterm1+nterm2+100))
      complex *16 fdjer(0:(nterm1+nterm2+100))
      complex *16 rz,ez,ztmp
      complex *16 ezs(-(nterm1+nterm2):(nterm1+nterm2))
      dimension rscales1(-(nterm1+nterm2):(nterm1+nterm2))
      dimension rscales2(0:(nterm1+nterm2))

      nterm = nterm1+nterm2 
    
      call zcart2polar(zc1,zc2,rz,ez)
      

      ezs(0) = 1 
      do m=1,nterm 
        ezs(m) = ezs(m-1)*ez 
      enddo

      do m=-nterm,-1 
        ezs(m) = 1/ezs(-m) 
      enddo

      ifder0 = 0 
      call jbessel2d(nterm,rz*zk,rscale1,fjs,ifder0,fjder)

      rscales1(0) = 1 
      do m=1,nterm1+nterm2 
        rscales1(m) = rscales1(m-1)*rscale1
      enddo 

      do m=-(nterm1+nterm2),-1 
        rscales1(m) = 1/rscales1(-m)
      enddo 


      rscales2(0) = 1 
      do m=1,nterm1+nterm2
        rscales2(m) = rscales2(m-1)*rscale2
      enddo 

      do m=-nterm1,nterm1 
        do n=-nterm2,nterm2
c          rfrac = (rscale2**abs(n))/(rscale1**(abs(m)-abs(m-n)))
          rfrac = rscales2(abs(n))*rscales1(abs(m-n))/rscales1(abs(m))

          if ((m-n).ge.0)then 
            ztmp = fjs(m-n)*rfrac
          else 
            ztmp = ((-1)**(n-m))*fjs(n-m)*rfrac
          endif 
          ztmp = ztmp*mp2(n)*(ezs(n-m))*((-1)**(m+n))
          mp1(m) = mp1(m)+ztmp
        enddo 
      enddo 

      end 