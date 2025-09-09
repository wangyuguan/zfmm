c
c   subroutines for translating the multipole and local 
c   expansions for 2d laplace
c
c
c
c
      subroutine zl2dlocloc(zc1,zc2,rscale1,rscale2,nterm1,nterm2,
     1 loc1,loc2,carray,ldc)  
c----------------------------------------------        
c     translate a local expansion loc2 at zc2 to a local expansion at zc1
c     and add to loc1      
c
c----------------------------------------------
c   INPUT PARAMETERS:
c   zc1(2)       : expansion center of target local expansion 
c   zc2(2)       : expansion center of original local expansion 
c   rscale1      : scaling parameter of target local expansion 
c   rscale2      : scaling parameter of original local expansion 
c   nterm1       : truncation limit of target local expansion 
c   nterm2       : truncation limit of original local expansion  
c   loc2(-nterm2:nterm2)   : coefficients of original local expansion 
c   carray(0:ldc,0:ldc)    : binomial coefficients 
c   ldc          : max order of binomial coefficients
c
c
c   OUTPUT PARAMETERS
c   loc1(-nterm1:nterm1)   : coefficients of target local expansion 
c----------------------------------------------
      implicit real *8 (a-h,o-z)
      complex *16 r0,e0
      complex *16 r0all(0:(nterm1+nterm2)),e0all(0:(nterm1+nterm2))
      complex *16 zc1(2),zc2(2)
      complex *16 loc1(-nterm1:nterm1),loc2(-nterm2:nterm2)
      dimension rscales1(0:nterm1),rscales2(0:nterm2)
      dimension carray(0:ldc,0:ldc)

      rscales1(0) = 1 
      do n=1,nterm1   
        rscales1(n) = rscales1(n-1)*rscale1 
      enddo 

      rscales2(0) = 1 
      do n=1,nterm2   
        rscales2(n) = rscales2(n-1)*rscale2 
      enddo 

      call zcart2polar(zc2,zc1,r0,e0)
      nterm = max(nterm1,nterm2)
      r0all = 1 
      e0all = 1 
      do n=1,nterm 
        r0all(n) = r0all(n-1)*r0 
        e0all(n) = e0all(n-1)*e0 
      enddo 

      loc1(0) = loc1(0)+loc2(0)
      do n=1,nterm2
        loc1(0) = loc1(0)+((-1)**n)*r0all(n)*e0all(n)*loc2(n)
     1            /(rscales2(n))
        loc1(0) = loc1(0)+((-1)**n)*r0all(n)/e0all(n)*loc2(-n)
     1            /(rscales2(n))
      enddo 

      do k=1,nterm1 
        do n=k,nterm2 
          loc1(k) = loc1(k)+((-1)**(n-k))*carray(n,k)*r0all(n-k)
     1  *e0all(n-k)*loc2(n)*rscales1(k)/rscales2(n)
          loc1(-k) = loc1(-k)+((-1)**(n-k))*carray(n,k)*r0all(n-k)
     1  /e0all(n-k)*loc2(-n)*rscales1(k)/rscales2(n)
        enddo 
      enddo 


      end
c
c
c
c
c
      subroutine zl2dmploc(zc1,zc2,rscale1,rscale2,nterm1,nterm2,loc,
     1 mp,carray,ldc)
c----------------------------------------------        
c     convert a multipole expansion mp at zc2 to a local expansion at zc1
c     and add to loc 
c----------------------------------------------
c   INPUT PARAMETERS:
c   zc1(2)       : expansion center of target local expansion 
c   zc2(2)       : expansion center of original multipole expansion 
c   rscale1      : scaling parameter of target local expansion
c   rscale2      : scaling parameter of original multipole expansion
c   nterm1        : truncation limit of target local expansions 
c   nterm2        : truncation limit of original multipole expansions 
c   mp(-nterm2:nterm2)   : coefficients of original multipole expansion 
c   carray(0:ldc,0:ldc)    : binomial coefficients 
c   ldc          : max order of binomial coefficients
c
c   OUTPUT PARAMETERS
c   loc(-nterm1:nterm1)   : coefficients of target local expansion 
c----------------------------------------------
      implicit real *8 (a-h,o-z)
      complex *16 r0,e0
      complex *16 r0all(0:(nterm1+nterm2)),e0all(0:(nterm1+nterm2))
      complex *16 zc1(2),zc2(2)
      complex *16 loc(-nterm1:nterm1),mp(-nterm2:nterm2)
      dimension rscales1(0:nterm1),rscales2(0:nterm2)
      dimension carray(0:ldc,0:ldc)


      rscales1(0) = 1 
      do n=1,nterm1   
        rscales1(n) = rscales1(n-1)*rscale1 
      enddo 

      rscales2(0) = 1 
      do n=1,nterm2   
        rscales2(n) = rscales2(n-1)*rscale2 
      enddo 

      done = 1 
      call zcart2polar(zc2,zc1,r0,e0)
 

      nterm = nterm1+nterm2 
      r0all(0) = 1 
      e0all(0) = 1
      do m=1,nterm 
        r0all(m) = r0all(m-1)*r0 
        e0all(m) = e0all(m-1)*e0 
      enddo

      loc(0) = loc(0)+mp(0)*log(r0)
      do n=1,nterm2 
        loc(0) = loc(0)+rscales2(n)*((-1)**n)*mp(n)*e0all(n)/r0all(n)
        loc(0) = loc(0)+rscales2(n)*((-1)**n)*mp(-n)/e0all(n)/r0all(n)
      enddo 

      do k=1,nterm1 
        loc(k)=loc(k)-rscales1(k)*mp(0)*((done/2)/k)/e0all(k)/r0all(k)
        loc(-k) = loc(-k)
     1             -rscales1(k)*mp(0)*((done/2)/k)*e0all(k)/r0all(k)
        do n=1,nterm2 
          loc(k) = loc(k)+((-1)**n)*mp(-n)
     1   *rscales1(k)*rscales2(n)*carray(n+k-1,k)/e0all(n+k)/r0all(n+k)
          loc(-k) = loc(-k)+((-1)**n)*mp(n)
     1   *rscales1(k)*rscales2(n)*carray(n+k-1,k)*e0all(n+k)/r0all(n+k)
        enddo 
      enddo 

      end
c
c
c
c
c
      subroutine zl2dmpmp(zc1,zc2,rscale1,rscale2,nterm1,nterm2,mp1,
     1 mp2,carray,ldc)
c----------------------------------------------        
c     translate a multipole expansion mp2 at zc2 to a multipole expansion at zc1
c     and add to mp1 
c----------------------------------------------
c   INPUT PARAMETERS:
c   mp1(2)       : expansion center of target multipole expansion 
c   mp2(2)       : expansion center of original multipole expansion 
c   rscale1      : scaling parameter of target multipole expansion 
c   rscale2      : scaling parameter of original multipole expansion 
c   nterm1       : truncation limit of target multipole expansion 
c   nterm2       : truncation limit of original multipole expansion  
c   mp2(-nterm2:nterm2)   : coefficients of original multipole expansion 
c   carray(0:ldc,0:ldc)    : binomial coefficients 
c   ldc          : max order of binomial coefficients
c
c   OUTPUT PARAMETERS
c   mp1(-nterm1:nterm1)   : coefficients of target multipole expansion 
c----------------------------------------------
      implicit real *8 (a-h,o-z)
      complex *16 r0,e0
      complex *16 r0all(0:(nterm1+nterm2)),e0all(0:(nterm1+nterm2))
      complex *16 zc1(2),zc2(2)
      complex *16 mp1(-nterm1:nterm1),mp2(-nterm2:nterm2)
      dimension rscales1(0:nterm1),rscales2(0:nterm2)
      dimension carray(0:ldc,0:ldc)

      done = 1 

      rscales1(0) = 1 
      do n=1,nterm1   
        rscales1(n) = rscales1(n-1)*rscale1 
      enddo 

      rscales2(0) = 1 
      do n=1,nterm2   
        rscales2(n) = rscales2(n-1)*rscale2 
      enddo 

      call zcart2polar(zc2,zc1,r0,e0)
      nterm = max(nterm1,nterm2)
      r0all(0) = 1 
      e0all(0) = 1 
      do n=1,nterm 
        r0all(n) = r0all(n-1)*r0 
        e0all(n) = e0all(n-1)*e0 
      enddo 

      mp1(0) = mp1(0)+mp2(0)


      do k=1,nterm1

        mp1(k) = mp1(k)-
     2            ((done/2)/k)*r0all(k)/e0all(k)*mp2(0)/rscales1(k) 
        mp1(-k) = mp1(-k)-
     1            ((done/2)/k)*r0all(k)*e0all(k)*mp2(0)/rscales1(k) 

        do n=1,k

          rfrac = rscales2(n)/rscales1(k) 

          if (n.le.nterm2) then 
            mp1(k) = mp1(k)
     1        +rfrac*carray(k-1,k-n)*r0all(k-n)/e0all(k-n)*mp2(n)
            mp1(-k) = mp1(-k)
     1             +rfrac*carray(k-1,k-n)*r0all(k-n)*e0all(k-n)*mp2(-n)
          endif   

        enddo 
      enddo 

      end 