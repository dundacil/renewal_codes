      implicit real*8 (a-h,p-z)
C Riscrittura di momdic per convincere Marco che il problema non e` come si calcola la funzione di correlazione  
C Code to compute the four times correlation function    
      parameter (np=12,nd=2**np,nr=1000,nf=18,nde=16)
      parameter (ncut=nint(0.005d0*nd))
      dimension xin(nd),xfr(nd,nf),xf2(nd,3)
      dimension rho(0:nr)
      common tgrade,amu
      character*120 cwd
      character*8 date
      character*10 ctime
      character*5 zone
      integer*4 values(8)

      external fhi,fhit
      f(x) = (amu-1)/(t0)/(1+x/t0)**amu

      ist = getcwd(cwd)
      j = len_trim(cwd)
      kdum = iachar(cwd(j:j))-48
      jdum = iachar(cwd(j-1:j-1))-48
      call date_and_time( date, ctime, zone, values )

      idum = -(jdum*1000+kdum*100)-values(3)-17
      print *,"IDUM",idum,values(3)

c      idum = -53
      amu = 3.50d0
      t0 = 1.d0
      tgrade = t0
!      ncut = 3000

      rho = 0.d0
      n_corr = 10
      alfa = -1.d0/(amu-1)


      xfr = 0.d0
      xf2 = 0.d0
      state = rnor(idum)
      nswitch = 1
      nave = 40000000

      do iave = 1,nave
         if (mod(iave,1000).eq.0) print *,iave
         io = 1
         dt = 0.1d0
         tmax = dt*nd
         x = power(idum,alfa,t0)
         tswitch = x
         j = int(x*10.d0)
         if (j.le.nr.and.j.ge.0) then 
            rho(j) = rho(j) + 1.d0
         endif
      
         time = 0.d0
         xin = 0.d0
         do while (time .le . tmax)
            xin(io) = state
            if(iave.eq.1) write(999,*) io*dt,xin(io),io
            time = time + dt
            io = io + 1

            if (time .gt. tswitch) then
               dummy = rnor(idum) ! call one for good measure, not really needed
               x = power(idum,alfa,t0)
               nswitch = nswitch + 1
               j = int(x*10.d0)
               if(j.le.nr.and.j.ge.0) then 
                  rho(j) = rho(j) + 1.d0
               endif
               tswitch = tswitch+x
               state = rnor(idum)
            endif
         enddo

         n4 = nd/8
         do k = 1,nf
            if (k.le.6) then
               i1 = 1
            else if (k.le.12) then
               i1 = n4
            else
               i1 = 1024
            endif
            if(mod(k,6).eq.1) then
               r2 = 0.25d0
               r3 = 0.75d0
            endif
            if(mod(k,6).eq.2) then
               r2 = 0.25d0
               r3 = 2.d0/3.d0
            endif
            if(mod(k,6).eq.3) then
               r2 = 0.25d0
               r3 = 0.5d0
            endif
            if(mod(k,6).eq.4) then
               r2 = 0.5d0
               r3 = 2.d0/3.d0
            endif
            if(mod(k,6).eq.5) then
               r2 = 1.d0/3.d0
               r3 = 2.d0/3.d0
            endif
            if(mod(k,6).eq.0) then
               r2 = 0.25d0
               r3 = 2.d0/3.d0
            endif
            do i4 = i1+ncut,nd
               i3 = nint(i1+r3*(i4-i1))
               i2 = nint(i1+r2*(i4-i1))
               xfr(i4,k) = xfr(i4,k) + xin(i1)*xin(i2)*xin(i3)*xin(i4)
               if(iave.eq.1) write(200+k,*) i1,i2,i3,i4 
               if(iave.eq.1) write(300+k,*) dt*(i1-1),dt*(i2-1),
     +                                   dt*(i3-1),dt*(i4-1) 
               if (mod(k,6).eq.1) then
                  ii1 = k/6+1
                  xf2(i4,ii1) = xf2(i4,ii1) + xin(i1)*xin(i4)
               endif
            enddo
            close(200+k)
         enddo
      enddo

      do k = 1,nf
         if (k.le.6) then
            i1 = 1
         else if (k.le.12) then
            i1 = n4
         else
            i1 = 2*n4
         endif
         do i4 = i1+ncut,nd
            write(100+k,*) (i4-1)*dt,xfr(i4,k)/nave
         enddo
      enddo
      do k = 1,nf,6
         if (k.le.6) then
            i1 = 1
         else if (k.le.12) then
            i1 = n4
         else
            i1 = 2*n4
         endif
         ii1 = k/6+1
         do i4 = i1+ncut,nd
            write(500+ii1,*)(i4-1)*dt,xf2(i4,ii1)/nave
         enddo
      enddo
      stop
      end

      function fhi(x)
      implicit real*8 (a-h,p-z)
      common t0,amu
      fhi = (amu-1)/(t0)/(1+x/t0)**amu
      return
      end

      function fhit(x)
      implicit real*8 (a-h,p-z)
      common t0,amu
      fhit = x*(amu-1)/(t0)/(1+x/t0)**amu
      return
      end
