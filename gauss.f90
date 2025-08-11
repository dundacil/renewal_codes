      program spagnolo
      use MPI
      implicit real*8 (a-h,o-z)
! Code to compute the time dependent probability distribution
      integer, parameter :: np=9,nd=2**np,nr=200,nf=18,nde=16,ninner=100
      parameter (ncut=nint(0.005d0*nd))
      integer :: ierr, rank, size, root
      dimension xin(nd),xfr(nd,nf),xf2(nd,3)
      DOUBLE PRECISION :: rho_local(0:nd,-nr:nr), rho_global(0:nd,-nr:nr)
      character*23 filename
!      character*120 cwd
!      character*8 date
!      character*10 ctime
!      character*5 zone
!      integer*4 values(8)
      integer :: status(MPI_STATUS_SIZE)
      common tgrade,amu
      
      external fhi,fhit
      f(x) = (amu-1)/(t0)/(1+x/t0)**amu


      call MPI_Init(ierr) 
      call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
      call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)

      do irun = 1,6
      idum = -17 - 1000*rank   ! Different seed for each rank
!      idum = -53
      if(irun==1) then
         amu = 1.50d0
         t0 = 1.d0
      endif
      if(irun==2) then
         amu = 1.50d0
         t0 = 20.d0
      endif
      if(irun==3) then
         amu = 2.50d0
         t0 = 1.d0
      endif
      if(irun==4) then
         amu = 2.50d0
         t0 = 20.d0
      endif
      if(irun==5) then
         amu = 3.50d0
         t0 = 1.d0
      endif
      if(irun==6) then
         amu = 3.50d0
         t0 = 20.d0
      endif    
      tgrade = t0
      root = 0
      idt0 = nint(t0/10)
      iut0 = nint(t0-idt0*10)
      ift0 = nint(10*(t0-iut0-10*idt0))
      ict0 = nint(100*(t0-iut0-10*idt0-0.1*ift0))
!      print *,idt0,iut0,ift0
      filename = "free-muxxxx-tz......out"
      write(filename(8:11),'(f4.2)') amu
      write(filename(15:15),'(I1)') idt0
      write(filename(16:16),'(I1)') iut0
      write(filename(18:18),'(I1)') ift0   
      write(filename(19:19),'(I1)') ict0    
!      print *,filename


      rho_local = 0.d0
      rho_global = 0.d0
      dr = 1.d0   ! rho bin size
      n_corr = 10
      alfa = -1.d0/(amu-1)    ! 1.0 < amu < 4.0 amu exponent in the powerlaw



      nave = 4000000 ! 0

      do iave = 1,nave
!         print *,"here2",iave
         state = rnor(idum)
!         if (ran2(idum).lt.0.5d0) then
!            state = 1.d0
!         else
!            state = -1.d0
!         endif
!          if (mod(iave,1000).eq.0) print *,iave
         dt = 0.1d0
         tmax = dt*nd*ninner
         y = 0.d0
         rho_local(0,0) = rho_local(0,0) + 1.d0
         
         time = 0.d0     
         x = power(idum,alfa,t0)
         tswitch = x

         do itime = 1,nd
            do inner = 1,ninner
!               if(iave.eq.1) write(999,*) time,y
               time = time + dt

               do while (time .gt. tswitch)
                  x = power(idum,alfa,t0)
                  tswitch = tswitch+x
                  state = rnor(idum)
!                  if (ran2(idum).lt.0.5d0) then
!                     state = 1.d0
!                    else
!                     state = -1.d0
!                  endif
!                  print *,'inside while ',time, tswitch, x
                  y = y + state
               enddo
            enddo
            ir = nint(y/dr)
            if(ir.ge.-nr.and.ir.le.nr) rho_local(itime,ir) = rho_local(itime,ir) + 1.d0
         enddo
      enddo
!      print *,"finito il loop!",rank,rho_local(0,0)
      call MPI_REDUCE(rho_local, rho_global,(2*nr+1)*(nd+1), MPI_DOUBLE_PRECISION, MPI_SUM, root, MPI_COMM_WORLD, ierr)
      if (rank==root) then
      open(10,file=filename)
      do itime = 0,nd
         do ir = -nr,nr
            write(10,*) itime*dt*ninner,dr*ir,rho_global(itime,ir)/nave/size/dr
         enddo
         write(10,*)
      enddo
      close(10)
      endif
      enddo ! Loop on irun
      call MPI_Finalize(ierr)
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
      
!      function power2(idum,alfa,t)
!      implicit real*8(a-h,p-z)
! Routine to generate according to the distribution
! (a-1)*T^(a-1)/(T+t)^a
! The argument of power is *defined*  alfa=1/(a-1)
            
!      ta = alfa*expdev(idum)
!      power2 = t*(dexp(ta) - 1)
      
!      return
!      end 

      function power(idum,alpha,t0)
      implicit real*8(a-h,p-z)
! Routine to generate according to the distribution
! (a-1)*T^(a-1)/(T+t)^a
! Set the argument alpha alpha=-1/(a-1) in the main
!         y = ran3(idum)
!         x = y**(-1.d0/(amu-1.d0))-1.d0
!         j = int(t0*x*10.d0)
!      y = ran3(idum)
!      y = rcarry(idum)
      y = ran2(idum)
      power = t0*(y**alpha-1.d0)
      return
      end 
      
      FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL*8 ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1.d0/IM1,IMM1=IM1-1,&
      IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,&
      IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2d-14,RNMX=1.d0-EPS)
!    Long period (> 2 × 1018) random number generator of L’Ecuyer with Bays-Durham shuffle
!    and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive
!    of the endpoint values). Call with idum a negative integer to initialize; thereafter, do not
!    alter idum between successive deviates in a sequence. RNMX should approximate the largest
!    floating value that is less than 1.
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then      !Initialize.
        idum=max(-idum,1)     ! Be sure to prevent idum = 0.
        idum2=idum
        do j=NTAB+8,1,-1 ! Load the shuffle table (after 8 warm-ups).
           k=idum/IQ1
           idum=IA1*(idum-k*IQ1)-k*IR1
           if (idum.lt.0) idum=idum+IM1
           if (j.le.NTAB) iv(j)=idum
        enddo
        iy=iv(1)
      endif
      k=idum/IQ1 ! Start here when not initializing.
      idum=IA1*(idum-k*IQ1)-k*IR1 ! Compute idum=mod(IA1*idum,IM1) without over-flows by Schrage’s method.
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2 ! Compute idum2=mod(IA2*idum2,IM2) likewise.
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV ! Will be in the range 1:NTAB.
      iy=iv(j)-idum2 ! Here idum is shuffled, idum and idum2 are combined to generate output.
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX) ! Because users don’t expect endpoint values.
      return
      END

      function rnor(idum)
      implicit real*8(a-h,p-z)
      parameter (twopi = 8.d0*datan(1.d0))
      SAVE iv,r,t,u1,u2
      data iv/0/

      if (iv.eq.0.or.idum.le.0) then
        u1 = ran2(idum)
        u2 = ran2(idum)
        r = dsqrt(-2.d0*dlog(u1))
        iv = 1
        t = twopi*u2
        rnor = r*dcos(t)
       else
        rnor = r*dsin(t)
        iv = 0
      endif
      return
      end


    
