!=======================================================================
! PROGRAMMA MOMGEN4_MPI (Modern Fortran with Kind Parameter and MPI)
!=======================================================================
      program momgen4_mpi
      
      use mpi
      use precision_kind_mod, ONLY: dp
      use dgauss_mod, ONLY: dgauss
      use random_numbers_mod, only: ran2, power !, dp ! Fix: Import ran2, power, and dp from the module
      
      ! --- IMPLICIT NONE MUST BE THE FIRST EXECUTABLE STATEMENT ---
      implicit none
      

! Code to compute the four times correlation function    
      integer, parameter :: np=12,nd=2**np,nr=1000,nf=18,nde=16
      integer, parameter :: ncut=nint(0.005_dp*nd)
      integer, parameter :: ng = 10000
      real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp),sq2=sqrt(2.0_dp)

      ! Dichiarazioni degli array e variabili 'dp'
      real(dp) :: xin(nd)
      real(dp) :: rho(0:nr)
      real(dp) :: pn(-ng:ng)
      
      ! Array locali (xfr_loc, xf2_loc) e globali (xfr_glo, xf2_glo)
      real(dp) :: xfr_loc(nd,nf), xfr_glo(nd,nf)
      real(dp) :: xf2_loc(nd,3), xf2_glo(nd,3)
      real(dp) :: xfr2_loc(nd,nf), xfr2_glo(nd,nf)
      
      real(dp) :: tgrade,amu,t0
      common /params/ tgrade,amu ! t0 is not in the original common block structure
      
      character(len=120) :: cwd
      character(len=8) :: date
      character(len=10) :: ctime
      character(len=5) :: zone
      integer :: values(8)
      
            
      ! Dichiarazione dei nomi delle statement functions
      real(dp) :: f, ff

      ! Statement functions (internal scope, access host variables)
      f(x) = (amu-1.0_dp)/(t0)/(1.0_dp+x/t0)**amu
      ff(x) = log((x*x+sq2*x+1.0_dp)/(x*x-sq2*x+1.0_dp))/(2.0_dp*pi)&
     +        atan((1.0_dp+sq2*x)/(sq2-x*x))/pi
     
      ! Variabili MPI e INTEGER
      integer :: ierr, rank, size, root, status(MPI_STATUS_SIZE)
      integer :: nave, nave_local, nave_offset, iave_start, iave_end
      integer :: i, j, k, ii1, i1, i2, i3, i4, io, im, il, iu, ish
      integer :: nswitch, n_corr
      integer :: kdum, jdum, ndum
      integer :: idum, ist 
      integer :: iave, n4 
      INTEGER :: host_rank_id


      ! Variabili di lavoro (con 'dp' kind)
      real(dp) :: x, xmin, xc, dpn, dt, tmax, tswitch, time, state
      real(dp) :: alfa, r2, r3
      real(dp) :: contribution 
      
      ! Inizializzazione MPI
      call MPI_Init(ierr)

      call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
      call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)
      
      
      root = 0
      
      host_rank_id = get_host_id() 
!      WRITE(*,*) 'This host ID is:', host_rank_id        
      call date_and_time( date, ctime, zone, values )
      
      idum = -17 - 100*values(3) - 2*host_rank_id 
!      print *,"cucu",rank,root,size
      
      if (rank .eq. root) then
         print *,"IDUM",idum,size,values(7)
      endif
      
      idum = idum - 2*rank

      amu = 1.50_dp
      t0 = 20.0_dp
      tgrade = t0

      ! Parametri per la generazione della PDF
      dpn = 1.0e-2_dp
      xmin = -dpn*ng
      
      ! NOTE: fpn is defined below CONTAINS, so it's visible here.
      pn(-ng) = dgauss(fpn,100.0_dp*xmin,xmin,1.0e-12_dp) 
      do i=-ng,ng-1
         x = i*dpn
         pn(i+1) = pn(i) + dgauss(fpn,x,x+dpn,1.0e-12_dp)
      enddo
      
      rho = 0.0_dp
      n_corr = 10
      alfa = -1.0_dp/(amu-1.0_dp)

      xfr_loc = 0.0_dp
      xfr2_loc = 0.0_dp
      xf2_loc = 0.0_dp
 
      ! --- Gestione della media (IAVE) con MPI ---
      nave = 50000 !0 ! 0    ! 25000 !0! 0 !0     40000 ! 40000000
      
      if (rank .eq. root) then
         print *, "Total runs (nave): ", nave
      endif
      
      ! Inizializzazione del generatore di numeri casuali
      call date_and_time( date, ctime, zone, values )
      do ish = 1, values(7) !+ rank
         x = ran2(idum) 
      enddo
      
      ! Loop di simulazione (IAVE)
      do iave = 1,nave
!         if (mod(iave,1000).eq.0 .and. rank .eq. root) print *,"IAVE",iave,"SIZE",size
         call date_and_time( date, ctime, zone, values )   
         do ish = 1, values(7) + rank
            x = ran2(idum) 
         enddo      
         ! --- Generazione dello stato iniziale (state) ---
         xc = ran2(idum)
         il = -ng
         iu = ng
         do j = 1,17
            im = (il+iu)/2
            if (pn(im).lt.xc) then
               il = im
            else
               iu = im
            endif
         enddo
         state = dpn*(real(il,dp)*(pn(iu)-xc)+real(iu,dp)*(xc-pn(il)))/(pn(iu)-pn(il)) 
         ! -----------------------------------------------
         
         io = 1
         dt = 0.1_dp
         tmax = dt*nd
         x = power(idum,alfa,t0)
         tswitch = x
         
         ! Aggiornamento rho (statistiche dei tempi di attesa)
         j = int(x*10.0_dp)
         if (j.le.nr.and.j.ge.0) then 
            rho(j) = rho(j) + 1.0_dp
         endif
      
         time = 0.0_dp
         xin = 0.0_dp
         nswitch = 0
       
         do while (time .le. tmax)
            xin(io) = state
            time = time + dt
            io = io + 1

            if (time .gt. tswitch) then
               x = power(idum,alfa,t0)
               nswitch = nswitch + 1
          
               ! Aggiornamento rho
               j = int(x*10.0_dp)
               if(j.le.nr.and.j.ge.0) then 
                  rho(j) = rho(j) + 1.0_dp
               endif
               tswitch = tswitch+x
           
               ! --- Generazione del nuovo stato (state) ---
               xc = ran2(idum)
               il = -ng
               iu = ng
               do j = 1,17
                  im = (il+iu)/2
                  if (pn(im).lt.xc) then
                     il = im
                  else
                     iu = im
                  endif
               enddo
               state = dpn*(real(il,dp)*(pn(iu)-xc)+real(iu,dp)*(xc-pn(il)))/(pn(iu)-pn(il))
               ! -------------------------------------------
            endif
         enddo

         ! Calcolo delle correlazioni locali per questa simulazione
         n4 = nd/8
         do k = 1,nf
            if (k.le.6) then
               i1 = 1
            else if (k.le.12) then
               i1 = n4
            else
               i1 = 2*n4
            endif
      
            if(mod(k,6).eq.1) then
               r2 = 0.25_dp
               r3 = 0.75_dp
            else if(mod(k,6).eq.2) then
               r2 = 0.25_dp
               r3 = 2.0_dp/3.0_dp
            else if(mod(k,6).eq.3) then
               r2 = 0.25_dp
               r3 = 0.5_dp
            else if(mod(k,6).eq.4) then
               r2 = 0.5_dp
               r3 = 2.0_dp/3.0_dp
            else if(mod(k,6).eq.5) then
               r2 = 1.0_dp/3.0_dp
               r3 = 2.0_dp/3.0_dp
            else if(mod(k,6).eq.0) then
               r2 = 0.25_dp
               r3 = 2.0_dp/3.0_dp
            end if
            
            do i4 = i1+ncut,nd
               i3 = nint(real(i1,dp)+r3*real(i4-i1,dp)) 
               i2 = nint(real(i1,dp)+r2*real(i4-i1,dp))
   
               ! Aggiornamento del contributo locale (xfr_loc)
               contribution = xin(i1)*xin(i2)*xin(i3)*xin(i4) 
               xfr_loc(i4,k) = xfr_loc(i4,k) + contribution
               xfr2_loc(i4,k) = xfr2_loc(i4,k) + contribution*contribution
               
               ! Aggiornamento del contributo locale (xf2_loc)
               if (mod(k,6).eq.1) then
                  ii1 = k/6+1
                  xf2_loc(i4,ii1) = xf2_loc(i4,ii1) + xin(i1)*xin(i4)
               endif
            enddo
         enddo
      enddo
      
      ! --- Riduzione MPI: Somma dei risultati locali ---
      
      ! Somma di xfr_loc in xfr_glo (solo sul root)
      call MPI_Reduce(xfr_loc, xfr_glo, nd*nf, MPI_DOUBLE_PRECISION,&
                      MPI_SUM, root, MPI_COMM_WORLD, ierr)

      ! Somma di xfr2_loc in xfr2_glo (solo sul root)
      call MPI_Reduce(xfr2_loc, xfr2_glo, nd*nf, MPI_DOUBLE_PRECISION,&
                      MPI_SUM, root, MPI_COMM_WORLD, ierr)

      ! Somma di xf2_loc in xf2_glo (solo sul root)
      call MPI_Reduce(xf2_loc, xf2_glo, nd*3, MPI_DOUBLE_PRECISION,&
                      MPI_SUM, root, MPI_COMM_WORLD, ierr)

      ! --- Output dei risultati (solo sul core root) ---
      if (rank .eq. root) then
         
         ! Media: dividiamo la somma globale per il numero totale di simulazioni (nave * size)
         
         do k=1,nf
            do i4=1,nd
               xfr_glo(i4,k) = xfr_glo(i4,k) / (real(nave, dp)*real(size,dp))
            enddo
         enddo
         do k=1,nf
            do i4=1,nd
               ! Calcolo della varianza: <x^2> - <x>^2
               xfr2_glo(i4,k) = xfr2_glo(i4,k) / (real(nave, dp)*real(size,dp)) - xfr_glo(i4,k)*xfr_glo(i4,k)
            enddo
         enddo         
         do k=1,3
            do i4=1,nd
               xf2_glo(i4,k) = xf2_glo(i4,k) / (real(nave, dp)*real(size,dp))
            enddo
         enddo

!         dt = 0.1_dp
         
         ! Scrittura di xfr (quattro tempi) con varianza (xfr2_glo)
         do k = 1,nf
            if (k.le.6) then
               i1 = 1
            else if (k.le.12) then
               i1 = n4
            else
               i1 = 2*n4
            endif
      
            do i4 = i1+ncut,nd
               write(100+k,*) real(i4-1,dp)*dt,xfr_glo(i4,k),&
                sqrt(xfr2_glo(i4,k)/(real(nave, dp)*real(size,dp)-1)),&
                sqrt(xfr2_glo(i4,k))
            enddo
            close(100+k)
         enddo
         
         ! Scrittura di xf2 (due tempi)
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
               write(500+ii1,*)real(i4-1,dp)*dt,xf2_glo(i4,ii1)
            enddo
            close(500+ii1)
          enddo
      endif

      ! Finalizzazione MPI
      call MPI_Finalize(ierr)
      
      ! =============================================================
      ! CONTAINS: Procedure Definitions
      ! =============================================================
      contains
      
      ! FUNZIONI INTERNE (ereditano 'dp', 'tgrade', 'amu', 't0' dal main program)

      PURE FUNCTION fhi(x) RESULT(res)
      USE precision_kind_mod, ONLY: dp
      implicit none
      real(dp), intent(in) :: x
      REAL(dp) :: res
      REAL(dp) :: t0,amu
      common /params/ t0,amu
      res = (amu-1.0_dp)/(t0)/(1.0_dp+x/t0)**amu
      return
      end function fhi

      PURE FUNCTION fhit(x) RESULT(res)
      USE precision_kind_mod, ONLY: dp
      implicit none
      real(dp), intent(in) :: x
      REAL(dp) :: res
      REAL(dp) :: t0,amu
      common /params/ t0,amu
      res = x*(amu-1.0_dp)/(t0)/(1.0_dp+x/t0)**amu
      return
      end function fhit

      PURE FUNCTION fpn(x) RESULT(res)
      USE precision_kind_mod, ONLY: dp
      IMPLICIT NONE
      REAL(dp), INTENT(IN) :: x
      REAL(dp) :: res
      real(dp), parameter :: pi=4.0_dp*atan(1.0_dp),sq2=sqrt(2.0_dp)
      real(dp), parameter :: cost=sq2/pi
      res = cost/(1.0_dp+x*x*x*x)
      return
      end function fpn
      
    INTEGER FUNCTION get_host_id()

    USE, INTRINSIC :: ISO_FORTRAN_ENV

    IMPLICIT NONE
    
    CHARACTER(LEN=128) :: hostname
    INTEGER :: status, i
    INTEGER(int32) :: host_hash ! Usiamo INT32 per garantire che il bit-masking sia corretto
    
    ! 1. Ottieni l'hostname
    CALL GET_ENVIRONMENT_VARIABLE('HOSTNAME', hostname, status=status)
    IF (status /= 0) THEN
      ! Fallback su 'USER' se HOSTNAME non è definito
      CALL GET_ENVIRONMENT_VARIABLE('USER', hostname, status=status)
    END IF

    ! 2. Algoritmo di Hashing Semplice
    host_hash = 0_INT32
    DO i = 1, LEN_TRIM(hostname)
      ! Hash = (Hash * 2) + ASCII_value (ISHFT è uno shift binario a sinistra)
      host_hash = host_hash + IACHAR(hostname(i:i))
      
      ! Questa linea è quella che generava l'errore, ora risolto con l'uso di ISO_FORTRAN_ENV
      ! 2147483647 è l'equivalente di 2**31 - 1, per mantenere il valore positivo.
      host_hash = iand(ISHFT(host_hash, 1), 2147483647_INT32) 
    END DO
    
    ! Restituisci un valore hash positivo
    get_host_id = ABS(host_hash) + 1
    
  END FUNCTION get_host_id

      end program momgen4_mpi
