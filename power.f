!      function power2(idum,alfa,t)
!      implicit real*8(a-h,p-z)
C Routine to generate according to the distribution
C (a-1)*T^(a-1)/(T+t)^a
C The argument of power is *defined*  alfa=1/(a-1)
            
!      ta = alfa*expdev(idum)
!      power2 = t*(dexp(ta) - 1)
      
!      return
!      end 

      function power(idum,alpha,t0)
      implicit real*8(a-h,p-z)
C Routine to generate according to the distribution
C (a-1)*T^(a-1)/(T+t)^a
C Set the argument alpha alpha=-1/(a-1) in the main
C         y = ran3(idum)
C         x = y**(-1.d0/(amu-1.d0))-1.d0
C         j = int(t0*x*10.d0)
C      y = ran3(idum)
c      y = rcarry(idum)
      y = ran2(idum)
      power = t0*(y**alpha-1.d0)
      return
      end 