SUBROUTINE DRP7(U,Ux,M,delta_x)
 IMPLICIT NONE
!****************introduction***********************
!this subroutine can get the partial(u)/partial(x)**
!using the DRP scheme*******************************
!***************************************************
 INTEGER::M,I,J
 REAL(KIND=8)::U(M),Ux,delta_x
 REAL(KIND=8)::a(-3:3),b(0:6),c(-1:5),d(-2:4)
!the coefficient of the DRP scheme
!a: coefficient of the central difference 
!b: coeeficient of the zero-point forward and -b stands for the coefficient of backward difference
!c: coefficient of the one-point forward and -c stands for the coefficient of backward difference
!d: coefficient of the two-point forward and -d stands for the coefficient of backward difference
!M=7
 a(-1)=-0.77088238051822552d0
 a(-2)=0.166705904414580469d0
 a(-3)=-0.02084314277031176d0
 a(0)=0.0d0
 a(1)=0.77088238051822552d0
 a(2)=-0.166705904414580469d0
 a(3)=0.02084314277031176d0
 !b(0)=-2.19228033900d0
 !b(1)=4.74861440100d0
 !b(2)=-5.10885191500d0
 !b(3)=4.46156710400d0
 !b(4)=-2.83349874100d0
 !b(5)=1.12832886100d0
 !b(6)=-0.20387637100d0
 !c(-1)=-0.20933762200d0
 !c(0)=-1.08487567600d0
 !c(1)=2.14777605000d0
 !c(2)=-1.38892832200d0
 !c(3)=0.76894976600d0
 !c(4)=-0.28181465000d0
 !c(5)=0.048230454000d0
 !d(-2)=0.049041958000d0
 !d(-1)=-0.4688403500d0
 !d(0)=-0.47476091400d0
 !d(1)=1.27327473700d0
 !d(2)=-0.51848452600d0
 !d(3)=0.16613853300d0
 !d(4)=-0.026369431000d0
!
 Ux=0.0d0
 do J=1,M
	Ux=Ux+a(J-(int(M/2)+1))*U(J)
 end do
 Ux=Ux/delta_x
END SUBROUTINE DRP7
