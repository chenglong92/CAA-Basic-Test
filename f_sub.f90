SUBROUTINE f(U,Ux,M,delta_x)
 IMPLICIT NONE
 integer::I,J,M
 real(kind=8)::U(M),Ux(M),U_temp(-2:M+3),DRP_input(7)
 real(kind=8)::delta_x
!***********************introduction******************************
!this program can help you get the right valueof time forwarding**
!method(LDDRK), and here use the periodic boundary condition
!*****************************************************************
!
!
!initialize the U_temp
!here use the periodic boundary conditions
 do I=1,M
 	U_temp(I)=U(I)
 end do
 do I=-2,0,1
 	U_temp(I)=U(M+I)
 	U_temp(M+I+3)=U(I+3)
 end do
!
 do I=1,M
 	DRP_input(1:7)=U_temp((I-3):(I+3))
 	call DRP7(DRP_input,Ux(I),7,delta_x)
 end do
END SUBROUTINE f
