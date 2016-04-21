PROGRAM EX01
!*************************introduction*************************
!this program can solve the one-dimensional migration equation*
!with one Gauss pulse wave,And the spatial discretization uses the DRP scheme and Lax-
!Wendroff Scheme, temporal discretization use the LDDRK scheme*
!********************************************************************* 
 IMPLICIT NONE
 INTEGER::I,J,k,M,N,t_step,X_range,y_range
 real(kind=8)::delta_x,delta_t
 parameter(X_range=-20,Y_range=450,delta_x=1,delta_t=0.8)
 parameter(M=(int((Y_range-X_range)/delta_x)+1),N=3,t_step=int(400/delta_t))
 REAL(KIND=8)::U(M),U1(M),U_temp(M),U_standard(M)
 REAL(KIND=8)::a(2,N)
 REAL(KIND=8)::t,error(2)
!here M stands for the number of discrete point
!N stands for the number of point stencil
!a(2,N) stands for the coefficient matrix
!a(1,:) stands for the standard 7-point scheme
!a(2,:) stands for the 7-point optimized scheme
!U(M): the LDDRK and DRP scheme
!U1(M): the Lax-Wendroff scheme
!U_standard(M): the theoritical solution
 t=0.0d0
 a(1,1)=0.79926643d0
 a(1,2)=-0.18941314d0
 a(1,3)=0.02651995d0
 a(2,1)=0.77088238051822552d0
 a(2,2)=-0.166705904414580469d0
 a(2,3)=0.02084314277031176d0
!initialize the velocity distribution
 DO I=1,M
 	U(I)=0.5d0*exp(-(log(2.0d0))*((-20+(I-1)*delta_x)*1.0d0/3.0d0)**2)
	U1(I)=0.5d0*exp(-(log(2.0d0))*((-20+(I-1)*delta_x)*1.0d0/3.0d0)**2)
 END DO
 DO I=1,M
 	U_standard(I)=0.5d0*exp(-(log(2.0d0))*((-20+(I-1)*delta_x-t)*1.0d0/3.0d0)**2)
 END DO
!
 OPEN(UNIT=10,file="velocity.dat")
 OPEN(UNIT=11,file="error.dat")
 DO I=1,M
 	write(unit=10,fmt="(5f20.10)")U(I),U1(I),U_standard(I),(U(I)-U_standard(I)),(U1(I)-U_standard(I))
! write(unit=10,fmt="(501f18.10)")U_standard(1,:)
 END DO
!use LDDRK scheme and DRP scheme to get the solution
 DO J=1,t_step
!
!the LDDRK scheme and DRP scheme
	do I=1,M
 		U_temp(I)=U(I)
	end do
	call LDDRK(U_temp,M,delta_t,delta_x)
	t=t+delta_t
	do I=1,M
 		U(I)=U(I)+U_temp(I)
	end do
!the Lax-Wendroff scheme
	call Lax_Wendroff(U1,M,delta_x,delta_t,1.0d0)
!
 	DO I=1,M
 		U_standard(I)=0.5d0*exp(-(log(2.0d0))*((-20+(I-1)*delta_x-t)*1.0d0/3.0d0)**2)
 	END DO
	if(mod(J,t_step/4)==0)then
		DO I=1,M
 			write(unit=10,fmt="(5f20.10)")U(I),U1(I),U_standard(I),(U(I)-U_standard(I)),(U1(I)-U_standard(I))
 		END DO
 	end if
	error=0.0d0
 	do I=1,M
 		error(1)=error(1)+(U(I)-U_standard(I))**2
 		error(2)=error(2)+(U1(I)-U_standard(I))**2
 	end do
 	write(unit=11,fmt="(2f20.10)")error(1),error(2)
 END DO
 CLOSE(UNIT=10)
 CLOSE(UNIT=11)
END PROGRAM EX01
