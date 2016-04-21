SUBROUTINE Lax_Wendroff(Un,M,delta_x,delta_t,a)
 IMPLICIT NONE
!*********************introduction********************
!this is the Lax_Wendroff scheme aiming at equation **
!partial_u/partial_t+a*partial_u/partial_x=0, and the*
!the solution in (n+1) time level can be got. in here*
!,in order to solve the migration equation and observe
!the movement of pulse wave over time, the periodic **
!conditions are taken*********************************
!*****************************************************
 integer::I,J,M
 real(kind=8)::Un(M),Un_temp(0:(M+1))
 real(kind=8)::delta_x,delta_t,a
!
 do I=1,M
 	Un_temp(I)=Un(I)
 end do
 Un_temp(0)=Un_temp(M)
 Un_temp(M+1)=Un_temp(1)
!
 do I=1,M
	Un(I)=Un_temp(I)-a*delta_t*(Un_temp(I+1)-Un_temp(I-1))/(2.0d0*delta_x)&
		&+delta_t**2*a**2*(Un_temp(I+1)-2.0d0*Un_temp(I)+Un_temp(I-1))/(2*delta_x**2)
 end do
END SUBROUTINE Lax_Wendroff
