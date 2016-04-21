SUBROUTINE LDDRK(Un,M,delta_t,delta_x)
 IMPLICIT NONE
!***********************introduction****************************
!consider the basic equation form such as f=partial_u/partial_x*
!this program can get the value for some section of LDDRK*******
!***************************************************************
 integer::I,J,M
 real(kind=8)::Un(M)
 real(kind=8)::delta_t,delta_x
 real(kind=8)::k1(M),k2(M),k3(M),k4(M),k_temp(M)
 real(kind=8)::w(4),beita(2:4)
!
!the coefficient of LDDRK
 w(1)=0.1630296d0
 w(4)=w(1)
 w(2)=0.348012d0
 w(3)=0.3259288d0
 beita(2)=0.5d0
 beita(3)=0.5d0
 beita(4)=1.0d0
!
!get the intermediate variates for LDDRK
 call f(Un,K1,M,delta_x)
 k1=-k1*delta_t
 k_temp=Un+beita(2)*K1
 call f(k_temp,k2,M,delta_x)
 k2=-k2*delta_t
 k_temp=(Un+beita(3)*k2)
 call f(k_temp,k3,M,delta_x)
 k3=-k3*delta_t
 k_temp=(Un+beita(4)*k3)
 call f(k_temp,k4,M,delta_x)
 k4=-k4*delta_t
!
 do I=1,M
 	Un(I)=w(1)*k1(I)+w(2)*k2(I)+w(3)*k3(I)+w(4)*k4(I)
 end do
END SUBROUTINE LDDRK
