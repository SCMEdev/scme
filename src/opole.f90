module opole
use printer_mod, only:printer
use data_types, only: dp, pi

implicit none

!call main()

contains !/////////////////



subroutine octupole_tensor(cec,cer2,charges,oct)
    integer, parameter   :: xyz=3,hho=3
    real(dp), intent(in) :: cec(xyz,hho),cer2(hho),charges(hho)
    real(dp), intent(out) :: oct(xyz,xyz,xyz)
    ! internal:
    real(dp) ch05
    integer i
    
    ! _linearly_indipendent_ octupole components added up from each charge site
    oct=0
    do i = 1,hho 
      ch05 = charges(i)*0.5_dp
      oct(1,1,1) = oct(1,1,1) + ch05*( 5_dp*cec(1,i)**3 - cer2(i)*(3*cec(1,i)) )
      oct(2,2,2) = oct(2,2,2) + ch05*( 5_dp*cec(2,i)**3 - cer2(i)*(3*cec(2,i)) )
      oct(3,3,3) = oct(3,3,3) + ch05*( 5_dp*cec(3,i)**3 - cer2(i)*(3*cec(3,i)) )
      oct(1,1,2) = oct(1,1,2) + ch05*( 5_dp*cec(1,i)**2*cec(2,i) - cer2(i)*(cec(2,i)) )
      oct(1,1,3) = oct(1,1,3) + ch05*( 5_dp*cec(1,i)**2*cec(3,i) - cer2(i)*(cec(3,i)) )
      oct(1,2,2) = oct(1,2,2) + ch05*( 5_dp*cec(1,i)*cec(2,i)**2 - cer2(i)*(cec(1,i)) )
      oct(1,2,3) = oct(1,2,3) + ch05*( 5_dp*cec(1,i)*cec(2,i)*cec(3,i) - 0.0_dp )
    enddo
    
    ! remaining _unique_ components from traceless condition
    oct(1,3,3) = -( oct(1,1,1) + oct(1,2,2) )  
    oct(2,3,3) = -( oct(1,1,2) + oct(2,2,2) )
    oct(2,2,3) = -( oct(1,1,3) + oct(3,3,3) )
    
    ! symmetrize 
    
    ! a,a,b:
    oct(1,2,1) = oct(1,1,2)
    oct(2,1,1) = oct(1,1,2)
    oct(1,3,1) = oct(1,1,3)
    oct(3,1,1) = oct(1,1,3)
    oct(2,1,2) = oct(1,2,2)
    oct(2,2,1) = oct(1,2,2)
    oct(3,1,3) = oct(1,3,3)
    oct(3,3,1) = oct(1,3,3)
    oct(3,2,3) = oct(2,3,3)
    oct(3,3,2) = oct(2,3,3)
    oct(2,3,2) = oct(2,2,3)
    oct(3,2,2) = oct(2,2,3)
    
    ! a,b,c:
    oct(1,3,2) = oct(1,2,3)
    oct(2,1,3) = oct(1,2,3)
    oct(2,3,1) = oct(1,2,3)
    oct(3,1,2) = oct(1,2,3)
    oct(3,2,1) = oct(1,2,3)
    

    
    ! o0(3,3,3) - oct(3,3,3) = ~1e-5 since o0 is not perfectrly traceless. 
    
end subroutine

end module
