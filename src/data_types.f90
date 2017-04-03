module data_types
  
  implicit none
  
  private
  public sp, dp, medint, longint, h2o
  
  integer, parameter :: sp = selected_real_kind(6,37) ! single precision, ~ 4 byte float
  integer, parameter :: dp = selected_real_kind(15, 307) ! double precision ~ 8 byte float
  integer, parameter :: medint = selected_int_kind(9) ! double precision ~ 4 byte float
  integer, parameter :: longint = selected_int_kind(15) ! 15 digits, ~ 8 byte integer

!position type
  type h2o
   real(dp), dimension(3) :: h1, h2, o
  end type

!force type  
!  type h2o_forces
!    real(dp), dimension(3) :: h1, h2, o
!  endtype
  
end module data_types
