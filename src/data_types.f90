module data_types
  
  implicit none
  
  private
  public sp, dp, medint, longint
  
  integer, parameter :: sp = selected_real_kind(6,37) ! single precision, ~ 4 byte float
  integer, parameter :: dp = selected_real_kind(15, 307) ! double precision ~ 8 byte float
  integer, parameter :: medint = selected_int_kind(9) ! double precision ~ 4 byte float
  integer, parameter :: longint = selected_int_kind(15) ! 15 digits, ~ 8 byte integer
  
end module data_types
