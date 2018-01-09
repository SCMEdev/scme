module compressed_foolin !>>>

use printer_mod, only: str, printer, printo
use compressed_arrays 
implicit none

contains !//////////////////////////////////////////////////////////////






subroutine print_trace_ind 
    integer i, pos2(2), k!tempa(1000), 
    pos2=1
    k = 7
    
    do i = 0, k, 2
        print*, len00(i), len00(i+1), pos2-1, pos00(i)+1, pos00(i+1)+1
        pos2(1) = pos2(1) + len00(i+0)
        pos2(2) = pos2(2) + len00(i+1)
        
    enddo
    print*, pos00
    print*, len00
    print*, (i,i=0,k)
end subroutine


subroutine main; ! MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN MAIN
call &
!print_trace_ind 
stupido
end



subroutine stupido
integer p1(10), p2(10), p3(10)
integer p10(0:10), p20(0:10), p30(0:10)
integer p11(11), p21(11), p31(11)
p1 =1; p2 =4; p3 = 0
p10=1; p20=4; p30 = 0
p11=1; p21=4; p31 = 0
call stupido2(p1,p2,p3) 
print*, p3
end


subroutine stupido2(p1,p2,p3)
    integer, parameter :: start=0, fin=10
    integer p1(start:fin),p2(start:fin),p3(start:fin), i
    do i = start,fin 
        p3(i) = p1(i)+p2(i)
    enddo
end


end module !<<< 
