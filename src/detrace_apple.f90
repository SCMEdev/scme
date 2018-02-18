module detrace_apple
! All routines herein are taken from Jon Applequist public repository. 
! I could see no license so I assume they are in the public domain. 
! the filenames coincide with the name of the routines. 
!   xtrace() removes the trace of a compressed tensor up to rank 10, and relies on the remaining routines. 
!   Jonatan Öström
private 
public detrace_a, ff, main, opdetr
contains !//

subroutine main
    implicit none
!    call io4
    print*, ff(0)
    
end subroutine

pure function ff(base)
    implicit none
    integer, intent(in) :: base
    integer i, iff
    real*8 ff
    iff = 1
    do i = base, 1, -2
      iff = iff*i
    enddo
    !(iff = iff*i, i = base, 1, -2) 
    ff = dble(iff)
endfunction

subroutine io4
    implicit none
    integer i
    integer, parameter :: rank = 4
    integer, parameter :: le=(rank+1)*(rank+2)/2
    real*8 inten(le),outen(le) 
    
    read(*,'(2f24.16)') (inten(i), i = 1,le)
    
    outen = detrace_a(inten,4)
    
    print'(2(f24.16))', (outen(i) / ff(2*rank-1), inten(i),i = 1,le)
    
end subroutine

function detrace_a(a,k) result(b) !wrapper function
    implicit none
    integer k
    real*8, dimension( ( k+1)*(k+2)/2 ) :: a, b
    call xtrace(a,b,k)
end function

function opdetr(a,k) result(b) !wrapper function
    implicit none
    integer k
    real*8, dimension( ( k+1)*(k+2)/2 ) :: a, b
    call xtrace(a,b,k)
    b = b/ff(2*k-1)
end function


! Original Jon Applequist Routines /////////////////////////////////////////////////

      SUBROUTINE XTRACE(A,B,K)
!
!  PURPOSE: CALCULATES THE TENSOR WHICH IS TRACELESS IN ALL INDEX
!    PAIRS FROM A GIVEN TOTALLY SYMMETRIC TENSOR.
!
!  DEFINITIONS:
!
!    *A    TOTALLY SYMMETRIC TENSOR IN COMPRESSED FORM.  DIMENSION
!          SPECIFIED BY CALLING PROGRAM MUST BE AT LEAST (K+1)*(K+2)/2.
!     B    TRACELESS TENSOR CALCULATED FROM A.  DIMENSION SPECIFIED
!          BY CALLING PROGRAM MUST BE AT LEAST (K+1)*(K+2)/2.
!    *K    RANK OF A AND B.  (K.LE.10)
!
!  SUBROUTINES REQUIRED: INDX4, INDX6, TRACE
!
!  PROGRAMMED BY J. APPLEQUIST, 6/30/83.
!
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(1), B(1), IA(10), FACT(6), L(5,2)
!
!         STORE COEFFICIENTS OF TRACES.
!
      MLIM=K/2 +1
      DO 15 M1=1,MLIM
      LIM=K-M1+1
      NFACT=1
      DO 10 LL=1,LIM
   10 NFACT=NFACT*(2*LL-1)
      IF (M1/2*2.EQ.M1) NFACT=-NFACT
   15 FACT(M1)=NFACT
      NLIM=(K+2)*(K+1)/2
!
!         SCAN COMPRESSED ARRAY OF COMPONENT INDICES.
!
      DO 100 N=1,NLIM
      CALL INDX4(K,N,IA,IER)
      CALL INDX6(IA,K,IN1,IN2,IN3)
      KOL1=IN3+1
      KOL2=KOL1+IN2
!
!         CALCULATE TERM WITH 0 DELTAS.
!
      B(N)=FACT(1)*A(N)
!
!       CALCULATE TERMS WITH 1 DELTA.
!
      DO 70 L21=2,K
      IF (L21.EQ.KOL1.OR.L21.EQ.KOL2) GO TO 70
      MAX1=L21-1
      MIN1=1
      L(1,2)=L21
      IF (L21.GT.KOL1) MIN1=KOL1
      IF (L21.GT.KOL2) MIN1=KOL2
      DO 69 L11=MIN1,MAX1
      L(1,1)=L11
      CALL TRACA(A,K,N,1,L(1,1),L(1,2),TR)
      B(N)=B(N)+FACT(2)*TR
!
!       CALCULATE TERMS WITH 2 DELTAS.
!
      DO 68 L22=L21,K
      IF (L22.EQ.KOL1.OR.L22.EQ.KOL2) GO TO 68
      MAX2=L22-1
      MIN2=1
      IF (L22.GT.KOL1) MIN2=KOL1
      IF (L22.GT.KOL2) MIN2=KOL2
      IF (L22.EQ.L21)  MIN2=L11+1
      IF (MIN2.GT.MAX2) GO TO 68
      L(2,2)=L22
      DO 67 L12=MIN2,MAX2
      L(2,1)=L12
      DO 29 J1=1,2
      DO 29 J2=1,2
   29 IF (L(2,J1).EQ.L(1,J2)) GO TO 67
      CALL TRACA(A,K,N,2,L(1,1),L(1,2),TR)
      B(N)=B(N)+FACT(3)*TR
!
!       CALCULATE TERMS WITH 3 DELTAS.
!
      DO 66 L23=L22,K
      IF (L23.EQ.KOL1.OR.L23.EQ.KOL2) GO TO 66
      MAX3=L23-1
      MIN3=1
      IF (L23.GT.KOL1) MIN3=KOL1
      IF (L23.GT.KOL2) MIN3=KOL2
      IF (L23.EQ.L22)  MIN3=L12+1
      IF (MIN3.GT.MAX3) GO TO 66
      L(3,2)=L23
      DO 65 L13=MIN3,MAX3
      L(3,1)=L13
      DO 36 J1=1,2
      DO 36 J2=1,2
      DO 36 JX=1,2
   36 IF (L(3,J1).EQ.L(JX,J2)) GO TO 65
      CALL TRACA(A,K,N,3,L(1,1),L(1,2),TR)
      B(N)=B(N)+FACT(4)*TR
!
!       CALCULATE TERMS WITH 4 DELTAS.
!
      DO 64 L24=L23,K
      IF (L24.EQ.KOL1.OR.L24.EQ.KOL2) GO TO 64
      MAX4=L24-1
      MIN4=1
      IF (L24.GT.KOL1) MIN4=KOL1
      IF (L24.GT.KOL2) MIN4=KOL2
      IF (L24.EQ.L23) MIN4=L13+1
      IF (MIN4.GT.MAX4) GO TO 64
      L(4,2)=L24
      DO 63 L14=MIN4,MAX4
      L(4,1)=L14
      DO 40 J1=1,2
      DO 40 J2=1,2
      DO 40 JX=1,3
   40 IF (L(4,J1).EQ.L(JX,J2)) GO TO 63
      CALL TRACA(A,K,N,4,L(1,1),L(1,2),TR)
      B(N)=B(N)+FACT(5)*TR
!
!       CALCULATE TERMS WITH 5 DELTAS.
!
      DO 62 L25=L24,K
      IF (L25.EQ.KOL1.OR.L25.EQ.KOL2) GO TO 62
      MAX5=L25-1
      MIN5=1
      IF (L25.GT.KOL1) MIN5=KOL1
      IF (L25.GT.KOL2) MIN5=KOL2
      IF (L25.EQ.L24)  MIN5=L14+1
      IF (MIN5.GT.MAX5) GO TO 62
      L(5,2)=L25
      DO 61 L15=MIN5,MAX5
      L(5,1)=L15
      DO 44 J1=1,2
      DO 44 J2=1,2
      DO 44 JX=1,4
   44 IF (L(5,J1).EQ.L(JX,J2)) GO TO 61
      CALL TRACA(A,K,N,5,L(1,1),L(1,2),TR)
      B(N)=B(N)+FACT(6)*TR
   61 CONTINUE
   62 CONTINUE
   63 CONTINUE
   64 CONTINUE
   65 CONTINUE
   66 CONTINUE
   67 CONTINUE
   68 CONTINUE
   69 CONTINUE
   70 CONTINUE
  100 CONTINUE
      RETURN
      END


      SUBROUTINE TRACA(A,K,N,M,I1,I2,TR)
!
!  PURPOSE: COMPUTES MULTIPLE TRACES OF A TOTALLY SYMMETRIC TENSOR.
!
!  DEFINITIONS:
!
!    *A   K-TH RANK TOTALLY SYMMETRIC TENSOR IN COMPRESSED FORM.
!         DIMENSION SPECIFIED BY CALLING PROGRAM MUST BE AT LEAST
!         (K+1)*(K+2)/2.
!    *K   RANK OF A (K.LE.10).
!    *N   NUMBER OF ELEMENT OF A ON WHICH TRACES ARE BASED
!         (N.LE.(K+1)*(K+2)/2).
!    *M   NUMBER OF INDEX PAIRS IN WHICH TRACE IS TO BE TAKEN (M.LE.5).
!    *I1  ARRAY OF FIRST INDICES OF PAIRS IN WHICH TRACES ARE TAKEN.
!    *I2  ARRAY OF SECOND INDICES OF PAIRS IN WHICH TRACES ARE TAKEN.
!     TR  COMPUTED TRACE.
!
!  SUBROUTINES REQUIRED:  INDX3, INDX4
!
!  PROGRAMMED BY J. APPLEQUIST, 6/30/83.
!
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(1), IA(10), I1(5), I2(5)
      CALL INDX4(K,N,IA,IER)
      TR=0.D0
      DO 30 J1=1,3
      IA(I1(1))=J1
      IA(I2(1))=J1
      IF (M.GT.1) GO TO 10
      CALL INDX3(K,IA,NX,IER)
      TR=TR+A(NX)
      GO TO 30
   10 DO 28 J2=1,3
      IA(I1(2))=J2
      IA(I2(2))=J2
      IF(M.GT.2) GO TO 12
      CALL INDX3(K,IA,NX,IER)
      TR=TR+A(NX)
      GO TO 28
   12 DO 26 J3=1,3
      IA(I1(3))=J3
      IA(I2(3))=J3
      IF (M.GT.3) GO TO 14
      CALL INDX3(K,IA,NX,IER)
      TR=TR+A(NX)
      GO TO 26
   14 DO 24 J4=1,3
      IA(I1(4))=J4
      IA(I2(4))=J4
      IF (M.GT.4) GO TO 16
      CALL INDX3(K,IA,NX,IER)
      TR=TR+A(NX)
      GO TO 24
   16 DO 22 J5=1,3
      IA(I1(5))=J5
      IA(I2(5))=J5
      CALL INDX3(K,IA,NX,IER)
      TR=TR+A(NX)
   22 CONTINUE
   24 CONTINUE
   26 CONTINUE
   28 CONTINUE
   30 CONTINUE
      RETURN
      END


      SUBROUTINE INDX4(K,N,IA,IER)
!
!  CONVERTS INDEX N OF AN ELEMENT OF A TOTALLY SYMMETRIC K-TH ORDER
!  TENSOR STORED IN A COMPRESSED 1-DIMENSIONAL ARRAY TO
!  COORDINATE INDICES IA(1)...IA(K).  ASSUMES 1ST INDICES VARY MOST
!  RAPIDLY IN THE ARRAY.
!  J. APPLEQUIST 7/9/81.
!
!  ERROR CODES:  IER=0   NO ERROR
!                IER=2   N IS GREATER THAN (K+2)*(K+1)/2.  INDICES NOT FOUND.
!
!  Revised 2/7/93 by J. Applequist to remove upper limit on K.
!
      DIMENSION IA(K)
      if (n .gt. (K+2)*(k+1)/2) then
        ier=2
        return
      else
        ier=0
      endif
      call indx13(n,k,n1,n2,n3)
      call indx7(ia,k,n1,n2,n3)
      return 
      end


      SUBROUTINE INDX6(IA,N,N1,N2,N3)
!
!  PURPOSE: DETERMINES THE NUMBERS (N1,N2,N3) OF 1'S, 2'3, AND 3'S
!           OCCURRING IN A GIVEN N-TH RANK TENSOR INDEX ARRAY IA.
!
!  J. APPLEQUIST  1/8/82
!
      DIMENSION IA(N)
      N1=0
      N2=0
      DO 10 I=1,N
      IF (IA(I).EQ.1) N1=N1+1
      IF (IA(I).EQ.2) N2=N2+1
   10 CONTINUE
      N3=N-N1-N2
      RETURN
      END


      SUBROUTINE INDX7(IA,N,N1,N2,N3)
!
!  PURPOSE: DETERMINES INDEX ARRAY IA FOR AN N-TH RANK COMPRESSED
!           TENSOR, GIVEN THE NUMBERS (N1,N2,N3) OF 1'S, 2'S, AND
!           3'S OCCURRING IN THE ARRAY.
!
!  J. APPLEQUIST  1/8/82
!
      DIMENSION IA(N)
      IF (N3.EQ.0) GO TO 15
      DO 10 I=1,N3
   10 IA(I)=3
   15 IF (N2.EQ.0) GO TO 25
      DO 20 I=1,N2
   20 IA(I+N3)=2
   25 IF (N1.EQ.0) GO TO 35
      DO 30 I=1,N1
   30 IA(N3+N2+I)=1
   35 RETURN
      END


      SUBROUTINE INDX13(K,N,N1,N2,N3)
!
!  PURPOSE: DETERMINES THE INDICES N1,N2,N3 OF AN ELEMENT OF A
!     COMPRESSED TENSOR OF RANK N, GIVEN N AND THE POSITION K IN THE
!     CANONICAL ARRAY.
!
!  PROGRAMMED BY J. APPLEQUIST, 10/21/88 (NOTES JA-TT-8)
!
      IROW=INT(0.5+0.5*SQRT(8.*K-7.))
      IPOS=K-IROW*(IROW-1)/2
      N1=N-IROW+1
      N2=IROW-IPOS
      N3=IPOS-1
      RETURN
      END


      SUBROUTINE INDX3(K,IA,N,IER)
!
!  CONVERTS COORDINATE INDICES IA(1)...IA(K) OF AN ELEMENT OF A
!  TOTALLY SYMMETRIC K-TH ORDER TENSOR TO A SINGLE INDEX
!  N IN A 1-DIMENSIONAL ARRAY, USING COMPRESSED STORAGE.  ASSUMES
!  1ST INDICES VARY MOST RAPIDLY IN THE ARRAY.  IF INDICES IA
!  ARE NOT IN ORDER OF DECREASING MAGNITUDE, N IS CALCULATED FOR THE
!  PERMUTATION OF IA WHICH IS SO ORDERED.
!  J. APPLEQUIST, 7/9/81.
!
!  ERROR CODES:   IER=0   NO ERROR
!
!  Revised 6/16/93 by J. Applequist to remove upper limit on K.
!
!  Subroutines required: indx6, indx 12
!
      DIMENSION IA(k)
      IER=0
      call indx6(ia,k,n1,n2,n3)
      call indx12(n,n1,n2,n3)
      return
      end


      SUBROUTINE INDX12(K,N1,N2,N3)
!
!  PURPOSE: DETERMINES THE POSITION K IN  A COMPRESSED CANONICAL
!     ARRAY, GIVEN INDICES N1,N2,N3 FOR AN ELEMENT OF A TENSOR
!     OF RANK N1+N2+N3.
!
!  PROGRAMMED BY J. APPLEQUIST, 10/21/88 (NOTES JA-TT-8)
!
      IROW=N2+N3+1
      IPOS=N3+1
      K=IROW*(IROW-1)/2+IPOS
      RETURN
      END
end module
