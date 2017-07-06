! cube_tools.f90 - library for accessing and manipulating (Gaussian) cube files
! Copyright (C) 2007  Jörg Meyer
! 
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
! 
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
! 
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

MODULE cube_tools

  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER, PUBLIC :: Cube_Stdin_Unit = 5
  INTEGER, PARAMETER, PUBLIC :: Cube_Stdout_Unit = 6
  INTEGER, PARAMETER, PUBLIC :: Cube_Stderr_Unit = 0

  INTEGER, PARAMETER, PUBLIC :: Cube_Comment_Length = 132
  INTEGER, PARAMETER, PUBLIC :: Cube_Integer_Length = 5
  INTEGER, PARAMETER, PUBLIC :: Cube_Coordinate_Length = 12
  INTEGER, PARAMETER, PUBLIC :: Cube_Coordinate_Digits = 6
  INTEGER, PARAMETER, PUBLIC :: Cube_MOs_per_Line = 10
  INTEGER, PARAMETER, PUBLIC :: Cube_Voxels_per_Line = 6

  CHARACTER* (*), PARAMETER :: Cube_Comment_Format = 'A'
  CHARACTER* (*), PARAMETER :: Cube_Integer_Format = 'I5'
  CHARACTER* (*), PARAMETER :: Cube_Coordinate_Format = 'F12.6'
  CHARACTER* (*), PARAMETER :: Cube_Voxel_Format = 'E13.5E2'

!!! do not use repetition counts in the following to circumvent gfortran bug

  CHARACTER* (*), PARAMETER :: Cube_Vector_Format = &
&	Cube_Coordinate_Format // ',' // Cube_Coordinate_Format // ',' // Cube_Coordinate_Format
! &     '3' // Cube_Coordinate_Format

  CHARACTER* (*), PARAMETER :: Cube_Comment_Line_Format = &
&       '(' // Cube_Comment_Format // ')'
  CHARACTER* (*), PARAMETER :: Cube_Origin_Line_Format = &
&	'(' // Cube_Integer_Format // ',' // Cube_Vector_Format // ')'
  CHARACTER* (*), PARAMETER :: Cube_Cell_Line_Format = &
&	'(' // Cube_Integer_Format // ',' // Cube_Vector_Format // ')'
  CHARACTER* (*), PARAMETER :: Cube_Atom_Line_Format = &
&	'(' // Cube_Integer_Format // ',' // Cube_Coordinate_Format // ',' // Cube_Vector_Format // ':)'
  CHARACTER* (*), PARAMETER :: Cube_MO_Line_Format = &
&	'(' // Cube_Integer_Format // ',' // Cube_Integer_Format // ',' // Cube_Integer_Format // ',' &
&           // Cube_Integer_Format // ',' // Cube_Integer_Format // ',' // Cube_Integer_Format // ',' &
&           // Cube_Integer_Format // ',' // Cube_Integer_Format // ',' // Cube_Integer_Format // ',' &
&           // Cube_Integer_Format // ':)'
! &     '(' // '10' // Cube_Voxel_Format // ':)'
  CHARACTER* (*), PARAMETER :: Cube_Voxel_Line_Format = &
&	'(' // Cube_Voxel_Format // ',' // Cube_Voxel_Format // ',' // Cube_Voxel_Format // ',' &
&           // Cube_Voxel_Format // ',' // Cube_Voxel_Format // ',' // Cube_Voxel_Format // ':)'
! &     '(' // '6' // Cube_Voxel_Format // ':)'

  ! taken from cube2xsf out of XCrySDen 1.4.1
  REAL, PARAMETER, PUBLIC :: Cube_Bohr_to_Angstrom = 0.529177
  REAL, PARAMETER, PUBLIC :: Cube_Bohr3_to_Angstrom3 = Cube_Bohr_to_Angstrom**3

  REAL, PARAMETER, PUBLIC :: Cube_Hartree_to_eV = 27.2113845

! --- end of parameter defintions ---

  TYPE, PUBLIC :: cube_atom_type
    INTEGER :: number
    REAL :: unknown, r(3)
  END TYPE cube_atom_type

  TYPE, PUBLIC :: cube_type
    CHARACTER(Cube_Comment_Length) :: comment1, comment2
    INTEGER :: N, na, nb, nc
    REAL, DIMENSION(3) :: r0, da, db, dc
    TYPE(cube_atom_type), DIMENSION(:), POINTER :: atoms
    INTEGER :: NMOs
    INTEGER, DIMENSION(:), POINTER :: MOs
    ! implement as one dimensional array for performance reasons
    REAL, DIMENSION(:,:,:,:), POINTER :: voxels
  END TYPE cube_type

  INTERFACE OPERATOR (+)
    MODULE PROCEDURE cube_plus_cube
  END INTERFACE

  INTERFACE OPERATOR (-)
    MODULE PROCEDURE cube_minus_cube
  END INTERFACE

  INTERFACE OPERATOR (*)
    MODULE PROCEDURE cube_times_scalar
  END INTERFACE

  PUBLIC :: cube_clear
  PUBLIC :: cube_read
  PUBLIC :: cube_write, cube_write_stdout, cube_write_stderr
  PUBLIC :: cube_write_xsf_atoms, cube_write_xsf_atoms_stdout, cube_write_xsf_atoms_stderr
  PUBLIC :: cube_write_xsf_slab, cube_write_xsf_slab_stdout, cube_write_xsf_slab_stderr
  PUBLIC :: cube_write_xsf_crystal, cube_write_xsf_crystal_stdout, cube_write_xsf_crystal_stderr
  PUBLIC :: cube_write_info, cube_write_info_stdout, cube_write_info_stderr
  PUBLIC :: cube_volume, cube_voxel_volume, cube_voxels_sum
  PUBLIC :: cube_plus, cube_minus
  PUBLIC :: OPERATOR(+), OPERATOR(-)


CONTAINS


  SUBROUTINE cube_clear(cube)
    TYPE(cube_type), INTENT(OUT) :: cube
    NULLIFY(cube%MOs)
    cube%NMOs = 0; 
    NULLIFY(cube%voxels)
    cube%na = 0; cube%nb = 0; cube%nc = 0; 
    cube%da = 0.0; cube%db = 0.0; cube%dc = 0.0;
    NULLIFY(cube%atoms)
    cube%N = 0; 
    cube%r0 = 0.0;
    cube%comment1 = ''; cube%comment2 = ''
  END SUBROUTINE cube_clear

  SUBROUTINE cube_read(cube,unit,ierror)
    TYPE(cube_type), INTENT(OUT) :: cube
    INTEGER, INTENT(IN) :: unit
    INTEGER, INTENT(OUT), OPTIONAL :: ierror
    INTEGER i,ie,n,na,nb,nc
    READ(UNIT=unit,FMT=Cube_Comment_Line_Format,IOSTAT=ie) cube%comment1
    READ(UNIT=unit,FMT=Cube_Comment_Line_Format,IOSTAT=ie) cube%comment2
    ! list directed input is used here (instead of the format specifiers) 
    ! in order to get be able to read in cube files that do not strictly obey the standard
    READ(UNIT=unit,FMT=*,IOSTAT=ie) cube%N, cube%r0
    READ(UNIT=unit,FMT=*,IOSTAT=ie) cube%na, cube%da
    READ(UNIT=unit,FMT=*,IOSTAT=ie) cube%nb, cube%db
    READ(UNIT=unit,FMT=*,IOSTAT=ie) cube%nc, cube%dc
    ALLOCATE(cube%atoms(ABS(cube%N)),STAT=ie)
    DO n=1,ABS(cube%N),1
       READ(UNIT=unit,FMT=*,IOSTAT=ie) cube%atoms(n)%number, cube%atoms(n)%unknown, (cube%atoms(n)%r(i), i=1,3)
    END DO
    IF (cube%N .LT. 0) THEN
       READ(UNIT=unit,FMT=*,IOSTAT=ie) cube%NMOs
       ALLOCATE(cube%MOs(cube%NMOs))
       BACKSPACE(UNIT=unit,IOSTAT=ie)
       READ(UNIT=unit,FMT=Cube_MO_Line_Format,IOSTAT=ie) cube%NMOs, (cube%MOs(n), n=1,cube%NMOs,1)
    ELSE
       cube%NMOs = 1
       NULLIFY(cube%MOs)
    END IF
    ALLOCATE(cube%voxels(cube%na,cube%nb,cube%nc,cube%NMOs),STAT=ie)
    DO na=1,ABS(cube%na),1
       DO nb=1,ABS(cube%nb),1
          READ(UNIT=unit,FMT=Cube_Voxel_Line_Format,IOSTAT=ie) &
&           ((cube%voxels(na,nb,nc,n), n=1,ABS(cube%NMOs)), nc=1,ABS(cube%nc),1)
       END DO
    END DO
    IF ( present(ierror) ) ierror=ie
  END SUBROUTINE cube_read

  SUBROUTINE cube_write(cube,unit,ierror)
    TYPE(cube_type), INTENT(IN) :: cube
    INTEGER, INTENT(IN) :: unit
    INTEGER, INTENT(OUT), OPTIONAL :: ierror
    INTEGER ie,n,na,nb,nc
    ! ormat specifiers are used here in order to get to strictly obey the standard
    WRITE(UNIT=unit,FMT=Cube_Comment_Line_Format,IOSTAT=ie) TRIM(cube%comment1)
    WRITE(UNIT=unit,FMT=Cube_Comment_Line_Format,IOSTAT=ie) TRIM(cube%comment2)
    WRITE(UNIT=unit,FMT=Cube_Origin_Line_Format,IOSTAT=ie) cube%N, cube%r0
    WRITE(UNIT=unit,FMT=Cube_Cell_Line_Format,IOSTAT=ie) cube%na, cube%da
    WRITE(UNIT=unit,FMT=Cube_Cell_Line_Format,IOSTAT=ie) cube%nb, cube%db
    WRITE(UNIT=unit,FMT=Cube_Cell_Line_Format,IOSTAT=ie) cube%nc, cube%dc
    DO n=1,ABS(cube%N),1
       WRITE(UNIT=unit,FMT=Cube_Atom_Line_Format,IOSTAT=ie) cube%atoms(n)%number, cube%atoms(n)%unknown, cube%atoms(n)%r
    END DO
    IF (cube%N .lt. 0) THEN
! sanity check?
!       IF (.NOT. ASSOCIATED(cube%MOs)) STOP 'N < 0 but no MOs present'
       WRITE(UNIT=unit,FMT=Cube_MO_Line_Format,IOSTAT=ie) cube%NMOs, (cube%MOs(n), n=1,cube%NMOs,1)
    END IF
    DO na=1,ABS(cube%na),1
       DO nb=1,ABS(cube%nb),1
          WRITE(UNIT=unit,FMT=Cube_Voxel_Line_Format,IOSTAT=ie) &
&           ((cube%voxels(na,nb,nc,n), n=1,ABS(cube%NMOs)), nc=1,ABS(cube%nc),1)
       END DO
    END DO
    IF ( present(ierror) ) ierror=ie
  END SUBROUTINE cube_write

  SUBROUTINE cube_write_stdout(cube)
    TYPE(cube_type), INTENT(IN) :: cube
    CALL cube_write(cube,Cube_Stdout_Unit)
  END SUBROUTINE cube_write_stdout

  SUBROUTINE cube_write_stderr(cube)
    TYPE(cube_type), INTENT(IN) :: cube
    CALL cube_write(cube,Cube_Stderr_Unit)
  END SUBROUTINE cube_write_stderr

  SUBROUTINE cube_write_xsf_atoms(cube,unit,ierror)
    TYPE(cube_type), INTENT(IN) :: cube
    INTEGER, INTENT(IN) :: unit
    INTEGER, INTENT(OUT), OPTIONAL :: ierror
    INTEGER ie
    WRITE(UNIT=unit,FMT='(A)',IOSTAT=ie) 'ATOMS'
    CALL cube_xsf_write_coordinates(cube,unit,ie)
    CALL cube_xsf_write_voxels(cube,unit,ie)
    IF ( present(ierror) ) ierror=ie
  END SUBROUTINE cube_write_xsf_atoms

  SUBROUTINE cube_write_xsf_atoms_stdout(cube)
    TYPE(cube_type), INTENT(IN) :: cube
    CALL cube_write_xsf_atoms(cube,Cube_Stdout_Unit)
  END SUBROUTINE cube_write_xsf_atoms_stdout

  SUBROUTINE cube_write_xsf_atoms_stderr(cube)
    TYPE(cube_type), INTENT(IN) :: cube
    CALL cube_write_xsf_atoms(cube,Cube_Stderr_Unit)
  END SUBROUTINE cube_write_xsf_atoms_stderr

  SUBROUTINE cube_write_xsf_slab(cube,unit,ierror)
    TYPE(cube_type), INTENT(IN) :: cube
    INTEGER, INTENT(IN) :: unit
    INTEGER, INTENT(OUT), OPTIONAL :: ierror
    INTEGER ie
    WRITE(UNIT=unit,FMT='(A)',IOSTAT=ie) 'SLAB'
    WRITE(UNIT=unit,FMT='(A)',IOSTAT=ie) 'PRIMVEC'
    CALL cube_xsf_write_cell(cube,unit,ie)
    WRITE(UNIT=unit,FMT='(A)',IOSTAT=ie) 'PRIMCOORD'
    WRITE(UNIT=unit,FMT='(I12,I12)',IOSTAT=ie) ABS(cube%N), 1
    CALL cube_xsf_write_coordinates(cube,unit,ie)
    CALL cube_xsf_write_voxels(cube,unit,ie)
    IF ( present(ierror) ) ierror=ie
  END SUBROUTINE cube_write_xsf_slab

  SUBROUTINE cube_write_xsf_slab_stdout(cube)
    TYPE(cube_type), INTENT(IN) :: cube
    CALL cube_write_xsf_slab(cube,Cube_Stdout_Unit)
  END SUBROUTINE cube_write_xsf_slab_stdout

  SUBROUTINE cube_write_xsf_slab_stderr(cube)
    TYPE(cube_type), INTENT(IN) :: cube
    CALL cube_write_xsf_slab(cube,Cube_Stderr_Unit)
  END SUBROUTINE cube_write_xsf_slab_stderr

  SUBROUTINE cube_write_xsf_crystal(cube,unit,ierror)
    TYPE(cube_type), INTENT(IN) :: cube
    INTEGER, INTENT(IN) :: unit
    INTEGER, INTENT(OUT), OPTIONAL :: ierror
    INTEGER ie
    WRITE(UNIT=unit,FMT='(A)',IOSTAT=ie) 'CRYSTAL'
    WRITE(UNIT=unit,FMT='(A)',IOSTAT=ie) 'PRIMVEC'
    CALL cube_xsf_write_cell(cube,unit,ie)
    WRITE(UNIT=unit,FMT='(A)',IOSTAT=ie) 'PRIMCOORD'
    WRITE(UNIT=unit,FMT='(I12,I12)',IOSTAT=ie) ABS(cube%N), 1
    CALL cube_xsf_write_coordinates(cube,unit,ie)
    CALL cube_xsf_write_voxels(cube,unit,ie)
    IF ( present(ierror) ) ierror=ie
  END SUBROUTINE cube_write_xsf_crystal

  SUBROUTINE cube_write_xsf_crystal_stdout(cube)
    TYPE(cube_type), INTENT(IN) :: cube
    CALL cube_write_xsf_crystal(cube,Cube_Stdout_Unit)
  END SUBROUTINE cube_write_xsf_crystal_stdout

  SUBROUTINE cube_write_xsf_crystal_stderr(cube)
    TYPE(cube_type), INTENT(IN) :: cube
    CALL cube_write_xsf_crystal(cube,Cube_Stderr_Unit)
  END SUBROUTINE cube_write_xsf_crystal_stderr

  SUBROUTINE cube_xsf_write_cell(cube,unit,ierror)
    TYPE(cube_type), INTENT(IN) :: cube
    INTEGER, INTENT(IN) :: unit
    INTEGER, INTENT(OUT), OPTIONAL :: ierror
    INTEGER ie
    WRITE(UNIT=unit,FMT=100,IOSTAT=ie) Cube_Bohr_to_Angstrom * (ABS(cube%na)-1) * cube%da
    WRITE(UNIT=unit,FMT=100,IOSTAT=ie) Cube_Bohr_to_Angstrom * (ABS(cube%nb)-1) * cube%db
    WRITE(UNIT=unit,FMT=100,IOSTAT=ie) Cube_Bohr_to_Angstrom * (ABS(cube%nc)-1) * cube%dc
100   FORMAT(1X,F18.14,6X,F18.14,6X,F18.14)
    IF ( present(ierror) ) ierror=ie
  END SUBROUTINE cube_xsf_write_cell

  SUBROUTINE cube_xsf_write_coordinates(cube,unit,ierror)
    TYPE(cube_type), INTENT(IN) :: cube
    INTEGER, INTENT(IN) :: unit
    INTEGER, INTENT(OUT), OPTIONAL :: ierror
    INTEGER ie,n
    DO n=1,ABS(cube%N)
       ! format specifier taken from cube2xsf
       WRITE(UNIT=unit,FMT='(I3,2X,F12.6,F12.6,F12.6)',IOSTAT=ie) cube%atoms(n)%number, Cube_Bohr_to_Angstrom * cube%atoms(n)%r
    END DO
    IF ( present(ierror) ) ierror=ie
  END SUBROUTINE cube_xsf_write_coordinates

  SUBROUTINE cube_xsf_write_voxels(cube,unit,ierror)
    TYPE(cube_type), INTENT(IN) :: cube
    INTEGER, INTENT(IN) :: unit
    INTEGER, INTENT(OUT), OPTIONAL :: ierror
    INTEGER ie,n,na,nb,nc
    WRITE(UNIT=unit,FMT='(A)',IOSTAT=ie) 'BEGIN_BLOCK_DATAGRID_3D'
    DO n=1,cube%NMOs
       IF (cube%N .lt. 0) THEN 
          WRITE(UNIT=unit,FMT='(A)',IOSTAT=ie) &
&		TRIM(cube%comment1) // ', ' // TRIM(cube%comment2) // ' [written by cube_tools]'
          WRITE(UNIT=unit,FMT='(A,I3.3)',IOSTAT=ie) 'DATAGRID_3D_g98Cube_molecular_orbital#', cube%MOs(n)
       ELSE
          WRITE(UNIT=unit,FMT='(A)',IOSTAT=ie) &
&		TRIM(cube%comment1) // ', ' // TRIM(cube%comment2) // ' [written by cube_tools]'
          WRITE(UNIT=unit,FMT='(A)',IOSTAT=ie) 'DATAGRID_3D_g98Cube'
       END IF
       WRITE(UNIT=unit,FMT='(I12,I12,I12)',IOSTAT=ie) ABS(cube%na), ABS(cube%nb), ABS(cube%nc)
       WRITE(UNIT=unit,FMT=100,IOSTAT=ie) Cube_Bohr_to_Angstrom * cube%r0
100   FORMAT(1X,F18.14,6X,F18.14,6X,F18.14)
       CALL cube_xsf_write_cell(cube,unit,ie)
       DO nc=1,ABS(cube%nc),1
          DO nb=1,ABS(cube%nb),1
             WRITE(UNIT=unit,FMT=Cube_Voxel_Line_Format,IOSTAT=ie) (cube%voxels(na,nb,nc,n), na=1,ABS(cube%na),1)
          END DO
       END DO
       WRITE(UNIT=unit,FMT='(A)',IOSTAT=ie) 'END_DATAGRID_3D'
    END DO
    WRITE(UNIT=unit,FMT='(A)',IOSTAT=ie) 'END_BLOCK_DATAGRID_3D'
    IF ( present(ierror) ) ierror=ie
  END SUBROUTINE cube_xsf_write_voxels

  SUBROUTINE cube_write_info(cube,unit,ierror)
    TYPE(cube_type), INTENT(IN) :: cube
    INTEGER, INTENT(IN) :: unit
    INTEGER, INTENT(OUT), OPTIONAL :: ierror
    INTEGER ie,n
    WRITE(UNIT=unit,FMT='(A,T40,A)',IOSTAT=ie) 'comment line 1 : ', TRIM(cube%comment1)
    WRITE(UNIT=unit,FMT='(A,T40,A)',IOSTAT=ie) 'comment line 2 : ', TRIM(cube%comment2)
    WRITE(UNIT=unit,FMT=*) ('-', n=1,78)
    WRITE(UNIT=unit,FMT='(A,T40,I20)',IOSTAT=ie) 'number of atoms : ', ABS(cube%N)
    WRITE(UNIT=unit,FMT=*) ('-', n=1,78)
    WRITE(UNIT=unit,FMT='(A,T40,I20)',IOSTAT=ie) 'number of voxels : ', cube_voxels(cube)
    WRITE(UNIT=unit,FMT='(A,T40,F20.6,F20.6)',IOSTAT=ie) 'volume of each voxel (Bohr^3, A^3) : ', &
&	cube_voxel_volume(cube), cube_voxel_volume(cube) * Cube_Bohr3_to_Angstrom3
    WRITE(UNIT=unit,FMT='(A,T40,F20.6,F20.6)',IOSTAT=ie) 'total volume (Bohr^3, A^3) : ', &
&	cube_volume(cube), cube_volume(cube) * Cube_Bohr3_to_Angstrom3
    WRITE(UNIT=unit,FMT=*) ('-', n=1,78)
    IF (cube%N .LT. 0) THEN
       WRITE(UNIT=unit,FMT='(A,T40,I20)',IOSTAT=ie) 'number of molecular orbitals : ', cube%NMOs
       DO n=1,ABS(cube%NMOs)
          WRITE(UNIT=unit,FMT='(1X,A,T40,I20)',IOSTAT=ie) '+ molecular orbital (MO) : ', n
          WRITE(UNIT=unit,FMT='(1X,A,T40,I20)',IOSTAT=ie) '| MO identification number : ', ABS(cube%MOs(n))
          WRITE(UNIT=unit,FMT='(1X,A,T40,F20.6)',IOSTAT=ie) '| sum of all voxels of this MO : ', cube_voxels_sum(cube,n)
       END DO
    ELSE
       WRITE(UNIT=unit,FMT='(A,T40)',IOSTAT=ie) 'no molecular orbitals present'
       WRITE(UNIT=unit,FMT='(A,T40,F20.6)',IOSTAT=ie) 'sum of all voxels : ', cube_voxels_sum(cube,cube%NMOs)
    END IF
    IF ( present(ierror) ) ierror=ie
  END SUBROUTINE

  SUBROUTINE cube_write_info_stdout(cube)
    TYPE(cube_type), INTENT(IN) :: cube
    CALL cube_write_info(cube,Cube_Stdout_Unit)
  END SUBROUTINE cube_write_info_stdout

  SUBROUTINE cube_write_info_stderr(cube)
    TYPE(cube_type), INTENT(IN) :: cube
    CALL cube_write_info(cube,Cube_Stderr_Unit)
  END SUBROUTINE cube_write_info_stderr


  FUNCTION cube_voxel_volume(cube)
    TYPE(cube_type), INTENT(IN) :: cube
    REAL :: cube_voxel_volume
    REAL :: cross_product(3)
    cross_product(1) = cube%db(2) * cube%dc(3) - cube%db(3) * cube%dc(2)
    cross_product(2) = cube%db(3) * cube%dc(1) - cube%db(1) * cube%dc(3)
    cross_product(3) = cube%db(1) * cube%dc(2) - cube%db(2) * cube%dc(1)
    cube_voxel_volume = ABS(cube%da(1) * cross_product(1) + cube%da(2) * cross_product(2) + cube%da(3) * cross_product(3))
  END FUNCTION cube_voxel_volume

  FUNCTION cube_volume(cube)
    TYPE(cube_type), INTENT(IN) :: cube
    REAL :: cube_volume
    cube_volume = cube_voxel_volume(cube) * (ABS(cube%na)-1) * (ABS(cube%nb)-1) * (ABS(cube%nc)-1)
  END FUNCTION cube_volume

  FUNCTION cube_voxels(cube)
    TYPE(cube_type), INTENT(IN) :: cube
    INTEGER :: cube_voxels
    cube_voxels = ABS(cube%na * cube%nb * cube%nc)
  END FUNCTION cube_voxels

  FUNCTION cube_voxels_sum(cube,mo)
    TYPE(cube_type), INTENT(IN) :: cube
    INTEGER, INTENT(IN) :: mo
    REAL :: cube_voxels_sum,s
    INTEGER :: na,nb,nc
    s = 0.0
    DO na=1,cube%na
       DO nb=1,cube%nb
          DO nc=1,cube%nc
             s = s + cube%voxels(na,nb,nc,mo)
          END DO
       END DO
    END DO
    cube_voxels_sum = s
  END FUNCTION cube_voxels_sum

  FUNCTION cube_same_grid(cube1,cube2)
    TYPE(cube_type), INTENT(IN) :: cube1, cube2
    LOGICAL cube_same_grid
    cube_same_grid = (cube1%na .EQ. cube2%na) .AND. (cube1%nb .EQ. cube2%nb) .AND. (cube1%nc .EQ. cube2%nc)
  END FUNCTION cube_same_grid

  FUNCTION cube_same_cell(cube1,cube2)
    TYPE(cube_type), INTENT(IN) :: cube1, cube2
    LOGICAL cube_same_cell
    cube_same_cell = &
&	( ALL((cube1%da - cube2%da) == 0) .AND. ALL((cube1%db - cube2%db) == 0) .AND. ALL((cube1%dc - cube2%dc) == 0) )
  END FUNCTION cube_same_cell

  FUNCTION cube_same_grid_and_cell(cube1,cube2)
    TYPE(cube_type), INTENT(IN) :: cube1, cube2
    LOGICAL cube_same_grid_and_cell
    cube_same_grid_and_cell = cube_same_grid(cube1,cube2) .AND. cube_same_cell(cube1,cube2)
  END FUNCTION cube_same_grid_and_cell

  FUNCTION cube_plus(cube1,cube2,ierror) RESULT (cube3)
    TYPE(cube_type), INTENT(IN) :: cube1, cube2
    INTEGER, INTENT(OUT), OPTIONAL :: ierror
    TYPE(cube_type) :: cube3
    INTEGER ie,n,na,nb,nc
    CALL cube_clear(cube3)
    IF (cube_same_grid_and_cell(cube1, cube2)) THEN
       cube3%comment1 = 'sum cube'
       cube3%comment2 = TRIM(cube1%comment1) // ' + ' // TRIM(cube2%comment1)
       cube3%r0 = cube1%r0
       cube3%da = cube1%da; cube3%db = cube1%db; cube3%dc = cube1%dc
       IF ( (cube1%N .GE. 0) .OR. (cube2%N .GE. 0) ) THEN
          cube3%N = ABS(cube1%N) + ABS(cube2%N)
       ELSE
          cube3%N = - ( ABS(cube1%N) + ABS(cube2%N) )
       END IF
       ALLOCATE(cube3%atoms(ABS(cube3%N)),STAT=ie)
       DO n=1,ABS(cube1%N),1
          cube3%atoms(n) = cube1%atoms(n)
       END DO
       DO n=ABS(cube1%N)+1,ABS(cube3%N),1
          cube3%atoms(n) = cube2%atoms(n-ABS(cube1%N))
       END DO
       IF (cube3%N .LT. 0) THEN
          cube3%NMOs = MAX(cube1%NMOs,cube2%NMOs)
          ALLOCATE(cube3%MOs(cube3%NMos),STAT=ie)
          DO n=1,cube3%NMOs,1
             cube3%MOs(n) = cube1%MOs(n)
          END DO
       ELSE
          cube3%NMOs = 1
          NULLIFY(cube3%MOs)
       END IF
       cube3%na = cube1%na; cube3%nb = cube1%nb; cube3%nc = cube1%nc;
       ALLOCATE(cube3%voxels(cube3%na,cube3%nb,cube3%nc,cube3%NMOs),STAT=ie)
       DO n=1,ABS(cube3%NMOs),1
          DO na=1,ABS(cube3%na),1
             DO nb=1,ABS(cube3%nb),1
                DO nc=1,ABS(cube3%nc),1
                   cube3%voxels(na,nb,nc,n) = cube1%voxels(na,nb,nc,n) + cube2%voxels(na,nb,nc,n)
                END DO
             END DO
          END DO
       END DO
    ELSE
       CALL cube_clear(cube3)
       ie=1
    END IF
    IF ( present(ierror) ) ierror=ie
  END FUNCTION cube_plus

  FUNCTION cube_plus_cube(cube1,cube2) RESULT (cube3)
    TYPE(cube_type), INTENT(IN) :: cube1, cube2
    TYPE(cube_type) :: cube3
    INTEGER ie
    cube3 = cube_plus(cube1,cube2,ie)
  END FUNCTION cube_plus_cube

  FUNCTION cube_minus(cube1,cube2,ierror) RESULT (cube3)
    TYPE(cube_type), INTENT(IN) :: cube1, cube2
    INTEGER, INTENT(OUT), OPTIONAL :: ierror
    TYPE(cube_type) :: cube3
    INTEGER ie,n,na,nb,nc
    IF (cube_same_grid_and_cell(cube1, cube2)) THEN
       cube3%comment1 = 'difference cube'
       cube3%comment2 = TRIM(cube1%comment1) // ' - ' // TRIM(cube2%comment1)
       cube3%r0 = cube1%r0
       cube3%da = cube1%da; cube3%db = cube1%db; cube3%dc = cube1%dc
       IF ( (cube1%N .GE. 0) .OR. (cube2%N .GE. 0) ) THEN
          cube3%N = ABS(cube1%N)
       ELSE
          cube3%N = - ABS(cube1%N)
       END IF
       ALLOCATE(cube3%atoms(ABS(cube3%N)),STAT=ie)
       DO n=1,ABS(cube1%N),1
          cube3%atoms(n) = cube1%atoms(n)
       END DO
       IF (cube3%N .LT. 0) THEN
          cube3%NMOs = MAX(cube1%NMOs,cube2%NMOs)
          ALLOCATE(cube3%MOs(cube3%NMos),STAT=ie)
          DO n=1,cube3%NMOs,1
             cube3%MOs(n) = cube1%MOs(n)
          END DO
       ELSE
          cube3%NMOs = 1
          NULLIFY(cube3%MOs)
       END IF
       cube3%na = cube1%na; cube3%nb = cube1%nb; cube3%nc = cube1%nc;
       ALLOCATE(cube3%voxels(cube3%na,cube3%nb,cube3%nc,cube3%NMOs),STAT=ie)
       DO n=1,ABS(cube3%NMOs),1
          DO na=1,ABS(cube3%na),1
             DO nb=1,ABS(cube3%nb),1
                DO nc=1,ABS(cube3%nc),1
                   cube3%voxels(na,nb,nc,n) = cube1%voxels(na,nb,nc,n) - cube2%voxels(na,nb,nc,n)
                END DO
             END DO
          END DO
       END DO
    ELSE
       CALL cube_clear(cube3)
       ie=1
    END IF
    IF ( present(ierror) ) ierror=ie
  END FUNCTION cube_minus

  FUNCTION cube_minus_cube(cube1,cube2) RESULT (cube3)
    TYPE(cube_type), INTENT(IN) :: cube1, cube2
    TYPE(cube_type) :: cube3
    INTEGER ie
    cube3 = cube_minus(cube1,cube2,ie)
  END FUNCTION cube_minus_cube

  FUNCTION cube_scale(cube_in,scalar,ierror) RESULT (cube_out)
    TYPE(cube_type), INTENT(IN) :: cube_in
    REAL, INTENT(IN) :: scalar
    INTEGER, INTENT(OUT), OPTIONAL :: ierror
    TYPE(cube_type) :: cube_out
    INTEGER ie,n,na,nb,nc
    CALL cube_clear(cube_out)
    cube_out%comment1 = 'rescaled cube'
    cube_out%comment2 = TRIM(cube_in%comment1)
    cube_out%r0 = cube_in%r0
    cube_out%da = cube_in%da; cube_out%db = cube_in%db; cube_out%dc = cube_in%dc
    cube_out%N = cube_in%N
    IF (cube_out%N .LT. 0) THEN
       cube_out%NMOs = cube_in%NMOs
       ALLOCATE(cube_out%MOs(cube_out%NMos),STAT=ie)
       DO n=1,cube_out%NMOs,1
          cube_out%MOs(n) = cube_in%MOs(n)
       END DO
    ELSE
       cube_out%NMOs = 1
       NULLIFY(cube_out%MOs)
    END IF
    cube_out%na = cube_in%na; cube_out%nb = cube_in%nb; cube_out%nc = cube_in%nc;
    DO n=1,ABS(cube_out%NMOs),1
       DO na=1,ABS(cube_out%na),1
          DO nb=1,ABS(cube_out%nb),1
             DO nc=1,ABS(cube_out%nc),1
                cube_out%voxels(na,nb,nc,n) = cube_in%voxels(na,nb,nc,n) * scalar
             END DO
          END DO
       END DO
    END DO
    IF ( present(ierror) ) ierror=ie
  END FUNCTION cube_scale

  FUNCTION cube_times_scalar(cube_in,scalar) RESULT (cube_out)
    TYPE(cube_type), INTENT(IN) :: cube_in
    REAL, INTENT(IN) :: scalar
    TYPE(cube_type) :: cube_out
    INTEGER ie
    cube_out = cube_scale(cube_in,scalar,ie)
  END FUNCTION cube_times_scalar

END MODULE cube_tools
