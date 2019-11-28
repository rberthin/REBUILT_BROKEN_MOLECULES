PROGRAM REBUILT

IMPLICIT NONE

INTEGER :: num_species, i, n, m, choice, io
INTEGER :: atom, mol
DOUBLE PRECISION :: boxx, boxy, boxz
DOUBLE PRECISION :: xdiff, ydiff, zdiff
INTEGER, ALLOCATABLE, DIMENSION(:) :: num_atoms, ref_atom, num_molecules
INTEGER, ALLOCATABLE, DIMENSION(:) :: xref, yref, zref
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: X, Y, Z
CHARACTER(LEN=20), ALLOCATABLE, DIMENSION(:,:) :: nom

WRITE(6,*) '*************************************'
WRITE(6,*) '*      REBUILT BROKEN MOLECULES     * '
WRITE(6,*) '*************************************'
WRITE(6,*)
WRITE(6,*) ' I read only xyz files !'
WRITE(6,*)
WRITE(6,*) ' SPECIES PARAMETER '
WRITE(6,*) '-------------------'
WRITE(6,*) ' How many species are there? '
READ(5,*) num_species

ALLOCATE(num_atoms(num_species))
ALLOCATE(ref_atom(num_species))
ALLOCATE(num_molecules(num_species))

! INPUT NUMBER OF ATOMS IN EACH SPECIE 
!------------------------------------------
DO i = 1, num_species
  WRITE(6,*) 'For the specie', i
  WRITE(6,*) ' How many molecules of this specie ? '
  READ(5,*) num_molecules(i)
  WRITE(6,*) ' How many atoms are there in one molecule? '
  READ(5,*) num_atoms(i)
  WRITE(6,*) ' Choose the reference atom  '
  WRITE(6,*) ' Number in the same order as in the position file '
  READ(5,*) ref_atom(i) 
ENDDO 
!------------------------------------------

! INPUT BOX LENGTH
!------------------------------------------
WRITE(6,*)
WRITE(6,*) ' BOX PARAMETER '
WRITE(6,*) '---------------'
WRITE(6,*)
WRITE(6,*) ' !! PUT IN THE SAME UNIT AS THE POSITIONS !!'
WRITE(6,*)
WRITE(6,*) ' Length of the box in the x direction ?'
READ(5,*) boxx
WRITE(6,*) ' Length of the box in the y direction ?'
READ(5,*) boxy
WRITE(6,*) ' Length of the box in the z direction ?'
READ(5,*) boxz
WRITE(6,*) '--------------------------'
WRITE(6,*) '**** END OF THE INPUT ****'
WRITE(6,*) '--------------------------'

!------------------------------------------------------------------
!*************************** REBUILDING ***************************
!------------------------------------------------------------------

OPEN(unit = 20, file = "box_rebuilt.dat", iostat=io)
OPEN(unit = 10, file = "trajectories.xyz", status='old', iostat=io)
READ(10,*)
READ(10,*)

DO i = 1, num_species
  
  atom = num_atoms(i)
  mol = num_molecules(i)
  
  ALLOCATE(X(atom, mol))
  ALLOCATE(Y(atom, mol))
  ALLOCATE(Z(atom, mol))
  ALLOCATE(nom(atom, mol))
  ALLOCATE(xref(mol))
  ALLOCATE(yref(mol))
  ALLOCATE(zref(mol))
  
  WRITE(6,*) '** Specie', i, '**'
 
  DO n = 1, atom
    DO m = 1, mol
      READ(10,*) nom(n,m), X(n,m), Y(n,m), Z(n,m)
      IF (n == ref_atom(i)) THEN
        xref(m) = X(n,m)
        yref(m) = Y(n,m)
        zref(m) = Z(n,m)
      ENDIF
    ENDDO
  ENDDO
  
  DO m = 1, mol
    DO n = 1, atom
      xdiff = X(n,m)-xref(m)
      ydiff = Y(n,m)-xref(m)
      zdiff = Z(n,m)-zref(m)
        
      IF(xdiff > boxx/2) THEN
        X(n,m) = X(n,m) - boxx
      ELSE IF (xdiff < -boxx/2) THEN
        X(n,m) = X(n,m) + boxx      
      ENDIF
        
      IF(ydiff > boxy/2) THEN
        Y(n,m) = Y(n,m) - boxy
      ELSE IF (ydiff < -boxy/2) THEN
        Y(n,m) = Y(n,m) + boxy
      ENDIF

      IF(zdiff > boxz/2) THEN
        Z(n,m) = Z(n,m) - boxz
      ELSE IF (zdiff < -boxz/2) THEN
        Z(n,m) = Z(n,m) + boxz
      ENDIF
    ENDDO
  ENDDO
  
  DO n = 1, atom
    DO m = 1, mol
      WRITE(20,*) nom(n,m), X(n,m), Y(n,m), Z(n,m)
    ENDDO
  ENDDO
  DEALLOCATE(X)
  DEALLOCATE(Y)
  DEALLOCATE(Z)
  DEALLOCATE(nom)  
ENDDO
WRITE(6,*) '--------------------------'
WRITE(6,*) '******* ALL DONE ! *******'
WRITE(6,*) '--------------------------'
END PROGRAM
