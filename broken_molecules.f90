PROGRAM REBUILT

IMPLICIT NONE

INTEGER :: num_species, i, n, m, io
INTEGER :: atom, mol, choice_input, num_xyz
DOUBLE PRECISION :: boxlen, boxx, boxy, boxz
DOUBLE PRECISION :: xdiff, ydiff, zdiff
INTEGER, ALLOCATABLE, DIMENSION(:) :: num_atoms, ref_atom, num_molecules
INTEGER, ALLOCATABLE, DIMENSION(:) :: xref, yref, zref
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: X, Y, Z
CHARACTER(LEN=20), ALLOCATABLE, DIMENSION(:,:) :: nom
CHARACTER(LEN=20) :: input_file
CHARACTER(LEN=1) :: inputfile_choice, box_choice

WRITE(6,*) '*************************************'
WRITE(6,*) '*      REBUILT BROKEN MOLECULES     * '
WRITE(6,*) '*************************************'
WRITE(6,*)
WRITE(6,*) ' I read only xyz files (for now..)!'
WRITE(6,*) 
WRITE(6,*) 'Read the input file parameters.inpt (1) or answer questions (2)?'
READ(5,*) choice_input
WRITE(6,*)
WRITE(6,*) ' SPECIES PARAMETER '
WRITE(6,*) '-------------------'

! *** READ PARAMETERS.INPT FILE ***
IF (choice_input == 1) THEN
  WRITE(6,*)'Parameters read from parameters.inpt'
  OPEN(unit = 11, file = "parameters.inpt", status='old', iostat=io)
  READ(11,*) input_file
  READ(11,*) num_species
  ALLOCATE(num_atoms(num_species))
  ALLOCATE(ref_atom(num_species))
  ALLOCATE(num_molecules(num_species))
  READ(11,*)    !read a comment 
  DO i = 1, num_species
    READ(11,*) num_molecules(i)
    READ(11,*) num_atoms(i)
    READ(11,*) ref_atom(i)
  ENDDO
  READ(11,*)   !read a comment
  READ(11,*) boxx
  READ(11,*) boxy
  READ(11,*) boxz
  ! *****     END     *****

ELSE IF (choice_input == 2) THEN
  WRITE(6,*) ' Is the xyz file named trajectories.xyz ? (y/n)'
  READ(5,*) inputfile_choice
  IF (inputfile_choice == 'n') THEN
    WRITE(6,*) ' Please enter the name of the xyz file '
    READ(5,*) input_file
  ELSE IF (inputfile_choice == 'y') THEN
    input_file = 'trajectories.xyz'
  ENDIF      
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
    WRITE(6,*) ' Choose the reference atom for this specie'
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
  WRITE(6,*) 'Are the 3 cell vectors the same size ?(y/n) '
  READ(5,*) box_choice
  IF (box_choice == 'y') THEN
    WRITE(6,*) 'Enter the length of the box '
    READ(5,*) boxlen
    boxx = boxlen
    boxy = boxlen
    boxz = boxlen
  ELSE IF (box_choice == 'n') THEN
    WRITE(6,*) ' Enter the length of the box in the x direction ?'
    READ(5,*) boxx
    WRITE(6,*) ' Enter the length of the box in the y direction ?'
    READ(5,*) boxy
    WRITE(6,*) ' Enter the length of the box in the z direction ?'
    READ(5,*) boxz
  ENDIF
ENDIF  
WRITE(6,*) '--------------------------'
WRITE(6,*) '**** END OF THE INPUT ****'
WRITE(6,*) '--------------------------'

!------------------------------------------------------------------
!*************************** REBUILDING ***************************
!------------------------------------------------------------------

OPEN(unit = 20, file = "box_rebuilt.dat", iostat=io)
OPEN(unit = 10, file = input_file, status='old', iostat=io)

READ(10,*) num_xyz
READ(10,*)

WRITE(20,*) num_xyz 
WRITE(20,*) 

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
      ydiff = Y(n,m)-yref(m)
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
