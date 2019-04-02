!*********************************************************
!**  SUPER_RDF 
!**  
!**  Christian J Burnham, University College Dublin, 2019
!**  @christianjburnham@gmail.com
!*********************************************************

      real(8) :: rcut
      integer :: itype,jtype,ibinmax
      character(len = 16) :: rdf_filename
      character(len = 32) :: input_filename

      open(14,file = 'super_rdf.control')
      read(14,*) input_filename
      read(14,*) rcut,ibinmax
      read(14,*) itype,jtype

      rdf_filename = 'rdf.dat'

      open(10,file = input_filename)
      open(77,file = rdf_filename)

      call initialise
      call read_coordinates

!     calculate water-molecule RDF between atoms of type itype and atoms of type jtype, 
!     where O = type 1, and H = type 2. 

      call find_rdf(itype,jtype,ibinmax,rcut)

      write(*,*) '***************************'
      write(*,*) '***   SUPER_RDF              '
      write(*,*) '***************************'
      write(*,*) 'Radial distrbution function' 
      write(*,*) 'between atoms of type ',itype
      write(*,*) 'and atoms of type ',jtype
      write(*,*) 'up to cut-off radius of ',rcut,'Angstroms'
      write(*,*) 'with ',ibinmax,'bins'
      write(*,*) 'printed to file ',rdf_filename
      end

      subroutine initialise
      implicit none
      include 'super_rdf.common'

      pi = 4.0d0 * datan(1.0d0)
      angfac = 180.0d0 / pi

      end subroutine initialise


      subroutine read_coordinates
      implicit none
      include 'super_rdf.common'
      real(8) :: amag,bmag,cmag,alpha,beta,gamma
      real(8) :: x,y,z
      character(len = 10) :: card
      character(len = 3) :: atname
      integer iatom,ios,imol,i,moltype

!     reads in the coordinates from input file
!     simplified code, which will only work with water molecules, but 
!     can be adapted to work with different types of molecules. 

      read(10,*) natoms

!     amag,bmag,alpha,beta,gamma define the unit cell vectors
      read(10,*) card,amag,bmag,cmag,alpha,beta,gamma

!     find the cartesian form of the cell vectors
      call get_cartesian_cvec(amag,bmag,cmag,alpha,beta,gamma,cvec)

!     assume we have water molecules, in order O,H,H,O,H,H, etc.

      moltype = 1
      mol_natoms(moltype) = 3

!     assign atom-types O = 1 and H = 2

      atom_type(1,moltype) = 1
      atom_type(2,moltype) = 2
      atom_type(3,moltype) = 2

      do i = 1,natoms
         read(10,*,iostat = ios) atname,x,y,z

!     imol is the molecular index
         imol = (i - 1)/3 + 1
         iatom = mod(i-1,3) + 1

!     pack the coordinates into lab_coord array

         lab_coord(1,iatom,imol) = x
         lab_coord(2,iatom,imol) = y
         lab_coord(3,iatom,imol) = z

!     mol_type is the molecular type. Here, we are only considering one type of molecule, so assign every 
!     molecule a mol_type of 1

         mol_type(imol) = 1
      end do 

      nmol = imol

!     again, assuming we have water molecules
      natoms_of_type(1) = imol 
      natoms_of_type(2) = 2 * imol 

      end subroutine read_coordinates


      subroutine find_rdf(i1,j1,ibinmax,rcut)
      implicit none
      include 'super_rdf.common'
      integer imol,jmol,imoltype,jmoltype,iatom,jatom
      integer ibin,ibinmax
      integer natoms_in_imol,natoms_in_jmol
      integer i1,j1
      real(8) :: s
      real(8), allocatable, dimension(:) :: rdf
      real(8) :: dr,rcut
      real(8) :: xi,yi,zi,xj,yj,zj
      real(8) :: rdis,rdis2
      real(8) :: xdif,ydif,zdif
      real(8) :: dcellx,dcelly,dcellz
      real(8) :: dcyz,dcxz,dcxy
      real(8) :: vol
      integer :: jcell1max,jcell2max,jcell3max
      integer :: jcell1,jcell2,jcell3
      integer :: error
      integer itype,jtype
      real(8), dimension(3,3) :: scvec,scveci
      real(8) :: astar_mag,bstar_mag,cstar_mag
      logical :: same_mol,okflag,same_cell,same_atom,loop

      allocate(rdf(0:ibinmax),stat=error)

!     iniitalise rdf histogram
      rdf = 0.0d0 

!     find the reciprocal lattice vectors
      call mat3inverse(cvec,cveci,okflag)

!     calculate the cell volume (assuming upper triangular cell vector matrix)
      vol = cvec(1,1) * cvec(2,2) * cvec(3,3)         
      
!     magnitude of the reciprocal lattice vectors
      astar_mag = dsqrt(cveci(1,1)**2 + cveci(1,2)**2 + cveci(1,3)**2)
      bstar_mag = dsqrt(cveci(2,1)**2 + cveci(2,2)**2 + cveci(2,3)**2)
      cstar_mag = dsqrt(cveci(3,1)**2 + cveci(3,2)**2 + cveci(3,3)**2)

!     calculate the size of the supercell

      jcell1max = int(2.0d0 * rcut * astar_mag) + 1
      jcell2max = int(2.0d0 * rcut * bstar_mag) + 1
      jcell3max = int(2.0d0 * rcut * cstar_mag) + 1

!     find the supercell lattice vectors
      scvec(:,1) = jcell1max * cvec(:,1)
      scvec(:,2) = jcell2max * cvec(:,2)
      scvec(:,3) = jcell3max * cvec(:,3)
      
!     find the supercell reciprocal lattice vectors
      call mat3inverse(scvec,scveci,okflag)      
      
!     loop over the i particles 
      do imol = 1,nmol
         imoltype = mol_type(imol)
         natoms_in_imol = mol_natoms(imoltype)
         do iatom = 1,natoms_in_imol

!     xyz of the i particle

            xi = lab_coord(1,iatom,imol)
            yi = lab_coord(2,iatom,imol)
            zi = lab_coord(3,iatom,imol)

!     loop over the j particles 
            do jmol = imol,nmol
               jmoltype = mol_type(jmol)
               natoms_in_jmol = mol_natoms(jmoltype)
               do jatom = 1,natoms_in_jmol
                  
!     loop over the unit cells in the supercell
!     dcellx,dcelly,dcellz are the lattice vectors of the unit cell origins wrt supercell origin

                  do jcell3 = 0,jcell3max-1
                     dcellz = jcell3 * cvec(3,3)
                     dcyz = jcell3 * cvec(2,3)
                     dcxz = jcell3 * cvec(1,3)
                     do jcell2 = 0,jcell2max-1
                        dcelly = jcell2 * cvec(2,2) + dcyz
                        dcxy = jcell2 * cvec(1,2)
                        do jcell1 = 0,jcell1max-1
                           dcellx = jcell1 * cvec(1,1) + dcxy + dcxz
                           
!     xyz of the j particle

                           xj = lab_coord(1,jatom,jmol)
                           yj = lab_coord(2,jatom,jmol)
                           zj = lab_coord(3,jatom,jmol)

!     check whether i and j are the same atom (in different cells)
                           same_atom = .false.
                           if(imol.eq.jmol.and.iatom.eq.jatom) same_atom = .true.
                           
!     decide whether to loop

                           loop = .false.

!     check to see if the i and j particles exist in the same cell, which 
!     happens when jcell1 = jcell2 = jcell3 = 0 

                           same_cell = .false.
                           if(jcell1.eq.0.and.jcell2.eq.0.and.jcell3.eq.0) same_cell = .true.
                           
!     loop if j > i
                           if((jmol.gt.imol).or.(imol.eq.jmol.and.jatom.gt.iatom)) loop = .true.
                           

!     for i and j describing the same atom in the same unit cell: don't loop
!     for i and j describing the same atom in different unit cells: allow loop

                           if(same_atom) then
                              loop = .false.
                              if(.not. same_cell) loop = .true.
                           endif
                           
                           if(loop) then
                              
                              xdif = xj - xi + dcellx
                              ydif = yj - yi + dcelly
                              zdif = zj - zi + dcellz
                              
                              call nearest(xdif,ydif,zdif,scvec,scveci)

!     If i and j describe the same atom in different unit cells, then, to avoid double counting,
!     allow loop to proceed for precisely half the cases.
!     Not the most efficient way to calculate the same-atom contributions, but number of these interactions 
!     scale only linearly, so not too much harm calculating them this way.

                              if(same_atom) then
                                 if((xdif.gt.0).or.(xdif.eq.0.and.ydif.gt.0)&
     &                                .or.(xdif.eq.0.and.ydif.eq.0.and.zdif.gt.0)) then 
                              else
                                 cycle
                              endif
                           endif
                           
                           rdis2 = xdif**2 + ydif**2 + zdif**2

!     decide whether this is an intramolecular or intermolecular interaction
                              same_mol = .false.
                              if(same_cell.and.(imol.eq.jmol)) same_mol = .true.

!     increment rdf histogram 

                           if(.not.same_mol.and.rdis2.le.rcut*rcut) then 
                              itype = atom_type(iatom,imoltype)
                              jtype = atom_type(jatom,jmoltype)

                              if(itype.eq.i1.and.jtype.eq.j1.or.itype.eq.j1.and.jtype.eq.i1) then 
                                 rdis = dsqrt(rdis2)
                                 ibin = anint((rdis / rcut) * ibinmax)
                                 rdf(ibin) = rdf(ibin) + 1.0d0 
                              endif
                           endif
                           
                        endif
                        
                     end do !jcell1 loop
                  end do !jcell2 looop
               end do !jcell3 loop
               
               
            end do !jatom loop 
         end do !jmol loop 
      end do !iatom loop 
      end do !imol loop

!     print out the rdf

      dr = rcut / dble(ibinmax)
      do ibin = 1,ibinmax-1
         rdis = rcut * dble(ibin) / dble(ibinmax)

!     if i and j are of different atom types, multiply by 1/2 to avoid 
!     double counting.
         s = 1.0d0 
         if(i1.ne.j1) s = 0.5d0 
         write(77,*) rdis,s*rdf(ibin)/(2.0d0 * pi * rdis*rdis * & 
     & dr * dble(natoms_of_type(i1)*natoms_of_type(j1))/vol)
      end do 
      write(77,*)
      call flush(77)

      end subroutine find_rdf


      subroutine nearest(xdif,ydif,zdif,cmat,cmati)
      implicit none
      real(8),dimension(3,3) :: cmat,cmati
      real(8) :: xdif,ydif,zdif
      real(8),dimension(3) :: vec1,vec2

!     applies the minimum image convention for triclinic cells
!     note, assumes that both cmat and cmati are upper triangular

      vec1(1) = anint(cmati(1,1) * xdif + cmati(1,2) * ydif + cmati(1,3) * zdif)
      vec1(2) = anint(cmati(2,2) * ydif + cmati(2,3) * zdif)
      vec1(3) = anint(cmati(3,3) * zdif)

      xdif = xdif - (cmat(1,1) * vec1(1) + cmat(1,2) * vec1(2) + cmat(1,3) * vec1(3))
      ydif = ydif - (cmat(2,2) * vec1(2) + cmat(2,3) * vec1(3))
      zdif = zdif - (cmat(3,3) * vec1(3))

      end subroutine nearest


      subroutine get_cartesian_cvec(amag,bmag,cmag,alpha,beta,gamma,cvec0)
      implicit none
      include 'super_rdf.common'
      real(8) :: alpha,beta,gamma,amag,bmag,cmag
      real(8) :: alpha_rad,beta_rad,gamma_rad
      real(8), dimension(6) :: param
      real(8), dimension(3,3) :: cvec0

!     converts from the a,b,c,alpha,beta,gamma unit cell form to cartesian 
!     cell vectors, in upper triangular form.

      alpha_rad = alpha / angfac
      beta_rad = beta / angfac
      gamma_rad = gamma / angfac

      cvec0 = 0.0d0

      cvec0(1,1) = amag
      cvec0(1,2) = bmag * dcos(gamma_rad)
      cvec0(2,2) = dsqrt(bmag**2 - cvec0(1,2)**2)
      cvec0(1,3) = cmag * dcos(beta_rad)
      cvec0(2,3) = (bmag * cmag * dcos(alpha_rad) - cvec0(1,2) * cvec0(1,3))/cvec0(2,2)
      cvec0(3,3) = dsqrt(cmag**2 - cvec0(1,3)**2 - cvec0(2,3)**2)

      end subroutine get_cartesian_cvec

      SUBROUTINE mat3inverse (A, AINV, OK_FLAG)

!***********************************************************************************************************************************
!  M33INV  -  Compute the inverse of a 3x3 matrix.
!
!  A       = input 3x3 matrix to be inverted
!  AINV    = output 3x3 inverse of matrix A
!  OK_FLAG = (output) .TRUE. if the input matrix could be inverted, and .FALSE. if the input matrix is singular.
!***********************************************************************************************************************************


      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN)  :: A
      DOUBLE PRECISION, DIMENSION(3,3), INTENT(OUT) :: AINV
      LOGICAL, INTENT(OUT) :: OK_FLAG

      DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-10
      DOUBLE PRECISION :: DET
      DOUBLE PRECISION, DIMENSION(3,3) :: COFACTOR


      DET =   A(1,1)*A(2,2)*A(3,3)  &
            - A(1,1)*A(2,3)*A(3,2)  &
            - A(1,2)*A(2,1)*A(3,3)  &
            + A(1,2)*A(2,3)*A(3,1)  &
            + A(1,3)*A(2,1)*A(3,2)  &
            - A(1,3)*A(2,2)*A(3,1)

      IF (ABS(DET) .LE. EPS) THEN
         AINV = 0.0D0
         OK_FLAG = .FALSE.
         RETURN
      END IF

      COFACTOR(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
      COFACTOR(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
      COFACTOR(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
      COFACTOR(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
      COFACTOR(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
      COFACTOR(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
      COFACTOR(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
      COFACTOR(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
      COFACTOR(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))

      AINV = TRANSPOSE(COFACTOR) / DET

      OK_FLAG = .TRUE.

      RETURN

      END SUBROUTINE mat3inverse
