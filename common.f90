!Sample program for solving Lid-driven cavity flow test using SIMPLE-algorithm
!Copyright (C) 2010  Michail Kiri
!
!This program is free software; you can redistribute it and/or
!modify it under the terms of the GNU General Public License
!as published by the Free Software Foundation; either version 2
!of the License, or (at your option) any later version.
!
!This program is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License
!along with this program; if not, write to the Free Software Foundation
!Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
!*********************************************************************
 module common
   implicit none
   private
    
   public :: init_workspace, init_all_cavity
   public :: init_grid, init_cells
   public :: init_exact, write_exact
   public :: out_array, out_paraview_vtk
   
   integer,  parameter, public :: dp = selected_real_kind(14,100) 
   real(dp), parameter    :: Pi =  3.14159265358979323844_dp
   integer, parameter, public :: Nvars = 4
   ! Reynolds number
   real(dp), public :: Re, Gam
   real(dp), public :: alfa, beta, omega
   integer, public :: Poiss_iter
   real(dp), public :: Poiss_acc
   logical, public :: rhiechow = .true.

   ! dimensionless sides
   real(dp) :: SLx = 1.0_dp
   real(dp) :: SLy = 1.0_dp

   ! Grid definition  
   integer, public ::  NXmax = 121   !numero nodi in x    
   integer, public ::  NYmax = 121   !numero nodi in y
   integer, public ::  NXmaxC   !numero dei dof in x
   integer, public ::  NYmaxC   !numero dei dof in y
   integer, public :: NX_center = 60
   integer, public :: NY_center = 26
   integer, public :: X_length = 2
   integer, public :: Y_length = 50

   ! common /var/
   real(dp), dimension(:,:), allocatable, public :: U, V 
   real(dp), dimension(:), allocatable, public :: X, Y  
   real(dp), dimension(:), allocatable, public :: Xc, Yc 
   real(dp), dimension(:,:,:), allocatable, public :: F

   ! common /var_Geom/
   logical, dimension(:,:), allocatable, public :: obstacle

   ! common /geomt/
   real(dp), dimension(:,:), allocatable, public :: DPx_c ,  DPy_c
             
   ! common /koeff/
   real(dp), dimension(:,:), allocatable, public :: Con_e, Con_n
   real(dp), dimension(:,:), allocatable, public :: Ap, Ae, As,  An, Aw
   real(dp), dimension(:,:), allocatable, public :: invAp
   real(dp), dimension(:,:,:), allocatable, public :: Sp

   ! exact solutions:
   real :: x_exct(17), v_exct(17)
   real :: y_exct(17), u_exct(17)
 
   logical, public :: readstart
   character(64), public :: filename
   character(1), public :: mode = 'H'

  contains 

  !****************************************************************************
  !****************************************************************************
  subroutine init_workspace(Nx,Ny)
    integer :: Nx, Ny
    integer :: err

   ! MEMORY: 18 arrays nx*ny
   !          3 arrays nx*ny*Nvars
    print*, Nx, Ny 
    print*,'Req Memory:',(17+3*Nvars)*Nx*Ny*dp/1e6,'Mb'
    
    allocate(U(Nx,Ny),stat=err)
    allocate(V(Nx,Ny),stat=err)
    allocate(Con_e(Nx,Ny),stat=err)
    allocate(Con_n(Nx,Ny),stat=err)
    allocate(invAp(Nx,Ny),stat=err)
    allocate(DPx_c(Nx,Ny),stat=err)
    allocate(DPy_c(Nx,Ny),stat=err)
    allocate(As(Nx,Ny),stat=err)
    allocate(An(Nx,Ny),stat=err)
    allocate(Ae(Nx,Ny),stat=err)
    allocate(Aw(Nx,Ny),stat=err)
    allocate(Ap(Nx,Ny),stat=err)
    allocate(obstacle(Nx,Ny),stat=err)
 
    allocate(F(Nx,Ny,Nvars),stat=err)
    allocate(Sp(Nx,Ny,Nvars),stat=err)

    if (err.ne.0) STOP 'ALLOCATION ERROR'


    U(:,:) = 0.0_dp
    V(:,:) = 0.0_dp 
 
    F(:,:,:) = 0.0_dp

    Con_e(:,:) = 0.0_dp
    Con_n(:,:) = 0.0_dp  
 
    invAp(:,:) = 0.0_dp

    Ap(:,:) = 0.0_dp
    As(:,:) = 0.0_dp
    An(:,:) = 0.0_dp
    Aw(:,:) = 0.0_dp
    Ae(:,:) = 0.0_dp
    Sp(:,:,:) = 0.0_dp
    obstacle(:,:) = .FALSE.

  end  subroutine init_workspace

  !****************************************************************************
  !****************************************************************************
  subroutine Init_all_cavity()

    integer :: i, j
    real(dp) :: tmp

    ! Gamma= (Reynolds)*(-1)
    Gam = 1.d0/Re

    ! Reading starting solution from file
    if (readstart) then
       open(44,file=filename)
       !read(44,*) 
       !read(44,*)
       do j = 1, NYmaxC
         do i = 1, NXmaxC
           read(44,*) tmp, tmp, F(i,j,1), F(i,j,2), F(i,j,4)
         end do
       end do  
       close(44)
    endif

    call init_forced_cavity()

    !call init_constant_flux()
    !call init_obstacle()

  end subroutine init_all_cavity

  ! condizione al contorno: velocità v positiva dalla parete inferiore
  ! slip-slit boundary condition u=1.0 on top surface
  subroutine init_forced_cavity()
    F(1     ,:     ,1) = 0.0_dp
    F(NXmaxC,:     ,1) = 0.0_dp
    F(:     ,1     ,1) = 0.0_dp
    F(:     ,NYmaxC,1) = 1.0_dp
    F(:     ,1     ,2) = 0.0_dp
    F(:     ,NYmaxC,2) = 0.0_dp
    F(1     ,:     ,2) = 0.0_dp
    F(NXmaxC,:     ,2) = 0.0_dp
  end subroutine init_forced_cavity

  ! Parabolic velocity profile on first and last grid points 
  subroutine init_constant_flux()
     integer :: j   
     do j = 1, NYmaxC
       F(1,      j  ,1) = 4.0_dp * Yc(j) * ( 1.0_dp - Yc(j))
       F(NXmaxC, j  ,1) = 4.0_dp * Yc(j) * ( 1.0_dp - Yc(j)) 
     end do
  end subroutine init_constant_flux

  subroutine init_obstacle()
    integer :: i,j   
    do i = NX_center - X_length/2 , NX_center + X_length/2
      do j = NY_center - Y_length/2, NY_center + Y_length/2
         if (i >= 1 .and. i <= NXMaxC .and. j >= 1 .and. j<= NYMaxC) then 
            obstacle(i,j) = .TRUE.
         end if  
         !if ((i-NX_center)*(i-NX_center) + (j-NY_center)*(j-NY_center) < X_length*X_length) then
         !   obstacle(i,j) = .TRUE.
         !end if   
      end do
    end do
  end subroutine init_obstacle
  !***************************************************************************
  !****************************************************************************
  subroutine init_grid()

    integer :: i, j
    real(dp) :: dx, dy
   
    allocate(X(NXmax))
    allocate(Y(NYmax))

    dx = SLx/(NXmax-1)
    dy = SLy/(NYmax-1)

    ! Definizione di una griglia di punti regolare (1..NXmax, 1..NYmax) 
    DO i=1, NXmax
       X(i) = dx * (i-1)
    ENDDO 

    DO j=1, NYmax
       Y(j) = dy * (j-1) 
    ENDDO
    ! ----------------------------------------------------------------    

    open (22,file='GRID.dat') !scrittura dei punti della griglia su file 
      WRITE(22,*)'VARIABLES = "X", "Y" ' 
      WRITE (22,*)' ZONE I=' ,NXmax, ', J=', NYmax, ', F=POINT'
      DO I=1, NXmax
         WRITE (22,*) X(I)  
      ENDDO
      DO J=1, NYmax
         WRITE (22,*) Y(J)  
      ENDDO 
    close(22)
    
  end subroutine init_grid


  subroutine init_cells()

    integer :: i, j, n

    !calculation Xc,Yc (calcolo delle posizioni dei volumetti) ------------------
    NXmaxC = NXmax + 1
    NYmaxC = NYmax + 1
    allocate(Xc(NXmaxC))
    allocate(Yc(NYmaxC))

    do i=2,NXmax
        Xc(i)=( X(i-1) + X(i) ) * 0.5_dp
    end do
    
    do j=2,NYmax
        Yc(j)=(  Y(j-1) + Y(j) ) * 0.5_dp 
    enddo
   
    Xc(1) = X(1)
    Xc(NXmaxC) = X(NXmax)
    
    Yc(1) = Y(1)
    Yc(NYmaxC) = Y(NYmax) 
  
  end subroutine init_cells

  !****************************************************************************
  ! co-located grid definition
  ! + represent original grid X(i,j), Y(i,j)
  ! o represent colocal grid Xc, Yc
  !                                         
  !  1  2     3                                NXmax+1 
  !  o--o--+--o--+--o--+--o--+--o--+--o--+--o--o
  !  o  o  |  o  |  o  |  o  |  o  |  o  |  o  o
  !  +-----+-----+-----+-----+-----+-----+-----+
  !  o  o  |  o  |  o  |  o  |  o  |  o  |  o  o
  !  +-----+-----+-----+-----+-----+-----+-----+
  !  o  o  |  o  |  o  |  o  |  o  |  o  |  o  o
  !  +-----+-----+-----+-----+-----+-----+-----+
  !  o  o  |  o  |  o  |  o  |  o  |  o  |  o  o
  !  +-----+-----+-----+-----+-----+-----+-----+
  !  o  o  |  o  |  o  |  o  |  o  |  o  |  o  o
  !  +-----+-----+==*=====*=====*==+-----+-----+
  !  o  o  |  o  *  X  |  X  |  X  *  o  |  o  o
  !  +-----+-----+-----+-----+-----+-----+-----+           flusso in presenza di un ostacolo
  !  o  o  |  o  *  X  |  X  |  X  *  o  |  o  o
  !  +-----+-----+-----+-----+-----+-----+-----+
  !  o  o  |  o  *  X  |  X  |  X  *  o  |  o  o
  !  +-----+-----+==*=====*=====*==+-----+-----+
  !  o  o  |  o  |  o  |  o  |  o  |  o  |  o  o
  !  +-----+-----+-----+-----+-----+-----+-----+
  !  o  o  |  o  |  o  |  o  |  o  |  o  |  o  o
  !  +-----+-----+-----+-----+-----+-----+-----+
  !  o  o  |  o  |  o  |  o  |  o  |  o  |  o  o
  !  +-----+-----+-----+-----+-----+-----+-----+
  !  o  o  |  o  |  o  |  o  |  o  |  o  |  o  o
  !  +-----+-----+-----+-----+-----+-----+-----+
  !  o  o  |  o  |  o  |  o  |  o  |  o  |  o  o
  !  +-----+-----+-----+-----+-----+-----+-----+
  !  o  o  |  o  |  o  |  o  |  o  |  o  |  o  o
  !  o--o--+--o--+--o--+--o--+--o--+--o--+--o--o
  !  1     2     3     4                       NXmax (nodi)
  !
  !  0.0                                       1.0
  !****************************************************************************
  !****************************************************************************
  
  subroutine init_exact()
    x_exct(1)=-0.500000+0.5; v_exct(1)= -0.000882
    x_exct(2)=-0.327052+0.5; v_exct(2)= -0.383699
    x_exct(3)=-0.397406+0.5; v_exct(3)= -0.297251
    x_exct(4)=-0.217948+0.5; v_exct(4)= -0.27788
    x_exct(5)=-0.428914+0.5; v_exct(5)= -0.222276
    x_exct(6)=-0.436290+0.5; v_exct(6)= -0.201989
    x_exct(7)=-0.444335+0.5; v_exct(7)= -0.181701
    x_exct(8)=-0.046595+0.5; v_exct(8)= -0.106804
    x_exct(9)= 0.001598+0.5; v_exct(9)= -0.060949
    x_exct(10)=0.118733+0.5; v_exct(10)= 0.057217
    x_exct(11)=0.235193+0.5; v_exct(11)= 0.186849
    x_exct(12)=0.352315+0.5; v_exct(12)= 0.333239
    x_exct(13)=0.45404 +0.5; v_exct(13)= 0.466401
    x_exct(14)=0.461386+0.5; v_exct(14)= 0.511382
    x_exct(15)=0.469392+0.5; v_exct(15)= 0.574884
    x_exct(16)=0.476719+0.5; v_exct(16)= 0.659554
    x_exct(17)=0.5     +0.5; v_exct(17)= 0.999118

    !x_exct(1)=0.00282901; v_exct(1)=-0.07028612
    !x_exct(2)=6.42882E-2; v_exct(2)=2.71771E+0
    !x_exct(3)=7.35542E-2; v_exct(3)=2.86215E+0
    !x_exct(4)=8.09753E-2; v_exct(4)=3.02728E+0
    !x_exct(5)=9.76486E-2; v_exct(5)=3.25421E+0
    !x_exct(6)=1.58722E-1; v_exct(6)=3.70750E+0
    !x_exct(7)=2.28897E-1; v_exct(7)=3.31348E+0
    !x_exct(8)=2.34432E-1; v_exct(8)=3.25139E+0
    !x_exct(9)=5.01948E-1; v_exct(9)=1.87997E-1
    !x_exct(10)=8.04508E-1; v_exct(10)=-3.33066E+0
    !x_exct(11)=8.59780E-1; v_exct(11)=-4.42685E+0
    !x_exct(12)=9.07688E-1; v_exct(12)=-5.33693E+0
    !x_exct(13)=9.44865E-1; v_exct(13)=-4.07737E+0
    !x_exct(14)=9.54199E-1; v_exct(14)=-3.51971E+0
    !x_exct(15)=9.61692E-1; v_exct(15)=-2.92069E+0
    !x_exct(16)=9.69195E-1; v_exct(16)=-2.25968E+0
    !x_exct(17)=1.00098E+0; v_exct(17)=-9.09091E-2
    y_exct(1)=-0.500000 +0.5; u_exct(1)=0.00069404
    y_exct(2)=-0.43768  +0.5; u_exct(2)=0.275621
    y_exct(3)=-0.429602 +0.5; u_exct(3)=0.290847
    y_exct(4)=-0.421523 +0.5; u_exct(4)=0.303994
    y_exct(5)=-0.406521 +0.5; u_exct(5)=0.326826
    y_exct(6)=-0.343624 +0.5; u_exct(6)=0.371038
    y_exct(7)=-0.273803 +0.5; u_exct(7)=0.330015
    y_exct(8)=-0.265724 +0.5; u_exct(8)=0.32307
    y_exct(9)=-0.000289 +0.5; u_exct(9)=0.0252893
    y_exct(10)= 0.304962+0.5; u_exct(10)=-0.318994
    y_exct(11)= 0.359781+0.5; u_exct(11)=-0.427191
    y_exct(12)= 0.406520+0.5; u_exct(12)=-0.515279
    y_exct(13)= 0.445182+0.5; u_exct(13)=-0.392034
    y_exct(14)= 0.453260+0.5; u_exct(14)=-0.336623
    y_exct(15)= 0.461339+0.5; u_exct(15)=-0.277749
    y_exct(16)= 0.46884 +0.5; u_exct(16)=-0.214023
    y_exct(17)= 0.5     +0.5; u_exct(17)=-6.20706e-17

    !y_exct(1)=-1.85185E-3; u_exct(1)=0.00000E+0
    !y_exct(2)=5.00000E-2; u_exct(2)=-1.84615E+0
    !y_exct(3)=5.92593E-2; u_exct(3)=-2.03077E+0
    !y_exct(4)=6.66667E-2; u_exct(4)=-2.21538E+0
    !y_exct(5)=1.00000E-1; u_exct(5)=-2.98462E+0
    !y_exct(6)=1.68519E-1; u_exct(6)=-3.81538E+0
    !y_exct(7)=2.77778E-1; u_exct(7)=-2.76923E+0
    !y_exct(8)=4.48148E-1; u_exct(8)=-1.07692E+0
    !y_exct(9)=4.96296E-1; u_exct(9)=-6.15385E-1
    !y_exct(10)=6.09259E-1; u_exct(10)=5.84615E-1
    !y_exct(11)=7.31481E-1; u_exct(11)=1.84615E+0
    !y_exct(12)=8.50000E-1; u_exct(12)=3.32308E+0
    !y_exct(13)=9.50000E-1; u_exct(13)=4.64615E+0
    !y_exct(14)=9.57407E-1; u_exct(14)=5.07692E+0
    !y_exct(15)=9.64815E-1; u_exct(15)=5.72308E+0
    !y_exct(16)=9.74074E-1; u_exct(16)=6.58462E+0
    !y_exct(17)=9.96296E-1; u_exct(17)=1.00000E+1
 
  end subroutine init_exact
  
  !****************************************************************************
  !****************************************************************************
  subroutine write_exact() 
    integer :: i, j

    open (23,file='Exact_V.dat') 
    !WRITE (23,*)' ZONE F=POINT, I=', 17

     	
    DO i=1,17     
      WRITE (23,*) x_exct(i), v_exct(i)
    ENDDO   
    close(23)

    ! ------------------------------------------- 
    open (23,file='Exact_U.dat') 
    !WRITE (23,*)' ZONE F=POINT, I=', 17
     DO j=1,17       
       WRITE (23,*) y_exct(j), u_exct(j)
    ENDDO
    close(23)
   
  end subroutine write_exact

  !****************************************************************************
  !****************************************************************************
  subroutine Out_array(F_out,NImax,NJmax,Filename)
     real(dp), dimension(:,:) :: F_out
     integer :: NIMAX, NJMAX
     character(10) :: Filename
 

     integer :: i, j

     Open(10,file=trim(Filename))

     do j=1,NJmax
       !write(10,'(I5)', advance='NO') j
       do i =1, NImax 
          write(10,'(F20.10)', advance='NO') F_out(i,j) 
       enddo
       write(10,*)
     enddo
     close(10)

   end subroutine out_array


   subroutine out_paraview_vtk(filename)
     character(*) :: filename
     

     integer :: nx, ny, nz, i,j
     character :: tab 

     tab = char(9) 
     nx = NXMAX-1
     ny = NYMAX-1
     nz = 1

     open(122,file=trim(filename))
     
     write(122,'(A)') "# vtk DataFile Version 2.0"
     write(122,'(A)') "U V velocities"
     write(122,'(A)') "ASCII"
     write(122,'(A)') "DATASET RECTILINEAR_GRID"
     write(122,'(A,3(I5))') "DIMENSIONS",nx,ny,nz
     write(122,'(A,I3,A)') "X_COORDINATES ",nx," float"
     do i = 2, NXMAX
       write(122,'(F15.7)') Xc(i) 
     end do
     write(122,'(A,I3,A)') "Y_COORDINATES ",ny," float"
     do j = 2, NYMAX
       write(122,'(F15.7)') Yc(j) 
     end do
     write(122,'(A,I3,A)') "Z_COORDINATES ",nz," float"
     write(122,'(F15.7)') 0.0 

     write(122,'(A,I6)') "POINT_DATA",nx*ny*nz
     write(122,'(A)') "VECTORS vectors float"
 
     do j = 2, NYMAX
       do i = 2, NXMAX
          write(122,'(F15.7,A,F15.7,A,F15.7)') F(i,j,1), tab, F(i,j,2), tab, 0.0 
       end do
     end do  

     close(122)

     open(122,file='p.vtk')

     write(122,'(A)') "# vtk DataFile Version 2.0"
     write(122,'(A)') "pressures"
     write(122,'(A)') "ASCII"
     write(122,'(A)') "DATASET RECTILINEAR_GRID"
     write(122,'(A,3(I5))') "DIMENSIONS",nx,ny,nz
     write(122,'(A,I3,A)') "X_COORDINATES ",nx," float"
     do i = 2, NXMAX
       write(122,'(F15.7)') Xc(i) 
     end do
     write(122,'(A,I3,A)') "Y_COORDINATES ",ny," float"
     do j = 2, NYMAX
       write(122,'(F15.7)') Yc(j) 
     end do
     write(122,'(A,I3,A)') "Z_COORDINATES ",nz," float"
     write(122,'(F15.7)') 0.0 

     write(122,'(A,I6)') "POINT_DATA",nx*ny*nz
     write(122,'(A)') "SCALARS pressure float 1"
     write(122,'(A)') "LOOKUP_TABLE default"    
     do j = 2, NYMAX
       do i = 2, NXMAX
          write(122,'(F15.7)') F(i,j,4) 
       end do
     end do  

     close(122)   

   end subroutine out_paraview_vtk


end module common
