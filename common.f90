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
 
   public :: init, init_all_cavity
   public :: geom, init_grid
   public :: init_exact, write_exact
   public :: out_array
   public :: out_paraview_vtk

   integer,  parameter, public :: dp = 8
   !integer,  parameter, public ::  nx=100, ny=100   
   integer, parameter, public :: Nvars = 4 
   ! Reynolds number
   real(dp), parameter :: Re = 1000.0_dp
 
   ! Grid definition  
   integer, public ::  NXmax = 121 
   integer, public ::  NYmax = 121 
   integer, public ::  NXmaxC = 122 
   integer, public ::  NYmaxC = 122 

   ! common /var/
   real(dp), dimension(:,:), allocatable, public :: U, V 
   real(dp), dimension(:,:), allocatable, public :: X, Y  
   real(dp), dimension(:,:), allocatable, public :: Xc, Yc 
   real(dp), dimension(:,:), allocatable, public :: Gam
   real(dp), dimension(:,:,:), allocatable, public :: F

   ! common /var_Geom/
   !real(dp), dimension(:,:), allocatable, public :: Dx_c, Dy_c 

   ! common /geomt/
   real(dp), dimension(:,:), allocatable, public :: DPx_c ,  DPy_c
             
   ! common /koeff/
   real(dp), dimension(:,:), allocatable, public :: Con_e, Con_n
   real(dp), dimension(:,:), allocatable, public :: DpU, DpV  
   real(dp), dimension(:,:), allocatable, public :: Ae, As,  An, Aw
   real(dp), dimension(:,:,:), allocatable, public :: Ap
   real(dp), dimension(:,:,:), allocatable, public :: Sp

   ! exact solutions:
   real :: x_exct(17), v_exct(17)
   real :: y_exct(17), u_exct(17)
 
   logical, public :: readstart
   character(64), public :: filename

  contains 

  !****************************************************************************
  !****************************************************************************
  Subroutine init(Nx,Ny)
    integer :: Nx, Ny
    integer :: err

   ! MEMORY: 18 arrays nx*ny
   !          3 arrays nx*ny*Nvars 
    print*,'Req Memory:',(17+3*Nvars)*Nx*Ny*dp/1e6,'Mb'
    
    allocate(U(Nx,Ny),stat=err)
    allocate(V(Nx,Ny),stat=err)
    allocate(Xc(Nx,Ny),stat=err)
    allocate(Yc(Nx,Ny),stat=err)
    allocate(X(Nx,Ny),stat=err)
    allocate(Y(Nx,Ny),stat=err)
    allocate(Gam(Nx,Ny),stat=err)
    !allocate(Ro(Nx,Ny),stat=err)
    allocate(Con_e(Nx,Ny),stat=err)
    allocate(Con_n(Nx,Ny),stat=err)
    allocate(DpU(Nx,Ny),stat=err)
    allocate(DpV(Nx,Ny),stat=err)
    allocate(DPx_c(Nx,Ny),stat=err)
    allocate(DPy_c(Nx,Ny),stat=err)
    !allocate(Dx_c(Nx,Ny),stat=err)
    !allocate(Dy_c(Nx,Ny),stat=err)
    allocate(As(Nx,Ny),stat=err)
    allocate(An(Nx,Ny),stat=err)
    allocate(Ae(Nx,Ny),stat=err)
    allocate(Aw(Nx,Ny),stat=err)
 
    allocate(F(Nx,Ny,Nvars),stat=err)
    allocate(Ap(Nx,Ny,Nvars),stat=err)
    allocate(Sp(Nx,Ny,Nvars),stat=err)

    if (err.ne.0) STOP 'ALLOCATION ERROR'

    U(:,:) = 0.0_dp
    V(:,:) = 0.0_dp 
 
    F(:,:,:) = 0.0_dp
    Xc(:,:) = 0.0_dp
    Yc(:,:) = 0.0_dp
    X(:,:) = 0.0_dp
    Y(:,:) = 0.0_dp
    Gam(:,:) = 0.0_dp
    !Ro(:,:) = 1.0_dp

    Con_e(:,:) = 0.0_dp
    Con_n(:,:) = 0.0_dp  
 
    DPx_c(:,:) = 0.0_dp
    DPy_c(:,:) = 0.0_dp

    DpU(:,:) = 0.0_dp
    DpV(:,:) = 0.0_dp  

    Ap(:,:,:) = 0.0_dp
    As(:,:) = 0.0_dp
    An(:,:) = 0.0_dp
    Aw(:,:) = 0.0_dp
    Ae(:,:) = 0.0_dp
    Sp(:,:,:) = 0.0_dp

  end  subroutine init

  !****************************************************************************
  !****************************************************************************
  subroutine Init_all_cavity()

     integer :: i,j
     real(dp) :: tmp

     F = 0.0_dp
     Gam = 1.d0/Re

     ! slip-slit boundary condition u=1.0 on top surface
     F(1     ,:     ,1) = 0.0_dp
     F(NXmaxC,:     ,1) = 0.0_dp
     F(:     ,1     ,1) = 0.0_dp
     F(:     ,NYmaxC,1) = 1.0_dp

     F(:     ,1     ,2) = 0.0_dp
     F(:     ,NYmaxC,2) = 0.0_dp
     F(1     ,:     ,2) = 0.0_dp
     F(NXmaxC,:     ,2) = 0.0_dp

     if (readstart) then
        open(44,file=filename)
        read(44,*) 
        read(44,*)
        do j = 1, NYmaxC
          do i = 1, NXmaxC
            read(44,*) tmp, tmp, F(i,j,1), F(i,j,2), F(i,j,4)
          end do
        end do  
        close(44)
     endif

  end subroutine init_all_cavity
 
  !***************************************************************************
  !****************************************************************************
  subroutine init_grid()

    real(dp) :: SLx, SLy, Xbeg, Ybeg
    integer :: i, j, err

    ! dimensionless sides
 
    SLx = 1.0_dp
    SLy = 1.0_dp
    ! offset ?
    Xbeg =  0.0_dp
    Ybeg =  0.0_dp

    ! Definition of regular grid (1..NXmax, 1..NYmax) 
    DO i=1, NXmax
      DO j=1, NYmax
         X(i,j) = Xbeg + (SLx /(NXmax-1))* (i-1)
      ENDDO
    ENDDO 

    DO j=1, NYmax
      DO i=1, NXmax
         Y(i,j) = Ybeg + (SLy /(NYmax-1))* (j-1) 
      ENDDO
    ENDDO

    ! grid spacing of regular grid
    !Do i=2,NXmaxC-1   
    !  Do j=2,NYmaxC-1
    !     Dx_c(i,j) = X(i,j) - X(i-1,j)
    !     Dy_c(i,j) = Y(i,j) - Y(i,j-1)
    !  enddo
    !enddo
    ! ----------------------------------------------------------------    

    open (22,file='GRID.dat') 
      WRITE(22,*)'VARIABLES = "X", "Y" ' 
      WRITE (22,*)' ZONE I=' ,NXmax, ', J=', NYmax, ', F=POINT'
      DO J=1, NYmax
      	DO I=1, NXmax
           WRITE (22,*) X(I,J), Y(I,J)  
        ENDDO
      ENDDO 
    close(22)
    
  end subroutine init_grid


  !****************************************************************************
  ! co-located grid definition
  ! + represent original grid X(i,j), Y(i,j)
  ! o represent colocal grid Xc, Yc
  !
  !
  !  1     2     3                             NXmax
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
  !  +-----+-----+-----+-----+-----+-----+-----+
  !  o  o  |  o  |  o  |  o  |  o  |  o  |  o  o
  !  +-----+-----+-----+-----+-----+-----+-----+
  !  o  o  |  o  |  o  |  o  |  o  |  o  |  o  o
  !  o--o--+--o--+--o--+--o--+--o--+--o--+--o--o
  !  1  2     3     4                          NXmax+1
  !****************************************************************************
  subroutine Geom()

    integer :: i, j

    !calculation Xc,Yc -----------------------------------------
    do j=2,NYmax
       do i=2,NXmax
            Xc(i,j)=(  X(i-1,j-1) + X(i-1,j  ) + &
                       X(i  ,j  ) + X(i  ,j-1)    ) * 0.25_dp 

            Yc(i,j)=(  Y(i-1,j-1) + Y(i-1,j  ) + &
                       Y(i  ,j  ) + Y(i  ,j-1)    ) * 0.25_dp 
       enddo
    enddo
  
    do I=2,NXmax
      Xc(i,1)       = ( X(i  ,1    ) + X(i-1,1    ) ) * 0.5_dp
      Xc(i,NYmaxC)  = ( X(i  ,NYmax) + X(i-1,NYmax) ) * 0.5_dp
      Yc(i,1      ) = ( Y(i  ,1    ) + Y(i-1,1    ) ) * 0.5_dp
      Yc(i,NYmaxC)  = ( Y(i  ,NYmax) + Y(i-1,NYmax) ) * 0.5_dp
    enddo 

    do j=2,NYmax
      Yc(1     ,j ) = ( Y(1     ,j) + Y(1    ,j-1) ) * 0.5_dp
      Yc(NXmaxC,j ) = ( Y(NXmax ,j) + Y(NXmax,j-1) ) * 0.5_dp
      Xc(1     ,j ) = ( X(1     ,j) + X(1    ,j-1) ) * 0.5_dp
      Xc(NXmaxC,j ) = ( X(NXmax ,j) + X(NXmax,j-1) ) * 0.5_dp
    enddo
   
    Xc(1     ,     1) = X(    1,    1)
    Xc(NXmaxC,     1) = X(NXmax,    1)
    Xc(     1,NYmaxC) = X(    1,NYmax)
    Xc(NXmaxC,NYmaxC) = X(NXmax,NYmax)

    Yc(1     ,     1) = Y(    1,    1)
    Yc(NXmaxC,     1) = Y(NXmax,    1)
    Yc(     1,NYmaxC) = Y(    1,NYmax)
    Yc(NXmaxC,NYmaxC) = Y(NXmax,NYmax)

  end subroutine geom

  !****************************************************************************
  !****************************************************************************
  subroutine init_exact()
    x_exct(1)=0.00282901; v_exct(1)=-0.07028612
    x_exct(2)=6.42882E-2; v_exct(2)=2.71771E+0
    x_exct(3)=7.35542E-2; v_exct(3)=2.86215E+0
    x_exct(4)=8.09753E-2; v_exct(4)=3.02728E+0
    x_exct(5)=9.76486E-2; v_exct(5)=3.25421E+0
    x_exct(6)=1.58722E-1; v_exct(6)=3.70750E+0
    x_exct(7)=2.28897E-1; v_exct(7)=3.31348E+0
    x_exct(8)=2.34432E-1; v_exct(8)=3.25139E+0
    x_exct(9)=5.01948E-1; v_exct(9)=1.87997E-1
    x_exct(10)=8.04508E-1; v_exct(10)=-3.33066E+0
    x_exct(11)=8.59780E-1; v_exct(11)=-4.42685E+0
    x_exct(12)=9.07688E-1; v_exct(12)=-5.33693E+0
    x_exct(13)=9.44865E-1; v_exct(13)=-4.07737E+0
    x_exct(14)=9.54199E-1; v_exct(14)=-3.51971E+0
    x_exct(15)=9.61692E-1; v_exct(15)=-2.92069E+0
    x_exct(16)=9.69195E-1; v_exct(16)=-2.25968E+0
    x_exct(17)=1.00098E+0; v_exct(17)=-9.09091E-2

    y_exct(1)=-1.85185E-3; u_exct(1)=0.00000E+0
    y_exct(2)=5.00000E-2; u_exct(2)=-1.84615E+0
    y_exct(3)=5.92593E-2; u_exct(3)=-2.03077E+0
    y_exct(4)=6.66667E-2; u_exct(4)=-2.21538E+0
    y_exct(5)=1.00000E-1; u_exct(5)=-2.98462E+0
    y_exct(6)=1.68519E-1; u_exct(6)=-3.81538E+0
    y_exct(7)=2.77778E-1; u_exct(7)=-2.76923E+0
    y_exct(8)=4.48148E-1; u_exct(8)=-1.07692E+0
    y_exct(9)=4.96296E-1; u_exct(9)=-6.15385E-1
    y_exct(10)=6.09259E-1; u_exct(10)=5.84615E-1
    y_exct(11)=7.31481E-1; u_exct(11)=1.84615E+0
    y_exct(12)=8.50000E-1; u_exct(12)=3.32308E+0
    y_exct(13)=9.50000E-1; u_exct(13)=4.64615E+0
    y_exct(14)=9.57407E-1; u_exct(14)=5.07692E+0
    y_exct(15)=9.64815E-1; u_exct(15)=5.72308E+0
    y_exct(16)=9.74074E-1; u_exct(16)=6.58462E+0
    y_exct(17)=9.96296E-1; u_exct(17)=1.00000E+1
 
  end subroutine init_exact
  
  !****************************************************************************
  !****************************************************************************
  subroutine write_exact() 
    integer :: i, j

    open (23,file='Exact_V.dat') 
    !WRITE (23,*)' ZONE F=POINT, I=', 17

     	
    DO i=1,17     
      WRITE (23,*) x_exct(i), v_exct(i)/10.
    ENDDO   
    close(23)

    ! ------------------------------------------- 
    open (23,file='Exact_U.dat') 
    !WRITE (23,*)' ZONE F=POINT, I=', 17
     DO j=1,17       
       WRITE (23,*) y_exct(j), u_exct(j)/10.
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
       write(122,'(F15.7)') Xc(i,1) 
     end do
     write(122,'(A,I3,A)') "Y_COORDINATES ",ny," float"
     do j = 2, NYMAX
       write(122,'(F15.7)') Yc(1,j) 
     end do
     write(122,'(A,I3,A)') "Z_COORDINATES ",nz," float"
     write(122,'(F15.7)') 0.0 

     write(122,'(A,I6)') "POINT_DATA",nx*ny*nz
     write(122,'(A)') "VECTORS vectors float"
 
     do j = 2, NXMAX
       do i = 2, NYMAX
          write(122,'(F15.7,A,F15.7,A,F15.7)') F(i,j,1), tab, F(i,j,2), tab, 0.0 
       end do
     end do  

     close(122)

   end subroutine out_paraview_vtk
end module common
