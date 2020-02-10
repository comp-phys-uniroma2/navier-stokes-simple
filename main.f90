!Sample program for solving Lid-Driven cavity test using SIMPLE-algorithm
! Main modul
!Copyright (C) 2010  Michail Kiriakov

!This program is free software; you can redistribute it and/or
!modify it under the terms of the GNU General Public License
!as published by the Free Software Foundation; either version 2
!of the License, or (at your option) any later version.

!This program is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.

!You should have received a copy of the GNU General Public License
!along with this program; if not, write to the Free Software
!Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

!-----------------------------------------------------------------------
Program Main
  use common
  use solvers 
  implicit none
   
  integer :: Niter, Niter_last, Niter_print, j, i, k
  character(64) :: arg
 
  if(iargc()<1) then
   print*,'ns Niter [filename]'
   stop  
  endif
 
  call getarg(1,arg)
  read(arg,*) Niter
 
  filename = ''
  if (iargc()>1) then
    call getarg(2,arg)
    read(arg,*) filename
  endif  
  if (trim(filename) /= '') then
    readstart = .true.
  else
    readstart=.false.
  endif  
  
  ! allocate all arrays with size NXmaxC=NXmax+1
  print*, 'Allocate workspace'
  call init(NXmaxC,NYmaxC)

  print*, 'Init Grid'
  call Init_grid() 

  print*, 'Init Geom'
  call Geom()

  print*, 'Init all'
  call Init_all_cavity()

  
  ! ------------------------------------------------
  ! Open(10,file='binary_dat_all.dat') 
  ! 	   read(10,*) F(:,:,1:Nvars)
  ! close(10)
  ! ------------------------------------------------
  !Niter=100
  !mNiter_print = 100 

  Do k = 1, Niter

    write(*,*) '***************',k,'*****************'
 
    call Solve_UV
    call Solve_Pressure_Correction
  
  ENDDO

  print*, 'Write output'

  ! write 'Exact_U.dat' and 'Exact_V.dat'
  !print*, 'Init exact'
  call init_exact
  call write_exact()

  !------------------------------------------------------
  open (23,file='Profiles_V.dat') 
    DO i=1,NXmaxC 
        WRITE (23,*) Xc(i,61), F(i,61,2)
    ENDDO
  close(23) 
  !--------------------------------------------------------
  open (24,file='Profiles_U.dat') 
    DO j=1,NYmaxC 
      WRITE (24,*) Yc(61,j), F(61,j,1)
    ENDDO
  close(24) 

  !********************************************************************************
  !********************************************************************************
  
  !call Out_array(F(:,:,1),NXmaxC,NYmaxC,'1_U_s.txt ')
  !call Out_array(F(:,:,2),NXmaxC,NYmaxC,'1_V_s.txt ')
  

  !-------------------------------------------------------------------
  ! Plot data for gnuplot
  !--------------------------------------------------------------------------
  open (23,file='Domain_all.dat') 
  WRITE(23,*)'VARIABLES = "Xp", "Yp" , "Up" , "Vp" , "Pp" ' 
  WRITE(23,*)' ZONE I=' ,NXmaxC, ', J=', NYmaxC, ', F=POINT'

  DO J=1, NYmaxC
    DO I=1, NXmaxC
       WRITE (23,*) Xc(I,J), Yc(I,J) , F(i,j,1) , F(i,j,2) , F(i,j,4)   
    ENDDO
  ENDDO

  close(23) 
  !--------------------------------------------------------------------------

  !-------------------------------------------------------------------
  ! Plot data for Paraview 
  !--------------------------------------------------------------------------
  call out_paraview_vtk('UV.vtk')

  !open (23,file='UV.raw') 
  !  WRITE (23,*) F(:,:,1)    
  !  WRITE (23,*) F(:,:,2)    
  !  WRITE (23,*) F(:,:,4)    
  !close(23)

  WRITE(*,*) 'all done'

end program main
