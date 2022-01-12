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
  character(1) :: ch

  if(iargc()<3) then
   print*,'ns Re Niter CDS|UPD|HLPA [filename=* alfa=* beta=* Poiss_acc=* Poiss_iter=* RhieChow=Y|N]'
   stop  
  endif
 
  call getarg(1,arg)
  read(arg,*) Re
  
  call getarg(2,arg)
  read(arg,*) Niter

  ! Select computation mode
  call getarg(3,arg)
  mode = 'U'

  if (trim(arg)=='CDS') mode = 'C'
  if (trim(arg)=='UPD') mode = 'U'
  if (trim(arg)=='HLPA') mode = 'H'

  readstart = .false.
  filename = ''
  alfa = 0.85_dp
  beta = 0.10_dp
  Poiss_iter = 500
  Poiss_acc = 1.d-7
  ch = 'Y'

  do i = 4, iargc()
    call getarg(i,arg)

    if (arg(1:9)=='filename=') then
      read(arg(10:),*) filename
      if (trim(filename) /= '') then
         readstart = .true.
      endif  
    elseif (arg(1:5)=='alfa=') then
      read(arg(6:),*) alfa
    elseif (arg(1:5)=='beta=') then
      read(arg(6:),*) beta 
    elseif (arg(1:11)=='Poiss_iter=') then
      read(arg(12:),*) Poiss_iter 
    elseif (arg(1:10)=='Poiss_acc=') then
      read(arg(11:),*) Poiss_acc
    elseif (arg(1:9)=='RhieChow=') then
      read(arg(10:),*) ch
    endif
  end do

  if (ch=='N') then 
    rhiechow = .false.
  elseif (ch=='Y') then
    rhiechow = .true.
  end if

  print*, filename
  print*, 'Re=',Re
  print*, 'alfa=',alfa
  print*, 'beta=',beta
  print*, 'Poiss_acc=',Poiss_acc 
  print*, 'Poiss_iter=',Poiss_iter 
  print*, 'Rhie&Chow=',ch 

  print*, 'Init Grid'
  call Init_grid() 

  print*, 'Init Cells'
  call init_cells()
  
  print*, 'Allocate workspace'
  call init_workspace(NXmaxC,NYmaxC)

  print*, 'Init all boundary conditions'
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

  call init_exact()
  call write_exact()

  !------------------------------------------------------
  open (23,file='Profiles_V.dat') 
    DO i=1,NXmaxC 
        WRITE (23,*) Xc(i), F(i,(NYmax+1)/2,2)
    ENDDO
  close(23) 
  
  !--------------------------------------------------------
  open (24,file='Profiles_U.dat') 
    DO j=1,NYmaxC 
      WRITE (24,*) Yc(j), F((NXmax+1)/2,j,1)
    ENDDO
  close(24) 

  !********************************************************************************
  !********************************************************************************
  !	Open(10,file='binary_dat_all.dat') 
  !	  write(10,*)F(:,:,1:Nvars)
  !	close(10)
  !----------------------------------------------------------------
  
  !call Out_array(F(:,:,1),NXmaxC,NYmaxC,'1_U_s.txt ')
  !call Out_array(F(:,:,2),NXmaxC,NYmaxC,'1_V_s.txt ')

  !-------------------------------------------------------------------
  !----------------------------------------------------------------
  ! NImax = NXmaxC
  ! NJmax = NYmaxC
  ! F_out     = F(:,:,5)
  ! Filename  ='1_T_s.txt' 
  ! call Out_array(F_out,NImax,NJmax,Filename)
  !-------------------------------------------------------------------
  !--------------------------------------------------------------------------
  ! This saves for gnuplot
  ! use plot 'flux.dat' using 1:2:3:4 with vectors !per visualizzare su gnuplot
  open (23,file='flux.dat')       

  DO J=1, NYmaxC
    DO I=1, NXmaxC
       WRITE (23,*) Xc(I), Yc(J), F(i,j,1), F(i,j,2), F(i,j,4)   
    ENDDO
  ENDDO
  close(23) 
  !--------------------------------------------------------------------------

  !-------------------------------------------------------------------
  ! Plot data for Paraview 
  !--------------------------------------------------------------------------
  call out_paraview_vtk('UV.vtk')


  WRITE(*,*) 'all done'

end program main
