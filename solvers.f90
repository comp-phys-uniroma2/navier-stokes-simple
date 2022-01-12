!Sample program for solving Lid-driven cavity flow using SIMPLE-algorithm
! solution of pressure-correction quation and correction U,V, and P
!Copyright (C) 2010  Michail Kiričkov
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
!You should have received a copy of the GNU General Public License
!along with this program; if not, write to the Free Software Foundation
!Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA

!**********************************************************************
module solvers
  use common
  implicit none
  private

  public :: TDMA
  public :: convergence_criteria
  public :: solve_pressure_correction   !correzioni sul campo di pressione
  public :: solve_uv   !calcolo delle componenti del campo di velocità


contains

!------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  !*********************************************************************
  ! Tri-diagonal Matrix Algorithm. Solve tri-diagonal equation system
  !-----------------------------------------------------------------------
  ! -Aw x(i-1) + Ap x(i) - Ae x(i+1) = S(i) + An x(i)[j+1] + As x(i)[j-1]
  !
  ! a(i) x(i-1) + b(i) x(i) + c(i) x(i+1) = d(i)
  !
  ! P(1) = -c(1)/b(1)
  ! Q(1) =  d(1)/b(1)
  !
  ! P(i) = -c(i) / [b(i)+a(i)P(i-1)]
  !
  ! Q(i) = [d(i)-a(i)Q(i-1)] / [b(i)+a(i)P(i-1)]
  !
  ! Back substitution:
  !
  ! x(n) = Q(n)
  ! x(i) = Q(i) + P(i) x(i+1) 
  !
  !-----------------------------------------------------------------------!
  subroutine TDMA(nf)
   integer :: nf

   real(dp) :: P(NXmaxC), Q(NXmaxC), temp, Spp
   real(dp) :: P1(NYmaxC), Q1(NYmaxC)
   integer :: i, j

   !-----------------------------------------------
   ! Sweep upward in y (increasing j-index)
   !-----------------------------------------------
   !$OMP PARALLEL DO PRIVATE(j,i)
   do j = 2, NYmax

      P(1) =  0.0_dp    
      Q(1) =  F(1,j,nf) 

      ! Forward Elimination
      do i = 2, NXmax
         temp = Ap(i,j) - Aw(i,j) * P(i-1)
         Spp  = Sp(i,j,nf) + As(i,j) * F(i,j-1,nf) + An(i,j)*F(i,j+1,nf)
         P(i) = Ae(i,j) / temp
         Q(i) = (Spp + Aw(i,j)*Q(i-1)) / temp
      enddo

      P(NXmaxC) = 0.0_dp
      Q(NXmaxC) = F(NXmaxC,j,nf)

      ! Back Substitution      
      do i = NXmax, 2, -1
         F(i,j,nf) = P(i)*F(i+1,j,nf) + Q(i)
      enddo

   enddo
   !$OMP END PARALLEL DO

   !-----------------------------------------------
   ! Sweep downward in y (decreasing j)
   !-----------------------------------------------
   !$OMP PARALLEL DO PRIVATE(j,i)
   do j =  NYmax, 2, -1
  
      P(1) =  0.0_dp
      Q(1) =  F(1,j,nf)   
      ! Forward Elimination
      do i = 2, NXmax         
         temp =  Ap(i,j) - Aw(i,j) * P(i-1)
         Spp = Sp(i,j,nf) + As(i,j) * F(i,j-1,nf) + An(i,j)*F(i,j+1,nf)
         P(i) = Ae(i,j) / temp
         Q(i) = (Spp + Aw(i,j)*Q(i-1)) / temp
      enddo

      P(NXmaxC) = 0.0_dp
      Q(NXmaxC) = F(NXmaxC,j,nf) 

      ! Back Substitution     
      do i = NXmax, 2, -1
         F(i,j,nf) = P(i)*F(i+1,j,nf) + Q(i)
      enddo

   enddo
   !$OMP END PARALLEL DO
 
  end subroutine TDMA
  
  !****************************************************************************
  !****************************************************************************
  subroutine Convergence_Criteria(nf,res_sum)
    integer :: nf
    real(dp) :: Res_sum

    
    real(dp) :: Res_vol
    integer  :: i,j

    Res_Sum = 0.0_dp

    Do j=2,NYmax
      Do i=2,NXmax
          Res_vol = Ap(i,j) * F(i  ,j  ,nf) - &
                  ( As(i,j) * F(i  ,j-1,nf) + &
                    An(i,j) * F(i  ,j+1,nf) + &
                    Aw(i,j) * F(i-1,j  ,nf) + &
                    Ae(i,j) * F(i+1,j  ,nf) ) - &
                    Sp(i,j,nf)

          Res_Sum = Res_Sum + abs(Res_vol)

       enddo
    enddo

  end subroutine convergence_criteria



  !************************************************************************
  !**********************************************************************
  ! solve_UV
  !**********************************************************************
  subroutine Solve_UV

    integer :: i, j, nf, niter 
    real(dp) :: Del_e,Del_w,Del_s,Del_n 
    real(dp) :: Area_w, Area_e, Area_s, Area_n
    real(dp) :: conv_w, conv_e, conv_s, conv_n
    real(dp) :: diff_w, diff_e, diff_s, diff_n
    real(dp) :: Res_sum_before, Res_sum_after
    real(dp) :: Pe, Pn, Pw, Ps
    real(dp) :: Vol, Dx, Dy, Urf, Delta_f


    !  calculation of fluxes
    !        N 
    !  j  +-----+----+
    !   W |  o  | E  
    ! j-1 +-----+----+
    !    i-1 S  i   i+1   
    !
    Sp = 0.0_dp

    Do i= 2, NXmax
       Do j= 2, NYmax

         Area_e = Y(j)-Y(j-1) 
         Area_w = Y(j)-Y(j-1) 
         Area_s = X(i)-X(i-1)  
         Area_n = X(i)-X(i-1) 
 
         Del_e  = Xc(i+1)-Xc(i) 
         Del_w  = Xc(i)-Xc(i-1)    
         Del_s  = Yc(j)-Yc(j-1)  
         Del_n  = Yc(j+1)-Yc(j)  
         !----------------------------------------------------


         ! upwind differencing (all other will be included into the source term) 

         ! %% Usare interpolazioni generali---------------------------------------------------
         ! Conv_e = Area_e * ((F(i+1,j,1)*(X(i)-Xc(i)) + F(i,j,1)*(Xc(i+1)-X(i)))/Del_e
         ! Conv_w = Area_w * ((F(i,j,1)*(X(i-1)-Xc(i-1)) +(Xc(i)-Xc(i-1))* F(i-1,j,1))/Del_w
         ! Conv_n = Area_n * ((F(i,j+1,2)* (Y(j)- Yc(J))) +(Yc(j+1)-Y(j))*F(i,j,2))/Del_n
         ! Conv_s = Area_s * ((F(i,j,2)*(Y(j-1)-Yc(j-1))) + (Yc(j)-Yc(j-1))*F(i,j-1,2))/Del_s

         Conv_e = Area_e *  ( F(i,j,1) + F(i+1,j  ,1) ) * 0.5_dp
         Conv_w = Area_w *  ( F(i,j,1) + F(i-1,j  ,1) ) * 0.5_dp
         Conv_n = Area_n *  ( F(i,j,2) + F(i  ,j+1,2) ) * 0.5_dp
         Conv_s = Area_s *  ( F(i,j,2) + F(i  ,j-1,2) ) * 0.5_dp

         ! Force 0 convection at the boundaries (ostacolo compreso)
         if (i.eq.2)     Conv_w = 0.0_dp
         if (i.eq.NXmax) Conv_e = 0.0_dp
         if (j.eq.2)     Conv_s = 0.0_dp
         if (j.eq.NYmax) Conv_n = 0.0_dp

         if (.not. obstacle(i,j)) then
           if ( obstacle(i-1,j) ) Conv_w = 0.0_dp
           if ( obstacle(i+1,j) ) Conv_e = 0.0_dp
           if ( obstacle(i,j-1) ) Conv_s = 0.0_dp
           if ( obstacle(i,j+1) ) Conv_n = 0.0_dp
         else 
           Conv_e = 0.0_dp
           Conv_w = 0.0_dp
           Conv_n = 0.0_dp
           Conv_s = 0.0_dp
         end if 

         Diff_e = Area_e * Gam / Del_e
         Diff_w = Area_w * Gam / Del_w 
         Diff_s = Area_s * Gam / Del_s
         Diff_n = Area_n * Gam / Del_n 

         ! Add Diffusive and Convective parts with standard Upwind
         select case(mode)
         case('C')  
           Aw(i,j) = Diff_w + Conv_w*0.5_dp  
           Ae(i,j) = Diff_e - Conv_e*0.5_dp  
           As(i,j) = Diff_s + Conv_s*0.5_dp  
           An(i,j) = Diff_n - Conv_n*0.5_dp  
         case('U') 
           Aw(i,j) = Diff_w + max( Conv_w,0.0_dp)
           Ae(i,j) = Diff_e + max(-Conv_e,0.0_dp)
           As(i,j) = Diff_s + max( Conv_s,0.0_dp)
           An(i,j) = Diff_n + max(-Conv_n,0.0_dp)
         case('H')
           Aw(i,j) = Diff_w + max( Conv_w,0.0_dp)
           Ae(i,j) = Diff_e + max(-Conv_e,0.0_dp)
           As(i,j) = Diff_s + max( Conv_s,0.0_dp)
           An(i,j) = Diff_n + max(-Conv_n,0.0_dp)
           if ( (i.GT.2).AND.(i.LT.NXmax) .and. &
                (j.GT.2).AND.(j.LT.NYmax)        ) then 
               ! NOTE: HLPA modifies Sp() !! 
               call HLPA_scheme(i,j,Conv_e,Conv_w,Conv_s,Conv_n)
           endif
         end select    

         Ap(i,j)= Aw(i,j) + Ae(i,j) + An(i,j) + As(i,j) 

     enddo
   enddo

   !------------------------------------------------------------------------
   !----------------------------- pressure gradient -----------------------
   do i = 2,NXmax
      do j = 2,NYmax

        DX = X(i) - X(i-1)   !passo dx della griglia
        DY = Y(j) - Y(j-1)   !passo dy della griglia

        VOL = DX * DY   !volume(superficie) di una cella


        ! %% Usare interpolazioni generali--------------------------------------
        ! Pe = ((F(i+1,j,4)*(X(i)-Xc(i)) + (Xc(i+1)-X(i))* F(i,j,4))/Del_e
        ! Pw = ((F(i,j,4)*(X(i-1)-Xc(i-1)) +(Xc(i)-Xc(i-1))* F(i-1,j,4))/Del_w
        ! Pn = ((F(i,j+1,4)* (Y(j)- Yc(J))) +(Yc(j+1)-Y(j))*F(i,j,4))/Del_n
        ! Ps = ((F(i,j,4)*(Y(j-1)-Yc(j-1))) + (Yc(j)-Yc(j-1))*F(i,j-1,4))/Del_s

        Pe = ( F(i,j,4) + F(i+1,j,4) ) * 0.5_dp 
        Pw = ( F(i,j,4) + F(i-1,j,4) ) * 0.5_dp
        Pn = ( F(i,j,4) + F(i,j+1,4) ) * 0.5_dp
        Ps = ( F(i,j,4) + F(i,j-1,4) ) * 0.5_dp

        if ( .not.obstacle(i,j)) then 
           if ( obstacle(i+1,j) ) Pe = F(i,j,4)
           if ( obstacle(i-1,j) ) Pw = F(i,j,4)
           if ( obstacle(i,j+1) ) Pn = F(i,j,4)
           if ( obstacle(i,j-1) ) Ps = F(i,j,4) 
        end if

        DPx_c(i,j) = (Pe-Pw)/DX
        DPy_c(i,j) = (Pn-Ps)/DY
        Sp(i,j,1) = Sp(i,j,1) - DPx_c(i,j) * VOL
        Sp(i,j,2) = Sp(i,j,2) - DPy_c(i,j) * VOL

      enddo
    end do

    !---------------------------- under-relaxation -----------------------------
    ! uij(new) = (1- alfa)* uij(old) + alfa* Ap^-1 [Sp + Te uij(new)] 
    ! Ap/alfa uij(new) = [Sp + (1-alfa)*Ap/alfa uij(old) ] + Te uij(new)
    !
    ! Note: This is a relaxation on the outer iteration loop 
    !
    Ap(:,:) = Ap(:,:) / Alfa

    Sp(:,:,1) = Sp(:,:,1) + (1.0_dp - Alfa )* Ap(:,:)*F(:,:,1) 
    Sp(:,:,2) = Sp(:,:,2) + (1.0_dp - Alfa )* Ap(:,:)*F(:,:,2) 
    !---------------------------------------------------------------------------
    write(*,*)'solve U'
    call apply_constrain(1)

    niter = 0
    call Convergence_Criteria(1,Res_sum_After)
       
    do while ( abs(Res_sum_After) > 1.d-12 .and. niter < 50)
       call TDMA(1)
       call apply_constrain(1)
       call Convergence_Criteria(1,Res_sum_After)
       niter= niter + 1
    enddo

    write(*,*) '|err|=',abs(Res_sum_After),niter


    write(*,*) 'solve V'
    call apply_constrain(2)

    niter = 0
    call Convergence_Criteria(2,Res_sum_After)

    do while ( abs(Res_sum_After) > 1.d-12 .and. niter < 50)
       call TDMA(2)
       call apply_constrain(2)
       call Convergence_Criteria(2,Res_sum_After)
       niter= niter + 1
    end do

    write(*,*) '|err|=',abs(Res_sum_After),niter
  
  end subroutine solve_UV


  !**********************************************************************
  !**********************************************************************
  subroutine Solve_Pressure_Correction()

    real(dp) :: Res_outer, Res_sum_after
    integer :: i, j, niter, relax_iter
    real(dp) :: PPe, PPs, PPw, PPn
    real(dp) :: DXc_e, S_e, Vol_e, APU_e, DPl_e, DPx_e, U_e
    real(dp) :: DYc_n, S_n, VOL_n, APV_n, DPl_n, DPy_n, V_n
    real(dp) :: Dx, Dy, Vol, Summ1, Summ2

    Do j=2,NYmax
      Do i=2,NXmax
        if (obstacle(i,j)) then
           invAp(i,j) = 0.d0
        else
           invAp(i,j) = 1.d0/Ap(i,j)
        endif
      enddo
    enddo   
   
    
    !----------------------------------------------------------------------
    !1./AP coefficient involved in the derivation of the pressure 
    !eq. and correction velocities. SIMPLE algorithm in dedicated books.
    !AP is the central coefficent there is one for U momentum equation, APU 
    !and one for V-momentum, APV.
    ! APU(i,j) = AP(i,j,1)
    ! APV(i,j) = AP(i,j,2)
    !---------------------------------------------------------------------
    ! Computing Ue Uw
    !
    ! Ue(i,j) is stored in Con_e(i,j)
    ! Uw(i,j) is stored in Con_e(i-1,j)
    ! 
    !---------------------------------------------------------------------
    do i=2,NXmax-1
       do j=2,NYmax

         ! distance between node P and node E:
         DXc_e = Xc(i+1) - Xc(i)
         S_e   =  Y(j) -  Y(j-1)
         VOL_e = DXc_e * S_e 

         !---RHIE AND CHOW ------------------------------------------
         ! Ul_e = ((F(i+1,j,1)*(X(i)-Xc(i)) + (Xc(i+1)-X(i))* F(i,j,1))/Del_e
         U_e  =       0.5_dp * (     F(i,j,1) +     F(i+1,j,1) )
         if (rhiechow) then
           APU_e =       0.5_dp * (   invAp(i,j)   +   invAp(i+1,j)   )
           DPl_e =       0.5_dp * ( DPx_c(i,j)   + DPx_c(i+1,j)   )
           DPx_e =       ( F(i+1,j,4) - F(i,j,4) ) / DXc_e
           U_e   = U_e  + APU_e * VOL_e * ( DPl_e - DPx_e)
         end if  
         !--------------------------------------------------------------
         Con_e(i,j) =   U_e  * S_e

         if (.not.obstacle(i,j)) then
           if (obstacle(i+1,j)) Con_e(i,j) = 0.0_dp
           if (obstacle(i-1,j)) Con_e(i-1,j) = 0.0_dp
         else
           Con_e(i,j) = 0.0_dp
           Con_e(i-1,j) = 0.0_dp
         end if

      enddo
    enddo

    Con_e(1,:) = F(1,:,1) * S_e
    Con_e(NXmax,:) = F(NXmax+1,:,1) * S_e

    !---------------------------------------------------------------------
    ! Computing Vn Vs
    !
    ! Vn(i,j) is stored in Con_n(i,j)
    ! Vs(i,j) is stores in Con_n(i,j-1)
    ! 
    !---------------------------------------------------------------------
    !-------------------------------------------------------------------
    ! horizontal  faces
    do i=2,NXmax
      do j=2,NYmax-1
 
        DYc_n = Yc(j+1) - Yc(j)
        S_n   =  X(i) -  X(i-1)
        
        VOL_n = DYc_n * S_n

        !---RHIE AND CHOW ------------------------------------------
        !V_n = ((F(i,j+1,2)* (Y(j)- Yc(J))) +(Yc(j+1)-Y(j))*F(i,j,2))/Del_n
        V_n  =  0.5_dp * (     F(i,j,2) +     F(i,j+1,2) )
        if (rhiechow) then
          APV_n =  0.5_dp * (   invAp(i,j)   +   invAp(i,j+1)   )
          DPl_n =  0.5_dp * ( DPy_c(i,j)   + DPy_c(i,j+1)   )
          DPy_n = ( F(i,j+1,4) - F(i,j,4) ) / DYc_n
          V_n   = V_n + APV_n * VOL_n * ( DPl_n - DPy_n ) 
        end if
        !----------------------------------------------------------- 
        Con_n(i,j) =   V_n  * S_n 
        
        if (.not.obstacle(i,j)) then
          if (obstacle(i+1,j)) Con_n(i,j) = 0.0_dp
          if (obstacle(i-1,j)) Con_n(i-1,j) = 0.0_dp
        else
          Con_n(i,j) = 0.0_dp
          Con_n(i-1,j) = 0.0_dp
        end if
      enddo
    enddo

    Con_n(:,1) =  F(:,1,2) * S_n
    Con_n(:,NYmax) = F(:,NYmax+1,2) * S_n

    !-----------------------------------------------------------------------
    ! Pressure Poisson Equation. The rhs is computed from u/v
    !
    ! A_p p'_C = A_e p'_E + A_n p'_N + A_w p'_W + A_s p'_S - D*
    !       D* = S_e (u_e - u_w) + S_n (v_n - v_s) 
    !
    !-----------------------------------------------------------------------
    do i=2,NXmax
      do j=2,NYmax

         S_e   =  Y(j) - Y(j-1)
         S_n   =  X(i) - X(i-1)
             
         ! %% Usare interpolazioni generali-----------------------------------------------------------
         ! Ae(i,j) = S_e * S_e* ((invAp(i+1,j)*(X(i)-Xc(i)) + (Xc(i+1)-X(i))*invAp(i,j))/Del_e
         ! Aw(i,j) = S_e * S_e* ((invAp(i,j)*(X(i-1)-Xc(i-1)) + (Xc(i)-Xc(i-1))*invAp(i-1,j))/Del_w
         ! An(i,j) = S_n * S_n* ((invAp(i,j+1)* (Y(j)- Yc(J))) + (Yc(j+1)-Y(j))*invAp(i,j))/Del_n
         ! As(i,j) = S_n * S_n* ((invAp(i,j)*(Y(j-1)-Yc(j-1))) + (Yc(j)-Yc(j-1))*invAp(i,j-1))/Del_s
         
         An(i,j) = 0.5_dp * (   invAp(i,j)   +   invAp(i,j+1)   )  * S_n * S_n
         As(i,j) = 0.5_dp * (   invAp(i,j)   +   invAp(i,j-1)   )  * S_n * S_n
         Ae(i,j) = 0.5_dp * (   invAp(i,j)   +   invAp(i+1,j)   )  * S_e * S_e
         Aw(i,j) = 0.5_dp * (   invAp(i,j)   +   invAp(i-1,j)   )  * S_e * S_e
 
         Ap(i,j) = (An(i,j) + As(i,j) + Ae(i,j) + Aw(i,j))
   
         ! Formula corretta sarebbe:
         ! - S_e (u_e - u_w) - S_n (v_n - v_s) 
         ! Con_e(i-1) = u_w(i) 
         Sp(i,j,3) = - ( Con_e(i,j) - Con_e(i-1,j) + &
                         Con_n(i,j) - Con_n(i,j-1) )


      enddo
    enddo
    !-------------------------------------------------------------------
    ! Under-relaxation (cannot work here!)
    !Ap(:,:) = Ap(:,:) / Omega
    !Sp(:,:,3) = Sp(:,:,4) + (1.0_dp - Omega )* Ap(:,:)*F(:,:,3) 
    !-------------------------------------------------------------------
    write(*,*) 'solve PP'

    call apply_constrain(3)
    Res_sum_After=1.0_dp
    niter = 0
   
    do while ( abs(Res_sum_After) > Poiss_acc .and. niter < Poiss_iter ) 
       call TDMA(3)
       call apply_constrain(3)
       call Convergence_Criteria(3,Res_sum_After)
       niter = niter + 1
    enddo

    write(*,*) '|err|=',abs(Res_sum_After), niter
 
    !--------------------------------------------------------------------
    !--------------------------------------------------------------------
    !  velocities and pressure correction

    Summ1 = 0.0_dp
    Summ2 = 0.0_dp
    Do i=2,NXmax
       Do j=2,NYmax  
     
        Dy = Y(j)-Y(j-1)
        Dx = X(i)-X(i-1)
    
        ! %% Usare interpolazioni generali-------------------------------------------
        ! PPe = ((F(i+1,j,3)*(X(i)-Xc(i)) + F(i,j,3)*(Xc(i+1)-X(i)))/Del_e
        ! PPw = ((F(i,j,3)*(X(i-1)-Xc(i-1)) +(Xc(i)-Xc(i-1))* F(i-1,j,3))/Del_w
        ! PPn = ((F(i,j+1,3)* (Y(j)- Yc(J))) +(Yc(j+1)-Y(j))*F(i,j,3))/Del_n
        ! PPs = ((F(i,j,3)*(Y(j-1)-Yc(j-1))) + (Yc(j)-Yc(j-1))*F(i,j-1,3))/Del_s
      
        PPe = 0.5_dp * ( F(i,j,3) + F(i+1,j,3) )
        PPw = 0.5_dp * ( F(i,j,3) + F(i-1,j,3) )
        PPn = 0.5_dp * ( F(i,j,3) + F(i,j+1,3) )
        PPs = 0.5_dp * ( F(i,j,3) + F(i,j-1,3) )

        ! override with a Neumann BC near obstacles
        if ( .not.obstacle(i,j) ) then 
           if (obstacle(i+1,j)) PPe = F(i,j,3)
           if (obstacle(i-1,j)) PPw = F(i,j,3)
           if (obstacle(i,j+1)) PPn = F(i,j,3)
           if (obstacle(i,j-1)) PPs = F(i,j,3) 
        end if

        ! Update velocities
        ! u(k+1) = u(*) - dP/dx  
        if (.not.obstacle(i,j)) then

           F(i,j,1) = F(i,j,1) - (PPe - PPw) * invAp(i,j) * Dy 
           F(i,j,2) = F(i,j,2) - (PPn - PPs) * invAp(i,j) * Dx
      
           ! Use a relaxation of beta this is also a outer loop relaxation 
           ! p(k+1) = p(k) + p' 
           ! p(k+1) = (1-beta) p(k) + beta (p(k) + p') = p(k) + beta p'
           !
           F(i,j,4) = F(i,j,4) + beta * F(i,j,3)
           
        end if
 
        Summ1 = Summ1 + F(i,j,3) 
        Summ2 = Summ2 + F(i,j,4) 

      enddo
    enddo

    write(*,*) 
    write(*,*) "p=",summ2,"p'=",summ1

    !----------------------------------------------------------------------
    ! Neumann BC for pressure
    F(:,1,4) = F(:,2,4)               
    F(:,NYmaxC,4) = F(:,NYmaxC-1,4)    
    F(1,:,4) = F(2,:,4)             
    F(NXmaxC,:,4) = F(NXmaxC-1,:,4)
    !do i = 1, NXMaxC
    !  do j = 1, NYMaxC 
    !    if (obstacle(i,j)) then
    !      if (obstacle(i+1,j)) F(i,j,4)=F(i+1,j,4)
    !      if (obstacle(i-1,j)) F(i,j,4)=F(i-1,j,4)
    !      if (obstacle(i,j+1)) F(i,j,4)=F(i,j+1,4)
    !      if (obstacle(i,j-1)) F(i,j,4)=F(i,j-1,4)
    !    end if
    !  end do
    !end do
    !----------------------------------------------------------------------   
    ! Pressure undefined to a constant. Fix at least 1 value of pressure
    F(:,:,4) = F(:,:,4) - F(3,4,4) 


  end subroutine solve_pressure_correction  



  subroutine apply_constrain(nf)
    integer :: nf


    integer :: i,j

    do i = 2, NXmax
       do j = 2, NYmax
          if (obstacle(i,j)) then
             Ap(i,j) = 1.0_dp
             Ae(i,j) = 0.0_dp
             Aw(i,j) = 0.0_dp
             An(i,j) = 0.0_dp
             As(i,j) = 0.0_dp
             F(i,j,nf) = 0.0_dp
             Sp(i,j,nf) = 0.0_dp !F(i,j,nf)
          end if   
       end do
    end do


  end subroutine apply_constrain


  !************************************************************************
  !**********************************************************************
  ! HLPA scheme (Hybrid Linear / Parabolic Approximation)
  ! Zhu J Rodi W Comp. Meth. Appl. Mat. & Eng. 92, 87-96 (1991)
  ! Approximate convection term 
  !**********************************************************************
  subroutine HLPA_scheme(i,j,Conv_e,Conv_w,Conv_s,Conv_n)
     integer :: i, j  
     real(dp) :: Conv_e,Conv_w,Conv_s,Conv_n
         
     real(dp) :: Fww, Fw, Fp, Fe, Delta_f
     integer :: nf

     DO nf = 1, 2
  
        !------------------ W face -------------------
        Fww = F(i-2,j,nf)
        Fw  = F(i-1,j,nf)
        Fp  = F(i  ,j,nf)
        Fe  = F(i+1,j,nf)
  
        call  HLPA(Conv_w,Fww,Fw,Fp,Fe,Delta_f)
  
        Sp(i,j,nf) = Sp(i,j,nf) + Conv_w * Delta_f
  
        !------------------ E face--------------------
        Fww = F(i-1,j,nf)
        Fw  = F(i  ,j,nf)
        Fp  = F(i+1,j,nf)
        Fe  = F(i+2,j,nf)
  
        call  HLPA(Conv_e,Fww,Fw,Fp,Fe,Delta_f)
  
        Sp(i,j,nf) = Sp(i,j,nf) - Conv_e * Delta_f 
  
        !------------------ S face--------------------
        Fww = F(i  ,j-2,nf)
        Fw  = F(i  ,j-1,nf)
        Fp  = F(i  ,j  ,nf)
        Fe  = F(i  ,j+1,nf)
  
        call  HLPA(Conv_s,Fww,Fw,Fp,Fe,Delta_f)
  
        Sp(i,j,nf) = Sp(i,j,nf) + Conv_s * Delta_f
   
        !------------------ N face--------------------
        Fww = F(i  ,j-1,nf)
        Fw  = F(i  ,j  ,nf)
        Fp  = F(i  ,j+1,nf)
        Fe  = F(i  ,j+2,nf)
  
        call  HLPA(Conv_n,Fww,Fw,Fp,Fe,Delta_f)
  
        Sp(i,j,nf) = Sp(i,j,nf) - Conv_n * Delta_f 
  
     end do 
  
  end subroutine HLPA_scheme
  !************************************************************************
  ! utility subroutine for HLPA
  subroutine HLPA(Uw,Fww,Fw,Fp,Fe,Delta_f)
    real(dp) :: Uw, Fww, Fw, Fp, Fe, Delta_f
   
    real(dp) :: Alpha

    Delta_f = 0.0_dp


    if (Uw.GE.0.0_dp) then

      if ( ABS( Fp - 2.0_dp * Fw + Fww ).LT.ABS( Fp - Fww ) ) then
         Alpha = 1.0_dp
      else
         Alpha = 0.0_dp
      end if

      if ( (Fp - Fww).NE.0.0_dp ) then 
         Delta_f = Alpha * (Fp - Fw)* (Fw - Fww) / (Fp - Fww)
      endif

    else

      if ( ABS( Fw - 2.0_dp* Fp + Fe ).LT.ABS( Fw - Fe ) ) then
         Alpha = 1.0_dp
      else
         Alpha = 0.0_dp
      endif
 
      if ( (Fw - Fe).NE.0.0_dp ) then   
         Delta_f = Alpha * (Fw - Fp)* (Fp - Fe) / (Fw - Fe)
      end if

    end if

  end subroutine HLPA
  
end module solvers
