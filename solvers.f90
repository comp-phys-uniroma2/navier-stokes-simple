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
  public :: solve_pressure_correction
  public :: solve_uv


  contains

  !*********************************************************************
  ! Tri-diagonal Matrix Algorithm. Solve tri-diagonal equation system
  !
  subroutine TDMA(nf)
   integer :: nf

   real(dp) :: P(NXmaxC), Q(NXmaxC), temp, Spp
   integer :: i, j

   !------------ Upward sweep ------------------------------------------
   Do j = 2, NYmax

      P(1) =  0.0_dp    
      Q(1) =  F(1,j,nf) 

      ! Forward Elimination
      Do i = 2, NXmax
         temp = Ap(i,j,nf) - Aw(i,j) * P(i-1)
         Spp  = Sp(i,j,nf) + As(i,j) * F(i,j-1,nf) + An(i,j) * F(i,j+1,nf) 
         P(i) = Ae(i,j) / temp
         Q(i) = (Spp + Aw(i,j)*Q(i-1)) / temp
      ENDDO

      P(NXmaxC) = 0.0_dp
      Q(NXmaxC) = F(NXmaxC,j,nf)

      ! Back Substitution      
      Do i = NXmax,2,-1
          F(i,j,nf) = P(i)*F(i+1,j,nf) + Q(i)
      ENDDO

   ENDDO 

   !------------- Backward sweep -------------------------------------------  
   DO j =  NYmax,2,-1
  
      P(1) =  0.0_dp
      Q(1) =  F(1,j,nf)   

      P(NXmaxC) = 0.0_dp
      Q(NXmaxC) = F(NXmaxC,j,nf) 

      ! Forward Elimination
      DO i = 2, NXmax         
         temp =  Ap(i,j,nf) - Aw(i,j) * P(i-1)
         Spp= Sp(i,j,nf) + As(i,j) * F(i,j-1,nf) + An(i,j) * F(i,j+1,nf)  
         P(i) = Ae(i,j) / temp
         Q(i) = (Spp + Aw(i,j)*Q(i-1)) / temp
      ENDDO

      ! Back Substitution     
      DO i = NXmax,2,-1
         F(i,j,nf) = P(i)*F(i+1,j,nf) + Q(i)
      ENDDO

    ENDDO
  
   end subroutine TDMA
  
  !****************************************************************************
  ! Function checks convergency on different variables depending on input nf
  ! The soubroutine checks fullfilling the equation:
  ! Ap uC - Ae uE  - An uN  - Aw uW  - As uS - Sp = 0 
  !****************************************************************************
  subroutine Convergence_Criteria(nf,res_sum)
    integer :: nf
    real(dp) :: Res_sum

    
    real(dp) :: Res_vol
    integer  :: i,j

    Res_Sum = 0.0_dp

    Do j=2,NYmax
      Do i=2,NXmax
          Res_vol = Ap(i,j,nf) * F(i  ,j  ,nf) - &
                  ( As(i,j) * F(i  ,j-1,nf) + &
                    An(i,j) * F(i  ,j+1,nf) + &
                    Aw(i,j) * F(i-1,j  ,nf) + &
                    Ae(i,j) * F(i+1,j  ,nf) ) - &
                    Sp(i,j,nf)

          Res_Sum = Res_Sum + Res_vol

       ENDDO
    ENDDO

  end subroutine convergence_criteria



  !************************************************************************
  !**********************************************************************
  ! solve_UV
  !**********************************************************************
  subroutine Solve_UV

    integer :: i, j, nf, niter 
    real(dp) :: Gam_e,Gam_w,Gam_s,Gam_n 
    real(dp) :: Del_e,Del_w,Del_s,Del_n 
    real(dp) :: Area_w, Area_e, Area_s, Area_n
    real(dp) :: conv_w, conv_e, conv_s, conv_n
    real(dp) :: diff_w, diff_e, diff_s, diff_n
    real(dp) :: Res_sum_before, Res_sum_after
    real(dp) :: Pe, Pn, Pw, Ps
    real(dp) :: Vol, Dx, Dy, Alfa, Urf, Delta_f


    ! Problem Assembly
    ! --------------------------------------------------------------------
    !  calculation of fluxes
    !        N 
    !  j  +-----+----+
    !   W |  o  | E  
    ! j-1 +-----+----+
    !    i-1 S  i   i+1   
    !
    Do i= 2,NXmax
      Do j= 2,NYmax
         Gam_e = ( Gam(i+1,j  ) + Gam(i  ,j  ) ) * 0.5_dp
         Gam_w = ( Gam(i-1,j  ) + Gam(i  ,j  ) ) * 0.5_dp
         Gam_s = ( Gam(i  ,j-1) + Gam(i  ,j  ) ) * 0.5_dp
         Gam_n = ( Gam(i  ,j+1) + Gam(i  ,j  ) ) * 0.5_dp

         Area_e = Y(i  ,j)-Y(i  ,j-1) 
         Area_w = Y(i-1,j)-Y(i-1,j-1) 
         Area_s = X(i,j-1)-X(i-1,j-1)  
         Area_n = X(i,j  )-X(i-1,j  ) 
 
         Del_e  = Xc(i+1,j)-Xc(i  ,j) 
         Del_w  = Xc(i  ,j)-Xc(i-1,j)    
         Del_s  = Yc(i,j  )-Yc(i,j-1)  
         Del_n  = Yc(i,j+1)-Yc(i,j  )  
         !----------------------------------------------------

         ! upwind differencing (all other will be included into the source term) 
         Conv_w = Area_w *  ( F(i,j,1) + F(i-1,j  ,1) ) * 0.5_dp
         Conv_e = Area_e *  ( F(i,j,1) + F(i+1,j  ,1) ) * 0.5_dp
         Conv_s = Area_s *  ( F(i,j,2) + F(i  ,j-1,2) ) * 0.5_dp
         Conv_n = Area_n *  ( F(i,j,2) + F(i  ,j+1,2) ) * 0.5_dp

         ! Force 0 convection at the boundaries 
         if (i.eq.2)     Conv_w = 0.0_dp
         if (i.eq.NXmax) Conv_e = 0.0_dp
         if (j.eq.2)     Conv_s = 0.0_dp
         if (j.eq.NYmax) Conv_n = 0.0_dp

         Diff_e = Area_e * Gam_e / Del_e
         Diff_w = Area_w * Gam_w / Del_w 
         Diff_s = Area_s * Gam_s / Del_s
         Diff_n = Area_n * Gam_n / Del_n 

         ! Add Diffusive and Convective parts with standard Upwind 
         Aw(i,j) = Diff_w + max( Conv_w,0.0_dp)
         Ae(i,j) = Diff_e + max(-Conv_e,0.0_dp)
         As(i,j) = Diff_s + max( Conv_s,0.0_dp)
         An(i,j) = Diff_n + max(-Conv_n,0.0_dp)

         Ap(i,j,1:2)= Aw(i,j) + Ae(i,j) + An(i,j) + As(i,j) 
         Sp(i,j,1:2) = 0.0_dp

         if ( (i.GT.2).AND.(i.LT.NXmax) .and. &
              (j.GT.2).AND.(j.LT.NYmax)        ) then 
            call HLPA_scheme(i,j,Conv_e,Conv_w,Conv_s,Conv_n)
         endif

     enddo
   enddo

   !------------------------------------------------------------------------
   !----------------------------- pressure gradient -----------------------
   do i = 2,NXmax
     do j = 2,NYmax

        DX = X(i,j) - X(i-1,j)
        DY = Y(i,j) - Y(i,j-1)

        VOL = DX * DY

        Pe = ( F(i,j,4) + F(i+1,j,4) ) * 0.5_dp
        Pw = ( F(i,j,4) + F(i-1,j,4) ) * 0.5_dp
        Pn = ( F(i,j,4) + F(i,j+1,4) ) * 0.5_dp
        Ps = ( F(i,j,4) + F(i,j-1,4) ) * 0.5_dp

        DPx_c(i,j) = (Pe-Pw)/DX
        DPy_c(i,j) = (Pn-Ps)/DY
        Sp(i,j,1) = Sp(i,j,1) - DPx_c(i,j) * VOL
        Sp(i,j,2) = Sp(i,j,2) - DPy_c(i,j) * VOL

      enddo
    end do

    !---------------------------- under-relaxation -----------------------------
    ! uij = (1-alfa)* uij(old) + alfa* Ap^-1 [Sp + Te uij] 
    !
    ! (Ap/alfa) uij = Sp + (1-alfa) * (Ap/alfa) uij(old)  + Te uij
    !  Ap' uij = Sp' + Te uij
    !--------------------------------------------------------------------------
    Alfa = 0.85_dp  

    ! Note: Rescaling Ap by 1/alfa 
    Ap(:,:,1:2) = Ap(:,:,1:2) / Alfa

    Sp(:,:,1) = Sp(:,:,1) + (1.0_dp - Alfa )* Ap(:,:,1)*F(:,:,1) 
    Sp(:,:,2) = Sp(:,:,2) + (1.0_dp - Alfa )* Ap(:,:,2)*F(:,:,2) 

    !---------------------------------------------------------------------------
    write(*,*)'solve U'

    Res_sum_After=1.d0
    
    niter = 0

    do while ( abs(Res_sum_After) > 1.d-12 .and. niter < 50)
       call TDMA(1)
       call Convergence_Criteria(1,Res_sum_After)
       niter= niter + 1
    enddo
    write(*,*) '|err|=',abs(Res_sum_After),niter

    write(*,*)'solve V' 

    Res_sum_After=1.d0

    niter = 0

    do while ( abs(Res_sum_After) > 1.d-12 .and. niter < 50)
       call TDMA(2)
       call Convergence_Criteria(2,Res_sum_After)
       niter= niter + 1
    end do              

    write(*,*) '|err|=',abs(Res_sum_After),niter
    
    ! Restore Ap ?!!??!?
    !Ap = alfa * Ap  

  end subroutine solve_UV


  !**********************************************************************
  !**********************************************************************
  subroutine Solve_Pressure_Correction()

    real(dp) :: Res_sum_before, Res_sum_after
    integer :: i, j, niter
    real(dp) :: PPe, PPs, PPw, PPn
    real(dp) :: DXc_e, S_e, Vol_e, APU_e, Ul_e, DPl_e, DPx_e, U_e
    real(dp) :: DYc_n, S_n, VOL_n, APV_n, Vl_n, DPl_n, DPy_n, V_n
    real(dp) :: Dx, Dy, Vol, Summ
    real(dp) :: alfa, beta

    DpU(2:NXmax,2:NYmax) = 1.d0/Ap(2:NXmax,2:NYmax,1)
    DpV(2:NXmax,2:NYmax) = 1.d0/Ap(2:NXmax,2:NYmax,2)

    !----------------------------------------------------------------------
    !1./AP coefficient involved in the derivation of the pressure 
    !eq. and correction velocities. SIMPLE algorithm in dedicated books.
    !AP is the central coefficent there is one for U momentum equation, APU 
    !and one for V-momentum, APV.
    ! APU(i,j) = AP(i,j,1)
    ! APV(i,j) = AP(i,j,2)
    !---------------------------------------------------------------------
  
    ! vertical faces are all stored on Con_e (East)
    Do i=2,NXmax-1
      Do j=2,NYmax

   
        !If ((i.ne.1).and.(i.ne.NXmax)) then

         ! distance between node P and node E:
         DXc_e = Xc(i+1,j) - Xc(i,j  ) 
         ! Est surface S_e
         S_e   =  Y(i  ,j) -  Y(i,j-1)
         VOL_e = DXc_e * S_e 

         !-------------------------------------------------------------
         ! Iterpolation of 1/APU at east location:
         APU_e =       0.5_dp * (   DpU(i,j)   +   DpU(i+1,j)   )
         ! Interpolation of velocity U on east location:
         Ul_e  =       0.5_dp * (     F(i,j,1) +     F(i+1,j,1) )
         ! the same for interpolation of DPx:
         DPl_e =       0.5_dp * ( DPx_c(i,j)   + DPx_c(i+1,j)   )
         ! Gradient of pressure at E location using central difference:
         DPx_e =       ( F(i+1,j,4) - F(i,j,4) ) / DXc_e
         ! Interpolated velocity at E location using Rhie&Chow interp: 
         U_e   = Ul_e  + APU_e * VOL_e * ( DPl_e - DPx_e)
         ! Convective flux.    
         Con_e(i,j) =   U_e  * S_e
         !--------------------------------------------------------------
         
        !End If
 

      enddo
    enddo
        
    Con_e(1,:) = 0.0_dp
    Con_e(NXmax,:) = 0.0_dp
    !-------------------------------------------------------------------
    !call Out_array(Con_e(:,:),NXmaxC,NYmaxC,'Con_e.dat ')

    ! horizontal  faces are all stored on Con_n 
    Do i=2,NXmax
     Do j=2,NYmax-1
       
        !If ((j.ne.1).and.(j.ne.NYmax)) then
           DYc_n = Yc(i,j+1) - Yc(i  ,j)
           S_n   =  X(i,j  ) -  X(i-1,j)

           VOL_n = DYc_n * S_n 

           !-----------------------------------------------------------
           APV_n =  0.5_dp * (   DpV(i,j)   +   DpV(i,j+1)   )
           Vl_n  =  0.5_dp * (     F(i,j,2) +     F(i,j+1,2) )
           DPl_n =  0.5_dp * ( DPy_c(i,j)   + DPy_c(i,j+1)   )
                 
           DPy_n = ( F(i,j+1,4) - F(i,j,4) ) / DYc_n
  
           V_n   = Vl_n  + APV_n * VOL_n * ( DPl_n - DPy_n )    
           !----------------------------------------------------------- 
       
           Con_n(i,j) =   V_n  * S_n 
        !End If
   
      ENDDO
    ENDDO
        
    Con_n(:,1) = 0.0_dp
    Con_n(:,NYmax) = 0.0_dp
    !call Out_array(Con_n(:,:),NXmaxC,NYmaxC,'Con_n.dat ')

    !-----------------------------------------------------------------------
    ! Assembly of PPE 
    !-----------------------------------------------------------------------
    do i=2,NXmax
      do j=2,NYmax

         S_e   =  Y(i,j) - Y(i,j-1)
         S_n   =  X(i,j) - X(i-1,j)
             
         An(i,j) = 0.5_dp * (   DpV(i,j)   +   DpV(i,j+1)   )  * S_n * S_n
         As(i,j) = 0.5_dp * (   DpV(i,j)   +   DpV(i,j-1)   )  * S_n * S_n
         Ae(i,j) = 0.5_dp * (   DpU(i,j)   +   DpU(i+1,j)   )  * S_e * S_e
         Aw(i,j) = 0.5_dp * (   DpU(i,j)   +   DpU(i-1,j)   )  * S_e * S_e
 
         Ap(i,j,3) = (An(i,j) + As(i,j) + Ae(i,j) + Aw(i,j))
   
         Sp(i,j,3) = - ( Con_n(i,j) - Con_n(i,j-1) + &
                         Con_e(i,j) - Con_e(i-1,j) )
                         
      enddo
    enddo
    !--------------------------------------------------------------------------
    !---------------------------- under-relaxation -----------------------------
    ! pij = (1-alfa)* pij(old) + alfa* Ap^-1 [Sp + Te pij] 
    !
    ! (Ap/alfa) pij = Sp + (1-alfa) * (Ap/alfa) pij(old)  + Te pij
    !  Ap' pij = Sp' + Te pij
    !--------------------------------------------------------------------------
    !Note: Rescaling Ap by 1/alfa 
    !   Alfa = 0.85_dp  
    !   Ap = Ap / Alfa
    !   Sp(:,:,3) = Sp(:,:,3) + (1.0_dp - Alfa )* Ap(:,:)*F(:,:,3)
    !---------------------------------------------------------------------------

    write(*,*) 'solve PPE'
    Res_sum_After=1.0_dp
    niter = 0 

    !Solve by a first iteration
    do while ( abs(Res_sum_After) > 1.d-7 .and. niter < 500)
       call TDMA(3)
       call Convergence_Criteria(3,Res_sum_After)
       niter= niter + 1
    end do      

    write(*,*) '|err|=',abs(Res_sum_After),niter

    !--------------------------------------------------------------------
    !--------------------------------------------------------------------
    !  velocities and pressure correction
    beta = 0.1_dp

    Summ = 0.0_dp
    Do i=2,NXmax  
      Do j=2,NYmax
        
        Dy = Y(i,j)-Y(i,j-1)
        Dx = X(i,j)-X(i-1,j)

        PPe = 0.5_dp * ( F(i,j,3) + F(i+1,j,3) ) 
        PPw = 0.5_dp * ( F(i,j,3) + F(i-1,j,3) ) 
        PPn = 0.5_dp * ( F(i,j,3) + F(i,j+1,3) ) 
        PPs = 0.5_dp * ( F(i,j,3) + F(i,j-1,3) ) 

        F(i,j,1) = F(i,j,1) + (PPw - PPe)  * DpU(i,j)  * Dy 
        F(i,j,2) = F(i,j,2) + (PPs - PPn)  * DpV(i,j)  * Dx 
        F(i,j,4) = F(i,j,4) + beta *  F(i,j,3) 
        Summ = Summ + F(i,j,4) 

      enddo
    enddo

    write(*,*)
    write(*,*) 'p=',summ

    !----------------------------------------------------------------------
    ! Fix pressure on the boundary nodes
    F(:,1,4) = F(:,2,4)              
    F(:,NYmaxC,4) = F(:,NYmaxC-1,4)   
    F(1,:,4) = F(2,:,4)             
    F(NXmaxC,:,4) = F(NXmaxC-1,:,4) 

    !fix at least 1 value of pressure  why!???
    F(:,:,4) = F(:,:,4) - F(3,4,4) 


  end subroutine solve_pressure_correction  

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
