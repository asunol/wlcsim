#include "../defines.inc"
!---------------------------------------------------------------!

!
!     This subroutine calculates the change in the self energy for
!     a small Monte Carlo move in the position.
!
!     Andrew Spakowitz
!     Written 6-29-04
!
!     Edited by Shifan
!
!     Edited by Quinn in 2016

subroutine MC_int_rep(wlc_p,wlc_d,I1,I2,forward)
use params
implicit none

!   iputs
TYPE(wlcsim_params), intent(in) :: wlc_p   ! <---- Contains output
TYPE(wlcsim_data), intent(inout) :: wlc_d
integer, intent(in) :: I1  ! Test bead position 1
integer, intent(in) :: I2  ! Test bead position 2

!   Internal variables
integer I                 ! For looping over bins
integer II                ! For looping over IB
integer IB                ! Bead index
integer rrdr ! -1 if r, 1 if r + dr
integer IX(2),IY(2),IZ(2)
real(dp) WX(2),WY(2),WZ(2)
real(dp) WTOT       ! total weight ascribed to bin
real(dp) RBin(3)    ! bead position
integer inDBin              ! index of bin
integer ISX,ISY,ISZ
LOGICAL isA   ! The bead is of type A
real(dp), dimension(-2:2) :: phi2
integer m_index  ! m from Ylm spherical harmonics
integer NBinX(3)
real(dp) temp    !for speeding up code
LOGICAL, intent(in) :: forward ! move forward
integer AminusB
NBinX = wlc_p%NBINX

wlc_d%NPHI = 0
! -------------------------------------------------------------
!
!  Calculate end beads
!
!--------------------------------------------------------------

do II = 1,2
  if (II.eq.1) then
      IB = I1
      if (forward) then
          rrdr = -1
      else
          rrdr = 1
      endif
  elseif (II.eq.2) then
      IB = I2
      if (forward) then
          rrdr = 1
      else
          rrdr = -1
      endif
  else
      print*, "Error in MC_int_rep, II = {1,2}"
      stop 1
  endif
   ! subract current and add new
   if (rrdr.eq.-1) then
       RBin(1) = wlc_d%R(1,IB)
       RBin(2) = wlc_d%R(2,IB)
       RBin(3) = wlc_d%R(3,IB)
   else
       RBin(1) = wlc_d%RP(1,IB)
       RBin(2) = wlc_d%RP(2,IB)
       RBin(3) = wlc_d%RP(3,IB)
   endif
   isA = wlc_d%AB(IB).eq.1
   if (WLC_P__TWO_TAIL) then
       print*, "The Reptation move is not currently set up for two Tail"
       stop 1
   endif
   if (wlc_p%CHI_L2_ON .and. isA) then
       if (rrdr == -1) then
           call Y2calc(wlc_d%U(:,IB),phi2)
       else
           call Y2calc(wlc_d%UP(:,IB),phi2)
       endif
   else
       ! You could give some MS parameter to B as well if you wanted
       phi2=0.0
   endif
   ! --------------------------------------------------
   !
   !  Interpolate beads into bins
   !
   ! --------------------------------------------------
   call interp(wlc_p,RBin,IX,IY,IZ,WX,WY,WZ)

   ! -------------------------------------------------------
   !
   ! Count beads in bins
   !
   ! ------------------------------------------------------
   !   Add or Subtract volume fraction with weighting from each bin
   !   I know that it looks bad to have this section of code twice but it
   !   makes it faster.
   if (isA) then
       do ISX = 1,2
          do ISY = 1,2
             do ISZ = 1,2
                WTOT = WX(ISX)*WY(ISY)*WZ(ISZ)
                inDBin = IX(ISX) + (IY(ISY)-1)*NBinX(1) + (IZ(ISZ)-1)*NBinX(1)*NBinX(2)
                ! Generate list of which phi's change and by how much
                I = wlc_d%NPHI
                do
                   if (I.eq.0) then
                      wlc_d%NPHI = wlc_d%NPHI + 1
                      wlc_d%inDPHI(wlc_d%NPHI) = inDBin
                      temp = rrdr*WTOT*WLC_P__BEADVOLUME/(WLC_P__DBIN**3)
                      wlc_d%DPHIA(wlc_d%NPHI) = temp
                      wlc_d%DPHIB(wlc_d%NPHI) = 0.0_dp
                      if(wlc_p%CHI_L2_ON) then
                          do m_index = -2,2
                              wlc_d%DPHI_l2(m_index,wlc_d%NPHI) = &
                                  + phi2(m_index)*temp
                          enddo
                      endif
                      exit
                   elseif (inDBin == wlc_d%inDPHI(I)) then
                      temp = rrdr*WTOT*WLC_P__BEADVOLUME/(WLC_P__DBIN**3)
                      wlc_d%DPHIA(I) = wlc_d%DPHIA(I) + temp
                      if(wlc_p%CHI_L2_ON) then
                          do m_index = -2,2
                              wlc_d%DPHI_l2(m_index,I) = wlc_d%DPHI_l2(m_index,I) &
                                  + phi2(m_index)*temp
                          enddo
                      endif
                      exit
                   else
                      I = I-1
                   endif
                enddo
             enddo
          enddo
       enddo
   else
       do ISX = 1,2
          do ISY = 1,2
             do ISZ = 1,2
                WTOT = WX(ISX)*WY(ISY)*WZ(ISZ)
                inDBin = IX(ISX) + (IY(ISY)-1)*NBinX(1) + (IZ(ISZ)-1)*NBinX(1)*NBinX(2)
                ! Generate list of which phi's change and by how much
                I = wlc_d%NPHI
                do
                   if (I.eq.0) then
                      wlc_d%NPHI = wlc_d%NPHI + 1
                      wlc_d%inDPHI(wlc_d%NPHI) = inDBin
                      wlc_d%DPHIA(wlc_d%NPHI) = 0.0_dp
                      wlc_d%DPHIB(wlc_d%NPHI) = rrdr*WTOT*WLC_P__BEADVOLUME/(WLC_P__DBIN**3)
                      if(wlc_p%CHI_L2_ON) then
                          do m_index = -2,2
                              ! This is somewhat wastefull, could eliminate for speedup by having another NPHI for L=2
                              wlc_d%DPHI_l2(m_index,wlc_d%NPHI) = 0.0
                          enddo
                      endif
                      exit
                   elseif (inDBin == wlc_d%inDPHI(I)) then
                      wlc_d%DPHIB(I) = wlc_d%DPHIB(I) + rrdr*WTOT*WLC_P__BEADVOLUME/(WLC_P__DBIN**3)
                      exit
                   else
                      I = I-1
                   endif
                enddo
             enddo !ISZ
          enddo !ISY
       enddo !ISX
   endif
enddo ! loop over IB  A.k.a. beads
! ---------------------------------------------------------------------
!



!----------------------------------------------------------
!
!  Now do intermediate Beads
!
!-----------------------------------------------------------
do IB = I1,I2-1
   if (wlc_d%AB(IB).eq.wlc_d%AB(IB + 1)) CYCLE
   if (forward) then
       RBin(1) = wlc_d%RP(1,IB)
       RBin(2) = wlc_d%RP(2,IB)
       RBin(3) = wlc_d%RP(3,IB)
       isA = wlc_d%AB(IB).eq.1
   else
       RBin(1) = wlc_d%R(1,IB)
       RBin(2) = wlc_d%R(2,IB)
       RBin(3) = wlc_d%R(3,IB)
       isA = wlc_d%AB(IB + 1).eq.1
   endif
   if(isA) then
       AminusB=1
   else
       AminusB=-1
   endif
   ! --------------------------------------------------
   !
   !  Interpolate beads into bins
   !
   ! --------------------------------------------------
   call interp(wlc_p,RBin,IX,IY,IZ,WX,WY,WZ)

   if (wlc_p%CHI_L2_ON) then
       if (forward) then
           call Y2calc(wlc_d%UP(:,IB),phi2)
       else
           call Y2calc(wlc_d%U(:,IB),phi2)
       endif
   else
       ! You could give some MS parameter to B as well if you wanted
       phi2=0.0
   endif
   ! -------------------------------------------------------
   !
   ! Count beads in bins
   !
   ! ------------------------------------------------------
   !   Add or Subtract volume fraction with weighting from each bin
   do ISX = 1,2
      do ISY = 1,2
         do ISZ = 1,2
            WTOT = WX(ISX)*WY(ISY)*WZ(ISZ)
            inDBin = IX(ISX) + (IY(ISY)-1)*NBinX(1) + (IZ(ISZ)-1)*NBinX(1)*NBinX(2)
            ! Generate list of which phi's change and by how much
            I = wlc_d%NPHI
            do
               if (I.eq.0) then
                  wlc_d%NPHI = wlc_d%NPHI + 1
                  wlc_d%inDPHI(wlc_d%NPHI) = inDBin
                  temp = AminusB*WTOT*WLC_P__BEADVOLUME/(WLC_P__DBIN**3)
                  wlc_d%DPHIA(wlc_d%NPHI) = temp
                  wlc_d%DPHIB(wlc_d%NPHI) = -temp
                  if(wlc_p%CHI_L2_ON) then
                      do m_index = -2,2
                          wlc_d%DPHI_l2(m_index,wlc_d%NPHI) = &
                              + phi2(m_index)*temp
                      enddo
                  endif
                  exit
               elseif (inDBin == wlc_d%inDPHI(I)) then
                  temp = AminusB*WTOT*WLC_P__BEADVOLUME/(WLC_P__DBIN**3)
                  wlc_d%DPHIA(I) = wlc_d%DPHIA(I) + temp
                  wlc_d%DPHIB(I) = wlc_d%DPHIB(I)-temp
                  if(wlc_p%CHI_L2_ON) then
                      do m_index = -2,2
                          wlc_d%DPHI_l2(m_index,I) = wlc_d%DPHI_l2(m_index,I) &
                              + phi2(m_index)*temp
                      enddo
                  endif
                  exit
               else
                  I = I-1
               endif
            enddo
         enddo
      enddo
   enddo
enddo ! loop over IB  A.k.a. beads
! ---------------------------------------------------------------------
!
! Calcualte change in energy
!
!---------------------------------------------------------------------
call hamiltonian(wlc_p,wlc_d,.false.)

RETURN
END

!---------------------------------------------------------------!
