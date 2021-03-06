#include "../defines.inc"
!-----------------------------------------------------------!
!
!         Calculate Energy dou to exteranl potential
!
!            Started by Quinn, Dec 2017
!     Sets: dEField and dx_Field based on R, RP, and, HA
!-----------------------------------------------------------

subroutine MC_external_field(wlc_p,wlc_d,IT1,IT2)
use params, only: dp,wlcsim_params,wlcsim_data
implicit none
type(wlcsim_params), intent(in) :: wlc_p
TYPE(wlcsim_data), intent(inout) :: wlc_d
integer, intent(in) :: IT1    ! Start test bead
integer, intent(in) :: IT2    ! Final test bead

integer ii
real(dp) vv(3)
real(dp), parameter :: center(3) = [WLC_P__LBOX_X/2.0_dp,&
                                    WLC_P__LBOX_Y/2.0_dp,&
                                    WLC_P__LBOX_Z/2.0_dp]
wlc_d%dx_Field = 0.0_dp
do ii = IT1, IT2
    vv = wlc_d%RP(:,ii)-center
    if (dot_product(vv,vv) < &
        (WLC_P__BINDING_R + WLC_P__CONFINEMENT_SPHERE_DIAMETER/2.0)**2) then
        wlc_d%dx_Field = wlc_d%dx_Field + 1.0_dp
    endif
    vv = wlc_d%R(:,ii)-center
    if (dot_product(vv,vv) < &
        (WLC_P__BINDING_R + WLC_P__CONFINEMENT_SPHERE_DIAMETER/2.0)**2) then
        wlc_d%dx_Field = wlc_d%dx_Field - 1.0_dp
    endif
enddo
wlc_d%DEField = wlc_p%HA * wlc_d%dx_Field

END
