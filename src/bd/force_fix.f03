module fixation
  use params, only : dp
  implicit none
  integer, allocatable :: fixed_pairs(:,:)
  logical, save :: is_allocated = .False.
contains

subroutine force_fix(FFIX,TFIX,R,U,NT,wlc_p)
  use params, only : wlcsim_params, dp

  implicit none

  integer, intent(in) :: NT
  real(dp), intent(in) :: R(3,NT), U(3,NT)
  type(wlcsim_params), intent(in) :: wlc_p
  real(dp), intent(out) :: FFIX(3,NT), TFIX(3,NT)

  integer, save :: num_fixed = 0

  if (.not. is_allocated) then
     allocate(fixed_pairs(2,NT))
     is_allocated = .True.
  endif

  TFIX = 0 
  ! adds new springs to fixed_pairs if necessary
  call update_fixed_pairs(R, wlc_p%NT, wlc_p%fix_r, fixed_pairs, num_fixed)

  call calculate_force(FFIX, R, NT, fixed_pairs, num_fixed, wlc_p%EPAR)
end subroutine

subroutine update_fixed_pairs(R, NT, fix_r, fixed_pairs, num_fixed)
  use params, only: dp
  implicit none
  integer, intent(in) :: NT
  integer, intent(inout) :: num_fixed
  real(dp), intent(in) :: R(3,NT)
  integer, intent(inout) :: fixed_pairs(2,NT)
  real(dp), intent(in) :: fix_r

  integer I,J
  real(dp) dr

  do I = 1,NT
     do J = I,NT
        if (abs(I-J) .le. 1) then
           cycle
        end if
        dr = norm2(R(:,I)-R(:,J))
        if ((dr .le. fix_r) .and. .not.contains_pair(fixed_pairs,NT,num_fixed,I,J)) then
           num_fixed = num_fixed + 1
           if (num_fixed .GE. NT) then
              print *, "Number of fixation connections exceeds number of beads, stopping arbitrarily!"
              stop 1
           endif
           fixed_pairs(1,num_fixed) = I
           fixed_pairs(2,num_fixed) = J
           print *, "The beads ", I, " and ", J, " are now connected."
        end if
     end do
  end do
end subroutine

function contains_pair(fixed_pairs,NT,num_fixed,I,J) result(is_contained)
  implicit none
  integer, intent(in) :: NT, num_fixed, I, J
  integer, intent(in) :: fixed_pairs(2,NT)
  logical is_contained

  integer K

  is_contained = .false.

  do K = 1,num_fixed
     if (fixed_pairs(1,K).eq.I .and. fixed_pairs(2,K).eq.J) then
        is_contained = .true.
        exit
     endif
  enddo
  return
end function
  

subroutine calculate_force(FFIX, R, NT, fixed_pairs, num_fixed, EPAR)
  use params, only : dp
  
  implicit none
  
  real(dp), intent(out) :: FFIX(3,NT)
  real(dp), intent(in) :: R(3,NT)
  real(dp), intent(in) :: EPAR
  integer, intent(in) :: NT
  integer, intent(in) :: fixed_pairs(2,NT)
  integer, intent(in) :: num_fixed

  integer I,J,K
  real(dp) DR(3)
  
  FFIX = 0
  do K = 1,num_fixed
     I = fixed_pairs(1,K)
     J = fixed_pairs(2,K)
     DR = R(:,I) - R(:,J)
     FFIX(:,I) = FFIX(:,I) - EPAR*DR
     FFIX(:,J) = FFIX(:,J) + EPAR*DR
  end do
end subroutine
end module
