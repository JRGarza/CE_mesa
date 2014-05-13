! ***********************************************************************
!
!   Copyright (C) 2012  Bill Paxton and Pablo Marchant
!
!   MESA is free software; you can use it and/or modify
!   it under the combined terms and restrictions of the MESA MANIFESTO
!   and the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,
!   or (at your option) any later version.
!
!   You should have received a copy of the MESA MANIFESTO along with
!   this software; if not, it is available at the mesa website:
!   http://mesa.sourceforge.net/
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************
 
      module mod_other_binary_jdot
          use binary_def, only : binary_info, binary_ptr

      ! NOTE: remember to set one of:
      ! use_other_jdot_mb = .true.
      ! use_other_jdot_gr = .true.
      ! use_other_jdot_ml = .true.
      ! use_other_jdot_tide = .true.
      ! use_other_extra_jdot = .true.
 
      
! you can add your own routine for use instead of the default ones

! here's how to do it.

! Before doing anything, let's make sure your working copy of run_binary_extras works.
! edit the extras_binary_controls routine 
!      subroutine extras_binary_controls(ierr)
!         use binary_data
!         integer, intent(out) :: ierr
!         ierr = 0
! 	       write(*,*) 'hello from extra_binary_controls'
!      end subroutine extras_controls

! then, in your work directory, do ./mk and ./rn to check that it is okay.
! assuming that worked, now edit extra_binary_controls to set the procedure pointer to other_jdot

!      subroutine extras_binary_controls(ierr)
!         integer, intent(out) :: ierr
!         ierr = 0
! 	       bhooks% other_jdot_mb => jdot_mb_routine
!      end subroutine extras_controls

!      subroutine jdot_mb_routine(other_jdot_mb, ierr)
!         use const_def, only: dp 
!         real(dp), intent(out) :: other_jdot_mb
!         integer, intent(out) :: ierr
!         ierr = 0
!         ! here is an (unrealistic) example
!         other_jdot_mb = 1
!      end subroutine jdot_mb_routine
         
      ! NOTE: if you'd like to have some inlist controls for your routine,
      ! you can use the x_ctrl array of real(dp) variables that is in &controls
      ! e.g., in the &controls inlist, you can set
      !     x_ctrl(1) = my_special_param
      ! then in your routine, you can access that by
      !     s% x_ctrl(1)
      ! of course before you can use s, you need to get it using the id argument.
      ! here's an example of how to do that -- add these lines at the start of your routine:
      !         use star_lib, only: star_ptr
      !         type (star_info), pointer :: s
      !         call star_ptr(id, s, ierr)
      !         if (ierr /= 0) then ! OOPS
      !            return
      !         end if
      ! 
      ! for integer control values, you can use x_integer_ctrl
      ! for logical control values, you can use x_logical_ctrl


      implicit none
      
            
      contains
      
      subroutine null_other_jdot_mb(ierr)
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         call binary_ptr(b,ierr)
         ierr = 0
         b% jdot_mb = 0 
      end subroutine null_other_jdot_mb

      subroutine null_other_jdot_gr(ierr)
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         call binary_ptr(b,ierr)
         ierr = 0
         b% jdot_gr = 0 
      end subroutine null_other_jdot_gr

      subroutine null_other_jdot_ml(ierr)
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         call binary_ptr(b,ierr)
         ierr = 0
         b% jdot_ml = 0 
      end subroutine null_other_jdot_ml

      subroutine null_other_jdot_tide(ierr)
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         call binary_ptr(b,ierr)
         ierr = 0
         b% jdot_tide = 0 
      end subroutine null_other_jdot_tide

      subroutine null_other_extra_jdot(ierr)
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         call binary_ptr(b,ierr)
         ierr = 0
         b% extra_jdot = 0 
      end subroutine null_other_extra_jdot

      subroutine null_other_jdot_ls(ierr)
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         call binary_ptr(b,ierr)
         ierr = 0
         b% jdot_ls = 0 
      end subroutine null_other_jdot_ls

      end module mod_other_binary_jdot

