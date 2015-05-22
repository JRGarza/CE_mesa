! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
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
 
      module other_energy

      ! NOTE: remember to set use_other_energy = .true. to enable this.
      
      
! you can add your own routine for use instead of the default_other_energy

! here's how to do it.

! in your working copy of run_star_extras, replace 
!      include 'standard_run_star_extras.inc'
! by the contents of the included file from star/job or mesa/include

! Don't make any edits to any of the files in star/job or mesa/include. 
! You do all of this in your private copy of run_star_extras.

! before doing anything else, let's make sure your working copy of run_star_extras works.
! edit the extras_controls routine 
!      subroutine extras_controls(s, ierr)
!         type (star_info), pointer :: s
!         integer, intent(out) :: ierr
!         ierr = 0
! 	       write(*,*) 'hello from extra_controls'
!      end subroutine extras_controls

! then, in your work directory, do ./mk and ./rn to check that it is okay.
! assuming that worked, now edit extra_controls to set the procedure pointer to other_energy

!      subroutine extras_controls(s, ierr)
!         type (star_info), pointer :: s
!         integer, intent(out) :: ierr
!         ierr = 0
! 	       s% other_energy => energy_routine
!      end subroutine extras_controls

!      subroutine energy_routine(id, ierr)
!         use const_def, only: Rsun
!         integer, intent(in) :: id
!         integer, intent(out) :: ierr
!         type (star_info), pointer :: s
!         integer :: k
!         ierr = 0
!         call star_ptr(id, s, ierr)
!         if (ierr /= 0) return
!         s% extra_heat(:) = 1 ! erg/g/sec
!         return
!         ! here is an example of calculating extra_heat for each cell.
!         do k = 1, s% nz
!            if (s% r(k) > 0.7*Rsun .and. s% r(k) < 0.9*Rsun) then
!               s% extra_heat(k) = 1d3*exp(-10*(s% r(k) - 0.8*Rsun)**2)
!            end if
!         end do
!      end subroutine energy_routine



! change your profile_columns.list to include
! 	extra_heat

! then do ./mk and ./rn
! check a profile plot of extra_heat to make sure it is nonzero.

! now you are ready to change your energy_routine to do what you want.

         
         
      ! NOTE: if you'd like to have some inlist controls for your routine,
      ! you can use the x_ctrl array of real(dp) variables that is in &controls
      ! e.g., in the &controls inlist, you can set
      !     x_ctrl(1) = <my_special_param>
      ! where <my_special_param> is a real value such as 0d0 or 3.59d0
      ! Then in your routine, you can access that by
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

      use star_def

      implicit none
      
            
      contains
      
      
      subroutine default_other_energy(id, ierr)
         use const_def, only: Rsun
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         s% extra_heat(:) = s% extra_power_source
         return
         ! here is an example of calculating extra_heat for each cell.
         do k = 1, s% nz
            if (s% r(k) > 0.7*Rsun .and. s% r(k) < 0.9*Rsun) then
               s% extra_heat(k) = 1d3*exp(-10*(s% r(k) - 0.8*Rsun)**2)
            end if
         end do
      end subroutine default_other_energy


      end module other_energy
      
      
      
      
