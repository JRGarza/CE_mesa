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
 
      module CE_torque


         
         
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
      use const_def

      implicit none


      
            
      contains
      

      ! Angular momentum prescription from CE
      ! Based off example from mesa/star/other/other_torque.f
              
      subroutine CE_inject_am(id, ierr)
      
         use const_def, only: Rsun
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         real(dp) :: a_tukey, mass_to_be_spun, ff
         real(dp) :: CE_companion_position, CE_companion_radius, CE_n_acc_radii, CE_ang_mom_transferred
         
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         
         ! Load angular momentum dissipated in the envelope from the orbit decrease
         CE_companion_position = s% xtra2
         CE_companion_radius = s% xtra3
         CE_n_acc_radii = s% xtra5
         CE_ang_mom_transferred = s% xtra6


         
         
         ! First calculate the mass in which the angular momentum will be deposited
         mass_to_be_spun = 0.0
         do k = 1, s% nz
            ff = TukeyWindow((s% r(k) - CE_companion_position*Rsun)/(CE_n_acc_radii * CE_companion_radius*Rsun), a_tukey)
            mass_to_be_spun = mass_to_be_spun + s% dm(k) * ff
         end do

         ! Now redo the loop and add the extra specific heat
         do k = 1, s% nz
            ff = TukeyWindow((s% r(k) - CE_companion_position*Rsun)/(CE_n_acc_radii * CE_companion_radius*Rsun), a_tukey)
            s% extra_jdot(k) = CE_ang_mom_transferred / s% dt / mass_to_be_spun * ff
         end do


         ! Need to add some controls here to make sure certain layers of the star are not 
         ! getting spun up past break up velocity


         
!         ! here is an (unrealistic) example of adding angular momentum to each cell
!         do k = 1, s% nz
!            ! note that can set extra_omegadot instead of extra_jdot if that is more convenient
!            ! set one or the other, not both.  set the one you are not using to 0 as in the following line.
!            !s% extra_jdot(k) = 0 ! rate at which specific angular momentum is changed
!            !s% extra_omegadot(k) = 0 ! rate at which omega is changed
!         end do
         
         
         contains
         
         real(dp) function TukeyWindow(x,a)
            use const_def, only: dp, pi
            real(dp), intent(in) :: x, a

            if ((x .ge. -0.5) .and. (x .le. 0.5) .and. (2.*x+a .ge. 0) .and. (-2.*x+a .ge. 0)) then
               TukeyWindow = 1.
            else if ((x .ge. -0.5) .and. (x .le. 0.5) .and. (2.*x+a .lt. 0)) then
               TukeyWindow = 0.5*(1.-sin(pi*x/a))
            else if ((x .ge. -0.5) .and. (x .le. 0.5) .and. (2.*x+a .gt. 0) .and. (-2.*x+a .lt. 0)) then
               TukeyWindow = 0.5*(1.+sin(pi*x/a))
            else
               TukeyWindow = 0.
            endif
            
         end function TukeyWindow
         
         
      end subroutine CE_inject_am


      end module CE_torque



