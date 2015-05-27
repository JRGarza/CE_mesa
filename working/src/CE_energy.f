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
 
      module CE_energy


         
         
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
      
      
      subroutine CE_inject_energy(id, ierr)
         use const_def, only: Rsun
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: CE_energy_rate, CE_companion_position, CE_companion_radius, CE_companion_mass
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         s% extra_heat(:) = 0.0d0
         !#CE: Reading values of parameters from the extra controls that we are using
         !#CE: Note that "extra_heat" is the specific energy added to the the  cell in units of erg/s/gr
         CE_energy_rate = s% x_ctrl(1)
         CE_companion_position = s% x_ctrl(2)
         CE_companion_radius = s% x_ctrl(3)
         CE_companion_mass = s% x_ctrl(4)
         CE_test_case = s% x_integer_ctrl(1)

         if (CE_test_case == 1) then
            do k = 1, s% nz
               if (s% r(k) > s% he_core_radius * Rsun) then
                  s% extra_heat(k) = CE_energy_rate / ((s% star_mass - s% he_core_mass) * Msun)
               end if
            end do
         else if (CE_test_case == 2) then
            do k = 1, s% nz
               if ((s% m(k) > s% he_core_mass * Msun) .and. (s% m(k) < (s% he_core_mass + 0.1d0) * Msun)) then
                  s% extra_heat(k) = CE_energy_rate / (0.1 * Msun)
               end if
            end do

         else
            return
         endif
      end subroutine CE_inject_energy


      end module CE_energy
      
      
      
      
