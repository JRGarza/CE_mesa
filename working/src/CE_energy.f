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
         real(dp) :: CE_energy_rate, CE_companion_position, CE_companion_radius, CE_companion_mass, CE_test_case
         real(dp) :: mass_to_be_heated, a_tukey = 0.1, ff, ff_integral, extra_heat_integral, m_bot
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

         mass_to_be_heated = 0.
         ff = 0.

         if (CE_test_case == 1) then
         
            ! mass (g) of the bottom of the (outer) convective envelope         
              ! Based on the inner edge of the convective envelope
!            m_bot = s% conv_mx1_bot * s% mstar
              ! Based on the helium core
            m_bot = s% he_core_mass * Msun

            !First calculate the mass in which the energy will be deposited
            do k = 1, s% nz
               ff = EnvelopeWindow(s% m(k), m_bot)
               mass_to_be_heated = mass_to_be_heated + s% dm(k) * ff
            end do
            !Now redo the loop and add the extra specific heat
            do k = 1, s% nz
               s% extra_heat(k) = CE_energy_rate / mass_to_be_heated * EnvelopeWindow(s% m(k), m_bot)
            end do
            
         else if (CE_test_case == 2) then

            !First calculate the mass in which the energy will be deposited
            do k = 1, s% nz
               ff = TukeyWindow(s% r(k)/(CE_companion_radius*Rsun) - CE_companion_position, a_tukey)
               mass_to_be_heated =  mass_to_be_heated + s% dm(k) * ff
            end do
            !Now redo the loop and add the extra specific heat
            extra_heat_integral = 0.
            do k = 1, s% nz
               ff = TukeyWindow(s% r(k)/(CE_companion_radius*Rsun) - CE_companion_position, a_tukey)
               s% extra_heat(k) = CE_energy_rate / mass_to_be_heated * ff
            end do
            
         else
            return
         endif
         
         contains

         real(dp) function EnvelopeWindow(m_interior, m_bot)
            use const_def, only: pi, Msun
            real(dp), intent(in) :: m_interior, m_bot

            ! arctan function causes a smooth transition from 0 to unity around m_bot.
            ! The 0.002 in the denominator sets the width of the transition, in this case
            ! calibrated so the transition region is roughly 0.01 Msun.
            EnvelopeWindow = 1./pi * atan((s% m(k) - m_bot) / (0.002*Msun)) + 0.5
         
         end function EnvelopeWindow

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

      end subroutine CE_inject_energy


      end module CE_energy
      
      
      
      
