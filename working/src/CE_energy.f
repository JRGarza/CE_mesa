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
      use CE_orbit, only: AtoP, TukeyWindow, calc_quantities_at_comp_position


      implicit none


      contains


      subroutine CE_inject_energy(id, ierr)
         use const_def, only: Rsun
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: CE_test_case, k
         real(dp) :: CE_companion_position

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         CE_test_case = s% x_integer_ctrl(1)
         CE_companion_position = s% xtra2

         ! If the star is in the initial relaxation phase, skip energy calculations
         if (s% doing_relax) return
         ! If companion is outside star, skip energy calculations
         if (CE_companion_position*Rsun > s% r(1)) return





         ! Call functions to calculate test cases
         if (CE_test_case == 1) then

            call CE_inject_case1(id, ierr)

         else if (CE_test_case == 2) then

            call CE_inject_case2(id, ierr)

         else if (CE_test_case == 3) then

            call CE_inject_case3(id, ierr)

          else if (CE_test_case == 4) then

             call CE_inject_case4(id, ierr)
         endif


      end subroutine CE_inject_energy




      subroutine CE_inject_case1(id, ierr)

         use const_def, only: Msun
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         real(dp) :: ff, mass_to_be_heated, m_bot
         real(dp) :: CE_energy_rate

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         CE_energy_rate = s% x_ctrl(1)

         ! mass (g) of the bottom of the (outer) convective envelope
           ! Based on the inner edge of the convective envelope
!          m_bot = s% conv_mx1_bot * s% mstar
           ! Based on the helium core
         m_bot = s% he_core_mass * Msun

         !First calculate the mass in which the energy will be deposited
         mass_to_be_heated = 0.0
         do k = 1, s% nz
            ff = EnvelopeWindow(s% m(k), m_bot)
            mass_to_be_heated = mass_to_be_heated + s% dm(k) * ff
         end do

         !Now redo the loop and add the extra specific heat
         do k = 1, s% nz
            s% extra_heat(k) = CE_energy_rate / mass_to_be_heated * EnvelopeWindow(s% m(k), m_bot)
         end do


         ! Save the total erg/second added in this time step
         s% xtra1 = CE_energy_rate



         contains

         real(dp) function EnvelopeWindow(m_interior, m_bot)
            use const_def, only: pi, Msun
            real(dp), intent(in) :: m_interior, m_bot

            ! arctan function causes a smooth transition from 0 to unity around m_bot.
            ! The 0.002 in the denominator sets the width of the transition, in this case
            ! calibrated so the transition region is roughly 0.01 Msun.
            EnvelopeWindow = 1./pi * atan((s% m(k) - m_bot) / (0.002*Msun)) + 0.5

         end function EnvelopeWindow

      end subroutine CE_inject_case1



      subroutine CE_inject_case2(id, ierr)

         use const_def, only: Msun
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         real(dp) :: ff, mass_to_be_heated, a_tukey
         real(dp) :: CE_energy_rate, CE_companion_position, CE_companion_radius, CE_companion_mass
         real(dp) :: CE_n_acc_radii

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! Get input controls
         CE_energy_rate = s% x_ctrl(1)
         CE_companion_position = s% xtra2
         CE_companion_radius = s% xtra3
         CE_companion_mass = s% xtra4
         CE_n_acc_radii = s% xtra5

         ! Tukey window scale
         a_tukey = 0.1
         !TODO  Change CE_companion_radius to R_acc
         ! First calculate the mass in which the energy will be deposited
         mass_to_be_heated = 0.0
         do k = 1, s% nz
            ff = TukeyWindow((s% r(k) - CE_companion_position*Rsun)/(CE_n_acc_radii * CE_companion_radius*Rsun), a_tukey)
            mass_to_be_heated = mass_to_be_heated + s% dm(k) * ff
         end do

         ! Now redo the loop and add the extra specific heat
         do k = 1, s% nz
            ff = TukeyWindow((s% r(k) - CE_companion_position*Rsun)/(CE_n_acc_radii * CE_companion_radius*Rsun), a_tukey)
            s% extra_heat(k) = CE_energy_rate / mass_to_be_heated * ff
         end do

         ! Save the total erg/second added in this time step
         s% xtra1 = CE_energy_rate

      end subroutine CE_inject_case2

      subroutine CE_inject_case3(id, ierr)

         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         real(dp) :: CE_energy_rate, CE_companion_position, CE_companion_radius, CE_companion_mass
         real(dp) :: CE_n_acc_radii
         real(dp) :: time, M2
         real(dp) :: I, F_drag, F_coef
         real(dp) :: a_tukey, mass_to_be_heated, ff
         real(dp) :: R_acc, R_acc_low, R_acc_high
         real(dp) :: v_rel, v_rel_div_csound, M_encl, rho_at_companion, scale_height_at_companion


         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! Alternative energy source here


         ! Get input controls
         CE_energy_rate = s% xtra1
         CE_companion_position = s% xtra2
         CE_companion_radius = s% xtra3
         CE_companion_mass = s% xtra4
         CE_n_acc_radii = s% xtra5


         call calc_quantities_at_comp_position(id, ierr)

         R_acc = s% xtra12
         R_acc_low = s% xtra13
         R_acc_high = s% xtra14
         M_encl = s% xtra15
         v_rel = s% xtra16
         v_rel_div_csound = s% xtra17
         rho_at_companion = s% xtra18
         scale_height_at_companion = s% xtra19

!         ! This is incorrect, but for now, not completely crazy
!         R_acc = (R_acc_low + R_acc_high) / 2.0


         ! Determine drag force
         ! Equations from Ostriker (1999) ApJ, 513, 252
         time = 1.0 ! FIX THIS: What is time here?
         if (v_rel_div_csound .lt. 1.0) then
            I = 0.5 * log((1.0+v_rel_div_csound)/(1.0-v_rel_div_csound)) - v_rel_div_csound
         else
            I = 0.5 * log((v_rel_div_csound+1.0)/(v_rel_div_csound-1.0)) + log(v_rel*time / R_acc)
         end if

         M2 = CE_companion_mass * Msun

         F_coef = 4.0 * pi * standard_cgrav * standard_cgrav * M2 * M2 * rho_at_companion / (v_rel*v_rel)
         F_drag = -F_coef * I

         ! Total energy rate= drag force * velocity
         CE_energy_rate = F_drag * v_rel


         ! Tukey window scale
         a_tukey = 0.1

         ! First calculate the mass in which the energy will be deposited
         mass_to_be_heated = 0.0
         do k = 1, s% nz
            ff = TukeyWindow((s% r(k) - CE_companion_position*Rsun)/(CE_n_acc_radii * R_acc), a_tukey)
            mass_to_be_heated = mass_to_be_heated + s% dm(k) * ff
         end do

         ! Now redo the loop and add the extra specific heat
         do k = 1, s% nz
            ff = TukeyWindow((s% r(k) - CE_companion_position*Rsun)/(CE_n_acc_radii * R_acc), a_tukey)
            s% extra_heat(k) = CE_energy_rate / mass_to_be_heated * ff
         end do


         ! Save the total erg/second added in this time step
         s% xtra1 = CE_energy_rate

      end subroutine CE_inject_case3



      subroutine CE_inject_case4(id, ierr)

         use const_def, only: Rsun, Msun, pi, standard_cgrav
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k, k_bottom
         real(dp) :: CE_energy_rate, CE_companion_position, CE_companion_radius, CE_companion_mass
         real(dp) :: CE_n_acc_radii
         real(dp) :: M2
         real(dp) :: I, F_drag
         real(dp) :: a_tukey, mass_to_be_heated, ff
         real(dp) :: F_DHL, f1, f2, f3, e_rho
         real(dp) :: R_acc, R_acc_low, R_acc_high
         real(dp) :: v_rel, v_rel_div_csound, M_encl, rho_at_companion, scale_height_at_companion
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! Alternative energy source here


         ! Get input controls
         CE_energy_rate = s% xtra1
         CE_companion_position = s% xtra2
         CE_companion_radius = s% xtra3
         CE_companion_mass = s% xtra4
         CE_n_acc_radii = s% xtra5


         call calc_quantities_at_comp_position(id, ierr)

         R_acc = s% xtra12
         R_acc_low = s% xtra13
         R_acc_high = s% xtra14
         M_encl = s% xtra15
         v_rel = s% xtra16
         v_rel_div_csound = s% xtra17
         rho_at_companion = s% xtra18
         scale_height_at_companion = s% xtra19


!         ! For a first approximation, let's use the average R_acc
!         R_acc = (R_acc_low + R_acc_high) / 2.0
         F_DHL = pi * R_acc**2 * rho_at_companion * v_rel**2


         f1 = 1.91791946d0
         f2 = -1.52814698d0
         f3 = 0.75992092
         e_rho = R_acc / scale_height_at_companion
         F_drag = F_DHL*(f1 + f2*e_rho +f3*e_rho**2)




         ! Total energy rate= drag force * velocity
         CE_energy_rate = F_drag * v_rel


         ! Tukey window scale
         a_tukey = 0.5

         ! First calculate the mass in which the energy will be deposited
         mass_to_be_heated = 0.0
         do k = 1, s% nz
            if (s% r(k) < CE_companion_position*Rsun) then
               R_acc = R_acc_low
            else
               R_acc = R_acc_high
            end if

            ff = TukeyWindow((s% r(k) - CE_companion_position*Rsun)/(CE_n_acc_radii * 2.0 * R_acc), a_tukey)
            mass_to_be_heated = mass_to_be_heated + s% dm(k) * ff
            !Energy should be deposited only on the envelope of the star and not in the core
            !When we reach the core boundary we exit the loop
            if (s% m(k) < s% he_core_mass * Msun) exit
         end do
         !this is the limit in k of the boundary between core and envelope
         k_bottom = k-1
         ! If companion is outside star, set mass_to_be_heated arbitrarily low
         if (mass_to_be_heated == 0.) mass_to_be_heated = 1.0

         ! Now redo the loop and add the extra specific heat
         do k = 1, k_bottom
            if (s% r(k) < CE_companion_position*Rsun) then
               R_acc = R_acc_low
            else
               R_acc = R_acc_high
            end if

            ff = TukeyWindow((s% r(k) - CE_companion_position*Rsun)/(CE_n_acc_radii * 2.0 * R_acc), a_tukey)
            s% extra_heat(k) = CE_energy_rate / mass_to_be_heated * ff
         end do


         ! Save the total erg/second added in this time step
         s% xtra1 = CE_energy_rate

      end subroutine CE_inject_case4







      end module CE_energy
