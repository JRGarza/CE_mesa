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

      module CE_orbit



      use star_def
      use const_def
      use CE_energy

      implicit none

      contains


      subroutine CE_orbit_adjust(id, ierr)
         use const_def, only: standard_cgrav, Msun, Rsun
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         real(dp) :: CE_energy_rate, CE_companion_position, CE_companion_mass
         real(dp) :: J_tmp, J_init, J_final
         real(dp) :: E_init, E_loss, E_final, E_tmp
         real(dp) :: M_inner, R_inner, M_outer, R_outer, M_final, R_final
         real(dp) :: M_slope, R_slope, M_int, R_int, M_encl
         real(dp) :: top, bottom, k_final
         real(dp) :: CE_ang_mom_transferred


         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return


         CE_energy_rate = s% xtra1
         CE_companion_position = s% xtra2
         CE_companion_mass = s% xtra4

         ! Mass for orbital energy calculation depends on mass contained within a radius
         ! Include all the enclosed cells
         ! Add to it the enclosed mass of the current cell
         k = 1
         do while (s% r(k) > CE_companion_position * Rsun)
            k = k + 1
         end do

         M_encl = s% m(k)
         M_encl = M_encl + s% dm(k-1) * (CE_companion_position*Rsun - s% r(k)) / (s% r(k-1) - s% r(k))


         ! Calculate the angular momentum
         J_tmp = (CE_companion_mass * Msun)**2 * M_encl**2 / (CE_companion_mass * Msun + M_encl)
         J_init = sqrt(standard_cgrav * J_tmp * CE_companion_position * Rsun)


         ! Calculate the energies
         E_init = -standard_cgrav * CE_companion_mass * Msun * M_encl / (2.0 * CE_companion_position * Rsun)
         E_loss = CE_energy_rate * s% dt
         E_final = E_init - E_loss


         ! Move from outside of star in to find cell containing companion
         E_tmp = 0d0
         k = 1
         do while (E_tmp > E_final)
            M_inner = s% m(k)
            R_inner = s% r(k)
            E_tmp = -standard_cgrav * CE_companion_mass * Msun * M_inner / (2.0 * R_inner)
!            write(*,*) "Cell number: ", k, " Radius: ", s% r(k)/Rsun, " Mass: ", s% m(k), " Energy: ", E_tmp
            k = k + 1
         end do

         ! save end points of cell containing companion
         M_inner = s% m(k-2)
         R_inner = s% r(k-2)
         M_outer = s% m(k-1)
         R_outer = s% r(k-1)

         ! We could choose to interpolate for R using M as the independent variable. Instead, here we
         ! linearly interpolate across cell (using k as the independent variable)
         M_slope = (M_outer - M_inner) / real((k-1) - (k-2))
         M_int = s% m(k-1) - M_slope * real(k-1)
         R_slope = (R_outer - R_inner) / real((k-1) - (k-2))
         R_int = s% r(k-1) - R_slope * real(k-1)

         ! Given the final energy, E_final, determine the resulting k that solves the equation
         top = 2.0 * E_final * R_int + standard_cgrav * CE_companion_mass * Msun * M_int
         bottom = 2.0 * E_final * R_slope + standard_cgrav * CE_companion_mass * Msun * M_slope
         k_final = -top / bottom

         ! Now use the interpolations and the derived k_final, determine the resulting separation and enclosed mass
         R_final = R_slope * k_final + R_int
         M_final = M_slope * k_final + M_int

         s% xtra2 = R_final/Rsun

         ! Calculate the angular momentum lost to the star's envelope
         J_tmp = (CE_companion_mass * Msun)**2 * M_final**2 / (CE_companion_mass * Msun + M_final)
         J_final = sqrt(standard_cgrav * J_tmp * R_final)

         !The angular momentum that is lost from the orbit of the companion
         ! is added to the envelope of the donor.
         CE_ang_mom_transferred = -(J_final - J_init)
         s% xtra6 = CE_ang_mom_transferred/s% dt


         ! For diagnostics

         write(*,*) "Final k: ", k_final
         write(*,*) "Previous Enclosed Mass: ", M_encl, " Final Enclosed Mass: ", M_final
         write(*,*) "Previous Separation = ", CE_companion_position, " Final Separation: ", R_final/Rsun
         write(*,*) "Previous Orbital Energy = ", E_init, " Final Orbital Energy: ", E_final
         write(*,*) "Previous Angular momentum = ", J_init, " Final Angular momentum: ", J_final
         write(*,*) "Dissipated Energy Rate: ", s% xtra1, " Dissipated Angular Momentum Rate: ", s% xtra6


      end subroutine CE_orbit_adjust

      end module CE_orbit
