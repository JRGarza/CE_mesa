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
         real(dp) :: M_slope, R_slope, M_int, R_int
         real(dp) :: top, bottom, k_final
         real(dp) :: orbital_ang_mom_lost
         real(dp) :: R_acc, v_rel, v_rel_div_csound, M_encl, rho_at_companion, scale_height_at_companion


         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return


         CE_energy_rate = s% xtra1
         CE_companion_position = s% xtra2
         CE_companion_mass = s% xtra4

         ! If companion is outside star, skip energy calculations
         if (CE_companion_position*Rsun > s% r(1)) return


         call calc_quantities_at_comp_position(id, ierr)

         R_acc = s% xtra12
         M_encl = s% xtra13
         v_rel = s% xtra14
         v_rel_div_csound = s% xtra15
         rho_at_companion = s% xtra16
         scale_height_at_companion = s% xtra17


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
            k = k + 1
         end do


         ! If companion is outside star, set k to 3
         if (k < 3) k=3


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
         !Saving as s% xtra9 the enclosed mass so that we output it in the history data
         s% xtra9 = M_final/Msun

         ! Calculate the angular momentum lost to the star's envelope
         J_tmp = (CE_companion_mass * Msun)**2 * M_final**2 / (CE_companion_mass * Msun + M_final)
         J_final = sqrt(standard_cgrav * J_tmp * R_final)

         !The angular momentum that is lost from the orbit of the companion
         ! is added to the envelope of the donor.
         orbital_ang_mom_lost = J_final - J_init
         !We save in s% xtra6 the total torque that will be applied to the Envelope
         s% xtra6 = -orbital_ang_mom_lost/s% dt

         ! Keep track of orbital energy and angular momentum
         s% xtra8 = E_final
         s% xtra10 = J_final

         ! For diagnostics

         write(*,*) "Final k: ", k_final
         write(*,*) "Previous Enclosed Mass: ", M_encl/Msun, " Final Enclosed Mass: ", M_final/Msun
         write(*,*) "Previous Separation = ", CE_companion_position, " Final Separation: ", R_final/Rsun
         write(*,*) "Previous Orbital Energy = ", E_init, " Final Orbital Energy: ", E_final
         write(*,*) "Total Stellar Energy = ", s% total_energy
         write(*,*) "Previous Angular momentum = ", J_init, " Final Angular momentum: ", J_final
         write(*,*) "Dissipated Energy Rate: ", s% xtra1, " Dissipated Angular Momentum Rate: ", s% xtra6


      end subroutine CE_orbit_adjust



      subroutine calc_quantities_at_comp_position(id, ierr)

         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         real(dp) :: CE_companion_position, CE_companion_mass, omega_at_companion, csound_at_companion, P, M2
         real(dp) :: r_acc, v_rel, v_rel_div_csound, M_encl, rho_at_companion, scale_height_at_companion
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return


         CE_companion_position = s% xtra2
         CE_companion_mass = s% xtra4
         ! Calculate quantities at the position of the companion

         ! Keplerian velocity calculation depends on mass contained within a radius
         ! Include all the enclosed cells
         ! Add to it the enclosed mass of the current cell
         k = 1
         do while (s% r(k) > CE_companion_position * Rsun)
            k = k + 1
         end do

         M_encl = s% m(k)
         M_encl = M_encl + s% dm(k-1) * (CE_companion_position*Rsun - s% r(k)) / (s% r(k-1) - s% r(k))

         M2 = CE_companion_mass * Msun

         ! Determine orbital period in seconds
         P = AtoP(M_encl, M2, CE_companion_position*Rsun)
         ! Determine Keplerian velocity. Then subtract the local rotation velocity
         v_rel = 2.0 * pi * CE_companion_position*Rsun / P
         omega_at_companion = s% omega(k) * (CE_companion_position*Rsun - s% r(k)) * &
              (s% omega(k-1)-s% omega(k)) / (s% r(k-1) - s% r(k))
         rho_at_companion = s% rho(k) * (CE_companion_position*Rsun - s% r(k)) * &
              (s% rho(k-1)-s% rho(k)) / (s% r(k-1) - s% r(k))
         v_rel = v_rel - omega_at_companion * CE_companion_position*Rsun ! local rotation velocity = omega * r

         ! Determine Mach number
         csound_at_companion = s% csound(k) + (CE_companion_position*Rsun - s% r(k)) * &
              (s% csound(k-1)-s% csound(k)) / (s% r(k-1) - s% r(k))
         v_rel_div_csound = v_rel / csound_at_companion

         ! Determine accretion radius
         R_acc = 2.0 * standard_cgrav * M2 / (v_rel*v_rel)

         scale_height_at_companion =  s% scale_height(k) + (CE_companion_position*Rsun - s% r(k)) * &
              (s% scale_height(k-1)-s% scale_height(k)) / (s% r(k-1) - s% r(k))


         !saving these values to xtra variable so that tey are used in different CE_inject cases,
         ! in the torque calculations, and saved in the history file
         s% xtra12 = R_acc
         s% xtra13 = M_encl
         s% xtra14 = v_rel
         s% xtra15 = v_rel_div_csound
         s% xtra16 = rho_at_companion
         s% xtra17 = scale_height_at_companion

      end subroutine calc_quantities_at_comp_position




      real (dp) function AtoP(M1, M2, A)
         real(dp), intent(in) :: M1 ! in g
         real(dp), intent(in) :: M2 ! in g
         real(dp), intent(in) :: A ! in cm

         ! Kepler's 3rd Law - return orbital period in seconds
         AtoP = 2.0*pi * sqrt(A*A*A / (standard_cgrav * (M1+M2)))

      end function AtoP


      real(dp) function TukeyWindow(x,a)
         use const_def, only: dp, pi
         real(dp), intent(in) :: x, a

         if ((x .le. -0.5) .or. (x .ge. 0.5)) then
            TukeyWindow = 0.
         else if (x .le. -0.5 + a) then
            TukeyWindow = 0.5 - 0.5*cos(pi*(x+0.5)/a)
         else if (x .ge. 0.5 - a) then
            TukeyWindow = 0.5 - 0.5*cos(-pi*(x-0.5)/a)
         else
            TukeyWindow = 1.
         endif

      end function TukeyWindow


      end module CE_orbit
