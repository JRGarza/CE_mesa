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
      use CE_energy, only: AtoP, TukeyWindow
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
         real(dp) :: CE_companion_position, CE_companion_mass, CE_n_acc_radii, CE_torque
         real(dp) :: time, R_acc, Mach, M_encl, M2, vel, A, P

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return


         ! Load angular momentum dissipated in the envelope from the orbit decrease
         CE_companion_position = s% xtra2
         CE_companion_mass = s% xtra4
         CE_n_acc_radii = s% xtra5
         CE_torque = s% xtra6

         write(*,*) "CE torque ", CE_torque

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
         vel = 2.0 * pi * CE_companion_position*Rsun / P
         vel = vel - s% omega(k) * s% rmid(k) ! local rotation velocity = omega * rmid

         ! Determine accretion radius
         R_acc = 2.0 * standard_cgrav * M2 / (vel*vel)

         ! First calculate the mass in which the angular momentum will be deposited
         mass_to_be_spun = 0.0
         do k = 1, s% nz
            ff = TukeyWindow((s% r(k) - CE_companion_position*Rsun)/(CE_n_acc_radii * R_acc), a_tukey)
            mass_to_be_spun = mass_to_be_spun + s% dm(k) * ff
         end do

         ! Now redo the loop and add the extra torque
         do k = 1, s% nz
            ff = TukeyWindow((s% r(k) - CE_companion_position*Rsun)/(CE_n_acc_radii * R_acc), a_tukey)
            s% extra_jdot(k) = CE_torque / mass_to_be_spun * ff
         end do



      end subroutine CE_inject_am


      end module CE_torque
