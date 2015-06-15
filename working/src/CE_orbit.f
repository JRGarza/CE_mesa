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
         real(dp) :: E_init, E_loss, E_final
         real(dp) :: M_inner

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return


         CE_energy_rate = s% xtra1
         CE_companion_position = s% x_ctrl(2)
         CE_companion_mass = s% x_ctrl(4)

         ! Mass for orbital energy calculation depends on mass contained within a radius
         ! This is a bit dirty, but should work
         k = 1
         M_inner = s% m(k)
         do while (s% r(k) > CE_companion_position * Rsun)
            M_inner = s% m(k) 
            k = k + 1
         end do

         ! Calculate the energies
         E_init = -standard_cgrav * CE_companion_mass * Msun * M_inner / (2.0 * CE_companion_position * Rsun)
         E_loss = CE_energy_rate * s% dt
         E_final = E_init - E_loss

         ! Calculate the new orbital separation
         CE_companion_position = -standard_cgrav * CE_companion_mass * Msun * M_inner / (2.0 * E_final * Rsun)

         s% x_ctrl(2) = CE_companion_position

         write(*,*) "Companion mass: ", CE_companion_mass*Msun, " Enclosed mass: ", M_inner
         write(*,*) "Separation: ", CE_companion_position*Rsun 
         write(*,*) "Current Orbital Separation = ", CE_companion_position 
         write(*,*) "Current Orbital Energy = ", E_init
 
      end subroutine CE_orbit_adjust

      end module CE_orbit
