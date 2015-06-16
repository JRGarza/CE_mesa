! ***********************************************************************
!
!   Copyright (C) 2012  Bill Paxton
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
 
      module CE_adjust_mdot

      use star_def
      use const_def

      implicit none
      
            
      contains
      
      ! set use_other_adjust_mdot = .true. to enable this.
      ! your routine will be called after winds and before mass adjustment
   
      subroutine CE_remove_unbound_envelope(id, ierr)

         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         real(dp) :: mass_to_remove

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         k=1
         mass_to_remove = 0.0d0
         do while ((k < s% nz) .and. (.not. is_bound(k)))
            mass_to_remove = mass_to_remove + s% dm(k)
            k=k+1
         enddo
         s% mstar_dot = s% mstar_dot - (mass_to_remove/Msun) / (s% dt / secyer) !In Msun/yr

         !write(*,*) "CE_adjust_mdot ", mass_to_remove/Msun, (mass_to_remove/Msun) / (s% dt / secyer)





         contains

            logical function is_bound(k)
               integer, intent(in) :: k
               real(dp) :: val, f_energy
               logical :: include_internal_energy


               include_internal_energy = s% x_logical_ctrl(1)
               f_energy = logic2dbl(include_internal_energy)

               if (k == 1) then
                  val = s% energy(k)
               else if (k == s% nz) then
                  val = (s% dm(k)*s% energy(k) + &
                              0.5d0*s% dm(k-1)*s% energy(k-1))/ &
                        (s% dm(k) + 0.5d0*s% dm(k-1))
               else
                  val = (s% dm(k)*s% energy(k) + &
                              s% dm(k-1)*s% energy(k-1))/ &
                        (s% dm(k) + s% dm(k-1))
               end if
               val = val * f_energy - s% cgrav(k)*s% m_grav(k)/s% r(k) + &
                           0.5d0*s% v(k)*s% v(k)

               if (val > 0.0d0) then
                  is_bound = .false.
               else
                  is_bound = .true.
               endif

            end function is_bound


            real(dp) function logic2dbl(a)
               logical, intent(in) :: a

               if (a) then
                  logic2dbl = 1.d0
               else
                  logic2dbl = 0.d0
               end if
            end function logic2dbl

      end subroutine CE_remove_unbound_envelope






      end module CE_adjust_mdot
      
      
      
      
