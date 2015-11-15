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

      module CE_timestep
      
      use const_def
      use star_def

      implicit none


      contains


      integer function CE_pick_next_timestep(s)
         use const_def, only: secyer
         type (star_info), pointer :: s
         
         real(dp) dt_E, dt_A, dt_J, dt_min
         integer why_Tlim
         

!         ierr = 0
!         call star_ptr(id, s, ierr)
!         if (ierr /= 0) return

         CE_pick_next_timestep = keep_going
         

         ! CE_energy_rate, CE_companion_position, and CE_ang_mom_transferred are all 
         ! the previous values. We want to limit the change in the next values which
         ! are saved to the s% xtra's

         dt_E = CE_check_energy(s)

         dt_A = CE_check_separation(s)

         dt_J = CE_check_ang_mom(s)

         ! Get minimum next timestep
         dt_min = min(dt_E, dt_A, dt_J)

         write(*,*) "Timestep: ", dt_E, dt_A, dt_J, dt_min 
         write(*,*) "Next Timestep: ", s% time_step, s% dt_next, dt_min

         ! Create flag, and throw warning if not working properly
         if (dt_min == dt_E) then
            why_Tlim = 1
         else if (dt_min == dt_A) then
            why_Tlim = 2
         else if (dt_min == dt_J) then
            why_Tlim = 3
         else
            stop 'Something wrong with CE timestep controls'
         end if
         
         if (dt_min <= s% min_timestep_limit) then
            CE_pick_next_timestep = retry
            return
         endif
         
         
         s% max_years_for_timestep = dt_min
         
!         ! Set next timestep
!         if (s% dt_next > dt_min) then
!            s% dt_next = dt_min
!            s% why_Tlim = why_Tlim
!         end if

!         write(*,*) "Final next timestep: ", log10(s% dt_next / secyer)


      end function CE_pick_next_timestep




      real(dp) function CE_check_energy(s)
         use const_def, only: secyer
         type (star_info), pointer :: s
         real(dp) :: dE_fraction, dE_limit


         ! Fractional energy change
         dE_fraction = abs((s% xtra1 - s% xtra1_old) / s% xtra1_old)

         ! Limit from inlist
         dE_limit = s% x_ctrl(8)

         ! Calculate next timestep
         if (dE_fraction > dE_limit) then
            CE_check_energy = s% time_step / (dE_fraction / dE_limit) * secyer + 1d-99
         else
            CE_check_energy = 1d99
         end if
               
         
      end function CE_check_energy




      real(dp) function CE_check_separation(s)
         use const_def, only: secyer
         type (star_info), pointer :: s
         real(dp) :: dA_fraction, dA_limit


         ! Fractional separation change
         dA_fraction = abs((s% xtra2 - s% xtra2_old) / s% xtra2_old)

         ! Limit from inlist
         dA_limit = s% x_ctrl(9)

         ! Calculate next timestep
         if (dA_fraction > dA_limit) then
            CE_check_separation = s% time_step / (dA_fraction / dA_limit) * secyer + 1d-99
         else
            CE_check_separation = 1d99
         end if

      end function CE_check_separation




      real(dp) function CE_check_ang_mom(s)
         use const_def, only: secyer
         type (star_info), pointer :: s
         real(dp) :: dJ_fraction, dJ_limit

         ! Fractional angular momentum change
         dJ_fraction = abs((s% xtra6 - s% xtra6_old) / s% xtra6_old)

         ! Limit from inlist
         dJ_limit = s% x_ctrl(10)

         ! Calculate next timestep
         if (dJ_fraction > dJ_limit) then
            CE_check_ang_mom = s% time_step / (dJ_fraction / dJ_limit) * secyer + 1d-99
         else
            CE_check_ang_mom = 1d99
         end if


      end function CE_check_ang_mom


      integer function worst_result(result1, result2)
         integer, intent(in) :: result1, result2
         
         if(result1 == terminate .or. result2 == terminate) then
            worst_result = terminate
            return
         end if

         if(result1 == backup .or. result2 == backup) then
            worst_result = backup
            return
         end if
         
         if(result1 == retry .or. result2 == retry) then
            worst_result = retry
            return
         end if
         
         if(result1 == redo .or. result2 == redo) then
            worst_result = redo
            return
         end if

         worst_result = keep_going
         return
                              
      end function worst_result


      end module CE_timestep







