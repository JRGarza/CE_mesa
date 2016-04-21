! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful,
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************

 ! copied the content from star/job/standard_run_satr_extras.inc here and made the necessary changes

      module run_star_extras

      use star_lib
      use run_star_support, only: failed
      use star_def
      use const_def
      ! Add here all the external modules for CE_mesa here
      use CE_orbit
      use CE_energy
      use CE_torque
      use CE_after_struct_burn_mix
      use CE_before_struct_burn_mix
      use CE_adjust_mdot
      use CE_timestep

      implicit none

      ! these routines are called by the standard run_star check_model
      contains


      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)

         ! Here we should point to the names of the "other_" functions to be used
         s% other_energy => CE_inject_energy
         s% other_torque => CE_inject_am ! NEEDS TESTING
         s% other_before_struct_burn_mix => calc_recombination_before_struct_burn_mix
         s% other_after_struct_burn_mix => CE_other_after_struct_burn_mix
         s% other_adjust_mdot => CE_other_adjust_mdot

         ! Reading values of parameters from the extra controls that we are using
         ! Note that "extra_heat" is the specific energy added to the the  cell in units of erg/s/gr

         !s% xtra1 -> CE_energy_rate. It is initially set to 0. It will be calculated when CE_energy is called
         s% xtra1 = 0.0d0
         !s% xtra3 -> CE_companion_radius
         s% xtra3 = s% x_ctrl(3)
         !s% xtra4 -> CE_companion_mass
         s% xtra4 = s% x_ctrl(4)
         !s% xtra5 -> CE_n_acc_radii
         s% xtra5 = s% x_ctrl(5)
         !s% xtra6 -> CE_torque. It is initially set to 0. It will be calculated when CE_torque is called
         s% xtra6 = 0.0d0
         !s% xtra7 -> CE_mdot. It is initially set to 0. It will be calculated when CE_adjust_mdot is called
         s% xtra7 = 0.0d0
         !s% xtra20 -> L_acc. Accretion luminosity, calculated in CE_energy
         s% xtra20 = 0.0

         !s% xtra7 -> CE_test_case
         s% ixtra1 = s% x_integer_ctrl(1)

         ! s% job% relax_omega = .true.
         ! s% job% new_omega = s% x_ctrl(15) * 2.*pi/AtoP(1.496112*Msun,s% xtra4*Msun,s% xtra2*Rsun)
         ! write(*,*) s% job% new_omega, s%x_ctrl(15), 1.496112*Msun,s% xtra4*Msun,s% xtra2*Rsun
         ! ! ! We set a very small timestep during the relaxation phase, so that the star does not evolve significantly
         ! ! s% job% relax_omega_max_yrs_dt = 1d-8
         !  s% job% set_initial_dt = .True.
         !  s% job% years_for_initial_dt = 1d-8


      end subroutine extras_controls








      integer function extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp) :: CE_companion_position, R_acc, CE_n_acc_radii
         integer :: CE_test_case
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_startup = 0
         if (.not. restart) then
            call alloc_extra_info(s)
         else ! it is a restart
            call unpack_extra_info(s)
         end if


         !If we are restarting from a photo, the rest of the synchronization and relaxing steps should be skipped
         if (restart) then

            s% job% set_initial_model_number = .false.

            s% job% change_v_flag = .true.
            s% job% change_initial_v_flag = .false.
            s% job% new_v_flag = .true.

            s% job% new_rotation_flag = .false.
            s% job% change_rotation_flag = .false.

            s% job% set_initial_age = .false.
            s% job% set_initial_model_number = .false.

            return
         endif


         ! Reading values of parameters from the extra controls that we are using
         ! Note that "extra_heat" is the specific energy added to the the  cell in units of erg/s/gr

         !s% xtra1 -> CE_energy_rate. It is initially set to 0. It will be calculated when CE_energy is called
         s% xtra1 = 0.0d0
         !s% xtra3 -> CE_companion_radius
         s% xtra3 = s% x_ctrl(3)
         !s% xtra4 -> CE_companion_mass
         s% xtra4 = s% x_ctrl(4)
         !s% xtra5 -> CE_n_acc_radii
         s% xtra5 = s% x_ctrl(5)
         !s% xtra6 -> CE_torque. It is initially set to 0. It will be calculated when CE_torque is called
         s% xtra6 = 0.0d0
         !s% xtra7 -> CE_mdot. It is initially set to 0. It will be calculated when CE_adjust_mdot is called
         s% xtra7 = 0.0d0

         !s% xtra7 -> CE_test_case
         s% ixtra1 = s% x_integer_ctrl(1)
         !s% xtra2 -> CE_companion_position = CE_companion_initial_position * Rsatr
         s% xtra2 = s% x_ctrl(2) * s% r(1) / Rsun


         s% job% relax_omega = .true.
         s% job% new_omega = s% x_ctrl(15) * 2.*pi/AtoP(s% m(1),s% xtra4*Msun,s% xtra2 * Rsun)
         ! We set a very small timestep during the relaxation phase, so that the star does not evolve significantly
         s% job% relax_omega_max_yrs_dt = 1d-8
         s% job% set_initial_dt = .True.
         s% job% years_for_initial_dt = 1d-8


         ! We are calling here the relax_omega, because we want to first have loaded the model so that we know its radius, and mass.
         if (s% rotation_flag .and. s% job% relax_omega) then
            write(*,*) 'new_omega =', s% job% new_omega
            call star_relax_uniform_omega( &
               id, 0, s% job% new_omega, s% job% num_steps_to_relax_rotation,&
               s% job% relax_omega_max_yrs_dt, ierr)
            if (failed('star_relax_uniform_omega',ierr)) return
            s% job% relax_omega = .false.
         else
            !call star_relax_num_steps(id, 100, 1d-8 * secyer, ierr)
         endif

         !After relaxation is done, the timestep automatically increases to a "large" timestep. Here we are tryying to make this
         !transition smoother
         s% dt_next = 1d-8 * secyer


         CE_companion_position = s% xtra2
         CE_test_case = s% ixtra1
         CE_n_acc_radii = s% x_ctrl(5)
         call calc_quantities_at_comp_position(id, ierr)
         R_acc = s% xtra12

         ! We need to increase the resolution around the area where the extra heat is deposited
         ! We will do this at the startup and also in the extra_check model, since the position
         ! of the companion will be changing
         if (CE_test_case == 2 .or. CE_test_case == 3 .or. CE_test_case == 4) then
            s% R_function2_param1 = CE_companion_position/(s%r(1)/Rsun) + 2.0* CE_n_acc_radii * R_acc/s%r(1)
            s% R_function2_param2 = CE_companion_position/(s%r(1)/Rsun) - 2.0* CE_n_acc_radii * R_acc/s%r(1)
         endif




      end function extras_startup


      ! returns either keep_going, retry, backup, or terminate.
      integer function extras_check_model(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr, result
         type (star_info), pointer :: s
         real(dp) :: CE_energy_rate, CE_companion_position, CE_companion_radius, CE_companion_mass
         real(dp) :: CE_ang_mom_transferred, R_acc, CE_n_acc_radii
         integer :: CE_test_case

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         result = keep_going

         ! Reading initial values of parameters from the extra controls that we are using
         ! Note that "extra_heat" is the specific energy added to the the  cell in units of erg/s/gr
         CE_energy_rate = s% xtra1
         CE_companion_position = s% xtra2
         CE_companion_radius = s% xtra3
         CE_companion_mass = s% xtra4
         CE_ang_mom_transferred = s% xtra6
         CE_test_case = s% ixtra1


         call calc_quantities_at_comp_position(id, ierr)
         R_acc = s% xtra12

         ! We need to increase the resolution around the area where the extra heat is deposited
         ! We will do this at the startup and also in the extra_check model, since the position
         ! of the companion will be changing
         CE_n_acc_radii = s% x_ctrl(5)
         if (CE_test_case == 2 .or. CE_test_case == 3 .or. CE_test_case == 4) then
            s% R_function2_param1 = CE_companion_position/(s%r(1)/Rsun) + 2.0* CE_n_acc_radii * R_acc/s%r(1)
            s% R_function2_param2 = CE_companion_position/(s%r(1)/Rsun) - 2.0* CE_n_acc_radii * R_acc/s%r(1)
         endif


         ! For test cases 1 and 2 (heating of the whole envelope and of the base of the envelope) the code below must be skipped
         if (s% x_integer_ctrl(1) .ne. 1 .and. s% x_integer_ctrl(1) .ne. 2) then
            ! Adjust orbital separation based on energy deposited
            call CE_orbit_adjust(id, ierr)
            ! Added timestep controls
            result = worst_result(result, CE_pick_next_timestep(s))
         endif




         ! if you want to check multiple conditions, it can be useful
         ! to set a different termination code depending on which
         ! condition was triggered.  MESA provides 9 customizeable
         ! termination codes, named t_xtra1 .. t_xtra9.  You can
         ! customize the messages that will be printed upon exit by
         ! setting the corresponding termination_code_str value.
         ! termination_code_str(t_xtra1) = 'my termination condition'

         ! by default, indicate where (in the code) MESA terminated
         if (result == terminate) s% termination_code = t_extras_check_model

         extras_check_model = result
      end function extras_check_model


      integer function how_many_extra_history_columns(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 13
      end function how_many_extra_history_columns


      subroutine data_for_extra_history_columns(id, id_extra, n, names, vals, ierr)
         integer, intent(in) :: id, id_extra, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         !note: do NOT add the extras names to history_columns.list
         ! the history_columns.list is only for the built-in log column options.
         ! it must not include the new column names you are adding here.

         names(1) = 'CE_energy_rate'
         vals(1) = s% xtra1
         names(2) = 'L_accretion'
         vals(2) = s% xtra20
         names(3) = 'CE_torque'
         vals(3) = s% xtra6
         names(4) = 'CE_companion_position_r'
         vals(4) = s% xtra2
         names(5) = 'CE_companion_position_m'
         vals(5) = s% xtra9
         names(6) = 'CE_ang_mom_transferred'
         vals(6) = s% xtra6
         names(7) = 'envelope_binding_energy'
         vals(7) = s% xtra11
         names(8) = 'R_acc'
         vals(8) = s% xtra12
         names(9) = 'R_acc_low'
         vals(9) = s% xtra13
         names(10) = 'R_acc_high'
         vals(10) = s% xtra14
         names(11) = 'v_rel'
         vals(11) = s% xtra16
         names(12) = 'v_over_c_sound'
         vals(12) = s% xtra17
         names(13) = 'eta_pulse_wind' ! From Yoon & Cantiello (2010)
         vals(13) = s% xtra21

      end subroutine data_for_extra_history_columns


      integer function how_many_extra_profile_columns(id, id_extra)
         use star_def, only: star_info
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 2
      end function how_many_extra_profile_columns


      subroutine data_for_extra_profile_columns(id, id_extra, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp
         integer, intent(in) :: id, id_extra, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         !note: do NOT add the extra names to profile_columns.list
         ! the profile_columns.list is only for the built-in profile column options.
         ! it must not include the new column names you are adding here.

         ! here is an example for adding a profile column
         !if (n /= 1) stop 'data_for_extra_profile_columns'
         !names(1) = 'beta'
         !do k = 1, nz
         !   vals(k,1) = s% Pgas(k)/s% P(k)
         !end do

         names(1) = 'ionization_energy'
         names(2) = 'eps_recombination'
         do k = 1, nz
           vals(k,1) = s% xtra1_array(k)
           vals(k,2) = s% xtra2_array(k)
         end do

      end subroutine data_for_extra_profile_columns


      ! returns either keep_going or terminate.
      ! note: cannot request retry or backup; extras_check_model can do that.
      integer function extras_finish_step(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
         call store_extra_info(s)

         ! to save a profile,
            ! s% need_to_save_profiles_now = .true.
         ! to update the star log,
            ! s% need_to_update_history_now = .true.

         ! see extras_check_model for information about custom termination codes
         ! by default, indicate where (in the code) MESA terminated
         if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
      end function extras_finish_step


      subroutine extras_after_evolve(id, id_extra, ierr)
         integer, intent(in) :: id, id_extra
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         !Saving those values back to xtra_controls so that restarts work.
         s% x_ctrl(2) = s% xtra2 / (s% r(1) / Rsun)
         s% x_ctrl(3) = s% xtra3
         s% x_ctrl(4) = s% xtra4


      end subroutine extras_after_evolve


      ! routines for saving and restoring extra data so can do restarts

         ! put these defs at the top and delete from the following routines
         !integer, parameter :: extra_info_alloc = 1
         !integer, parameter :: extra_info_get = 2
         !integer, parameter :: extra_info_put = 3


      subroutine alloc_extra_info(s)
         integer, parameter :: extra_info_alloc = 1
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_alloc)
      end subroutine alloc_extra_info


      subroutine unpack_extra_info(s)
         integer, parameter :: extra_info_get = 2
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_get)
      end subroutine unpack_extra_info


      subroutine store_extra_info(s)
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_put)
      end subroutine store_extra_info


      subroutine move_extra_info(s,op)
         integer, parameter :: extra_info_alloc = 1
         integer, parameter :: extra_info_get = 2
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         integer, intent(in) :: op

         integer :: i, j, num_ints, num_dbls, ierr

         i = 0
         ! call move_int or move_flg
         num_ints = i

         i = 0
         ! call move_dbl

         num_dbls = i

         if (op /= extra_info_alloc) return
         if (num_ints == 0 .and. num_dbls == 0) return

         ierr = 0
         call star_alloc_extras(s% id, num_ints, num_dbls, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in star_alloc_extras'
            write(*,*) 'alloc_extras num_ints', num_ints
            write(*,*) 'alloc_extras num_dbls', num_dbls
            stop 1
         end if

         contains

         subroutine move_dbl(dbl)
            real(dp) :: dbl
            i = i+1
            select case (op)
            case (extra_info_get)
               dbl = s% extra_work(i)
            case (extra_info_put)
               s% extra_work(i) = dbl
            end select
         end subroutine move_dbl

         subroutine move_int(int)
            integer :: int
            i = i+1
            select case (op)
            case (extra_info_get)
               int = s% extra_iwork(i)
            case (extra_info_put)
               s% extra_iwork(i) = int
            end select
         end subroutine move_int

         subroutine move_flg(flg)
            logical :: flg
            i = i+1
            select case (op)
            case (extra_info_get)
               flg = (s% extra_iwork(i) /= 0)
            case (extra_info_put)
               if (flg) then
                  s% extra_iwork(i) = 1
               else
                  s% extra_iwork(i) = 0
               end if
            end select
         end subroutine move_flg

      end subroutine move_extra_info

      end module run_star_extras
