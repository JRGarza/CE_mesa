 
      module run_binary
      implicit none
      
      contains
      
      subroutine do_run_binary(tst)
         use binary_lib, only: run1_binary
         use run_star_extras
         use run_binary_extras
         
         logical, intent(in) :: tst
         
         integer :: ierr
         
         call run1_binary(tst, &
            ! star extras
            extras_controls, &
            extras_startup, &
            extras_check_model, &
            how_many_extra_history_columns, &
            data_for_extra_history_columns, &
            how_many_extra_profile_columns, &
            data_for_extra_profile_columns, &
            extras_finish_step, &
            extras_after_evolve, &
            ! binary extras
            extras_binary_controls, &
            how_many_extra_binary_history_columns, &
            data_for_extra_binary_history_columns, &

            ierr)

      end subroutine do_run_binary

      end module run_binary
      
