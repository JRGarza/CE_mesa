!   This is a test program to attempt to externalize the routines in /mesa/binary
!   for the calculation of orbital separation evolution. This is not necessary, but 
!   done to make the code more modular, specifically to ease the inclusion of
!   common envelope evolution into MESA. Woo Woo!


      module binary_separation
            
      use star_lib
      use crlibm_lib
      use star_def
      use const_def
      use binary_def
      
      implicit none

      contains
      
      
      logical function exit_CE(b)
         use binary_def, only: binary_info
         use binary_mdot
         type (binary_info), pointer :: b

         real(dp) :: temp_mdot

         include 'formats.inc'

         
         ! calculate the Kolb explicit MT rate
         call get_info_for_ritter(b)
         call get_info_for_kolb(b)
         temp_mdot = b% mdot_thin + b% mdot_thick
         
         ! Exit CE if Mdot is less than some small mass transfer rate
         if(temp_mdot .lt. 1.0d-8) then
            exit_CE = .true.
         else
            exit_CE = .false.
         end if

      end function exit_CE
      
      logical function check_CE(b)
         use binary_def, only: binary_info
         type (binary_info), pointer :: b
         type (star_info), pointer :: s
         real(dp) :: mdot, rl2, R_accretor
         CHARACTER(LEN=40) :: FMT_1
         CHARACTER(LEN=40) :: FMT_2
         
         integer :: i
         real(dp) :: temp_M, temp_R, temp_t
         real(dp) :: sound_xing_time
         real(dp) :: q_rl, q_rl_old, zeta_star, zeta_rl
         
         include 'formats.inc'
         
         b => binary
         s => b% s_donor

         mdot = b% mtransfer_rate_old

         R_accretor = accretor_mass_radius_relation(b% m2 * msol)
         rl2 = b% rl(b% a_i)
         
         FMT_1 = '(A,1pe16.9)'
         FMT_2 = '(A,1pe16.9,A,1pe16.9)'
         write(*,FMT_1) "Old mass transfer rate (msun/yr): ", b% mtransfer_rate_old * secyer / msol
         write(*,FMT_1) " Current accretion rate (msun/yr): ", b% mtransfer_rate * secyer / msol
!         write(*,FMT_1) "Mass change: ", s% mstar_old - s% mstar
         write(*,FMT_1) "Time step (yr): ",b% s_donor% dt / secyer
         write(*,FMT_2) "mstar = ", s% mstar, " dynamical timescale = ", 1.0/dsqrt(standard_cgrav * s% rho(1))
     !    write(*,FMT_1) "MESA calculated dynamical timescale = ", b% s_donor% dynamic_timescale
     !    write(*,FMT_1) "Textbook dynamical timescale - sqrt(R^3/GM) = ", &
     !          dsqrt(b% r(b% d_i)*b% r(b% d_i)*b% r(b% d_i)/standard_cgrav/b% m(b% d_i))
     !    write(*,FMT_1) "Outer layer dynamical timescale = ", 1.0/dsqrt(standard_cgrav * s% rho(1))\
         
         write(*,'(3x(A,1pe16.9))') 'R_1 = ', b% r(b% d_i), ' R_2 = ', R_accretor, &
               'separation = ', b% separation
         write(*,'(2x(A,1pe16.9))') "Donor RL: ", b% rl(b% d_i), " Accretor RL:", b% rl(b% a_i)

 !        i = 0
 !        temp_M = 0.0
 !        sound_xing_time = 0.0
 !        do while(temp_M .lt. b% mtransfer_rate * b% dt)
 !           i = i + 1
 !           temp_M = temp_M + s% dm(i)
 !           sound_xing_time = sound_xing_time + s% dr_div_csound(i)
 !        enddo

!         write(*,FMT_1) "Sound crossing time of star = ", sound_xing_time 
      !   write(*,FMT_1) "Outer layer mass / Sound crossing time = ", s% dm(1) / sound_xing_time
 
         temp_M = 0.0;
         i = 0
         temp_t = 1.0/dsqrt(standard_cgrav * s% rho(1))
         do while((temp_M < s% dt * s% mstar/ temp_t) .and. (i .lt. s% nz))
            i = i + 1
            temp_M = temp_M + s% dm(i)
         enddo
         
         
!         write(*,*) "# of cells removed for instability to ensue = ", i
!         write(*,FMT_1) "Percent of star radius required to reach this = ", (s% r(1) - s% r(i)) / s% r(1)
         write(*,FMT_1) "Temporary instability Mdot = ", 1.0d-2
         write(*,FMT_1) " Mass Transfer Rate (msun/yr): ", dabs(mdot) * secyer / msol
         
         
!         if( dabs(mdot) > 1.0d0*(s% mstar)*dsqrt(standard_cgrav * s% rho(1))) then
         ! Dynamically Unstable RLOF
         if( dabs(mdot) > 1.0d-2 * msol / secyer ) then
            check_CE = .true.
            write(*,*) "Dynamically unstable RLOF"

         else if ((b% r(b% d_i) > b% rl(b% d_i)) .and. (R_accretor > b% rl(b% a_i))) then
            ! Contact Binary
            check_CE = .true.
            write(*,*) "Contact binary"

         ! Double Common Envelope
         else if (b% rl(b% a_i) < b% r(b% a_i)) then
            check_CE = .true.
            write(*,*) "Double Common Envelope"
         else
            check_CE = .false.
         endif
         
         if (check_CE) write(*,*) "Entered Common Envelope"
     
      end function check_CE

      logical function check_merger(b)
         use binary_def, only: binary_info
         type (binary_info), pointer :: b
         type (star_info), pointer :: s
         real(dp) :: R_donor, R_accretor
         include 'formats.inc'
         
         b => binary
         s => b% s_donor

         R_accretor = accretor_mass_radius_relation(b% m(b% a_i))

         if (s% center_he4 .lt. 1d-3) then      ! AGB star
            R_donor = s% c_core_radius*Rsun
         else if (s% center_h1 .lt. 1d-3) then  ! RGB star
            R_donor = s% he_core_radius*Rsun
         else                                   ! MS star
            R_donor = s% r(1)
         endif

         if (R_accretor + R_donor .gt. b% separation) then
            check_merger = .true.
         else
            check_merger = .false.
         endif
      end function check_merger

      real(dp) function accretor_mass_radius_relation(mass)
         real(dp), intent(in) :: mass

         ! Basic mass-radius relation for stars less than 2 Msun
         accretor_mass_radius_relation = Rsun * (mass/Msun) ** (0.9)

      end function accretor_mass_radius_relation

      real(dp) function calc_naive_A_f(s)
         type (binary_info), pointer :: b
         type (star_info), pointer :: s
         real(dp) :: E_bind_ml, E_orb_init, af_E_orb_f
         real(dp) :: E_int_spec, E_bind_spec, E_kin_spec, E_enth_spec
         real(dp) :: Delm2
         integer :: i
         include 'formats.inc'

         b => binary
         s => b% s_donor

         ! ------------ A naive calculation of A_f --------------- !

         ! Set all the specific energies to zero
         E_int_spec = 0.0d0
         E_bind_spec = 0.0d0
         E_kin_spec = 0.0d0
         E_enth_spec = 0.0d0         
         E_bind_ml = 0.0d0

         ! Loop over grid cells, outside going in
         i = 1
         Delm2 = s% mstar - s% he_core_mass*Msun
         do while ((s% mstar * (1.0d0 - s% q(i)) .le. Delm2) .and. (i < s% nz))
            E_int_spec = dexp(s% lnE(i))
            E_bind_spec = - standard_cgrav * s% mstar / s% r(i)
            E_kin_spec = s% v(i) * s% v(i) / 2.0
            E_enth_spec = s% P(i)/s% rho(i)

            E_bind_ml = E_bind_ml - (E_int_spec + E_bind_spec + E_kin_spec + E_enth_spec) &
                        * s% mstar * s% dq(i)

            i = i +1
         enddo

         if (i > 1) then
            E_int_spec =  dexp(s% lnE(i))
            E_bind_spec = - standard_cgrav * s% mstar / s% r(i)
            E_kin_spec = s% v(i) * s% v(i) / 2.0
            E_enth_spec = s% P(i)/s% rho(i)

            E_bind_ml = E_bind_ml - (E_int_spec + E_bind_spec + E_kin_spec + E_enth_spec) &
                         * (Delm2 - s% mstar * (1.0d0 - s% q(i-1)))
         endif

         E_orb_init = b% alpha_CE * standard_cgrav * b% m(b% a_i)*Msun * s% mstar / (2.0d0 * b% separation)
         af_E_orb_f = b% alpha_CE * standard_cgrav * b% m(b% a_i)*Msun * s% he_core_mass*Msun / 2.0d0
         calc_naive_A_f = af_E_orb_f / (E_orb_init + E_bind_ml)
 

      end function calc_naive_A_f

      subroutine new_separation_CE(b)
         use binary_def, only: binary_info
         type (binary_info), pointer :: b
         type (star_info), pointer :: s
         include 'formats.inc'
         
         real(dp) :: Delm2, E_orb_init, E_orb_fin, af_E_orb_f
         real(dp) :: E_bind_ml, E_int_spec, E_bind_spec, E_kin_spec, E_enth_spec
         
         integer :: i
         
         s => b% s_donor

         Delm2 = s% mstar_old - s% mstar

         ! Set all the specific energies to zero
         E_int_spec = 0.0d0
         E_bind_spec = 0.0d0
         E_kin_spec = 0.0d0
         E_enth_spec = 0.0d0         
         E_bind_ml = 0.0d0
                  
         ! Loop over grid cells, outside going in
         i = 1
         do while ((s% mstar_old * (1.0d0-s% q_old(i)) .le. Delm2) .and. (i < s% nz_old)) 


            E_int_spec =  dexp(b% CE_lnE_old(i))
            E_bind_spec = - standard_cgrav * s% mstar_old / dexp(s% xh_old(s% i_lnR,i))  
            E_kin_spec = b% CE_vel_old(i) * b% CE_vel_old(i) / 2.0
            E_enth_spec = b% CE_P_old(i) / b% CE_rho_old(i)
            E_bind_ml = E_bind_ml + (E_int_spec + E_bind_spec + E_kin_spec + E_enth_spec) &
                * s% mstar * s% dq(i)

!            write(*,*) Delm2, " ", s% mstar_old * (1.0d0-s% q_old(i)), " ", &
!                       E_int_spec, " ", E_bind_spec, " ", E_kin_spec, " ", E_enth_spec

            i = i + 1
         enddo

         if (i > 1) then
            E_int_spec =  dexp(b% CE_lnE_old(i))
            E_bind_spec = - standard_cgrav * s% mstar_old / dexp(s% xh_old(s% i_lnR,i))  
            E_kin_spec = b% CE_vel_old(i) * b% CE_vel_old(i) / 2.0
            E_enth_spec = b% CE_P_old(i) / b% CE_rho_old(i)
            E_bind_ml = E_bind_ml + (E_int_spec + E_bind_spec + E_kin_spec + E_enth_spec) &
                * (Delm2 - s% mstar_old * (1.0d0 - s% q_old(i-1)))               
         endif
               
         ! Calculate the initial orbital energy
         E_orb_init = - standard_cgrav * b% m_old(b% a_i) * s% mstar_old / (2.0d0 * b% separation)

         ! Calculate the final orbital energy
         E_orb_fin = E_bind_ml / b% alpha_CE + E_orb_init

         write(*,*) "Initial Orb energy: ", E_orb_init, " Final Orb Energy: ", E_orb_fin

         ! Calculate the final orbital separation
         af_E_orb_f = - standard_cgrav * b% m(b% a_i) * s% mstar / (2.0d0 * E_orb_fin)

         write(*,*) "Initial Separation: ", b% separation, " Final separation: ", af_E_orb_f

         ! Set orbital separation to calculated value
         b% separation = af_E_orb_f

      end subroutine new_separation_CE


      subroutine new_separation_jdot(b)
         use binary_jdot, only: get_jdot
         type(binary_info), pointer :: b
         include 'formats.inc'



         ! solve the winds in the system for jdot calculation,
         ! these don't include mass lost due to mass_transfer_efficiency < 1.0
         b% mdot_system_wind(b% d_i) = b% s_donor% mstar_dot - b% mtransfer_rate
         if (b% evolve_both_stars) then
            b% mdot_system_wind(b% a_i) = b% s_accretor% mstar_dot + &
                b% mtransfer_rate * b% xfer_fraction
         else
            b% mdot_system_wind(b% a_i) = 0.0d0
         end if

         ! get jdot and update orbital J
         b% jdot = get_jdot(b% mtransfer_rate, b% xfer_fraction)
         b% angular_momentum_j = b% angular_momentum_j + b% jdot*b% s1% time_step*secyer

         if (b% angular_momentum_j <= 0) then
            stop 'bad angular_momentum_j'
         end if
         
         ! use the new j to calculate new separation
         
         b% separation = ((b% angular_momentum_j/(b% m(1)*b% m(2)))**2) *&
             (b% m(1)+b% m(2)) / b% s1% cgrav(1)
         if (b% separation < b% min_binary_separation) &
            b% min_binary_separation = b% separation





      end subroutine new_separation_jdot


      end module binary_separation