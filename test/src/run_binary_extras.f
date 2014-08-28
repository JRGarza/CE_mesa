! ***********************************************************************
!
!   Copyright (C) 2012  Bill Paxton
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
 
      module run_binary_extras 

      use star_lib
      use star_def
      use const_def
      use const_def
      use chem_def
      use binary_def
      
      implicit none
      
      contains
      
      subroutine extras_binary_controls(ierr)
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         ierr = 0
         call binary_ptr(b,ierr)
      end subroutine extras_binary_controls

      integer function how_many_extra_binary_history_columns(b)
         use binary_def, only: binary_info
         type (binary_info), pointer :: b
!         how_many_extra_binary_history_columns = 0
         how_many_extra_binary_history_columns = 2

      end function how_many_extra_binary_history_columns
      
      subroutine data_for_extra_binary_history_columns(b, n, names, vals, ierr)
!      subroutine data_for_extra_binary_history_columns(b, s, n, names, vals, ierr)
         use const_def, only: dp
         use star_def, only: star_info
         type (binary_info), pointer :: b
!         type (star_info), pointer :: s
         integer, intent(in) :: n
         character (len=maxlen_binary_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         
         include 'formats'
         
         ierr = 0
       
         names(1) = "surface_dynamic_timescale"
         vals(1) = 1.0/dsqrt(standard_cgrav * b% s_donor% rho(1))
!         vals(1) = 1.0/dsqrt(standard_cgrav * s% rho(1))
         names(2) = "donor_radius"
         vals(2) = b% s_donor% r(1)
!         vals(2) = s% r(1)

      end subroutine data_for_extra_binary_history_columns

      end module run_binary_extras
      
