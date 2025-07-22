! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
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

module run_star_extras

   use star_lib
   use star_def
   use const_def
   use math_lib
   use rates_def
   use chem_def

   implicit none

   ! s% xtra
   integer, parameter :: i_min_D_mix = 1

   ! s% lxtra
   integer, parameter :: i_use_dedt = 1
   integer, parameter :: i_use_conv_premix = 2
   integer, parameter :: i_use_gold_tolerances = 3

   ! Extra controls options

   ! s% x_ctrl
   integer, parameter :: i_aovhe_mod_pen_conv = 1
   integer, parameter :: i_aovhe_mod_pen_conv_exp = 2

   integer, parameter :: i_save_dT = 5
   integer, parameter :: i_save_dlogL = 6
   integer, parameter :: i_save_dHc = 7
   integer, parameter :: i_save_dlogHc = 8
   integer, parameter :: i_save_min_logHc = 9
   integer, parameter :: i_save_dHec = 10
   integer, parameter :: i_save_dlogHec = 11
   integer, parameter :: i_save_min_logHec = 12

   integer, parameter :: i_IGW_exponent = 15
   integer, parameter :: i_IGW_D_ext = 16
   integer, parameter :: i_IGW_D_ext_postMS = 17

   integer, parameter :: i_turb_constant = 18
   integer, parameter :: i_turb_exponent = 19
   integer, parameter :: i_turb_reference = 20

   ! s% x_integer_ctrl
   integer, parameter :: i_num_deltanu_for_q = 1

   ! s% x_logical_ctrl
   integer, parameter :: i_save_mesh = 1
   integer, parameter :: i_verbose_coupling = 2
   integer, parameter :: i_grow_aovhe = 5

   ! For saving routine
   real(dp) :: prev_Teff, prev_L, prev_Hc, prev_Hec
   logical :: first_step = .true.

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


      ! the extras functions in this file will not be called
      ! unless you set their function pointers as done below.
      ! otherwise we use a null_ version which does nothing (except warn).

      s% extras_startup => extras_startup
      s% extras_start_step => extras_start_step
      s% extras_check_model => extras_check_model
      s% extras_finish_step => extras_finish_step
      s% extras_after_evolve => extras_after_evolve
      s% how_many_extra_history_columns => how_many_extra_history_columns
      s% data_for_extra_history_columns => data_for_extra_history_columns
      s% how_many_extra_profile_columns => how_many_extra_profile_columns
      s% data_for_extra_profile_columns => data_for_extra_profile_columns

      s% how_many_extra_history_header_items => how_many_extra_history_header_items
      s% data_for_extra_history_header_items => data_for_extra_history_header_items
      s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
      s% data_for_extra_profile_header_items => data_for_extra_profile_header_items

      s% other_D_mix => turbulent_mixing
      s% other_adjust_mlt_gradT_fraction => other_adjust_mlt_gradT_fraction_Peclet
      s% other_overshooting_scheme => extended_convective_penetration

      ! To get through the helium flash turn off these options if they are used.
      s% lxtra(i_use_dedt) = (s% energy_eqn_option == 'dedt')
      s% lxtra(i_use_conv_premix) = s% do_conv_premix
      s% lxtra(i_use_gold_tolerances) = s% use_gold_tolerances

      s% xtra(i_min_D_mix) = s% min_D_mix
   end subroutine extras_controls


   subroutine extras_startup(id, restart, ierr)
      integer, intent(in) :: id
      logical, intent(in) :: restart
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
   end subroutine extras_startup


   integer function extras_start_step(id)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      extras_start_step = 0

      if (s% center_h1 < 1d-9) then  ! Switch to lower minDmix after MS
         s% min_D_mix = s%xtra(i_min_D_mix) * 0.2d0
      else
         s% min_D_mix = s%xtra(i_min_D_mix)
      end if
   end function extras_start_step


   ! returns either keep_going, retry, or terminate.
   integer function extras_check_model(id)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      extras_check_model = keep_going

      ! To get through the helium flash if using the dedt energy eqn,
      ! turn off dedt form of energy eqn if the time step becomes too small

      if (s% power_he_burn > 3d8) then
         if (s% lxtra(i_use_dedt)) then
            write(*, *) 'dedt off', s% model_number
            s% energy_eqn_option = 'eps_grav'
         end if

         if (s% lxtra(i_use_conv_premix)) then
            write(*, *) 'cpm off', s% model_number
            s% do_conv_premix = .false.
         end if

         if (s% lxtra(i_use_gold_tolerances)) then
            write(*, *) 'gold off', s% model_number
            s% use_gold_tolerances = .false.
         end if
      end if

      if (s% power_he_burn < 2.5d8) then
         if (s% lxtra(i_use_dedt)) then
            s% energy_eqn_option = 'dedt'
         end if

         if (s% lxtra(i_use_conv_premix)) then
            s% do_conv_premix = .true.
         end if

         if (s% lxtra(i_use_gold_tolerances)) then
            s% use_gold_tolerances = .true.
         end if
      end if

      ! if you want to check multiple conditions, it can be useful
      ! to set a different termination code depending on which
      ! condition was triggered.  MESA provides 9 customizeable
      ! termination codes, named t_xtra1 .. t_xtra9.  You can
      ! customize the messages that will be printed upon exit by
      ! setting the corresponding termination_code_str value.
      ! termination_code_str(t_xtra1) = 'my termination condition'

      ! by default, indicate where (in the code) MESA terminated
      if (extras_check_model == terminate) s% termination_code = t_extras_check_model
   end function extras_check_model


   integer function how_many_extra_history_columns(id)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      if (s% x_integer_ctrl(i_num_deltanu_for_q) <= 0) then
         how_many_extra_history_columns = 32 + 2 + 14
      else
         how_many_extra_history_columns = 32 + s% x_integer_ctrl(i_num_deltanu_for_q) * 2 * 6 + 2 + 14
      end if
   end function how_many_extra_history_columns


   subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
      integer, intent(in) :: id, n
      character (len = maxlen_history_column_name) :: names(n)
      real(dp) :: vals(n)
      integer, intent(out) :: ierr
      type (star_info), pointer :: s

      integer :: k, k2, k3, max_eps_h_k
      character (len = 1) :: pm
      real(dp) :: nu_q
      real(dp) :: out(30, 2)
      integer :: out_int(2, 2)

      integer :: i, k_l
      real(dp) :: r_bCZ, m_bCZ, alfa, beta

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      k = 1
      ! Calculate bottom of conv envelope
      max_eps_h_k = maxloc(s% eps_nuc_categories(ipp, 1:s% nz) + s% eps_nuc_categories(icno, 1:s% nz), 1)
      r_bCZ = -1d99
      m_bCZ = -1d99
      do i = max_eps_h_k, 2, -1
         if ((s% gradr(i) < s% grada(i)) .and. (s% gradr(i - 1) > s% grada(i - 1)) .and. (s% xa(s% net_iso(ih1),i) > 0.5)) then  ! botCZ between i and i-1
            alfa = (s% gradr(i) - s% grada(i)) / ((s% gradr(i) - s% grada(i)) - (s% gradr(i - 1) - s% grada(i - 1)))
            beta = 1 - alfa
            r_bCZ = (beta * s% r(i)**3 + alfa * s% r(i - 1)**3)**(1d0/3d0) / rsun
            m_bCZ = (beta * s% m(i) + alfa * s% m(i - 1)) / msun
            exit
         end if
      end do
      names(k) = 'r_botCZ'
      vals(k) = r_bCZ
      k = k + 1
      names(k) = 'm_botCZ'
      vals(k) = m_bCZ
      k = k + 1

      ! Calc stuff for dq/dnu
      if (s% x_integer_ctrl(i_num_deltanu_for_q) >= 0) then
         do k2 = 1, s% x_integer_ctrl(i_num_deltanu_for_q)
            do k3 = 1, 2
               if (k3 == 1) then
                  pm = 'p'
                  nu_q = s%nu_max + k2 * s%delta_nu
               else
                  pm = 'm'
                  nu_q = s%nu_max - k2 * s%delta_nu
               end if

               call do_strong_coupling(id, nu_q, out, out_int)
               write(names(k), '(A, I2.2)') 'r1_nu_' // pm, k2
               vals(k) = out(1, 1)
               k = k + 1
               write(names(k), '(A, I2.2)') 'r2_nu_' // pm, k2
               vals(k) = out(2, 1)
               k = k + 1
               write(names(k), '(A, I2.2)') 'm1_nu_' // pm, k2
               vals(k) = out(3, 1)
               k = k + 1
               write(names(k), '(A, I2.2)') 'm2_nu_' // pm, k2
               vals(k) = out(4, 1)
               k = k + 1
               write(names(k), '(A, I2.2)') 'q_nu_' // pm, k2
               vals(k) = out(5, 1)
               k = k + 1
               write(names(k), '(A, I2.2)') 'X_int_' // pm, k2
               vals(k) = out(6, 1)
               k = k + 1
            end do
         end do
      end if

      call do_strong_coupling(id, s% nu_max, out, out_int)

      names(k) = 'r_1'
      vals(k) = out(1, 1)  ! r_1/rsun
      k = k + 1

      names(k) = 'r_2'
      vals(k) = out(2, 1)  ! r_2/rsun
      k = k + 1

      names(k) = 'm_1'
      vals(k) = out(3, 1)  ! m_1
      k = k + 1

      names(k) = 'm_2'
      vals(k) = out(4, 1)  ! m_2
      k = k + 1

      names(k) = 'coupling_strong'
      vals(k) = out(5, 1)  ! q
      k = k + 1

      names(k) = 'X_integral_part'
      vals(k) = out(6, 1)  ! PQ_integral / pi
      k = k + 1

      names(k) = 'X_gradient_part'
      vals(k) = out(7, 1)  ! G2_div_2kappa_s0
      k = k + 1

      names(k) = 'dlnc_ds_s0_part_ip'
      vals(k) = out(8, 1)  ! dlnc_ds_s0_part_ip
      k = k + 1

      names(k) = 'dlnc_ds_s0_part_Sl'
      vals(k) = out(9, 1)  ! dlnc_ds_s0_part_Sl
      k = k + 1

      names(k) = 'dlnc_ds_s0_part_N'
      vals(k) = out(10, 1)  ! dlnc_ds_s0_part_N
      k = k + 1

      names(k) = 'kappa_s0'
      vals(k) = out(11, 1)  ! kappa_s0
      k = k + 1

      names(k) = 'NuAJ_s0'
      vals(k) = out(12, 1)  ! NuAJ_s0
      k = k + 1

      names(k) = 'dlnPds_s0'
      vals(k) = out(13, 1)  ! dlnPds_s0
      k = k + 1

      names(k) = 'dlnQds_s0'
      vals(k) = out(14, 1)  ! dlnQds_s0
      k = k + 1

      names(k) = 'k_P'
      vals(k) = out_int(1, 1)  ! k_P
      k_l = out_int(1, 1)
      k = k + 1

      names(k) = 'k_Q'
      vals(k) = out_int(2, 1)  ! k_Q
      k_l = max(k_l, out_int(2, 1))
      k = k + 1

      ! For second evan region
      names(k) = 'r_1b'
      vals(k) = out(1, 2)  ! r_1/rsun
      k = k + 1

      names(k) = 'r_2b'
      vals(k) = out(2, 2)  ! r_2/rsun
      k = k + 1

      names(k) = 'm_1b'
      vals(k) = out(3, 2)  ! m_1
      k = k + 1

      names(k) = 'm_2b'
      vals(k) = out(4, 2)  ! m_2
      k = k + 1

      names(k) = 'coupling_strong_b'
      vals(k) = out(5, 2)  ! q
      k = k + 1

      names(k) = 'X_integral_part_b'
      vals(k) = out(6, 2)  ! PQ_integral / pi
      k = k + 1

      names(k) = 'X_gradient_part_b'
      vals(k) = out(7, 2)  ! G2_div_2kappa_s0
      k = k + 1

      names(k) = 'dlnc_ds_s0_part_ip_b'
      vals(k) = out(8, 2)  ! dlnc_ds_s0_part_ip
      k = k + 1

      names(k) = 'dlnc_ds_s0_part_Sl_b'
      vals(k) = out(9, 2)  ! dlnc_ds_s0_part_Sl
      k = k + 1

      names(k) = 'dlnc_ds_s0_part_N_b'
      vals(k) = out(10, 2)  ! dlnc_ds_s0_part_N
      k = k + 1

      names(k) = 'kappa_s0b'
      vals(k) = out(11, 2)  ! kappa_s0
      k = k + 1

      names(k) = 'NuAJ_s0b'
      vals(k) = out(12, 2)  ! NuAJ_s0
      k = k + 1

      names(k) = 'dlnPds_s0b'
      vals(k) = out(13, 2)  ! dlnPds_s0
      k = k + 1

      names(k) = 'dlnQds_s0b'
      vals(k) = out(14, 2)  ! dlnQds_s0
      k = k + 1

      names(k) = 'k_u2b'
      vals(k) = out_int(1, 2)  ! k_u2
      k = k + 1

      names(k) = 'k_l2b'
      vals(k) = out_int(2, 2)  ! k_l2
      k = k + 1

      k2 = 17
      names(k) = 'a_0_grd_P'
      vals(k) = out(k2, 1)
      k = k + 1
      k2 = k2 + 1
      names(k) = 'a_1_grd_P'
      vals(k) = out(k2, 1)
      k = k + 1
      k2 = k2 + 1
      names(k) = 'a_2_grd_P'
      vals(k) = out(k2, 1)
      k = k + 1
      k2 = k2 + 1
      names(k) = 'a_3_grd_P'
      vals(k) = out(k2, 1)
      k = k + 1
      k2 = k2 + 1

      names(k) = 'a_0_grd_Q'
      vals(k) = out(k2, 1)
      k = k + 1
      k2 = k2 + 1
      names(k) = 'a_1_grd_Q'
      vals(k) = out(k2, 1)
      k = k + 1
      k2 = k2 + 1
      names(k) = 'a_2_grd_Q'
      vals(k) = out(k2, 1)
      k = k + 1
      k2 = k2 + 1
      names(k) = 'a_3_grd_Q'
      vals(k) = out(k2, 1)
      k = k + 1
      k2 = k2 + 1

      names(k) = 'a_0_int_P'
      vals(k) = out(k2, 1)
      k = k + 1
      k2 = k2 + 1
      names(k) = 'a_1_int_P'
      vals(k) = out(k2, 1)
      k = k + 1
      k2 = k2 + 1
      names(k) = 'a_2_int_P'
      vals(k) = out(k2, 1)
      k = k + 1
      k2 = k2 + 1

      names(k) = 'a_0_int_Q'
      vals(k) = out(k2, 1)
      k = k + 1
      k2 = k2 + 1
      names(k) = 'a_1_int_Q'
      vals(k) = out(k2, 1)
      k = k + 1
      k2 = k2 + 1
      names(k) = 'a_2_int_Q'
      vals(k) = out(k2, 1)
      k = k + 1
      k2 = k2 + 1

   end subroutine data_for_extra_history_columns


   integer function how_many_extra_profile_columns(id)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      how_many_extra_profile_columns = 5
   end function how_many_extra_profile_columns


   subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
      integer, intent(in) :: id, n, nz
      character (len = maxlen_profile_column_name) :: names(n)
      real(dp) :: vals(nz, n)
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      integer :: k
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      ! note: do NOT add the extra names to profile_columns.list
      ! the profile_columns.list is only for the built-in profile column options.
      ! it must not include the new column names you are adding here.

      names(1) = 'q_J'
      names(2) = 'q_A'
      names(3) = 'q_Nu'
      names(4) = 'q_P'
      names(5) = 'q_Q'

      !Have to repopulate xtra arrays
      call calc_JANu(id)
      call calc_PQ(id, (2d0 * pi / 1d6) * s% nu_max)

      do k = 1, nz
         vals(k, 1) = s% xtra1_array(k)
         vals(k, 2) = s% xtra2_array(k)
         vals(k, 3) = s% xtra3_array(k)
         vals(k, 4) = s% xtra4_array(k)
         vals(k, 5) = s% xtra5_array(k)
      end do

   end subroutine data_for_extra_profile_columns


   integer function how_many_extra_history_header_items(id)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      how_many_extra_history_header_items = 0
   end function how_many_extra_history_header_items


   subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
      integer, intent(in) :: id, n
      character (len = maxlen_history_column_name) :: names(n)
      real(dp) :: vals(n)
      type(star_info), pointer :: s
      integer, intent(out) :: ierr
      ierr = 0
      call star_ptr(id, s, ierr)
      if(ierr/=0) return

      ! here is an example for adding an extra history header item
      ! also set how_many_extra_history_header_items
      ! names(1) = 'mixing_length_alpha'
      ! vals(1) = s% mixing_length_alpha

   end subroutine data_for_extra_history_header_items


   integer function how_many_extra_profile_header_items(id)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      how_many_extra_profile_header_items = 0
   end function how_many_extra_profile_header_items


   subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
      integer, intent(in) :: id, n
      character (len = maxlen_profile_column_name) :: names(n)
      real(dp) :: vals(n)
      type(star_info), pointer :: s
      integer, intent(out) :: ierr
      ierr = 0
      call star_ptr(id, s, ierr)
      if(ierr/=0) return

      ! here is an example for adding an extra profile header item
      ! also set how_many_extra_profile_header_items
      ! names(1) = 'mixing_length_alpha'
      ! vals(1) = s% mixing_length_alpha

   end subroutine data_for_extra_profile_header_items


   ! returns either keep_going or terminate.
   ! note: cannot request retry; extras_check_model can do that.
   integer function extras_finish_step(id)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s
      logical :: save_now

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      extras_finish_step = keep_going

      ! to save a profile,
      ! s% need_to_save_profiles_now = .true.
      ! to update the star log,
      ! s% need_to_update_history_now = .true.

      ! see extras_check_model for information about custom termination codes
      ! by default, indicate where (in the code) MESA terminated
      if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step

      call saving_routine(id, save_now, extras_finish_step)
      if (save_now) then
         s% need_to_save_profiles_now = .true.
      end if
   end function extras_finish_step


   subroutine extras_after_evolve(id, ierr)
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
   end subroutine extras_after_evolve

   subroutine calc_JANu(id)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s

      integer :: nz, k
      real(dp) :: J, S_red, A, Nu

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      nz = s% nz
      !$OMP DO
      do k = 1, nz
         J = 1d0 - 4d0 * pi / 3d0 * s% rho_face(k) * (s% r(k)**3 / s%m(k))
         S_red = J * (sqrt(2d0) * s% csound_face(k) / s% r(k))

         A = s% brunt_N2(k) * s% r(k) / s% grav(k) / J
         Nu = 2d0 * J * s% grav(k) / (S_red**2 * s% r(k))

         s% xtra1_array(k) = J
         s% xtra2_array(k) = A
         s% xtra3_array(k) = Nu
      end do
      !$OMP END DO

   end subroutine calc_JANu


   subroutine calc_PQ(id, sigma)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s

      real(dp), intent(in) :: sigma
      integer :: nz, k, i, l, nz_skip, k_P, k_Q, k_l2, k_u2
      integer, parameter :: npass = 10, extra_pts = 40
      real(dp) :: J, S_red, P, Q, x1, x2, x3, y1, y2, y3, a, b, c, new_y
      real(dp) :: P_km1, Q_km1, w1, w2

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      nz = s% nz
      !$OMP DO
      do k = 1, nz
         J = s% xtra1_array(k)
         S_red = J * (sqrt(2d0) * s% csound_face(k) / s% r(k))
         P = 2d0 * J * (1d0 - (sigma / S_red)**2)
         Q = J * (1d0 - s% brunt_N2(k) / (sigma * J)**2)

         s% xtra4_array(k) = P
         s% xtra5_array(k) = Q
      end do
      !$OMP END DO

      ! only smooth in evan zone
      call get_evanescent_zones(id, sigma, k_P, k_Q, k_u2, k_l2)
      if ((k_P /= 0) .and. (k_Q /= 0)) then
         nz = max(k_P, k_Q, k_u2, k_l2)
         if ((k_u2 /= 0) .and. (k_l2 /= 0)) then
            nz_skip = min(k_P, k_Q, k_u2, k_l2)
         else
            nz_skip = min(k_P, k_Q)
         endif
         nz = min(s% nz - extra_pts, nz + extra_pts)
         nz_skip = max(1, nz_skip - extra_pts)
         do k = nz_skip, nz
            if (s% brunt_N2(k) > 0d0) then
               nz_skip = k + 5  ! Don't smooth conv boundary
               exit
            end if
         end do

         !first a moving avg pass
         P_km1 = s% xtra4_array(nz_skip - 1)
         Q_km1 = s% xtra5_array(nz_skip - 1)
         P = s% xtra4_array(nz_skip)
         Q = s% xtra5_array(nz_skip)
         do k = nz_skip, nz
            x1 = s% lnr(k - 1)
            x2 = s% lnr(k)
            x3 = s% lnr(k + 1)

            w1 = (abs((x1 - x2) / (x1 - x3)))
            w2 = (abs((x2 - x3) / (x1 - x3)))

            P = s%xtra4_array(k)
            s% xtra4_array(k) = (w2 * P_km1 + P + w1 * s% xtra4_array(k + 1)) / (1 + w1 + w2)
            P_km1 = P
            Q = s%xtra5_array(k)
            s% xtra5_array(k) = (w2 * Q_km1 + Q + w1 * s% xtra5_array(k + 1)) / (1 + w1 + w2)
            Q_km1 = Q
         end do

         ! Smooth P & Q
         do l = 1, 2
            do i = 1, npass
               do k = nz_skip, nz
                  ! use a bigger window for second half of passes to smooth gradient better
                  !                     x1 = s% lnr(k-2 - i/(npass/2))
                  ! swap back and forth between 2 points ahead or behind target point
                  x1 = s% lnr(k - 2 + mod(i - 1, 2) * 4)
                  x2 = s% lnr(k - 1)
                  x3 = s% lnr(k + 1 + i / (npass / 2))
                  if (l == 1) then
                     y1 = s% xtra4_array(k - 2 + mod(i - 1, 2) * 4)
                     y2 = s% xtra4_array(k - 1)
                     y3 = s% xtra4_array(k + 1 + i / (npass / 2))
                  elseif (l == 2) then
                     y1 = s% xtra5_array(k - 2 + mod(i - 1, 2) * 4)
                     y2 = s% xtra5_array(k - 1)
                     y3 = s% xtra5_array(k + 1 + i / (npass / 2))
                  end if
                  call quadr_fit_3pt(x1, x2, x3, y1, y2, y3, c, b, a)
                  new_y = a * s% lnr(k)**2 + b * s% lnr(k) + c
                  if ((new_y > min(y1, y2, y3)) .and. (new_y < max(y1, y2, y3))) then
                     if (l == 1) then
                        s% xtra4_array(k) = new_y
                     else if (l == 2) then
                        s% xtra5_array(k) = min(1d0, new_y)
                     end if
                  end if
               end do
            end do
         end do
      end if

   end subroutine calc_PQ


   subroutine get_evanescent_zones(id, sigma, k_P, k_Q, k_u2, k_l2)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s

      real(dp), intent(in) :: sigma
      integer, intent(out) :: k_P, k_Q
      integer, intent(out) :: k_u2, k_l2

      integer :: nz, k, k_spike_lo, k_spike_hi, nz_skip

      real(dp) :: P, alfa, beta, J, Nred2, Sred2
      real(dp) :: r_1b, m_1b, r_2b, m_2b
      integer k_s, k_e, max_eps_h_k

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      max_eps_h_k = maxloc(s% eps_nuc_categories(ipp, 1:s% nz) + s% eps_nuc_categories(icno, 1:s% nz), 1)

      k_P = 0
      k_Q = 0
      k_u2 = 0
      k_l2 = 0
      k_spike_lo = 0
      k_spike_hi = 0
      nz = s% nz
      nz_skip = 0

      ! Skip outer layers
      do k = 1, nz
         if (0.9d0 * s% r(1) > s% r(k)) then
            nz_skip = k
            exit
         end if
      end do

      do k = nz_skip, nz
         if (k_P == 0) then
            if (s% xtra4_array(k) >= 0 .and. s% xtra4_array(k - 1) < 0) then
               k_P = k  ! zero is between k and k-1
               exit
            end if
         end if
      end do

      do k = nz_skip, max_eps_h_k
         if (k_spike_lo == 0) then
            if (s% xtra5_array(k) >= 0 .and. s% xtra5_array(k - 1) < 0) then
               k_spike_lo = k
               exit
            end if
         end if
      end do

      ! Search upward from spike
      if (k_spike_lo /= 0) then
         do k = max(nz_skip, k_spike_lo), 2, -1
            if (k_spike_hi == 0) then
               if (s% xtra5_array(k) <= 0 .and. s% xtra5_array(k - 1) > 0) then
                  k_spike_hi = k
                  exit
               end if
            end if
         end do
      end if

      ! Search downward for 'normal' edge from spike if there or from nz_skip
      do k = max(nz_skip, k_spike_lo), max_eps_h_k
         if (k_Q == 0) then
            if (s% xtra5_array(k) <= 0 .and. s% xtra5_array(k - 1) > 0) then
               k_Q = k
               exit
            end if
         end if
      end do

      if (k_spike_lo < min(k_P, k_Q)) then
         ! Spike above 'normal' evanescent zone or no spike
         k_P = k_P
         k_Q = k_Q
         k_u2 = k_spike_hi
         k_l2 = k_spike_lo
      else if (k_spike_lo > min(k_P, k_Q)) then
         ! Spike inside 'normal' evanescent zone
         k = k_Q  ! temporarily hold k_Q as values are switched around
         k_P = k_P
         k_Q = k_spike_hi
         k_u2 = k_spike_lo
         k_l2 = k  ! = k_Q
      end if

      ! Flash evanescent zone
      if (s% power_he_burn >= 1) then
         ! start from h-shell burning
         do k = max_eps_h_k, s% nz - 1  ! ignore central cell
            if ((s% xtra5_array(k) >= 0) .and. (s% xtra5_array(k - 1) < 0) .and. (k_u2 == 0)) then
               k_u2 = k
            end if

            if ((s% xtra5_array(k) <= 0) .and. (s% xtra5_array(k - 1) > 0) .and. (k_l2 == 0)) then
               k_l2 = k
            end if

            if ((k_u2 /= 0) .and. (k_l2 /= 0)) then
               call calc_evan_radii(id, k_u2, k_l2, .false., r_1b, r_2b, m_1b, m_2b)
               if (abs(m_1b - m_2b) < 2d-3) then
                  ! Zone too small
                  k_u2 = 0
                  k_l2 = 0
               else
                  exit
               end if
            end if
         end do
      end if

      ! if too close to center can't calculate gradients correctly
      if (k_l2 >= s% nz - 10 ) then
         k_l2 = 0
      end if

   end subroutine get_evanescent_zones


   subroutine calc_evan_radii(id, k_P, k_Q, main_evn_zone, r_1, r_2, m_1, m_2)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s

      integer, intent(in) :: k_P, k_Q
      logical, intent(in) :: main_evn_zone
      real(dp), intent(out) :: r_1, r_2, m_1, m_2
      real(dp) :: alfa, beta

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      r_1 = -1d99
      r_2 = -1d99
      m_1 = -1d99
      m_2 = -1d99

      if (k_P > 1) then
         if (main_evn_zone) then
            alfa = s% xtra4_array(k_P) / (s% xtra4_array(k_P) - s% xtra4_array(k_P - 1))
         else
            alfa = s% xtra5_array(k_P) / (s% xtra5_array(k_P) - s% xtra5_array(k_P - 1))
         end if
         beta = 1d0 - alfa
         ! r1 and r2 are reversed compared to Pincon 2019 (using Takata r1 and r2 here).
         r_1 = (beta * s% r(k_P)**3 + alfa * s% r(k_P - 1)**3)**(1d0/3d0)
         m_1 = beta * s% m(k_P) + alfa * s% m(k_P - 1)
      end if

      if (k_Q > 1) then
         alfa = s% xtra5_array(k_Q) / (s% xtra5_array(k_Q) - s% xtra5_array(k_Q - 1))
         beta = 1d0 - alfa
         ! r1 and r2 are reversed compared to Pincon 2019 (using Takata r1 and r2 here).
         r_2 = (beta * s% r(k_Q)**3 + alfa * s% r(k_Q - 1)**3)**(1d0/3d0)
         m_2 = beta * s% m(k_Q) + alfa * s% m(k_Q - 1)
      end if

      m_1 = m_1 / msun
      m_2 = m_2 / msun

   end subroutine calc_evan_radii


   function quadr_PQ(s, a_int)
      real(dp), intent(in) :: s
      real(dp), intent(in) :: a_int(2, 4)
      real(dp) :: quadr_PQ

      quadr_PQ = (s**2 * a_int(1, 1) + s * a_int(1, 2) + a_int(1, 3)) * (s**2 * a_int(2, 1) + s * a_int(2, 2) + a_int(2, 3))
      quadr_PQ = sqrt(quadr_PQ)
   end function quadr_PQ


   function cubic_PQ(s, a_int)
      real(dp), intent(in) :: s
      real(dp), intent(in) :: a_int(2, 4)
      real(dp) :: cubic_PQ
      real(dp) :: P_part, Q_part
      integer :: k

      P_part = a_int(1, 1)
      Q_part = a_int(2, 1)
      do k = 2, 4
         P_part = P_part + s**(k - 1) * a_int(1, k)
         Q_part = Q_part + s**(k - 1) * a_int(2, k)
      end do
      cubic_PQ = sqrt(P_part * Q_part)
   end function cubic_PQ


   subroutine simpne (ntab, x, y, result)

      !*****************************************************************************80
      !
      !! SIMPNE approximates the integral of unevenly spaced data.
      !
      !  Discussion:
      !
      !    The routine repeatedly interpolates a 3-point Lagrangian polynomial
      !    to the data and integrates that exactly.
      !
      !  Modified:
      !
      !    10 February 2006
      !
      !  Reference:
      !
      !    Philip Davis, Philip Rabinowitz,
      !    Methods of Numerical Integration,
      !    Second Edition,
      !    Dover, 2007,
      !    ISBN: 0486453391,
      !    LC: QA299.3.D28.
      !
      !  Parameters:
      !
      !    Input, integer ( kind = 4 ) NTAB, number of data points.
      !    NTAB must be at least 3.
      !
      !    Input, real ( kind = 8 ) X(NTAB), contains the X values of the data,
      !    in order.
      !
      !    Input, real ( kind = 8 ) Y(NTAB), contains the Y values of the data.
      !
      !    Output, real ( kind = 8 ) RESULT.
      !    RESULT is the approximate value of the integral.
      !
      implicit none

      integer (kind = 4) ntab

      real(dp) :: del(3)
      real(dp) :: e
      real(dp) :: f
      real(dp) :: feints
      real(dp) :: g(3)
      integer (kind = 4) i
      integer (kind = 4) n
      real(dp) :: pi(3)
      real(dp) :: result
      real(dp) :: sum1
      real(dp) :: x(ntab)
      real(dp) :: x1
      real(dp) :: x2
      real(dp) :: x3
      real(dp) :: y(ntab)

      result = 0.0d0

      if (ntab <= 2) then
         write (*, '(a)') ' '
         write (*, '(a)') 'SIMPNE - Fatal error!'
         write (*, '(a)') '  NTAB <= 2.'
         stop 1
      end if

      n = 1

      do

         x1 = x(n)
         x2 = x(n + 1)
         x3 = x(n + 2)
         e = x3 * x3 - x1 * x1
         f = x3 * x3 * x3 - x1 * x1 * x1
         feints = x3 - x1

         del(1) = x3 - x2
         del(2) = x1 - x3
         del(3) = x2 - x1

         g(1) = x2 + x3
         g(2) = x1 + x3
         g(3) = x1 + x2

         pi(1) = x2 * x3
         pi(2) = x1 * x3
         pi(3) = x1 * x2

         sum1 = 0.0D+00
         do i = 1, 3
            sum1 = sum1 + y(n - 1 + i) * del(i) &
                  * (f / 3.0D+00 - g(i) * 0.5D+00 * e + pi(i) * feints)
         end do
         result = result - sum1 / (del(1) * del(2) * del(3))

         n = n + 2

         if (ntab <= n + 1) then
            exit
         end if

      end do

      if (mod (ntab, 2) /= 0) then
         return
      end if

      n = ntab - 2
      x3 = x(ntab)
      x2 = x(ntab - 1)
      x1 = x(ntab - 2)
      e = x3 * x3 - x2 * x2
      f = x3 * x3 * x3 - x2 * x2 * x2
      feints = x3 - x2

      del(1) = x3 - x2
      del(2) = x1 - x3
      del(3) = x2 - x1

      g(1) = x2 + x3
      g(2) = x1 + x3
      g(3) = x1 + x2

      pi(1) = x2 * x3
      pi(2) = x1 * x3
      pi(3) = x1 * x2

      sum1 = 0.0D+00
      do i = 1, 3
         sum1 = sum1 + y(n - 1 + i) * del(i) * &
               (f / 3.0D+00 - g(i) * 0.5D+00 * e + pi(i) * feints)
      end do

      result = result - sum1 / (del(1) * del(2) * del(3))

      return
   end


   subroutine calc_PQ_integral_v2(id, r_1, r_2, k_u, k_l, PQ_integral, a_int)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s

      real(dp), intent(in) :: r_1, r_2
      integer, intent(in) :: k_u, k_l
      integer :: n, k, i

      real(dp), intent(out) :: PQ_integral
      real(dp) :: P, Q, r, lnr0, s_0, s_j, h, integ_approx
      real(dp), intent(in) :: a_int(2, 4)
      real(dp), allocatable :: x(:), y(:)
      logical :: conv_close
      integer, parameter :: min_n = 6

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      s_0 = 0.5d0 * (log(r_1) - log(r_2))
      n = k_l - k_u + 2
      if (n >= min_n) then
         allocate(x(n))
         allocate(y(n))
         x(1) = -abs(s_0)
         x(n) = abs(s_0)
         y(1) = 0d0
         y(n) = 0d0
         lnr0 = log(r_1 * r_2) / 2d0
         do k = 2, n - 1
            i = k_l - (k - 1)
            x(k) = s% lnr(i) - lnr0
            P = s% xtra4_array(i)
            Q = s% xtra5_array(i)
            y(k) = sqrt(P * Q)
         end do
         call simpne(n, x, y, PQ_integral)
         deallocate(x)
         deallocate(y)
      else
         PQ_integral = -1d0
      end if

      n = 51
      allocate(x(n))
      allocate(y(n))

      x(1) = -abs(s_0)
      x(n) = abs(s_0)
      y(1) = 0d0
      y(n) = 0d0

      do k = 2, n - 1
         x(k) = -abs(s_0) * cos((k-1d0)/(n-1d0) * pi)  ! Chebyshev nodes of the second kind
         y(k) = cubic_PQ(x(k), a_int)
      end do
      call simpne(n, x, y, integ_approx)
      if (s% x_logical_ctrl(i_verbose_coupling)) then
         write(*, *) PQ_integral, integ_approx
      end if

      ! check if approx is close to actual values
      x(1:2) = 0d0
      y(1:5) = 0d0
      do k = k_u, k_l
         x(1) = s % lnr(k) - lnr0  ! s

         ! calc P and Q from polynomial
         y(1) = a_int(1, 1)  ! P
         y(2) = a_int(2, 1)  ! Q
         do i = 2, 4
            y(1) = y(1) + x(1)**(i - 1) * a_int(1, i)
            y(2) = y(2) + x(1)**(i - 1) * a_int(2, i)
         end do

         y(3) = y(3) + (y(1) / s% xtra4_array(k) - 1d0)**2
         y(4) = y(4) + (y(2) / s% xtra5_array(k) - 1d0)**2
         x(2) = x(2) + 1d0
      end do
      y(3) = sqrt(y(3)) / x(2)
      y(4) = sqrt(y(4)) / x(2)
      if (s% x_logical_ctrl(i_verbose_coupling)) then
         write(*, *) 'mean abs% difference from poly', y(3), y(4)
      end if

      ! try to always use numerical integral
      if (((s% brunt_N2(k_u) <= 0d0) .or. (s% brunt_N2(k_l) <= 0d0)) .and. (n >= min_n)) then  ! part convective
         PQ_integral = PQ_integral
      else if ((n >= min_n) .and. max(y(3), y(4)) > 5d-2) then  ! bad fit
         PQ_integral = PQ_integral
      else if (.not. (integ_approx == integ_approx)) then  ! check for nan
         PQ_integral = PQ_integral
      else
         PQ_integral = integ_approx
      end if
      deallocate(x)
      deallocate(y)
   end subroutine calc_PQ_integral_v2


   subroutine populate_weights(z, x, nd, m, c)
      !
      !  Input Parameters
      !    z            -  location where approximations are to be
      !                    accurate
      !    x(0:nd)      -  grid point locations, found in x(0:n)
      !    nd           -  dimension of x- and c-arrays in calling
      !                    program x(0:nd) and c(0:nd, 0:m), respectively
      !    m            -  highest derivative for which weights are
      !                    sought
      !
      !  Output Parameter
      !    c(0:nd,0:m)  -  weights at grid locations x(0:n) for
      !                    derivatives of order 0:m, found in c(0:nd, 0:m)
      !
      !  References:
      !      Generation of Finite Difference Formulas on Arbitrarily
      !          Spaced Grids, Bengt Fornberg,
      !          Mathematics of compuation, 51, 184, 1988, 699--706,
      !          doi: 10.1090/S0025-5718-1988-0935077-0
      !
      !      Classroom note: Calculation of weights in finite
      !          difference formulas, Bengt Fornberg,
      !          SIAM review, 40, 3, 1998, 685--691,
      !          doi: 10.1137/S0036144596322507
      ! https://github.com/bjodah/finitediff/

      real(dp), intent(in) :: z
      integer, intent(in) :: nd, m
      real(dp), intent(in) :: x(0:nd)
      real(dp), intent(out) :: c(0:nd, 0:m)

      real(dp) :: c1, c2, c3, c4, c5
      integer :: i, j, k, mn

      c1 = 1
      c4 = x(0) - z
      c = 0
      c(0, 0) = 1
      do i = 1, nd
         mn = min(i, m)
         c2 = 1
         c5 = c4
         c4 = x(i) - z
         do j = 0, i - 1
            c3 = x(i) - x(j)
            c2 = c2 * c3
            if (j == i - 1) then
               do k = mn, 1, -1
                  c(i, k) = c1 * (k * c(i - 1, k - 1) - c5 * c(i - 1, k)) / c2
               end do
               c(i, 0) = -c1 * c5 * c(i - 1, 0) / c2
            endif
            do k = mn, 1, -1
               c(j, k) = (c4 * c(j, k) - k * c(j, k - 1)) / c3
            end do
            c(j, 0) = c4 * c(j, 0) / c3
         end do
         c1 = c2
      end do
   end subroutine populate_weights


   subroutine quadr_fit_3pt(x1, x2, x3, y1, y2, y3, a0, a1, a2)
      real(dp), intent(in) :: x1, x2, x3, y1, y2, y3
      real(dp), intent(out) :: a0, a1, a2

      a2 = (x1 * (y3 - y2) + x2 * (y1 - y3) + x3 * (y2 - y1)) / ((x1 - x2) * (x1 - x3) * (x2 - x3))
      a1 = (y2 - y1) / (x2 - x1) - a2 * (x1 + x2)
      a0 = (y1 + y2 + y3 - a2 * (x1**2 + x2**2 + x3**2) - a1 * (x1 + x2 + x3)) / 3d0

   end subroutine quadr_fit_3pt

   ! Calculate polynomial coefficients for P and Q for the integral
   subroutine calc_quadr_coeffs_int(id, r_1, r_2, k_u, k_l, main_evn_zone, a_int)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s

      logical, intent(in) :: main_evn_zone
      integer, intent(in) :: k_u, k_l
      real(dp), intent(in) :: r_1, r_2
      real(dp), intent(out) :: a_int(2, 4)

      integer, parameter :: n_evan = 6
      integer :: i, j, k, k_mid, k_P, k_Q, k_star
      real(dp) :: r_0, lnr_0, s_0, alfa, beta, Q_0, Q_1, P_0, P_1
      real(dp) :: x1, x2, x3, y1, y2, y3
      real(dp) :: c(n_evan, 0:1)
      real(dp) :: s_arr(n_evan)

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      a_int(:, :) = 0d0

      r_0 = sqrt(r_1 * r_2)
      lnr_0 = log(r_0)
      s_0 = 0.5d0 * (log(r_1) - log(r_2))

      if (r_1 > r_2) then
         k_P = k_u
         k_Q = k_l
      else
         k_P = k_l
         k_Q = k_u
      end if

      if (main_evn_zone) then
         alfa = s% xtra4_array(k_P) / (s% xtra4_array(k_P) - s% xtra4_array(k_P - 1))
      else
         alfa = s% xtra5_array(k_P) / (s% xtra5_array(k_P) - s% xtra5_array(k_P - 1))
      end if
      beta = 1 - alfa
      Q_1 = beta * s% xtra5_array(k_P) + alfa * s% xtra5_array(k_P - 1)  ! Q at P=0

      alfa = s% xtra5_array(k_Q) / (s% xtra5_array(k_Q) - s% xtra5_array(k_Q - 1))
      beta = 1 - alfa
      P_1 = beta * s% xtra4_array(k_Q) + alfa * s% xtra4_array(k_Q - 1)  ! P at Q=0

      k_mid = 0
      do k = k_u, s% nz
         if (s% r(k) < r_0) then
            ! r_0 between k and k-1
            k_mid = k
            exit
         end if
      end do

      if (abs(k_P - k_Q) > 3) then
         do k = 1, size(s_arr)
            k_star = (k - 1) + k_Q - n_evan / 2
            s_arr(k) = s% lnr(k_star) - lnr_0
         end do
         call populate_weights(-s_0, s_arr, size(s_arr) - 1, 1, c)
         P_1 = 0d0
         do k = 1, n_evan
            k_star = (k - 1) + k_Q - n_evan / 2
            P_1 = P_1 + c(k, 0) * s% xtra4_array(k_star)
         end do

         do k = 1, size(s_arr)
            k_star = (k - 1) + k_P - n_evan / 2
            s_arr(k) = s% lnr(k_star) - lnr_0
         end do
         call populate_weights(s_0, s_arr, size(s_arr) - 1, 1, c)
         Q_1 = 0d0
         do k = 1, n_evan
            k_star = (k - 1) + k_P - n_evan / 2
            Q_1 = Q_1 + c(k, 0) * s% xtra5_array(k_star)
         end do

         do k = 1, size(s_arr)
            k_star = (k - 1) + k_mid - n_evan / 2
            s_arr(k) = s% lnr(k_star) - lnr_0
         end do
         call populate_weights(0d0, s_arr, size(s_arr) - 1, 1, c)
         P_0 = 0d0
         Q_0 = 0d0
         do k = 1, n_evan
            k_star = (k - 1) + k_mid - n_evan / 2
            P_0 = P_0 + c(k, 0) * s% xtra4_array(k_star)
            Q_0 = Q_0 + c(k, 0) * s% xtra5_array(k_star)
         end do
         a_int(1, 3) = (P_1 - 2 * P_0) / (2 * s_0**2)
         a_int(1, 2) = -P_1 / (2 * s_0)
         a_int(1, 1) = P_0

         a_int(2, 3) = (Q_1 - 2 * Q_0) / (2 * s_0**2)
         a_int(2, 2) = Q_1 / (2 * s_0)
         a_int(2, 1) = Q_0
      else
         x1 = s% lnr(k_P - 2) - lnr_0
         x2 = s_0
         x3 = s% lnr(k_P + 2) - lnr_0
         y1 = s% xtra4_array(k_P - 2)
         y2 = 0d0
         y3 = s% xtra4_array(k_P + 2)
         call quadr_fit_3pt(x1, x2, x3, y1, y2, y3, a_int(1, 1), a_int(1, 2), a_int(1, 3))

         x1 = s% lnr(k_Q - 2) - lnr_0
         x2 = -s_0
         x3 = s% lnr(k_Q + 2) - lnr_0
         y1 = s% xtra5_array(k_Q - 2)
         y2 = 0d0
         y3 = s% xtra5_array(k_Q + 2)
         call quadr_fit_3pt(x1, x2, x3, y1, y2, y3, a_int(2, 1), a_int(2, 2), a_int(2, 3))
      endif
   end subroutine calc_quadr_coeffs_int

   ! Calculate polynomial coefficients for P and Q for the gradient
   subroutine calc_quadr_coeffs_grd(id, r_1, r_2, k_u, k_l, a_grd, save_mesh)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s

      integer, parameter :: deg = 3, npoints = 11, nsmooth = 0, npre = -9, npost = 9
      integer, intent(in) :: k_u, k_l
      real(dp), intent(in) :: r_1, r_2
      logical, intent(in) :: save_mesh
      real(dp), intent(out) :: a_grd(2, 4)

      integer :: i, j, k, l, k_mid, i_arr(npre - 1:npoints + npost + 1), k_P, k_Q, nstart, nend, step
      real(dp) :: r_0, lnr_0, s_0, x1, x2, x3, x4, y1, y2, y3, y4, alfa, beta, a2, a1, a0
      real(dp) :: s11, s12, s22, sy1, sy2, gamma, avg_x, avg_y, w(npoints), w_smooth(-nsmooth:nsmooth)
      real(dp) :: x(npre:npoints + npost), y(npre:npoints + npost)
      real(dp) :: A(npoints, deg + 1), b(npoints, 1), work(max(64, 16 * (deg + 1)))
      integer :: m, n, nrhs, lda, ldb, lwork, info
      real(dp) :: c(npoints, 0:deg)

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      a_grd(:, :) = 0d0

      r_0 = sqrt(r_1 * r_2)
      lnr_0 = log(r_0)
      s_0 = (log(r_1) - log(r_2)) / 2

      k_mid = 0
      do k = k_u, s% nz
         if (s% r(k) < r_0) then
            ! r_0 between k and k-1
            k_mid = k
            exit
         end if
      end do

      if ((r_1 >= s% r(k_u)) .and. (r_1 < s% r(k_u - 1)))  then
         k_P = k_u
         k_Q = k_l
      else
         k_P = k_l
         k_Q = k_u
      end if
      if (s% x_logical_ctrl(i_verbose_coupling)) then
         write(*, *) 'model', s% model_number
      end if
      do j = 1, 2

         ! Avoid r_1 and r_2 as the couple points around them are very noisy due to approaching 0/0
         i = 0

         if (((r_1 > r_2) .and. (j == 2)) .or. ((r_1 < r_2) .and. (j == 1))) then
            nstart = npre - 1
            nend = npoints + npost + 1
            step = 1
         else
            nend = npre - 1
            nstart = npoints + npost + 1
            step = -1
         end if
         n = 2
         do k = nstart, nend, step
            i_arr(k) = k + k_mid - (npoints / 2) + i
            if (j == 1) then
               if ((i_arr(k) == k_P - n * step) .or. i_arr(k) == k_P - 1 - n * step) then
                  i = i + 2 * n * step
                  i_arr(k) = k_P
               end if
            else if (j == 2) then
               if ((i_arr(k) == k_Q - step) .or. i_arr(k) == k_Q - 1 - step) then
                  i = i + 2 * n * step
                  i_arr(k) = k_Q
               end if
            end if
         end do
         do k = npre, npoints + npost
            i = i_arr(k)

            x1 = s% lnr(i_arr(k - 1)) - lnr_0
            x2 = s% lnr(i_arr(k + 1)) - lnr_0

            if (j == 1) then  ! P
!               y1 = s% xtra4_array(i-1)
!               y2 = s% xtra4_array(i+1)
!               y1 = y1 + s% xtra4_array(i-1+l) * (x1 + s_0)
!               y2 = y2 + s% xtra4_array(i  +l) * (x2 + s_0)
               if (i_arr(k - 1) == k_P) then
                  ! use l'Hopital to get value at s=s_0
                  call populate_weights(s_0 + lnr_0, s% lnr(k_P - npoints / 2:k_P + npoints / 2), npoints - 1, 1, c)
                  y1 = 0d0
                  do n = 1, npoints
                     y1 = y1 + c(n, 1) * s% xtra4_array(k_P - npoints / 2 + (n - 1))
                  end do
                  y1 = y1 / (-1d0)
!                  y1 = (s% xtra4_array(k_P-1) - s% xtra4_array(k_P))/ &
!                         (s% lnr(k_P-1) - s% lnr(k_P)) / (-1d0)
                  x1 = s_0
               else
                  y1 = s% xtra4_array(i_arr(k - 1)) / (s_0 - x1)
               end if

               if (i_arr(k + 1) == k_P) then
                  ! use l'Hopital to get value at s=s_0
                  call populate_weights(s_0 + lnr_0, s% lnr(k_P - npoints / 2:k_P + npoints / 2), npoints - 1, 1, c)
                  y2 = 0d0
                  do n = 1, npoints
                     y2 = y2 + c(n, 1) * s% xtra4_array(k_P - npoints / 2 + (n - 1))
                  end do
                  y2 = y2 / (-1d0)
!                  y2 = (s% xtra4_array(k_P-1) - s% xtra4_array(k_P))/ &
!                         (s% lnr(k_P-1) - s% lnr(k_P)) / (-1d0)
                  x2 = s_0
               else
                  y2 = s% xtra4_array(i_arr(k + 1)) / (s_0 - x2)
               end if
            elseif (j == 2) then  ! Q
!               y1 = s% xtra5_array(i-1)
!               y2 = s% xtra5_array(i+1)
!               y1 = y1 + s% xtra5_array(i-1+l) * (s_0 - x1)
!               y2 = y2 + s% xtra5_array(i  +l) * (s_0 - x2)
               y1 = s% xtra5_array(i_arr(k - 1)) / (x1 + s_0)
               y2 = s% xtra5_array(i_arr(k + 1)) / (x2 + s_0)
               if (i_arr(k - 1) == k_Q) then
                  ! use l'Hopital to get value at s=-s_0
                  call populate_weights(-s_0 + lnr_0, s% lnr(k_Q - npoints / 2:k_Q + npoints / 2), npoints - 1, 1, c)
                  y1 = 0d0
                  do n = 1, npoints
                     y1 = y1 + c(n, 1) * s% xtra5_array(k_Q - npoints / 2 + (n - 1))
                  end do
                  y1 = y1 / (1d0)
!                  y1 = (s% xtra5_array(k_Q-1) - s% xtra5_array(k_Q+1))/ &
!                        (s% lnr(k_Q-1) - s% lnr(k_Q+1)) / 1d0
                  x1 = -s_0
               else
                  y1 = s% xtra5_array(i_arr(k - 1)) / (x1 + s_0)
               end if

               if (i_arr(k + 1) == k_Q) then
                  ! use l'Hopital to get value at s=-s_0
                  call populate_weights(-s_0 + lnr_0, s% lnr(k_Q - npoints / 2:k_Q + npoints / 2), npoints - 1, 1, c)
                  y2 = 0d0
                  do n = 1, npoints
                     y2 = y2 + c(n, 1) * s% xtra5_array(k_Q - npoints / 2 + (n - 1))
                  end do
                  y2 = y2 / (1d0)
!                  y2 = (s% xtra5_array(k_Q-1) - s% xtra5_array(k_Q+1))/ &
!                         (s% lnr(k_Q-1) - s% lnr(k_Q)) / 1d0
                  x2 = -s_0
               else
                  y2 = s% xtra5_array(i_arr(k + 1)) / (x2 + s_0)
               end if
            endif
            x(k) = (x1 + x2) / 2d0
            y(k) = (y1 - y2) / (x1 - x2)
            ! Fit to P and Q
!            x(k) = x2
!            y(k) = y2
         end do
         if (save_mesh) then
            open(2, FILE = 'LOGS/unsmoothed.dat', status = 'unknown', position = 'append')
            write(2, *) s% model_number, j, x, y
            close(2)
         end if
         do i = 1, 8
            exit  ! dont do extra smoothing
            ! swap smoothing direction for P and Q to smooth towards the 0/0 boundary
            if (i == 1) then
               nstart = npre + 2
               nend = npoints + npost - 1
               step = 1
            else
               nend = npre + 2
               nstart = npoints + npost - 1
               step = -1
            end if
            do k = npre + 2, npoints + npost - 1, step
               x1 = x(k - 2)
               x2 = x(k - 1)
               x3 = x(k + 1)

               y1 = y(k - 2)
               y2 = y(k - 1)
               y3 = y(k + 1)

               call quadr_fit_3pt(x1, x2, x3, y1, y2, y3, a0, a1, a2)
               y4 = a2 * x(k)**2 + a1 * x(k) + a0
               y(k) = y4
            end do
         end do
         if (save_mesh) then
            open(2, FILE = 'LOGS/smoothed.dat', status = 'unknown', position = 'append')
            write(2, *) s% model_number, j, x, y
            close(2)
         end if
         ! remesh
         if (.true.) then
            x1 = -max(abs(s_0), 0.05d0) / npoints  ! new mesh between 0.5s0 and -0.5s0
            ! Force mesh to be within available points and not near
            x2 = abs(x1 * npoints / 2)
            if (x2 > min(abs(x(npre + 2)), abs(x(npoints + npost - 1)))) then
               x1 = -min(abs(x(npre + 2)), abs(x(npoints + npost - 1))) / (npoints / 2)
            end if

            do k = 1, npoints
               A(k, 1) = x1 * (k - 1 - npoints / 2) * exp(-(npoints / 2 - abs(k - 1 - npoints / 2)) / (0.5d0 * npoints / 2))
            end do

            l = npre + 2
            do k = 1, npoints
               x1 = A(k, 1)
               A(k, 2) = -1d99
               do i = l, npoints + npost - 1
                  if ((x1 >= x(i)) .and. (x1 < x(i - 1))) then  ! point between i and i-1
                     x2 = x(i)
                     x3 = x(i - 1)
                     x4 = x(i - 2)
                     y2 = y(i)
                     y3 = y(i - 1)
                     y4 = y(i - 2)
                     call quadr_fit_3pt(x2, x3, x4, y2, y3, y4, a0, a1, a2)
                     A(k, 2) = a2 * x1**2 + a1 * x1 + a0

                     x4 = x(i + 1)
                     y4 = y(i + 1)
                     call quadr_fit_3pt(x2, x3, x4, y2, y3, y4, a0, a1, a2)
                     A(k, 2) = A(k, 2) + a2 * x1**2 + a1 * x1 + a0
                     A(k, 2) = A(k, 2) / 2d0
                     l = i
                     exit
                  end if
               end do
            end do
            x(1:npoints) = A(1:npoints, 1)
            y(1:npoints) = A(1:npoints, 2)
         end if
         if (save_mesh) then
            open(2, FILE = 'LOGS/meshed.dat', status = 'unknown', position = 'append')
            write(2, *) s% model_number, j, x(1:npoints), y(1:npoints)
            close(2)
         end if

         ! Set up vandermonde matrix
         do k = 1, npoints
            A(k, 1) = 1d0
            do i = 1, deg
               A(k, i + 1) = x(k) ** i
            end do
         end do
         b(:, 1) = y(1:npoints)
         m = npoints
         n = deg + 1
         nrhs = 1
         lda = m
         ldb = m
         lwork = 64
         call dgels('N', m, n, nrhs, A, lda, b, ldb, work, lwork, info)
         if (info /= 0) stop 'lapack error'

         a_grd(j, 1:deg + 1) = b(1:deg + 1, 1)
         if (s% x_logical_ctrl(i_verbose_coupling)) then
            write(*, *) 'dgels'
            write(*, *) b(1:deg + 1, 1)

            write(*, *) k_u, k_mid, k_l
            write(*, *) s% model_number
         end if
      end do

   end subroutine calc_quadr_coeffs_grd


   subroutine calc_G2_div_2kappa_s0_v2(id, r_1, r_2, k_u, k_l, G2_div_2kappa_s0, kappa_s0, NuAJ_s0, &
         dlnPds_s0, dlnQds_s0, dlnc_ds_s0_part_grd, dlnc_ds_s0_part_Sl, dlnc_ds_s0_part_N, a_grd, a_int)

      integer, intent(in) :: id
      integer :: ierr

      type (star_info), pointer :: s

      integer, intent(in) :: k_u, k_l
      real(dp), intent(in) :: r_1, r_2
      real(dp), intent(out) :: G2_div_2kappa_s0, kappa_s0, NuAJ_s0, dlnPds_s0, dlnQds_s0
      real(dp), intent(out) :: dlnc_ds_s0_part_grd, dlnc_ds_s0_part_Sl, dlnc_ds_s0_part_N
      integer :: i, n_evan, k, k_mid, k_star, k_interp, k_P, k_Q
      real(dp) :: s_0, r_0, ss_s0, w, total_w
      real(dp) :: dlnc_ds_s0_part, dlncds_s0, dssds_s0, lnP_s0, lnQ_s0, lnPQ_s0, dlnPQds_s0
      real(dp) :: dPds_s0, dQds_s0, P_s0, Q_s0, alfa, beta
      real(dp), dimension(15) :: x, y
      real(dp) :: sx, sy, sxx, syy, sxy, y1, y2
      real(dp), intent(in) :: a_grd(2, 4), a_int(2, 4)

      logical :: use_approx
      ! For appyling fd
      real(dp) :: out(2)  ! Only need 0 and/or 1st order derivs/interps
      integer :: nin, j
      real(dp), dimension(:), allocatable :: s_arr
      real(dp), dimension(:, :), allocatable :: c

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      r_0 = sqrt(r_1 * r_2)
      s_0 = 0.5d0 * (log(r_1) - log(r_2))

      k_mid = 0
      do k = k_u, s% nz
         if (s% r(k) < r_0) then
            ! r_0 between k and k-1
            k_mid = k
            exit
         end if
      end do

      n_evan = max(4, min(k_l - k_u, 10))
      n_evan = 11

      P_s0 = 0d0
      Q_s0 = 0d0
      dPds_s0 = 0d0
      dQds_s0 = 0d0
      NuAJ_s0 = 0d0
      kappa_s0 = 0d0
      allocate(s_arr(n_evan))
      allocate(c(n_evan, 0:1))

      do k = 1, n_evan
         k_star = (k - 1) + k_mid - n_evan / 2
         s_arr(k) = log(s% r(k_star) / r_0)
      end do

      ! populate_weights(xtgt, xdata, nin-1, maxorder, c)
      call populate_weights(0d0, s_arr, n_evan - 1, 1, c)
      do k = 1, n_evan
!         k_star = (k - 1) + max(k_u, (k_mid - 1) - 5)
         k_star = (k - 1) + k_mid - n_evan / 2
!         P_s0 = P_s0 + c(k, 0)*s% xtra4_array(k_star)
!         Q_s0 = Q_s0 + c(k, 0)*s% xtra5_array(k_star)
!         dPds_s0 = dPds_s0 + c(k, 1)*s% xtra4_array(k_star)
!         dQds_s0 = dQds_s0 + c(k, 1)*s% xtra5_array(k_star)

         NuAJ_s0 = NuAJ_s0 + c(k, 0) * (s% xtra3_array(k_star) - s% xtra2_array(k_star) - s% xtra1_array(k_star))
!        kappa_s0 = kappa_s0 + c(k, 0)*(sqrt(s% xtra4_array(k_star) * s% xtra5_array(k_star) / (s_0**2 - s_arr(k)**2)))
      end do
      deallocate(c)
      deallocate(s_arr)

!     dlnPds_s0 = dPds_s0 / P_s0
!     dlnQds_s0 = dQds_s0 / Q_s0
!     if (s% x_logical_ctrl(i_verbose_coupling)) then
!        write(*,*) P_s0, a_int(1,1)
!        write(*,*) Q_s0, a_int(2,1)
!     end if
      ! if fitting to dP and dQ in gradient fit
      ! appears to be a little better than interpolation (for P_s0 and Q_s0)
!      dssds_s0 = 2d0/s_0  ! d/ds ln (s+s0)/(s0-s) at s=0
!      P_s0 = a_int(1,1)
!      Q_s0 = a_int(2,1)
!      dlnPds_s0 = a_grd(1, 1)/P_s0
!      dlnQds_s0 = a_grd(2, 1)/Q_s0

      ! if fitting to P and Q direnctly in gradient fit
      !         dssds_s0 = 2d0/s_0  ! d/ds ln (s+s0)/(s0-s) at s=0
      !         P_s0 = a_grd(1,1)
      !         Q_s0 = a_grd(2,1)
      !         dlnPds_s0 = a_grd(1, 2)/P_s0
      !         dlnQds_s0 = a_grd(2, 2)/Q_s0

      ! if fitting P/(s0-s) and Q/(s+s0)
      !         dssds_s0 = 0d0 ! incorporated in dlnPds and dlnQds
      !         dlnPds_s0 = a_grd(1, 2)/a_grd(1, 1)
      !         dlnQds_s0 = a_grd(2, 2)/a_grd(2, 1)

      ! if fitting to d/ds P/(s0-s) and d/ds Q/(s+s0)
      dssds_s0 = 0d0 ! incorporated in dlnPds and dlnQds
      P_s0 = a_int(1, 1)
      Q_s0 = a_int(2, 1)
      dlnPds_s0 = a_grd(1, 1) / (P_s0 / s_0)
      dlnQds_s0 = a_grd(2, 1) / (Q_s0 / s_0)

      kappa_s0 = sqrt(P_s0 * Q_s0) / abs(s_0)
      dlnc_ds_s0_part_grd = 0.25d0 * (dlnPds_s0 - dlnQds_s0 + dssds_s0)
      dlncds_s0 = dlnc_ds_s0_part_grd - 0.5d0 * NuAJ_s0
      G2_div_2kappa_s0 = dlncds_s0**2 / (2d0 * kappa_s0)

      call calc_dlnc_ds_s0_part_ap(id, r_1, r_2, k_mid, dlnc_ds_s0_part_Sl, dlnc_ds_s0_part_N)
   end subroutine calc_G2_div_2kappa_s0_v2


   subroutine calc_G2_div_2kappa_s0(id, r_1, r_2, k_u, k_l, G2_div_2kappa_s0, kappa_s0, NuAJ_s0, &
         dlnPds_s0, dlnQds_s0, dlnc_ds_s0_part_ip, dlnc_ds_s0_part_Sl, dlnc_ds_s0_part_N, a_int)

      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s

      integer, intent(in) :: k_u, k_l
      real(dp), intent(in) :: r_1, r_2
      real(dp), intent(in) :: a_int(2, 4)
      real(dp), intent(out) :: G2_div_2kappa_s0, kappa_s0, NuAJ_s0, dlnPds_s0, dlnQds_s0
      real(dp), intent(out) :: dlnc_ds_s0_part_ip, dlnc_ds_s0_part_Sl, dlnc_ds_s0_part_N
      integer :: n_evan, k, k_mid, k_star, k_interp
      real(dp) :: s_0, r_0, ss_s0, yP, yQ, yP_s0, yQ_s0, P_s0, Q_s0
      real(dp) :: yPQ_s0, dyPds_s0, dyQds_s0, dyPQds_s0, dlnyPQds_s0
      real(dp), dimension(:), allocatable :: s_arr
      real(dp) :: dlnc_ds_s0_part, dlncds_s0, dssds_s0
      logical :: use_approx

      ! For appyling fd
      real(dp), dimension(:, :), allocatable :: c

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      r_0 = sqrt(r_1 * r_2)
      s_0 = 0.5d0 * (log(r_1) - log(r_2))

      n_evan = 11

      if (.false.) then
         dlnPds_s0 = 0d0
         dlnPds_s0 = 0d0 / dlnPds_s0
         dlnQds_s0 = 0d0 / dlnPds_s0
         NuAJ_s0 = 0d0 / dlnPds_s0
         kappa_s0 = 0d0 / dlnPds_s0
         dlncds_s0 = 0d0 / dlnPds_s0
         G2_div_2kappa_s0 = 0d0 / dlnPds_s0
         dlnc_ds_s0_part_Sl = 0d0 / dlnPds_s0
         dlnc_ds_s0_part_ip = 0d0 / dlnPds_s0
         return
      end if

      k_mid = 0
      do k = min(k_u, k_l) - 1, max(k_l, k_u) + 1
         if (s% r(k) < r_0) then
            ! r_0 between k and k-1
            k_mid = k
            exit
         end if
      end do

      allocate(s_arr(n_evan))
      allocate(c(n_evan, 0:1))

      do k = 1, n_evan
         k_star = (k - 1) + max(k_u, (k_mid - 1) - 5)
         s_arr(k) = log(s% r(k_star) / r_0)
      end do

      ! populate_weights(xtgt, xdata, nin-1, maxorder, c)
      call populate_weights(0d0, s_arr, n_evan - 1, 1, c)

      P_s0 = 0d0
      Q_s0 = 0d0
      yP_s0 = 0d0
      yQ_s0 = 0d0
      yPQ_s0 = 0d0
      dyPds_s0 = 0d0
      dyQds_s0 = 0d0
      dyPQds_s0 = 0d0
      NuAJ_s0 = 0d0
      kappa_s0 = 0d0

      do k = 1, n_evan
         !             k_star = (k - 1) + max(k_u, (k_mid - 1) - 5)
         k_star = (k - 1 - n_evan / 2) + k_mid - 1

         !             yP = s% xtra4_array(k_star)
         !             yQ = s% xtra5_array(k_star)
         !             yP = log(s% xtra4_array(k_star))
         !             yQ = log(s% xtra5_array(k_star))

         yP = s% xtra4_array(k_star) / (s_0 - s_arr(k))
         yQ = s% xtra5_array(k_star) / (s_arr(k) + s_0)
         !             yP = s% xtra4_array(k_star) * (s_arr(k) + s_0)
         !             yQ = s% xtra5_array(k_star) * (s_0 - s_arr(k))

         yP_s0 = yP_s0 + c(k, 0) * yP
         yQ_s0 = yQ_s0 + c(k, 0) * yQ
         yPQ_s0 = yPQ_s0 + c(k, 0) * (yP / yQ)

         dyPds_s0 = dyPds_s0 + c(k, 1) * yP
         dyQds_s0 = dyQds_s0 + c(k, 1) * yQ
         dyPQds_s0 = dyPQds_s0 + c(k, 1) * (yP / yQ)

         NuAJ_s0 = NuAJ_s0 + c(k, 0) * (s% xtra3_array(k_star) - s% xtra2_array(k_star) - s% xtra1_array(k_star))

         !             ! for kappa_s0 if not using P/(s0-s) and Q/(s-s0)
         !             P_s0 = P_s0 + c(k, 0) * s% xtra4_array(k_star)
         !             Q_s0 = Q_s0 + c(k, 0) * s% xtra5_array(k_star)
      end do
      !         kappa_s0 = sqrt(yP_s0*yQ_s0)
      ! use P_s0 and Q_s0 from integral fit works much better than the fits tried above
      P_s0 = a_int(1, 1)
      Q_s0 = a_int(2, 1)
      kappa_s0 = sqrt(P_s0 * Q_s0) / abs(s_0)

      deallocate(c)
      deallocate(s_arr)

      dlnPds_s0 = dyPds_s0 / yP_s0
      dlnQds_s0 = dyQds_s0 / yQ_s0
      dlnyPQds_s0 = dyPQds_s0 / yPQ_s0
      if (yP == s% xtra4_array(k_star)) then
         dssds_s0 = 2d0 / s_0  ! d/ds ln (s+s0)/(s0-s) at s=0
      elseif (yP == log(s% xtra4_array(k_star))) then
         dssds_s0 = 2d0 / s_0  ! d/ds ln (s+s0)/(s0-s) at s=0
         dlnPds_s0 = dyPds_s0
         dlnQds_s0 = dyQds_s0
         dlnyPQds_s0 = dyPQds_s0
      else
         dssds_s0 = 0d0
      end if

      dlnc_ds_s0_part_ip = 0.25d0 * (dlnPds_s0 - dlnQds_s0 + dssds_s0)
      if (dlnc_ds_s0_part_ip /= dlnc_ds_s0_part_ip) then  ! isnan check
         dlnc_ds_s0_part_ip = 0.25d0 * (dlnyPQds_s0 + dssds_s0)
      end if
      if (s% x_logical_ctrl(i_verbose_coupling)) then
         write(*, *) 'dlncds_s0_part_ip = ', 0.25d0 * (dlnPds_s0 - dlnQds_s0 + dssds_s0), 0.25d0 * (dlnyPQds_s0 + dssds_s0)
      end if
      call calc_dlnc_ds_s0_part_ap(id, r_1, r_2, k_mid, dlnc_ds_s0_part_Sl, dlnc_ds_s0_part_N)

      !         use_approx = (s% power_he_burn <= 1d-9 .and. s% center_he4 >=0.95)
      use_approx = .false.
      if (use_approx) then
         dlnc_ds_s0_part = dlnc_ds_s0_part_Sl
      else
         dlnc_ds_s0_part = dlnc_ds_s0_part_ip
      end if

      dlncds_s0 = dlnc_ds_s0_part - 0.5d0 * NuAJ_s0
      G2_div_2kappa_s0 = dlncds_s0**2 / (2d0 * kappa_s0)

   end subroutine calc_G2_div_2kappa_s0

   subroutine calc_dlnc_ds_s0_part_ap(id, r_1, r_2, k_mid, dlnc_ds_s0_part_Sl, dlnc_ds_s0_part_N)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s

      integer, intent(in) :: k_mid
      real(dp), intent(in) :: r_1, r_2
      real(dp), intent(out) :: dlnc_ds_s0_part_Sl, dlnc_ds_s0_part_N

      real(dp) :: Sred_km1, Sred_k, Nred_km1, Nred_k
      real(dp) :: beta_Sl, beta_N, r_0, J, alfa_Sl, alfa_N
      integer :: k, km1, kp1

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      km1 = k_mid - 1
      J = s% xtra1_array(km1)
      !         Sred_km1 = J * (sqrt(2d0)*s% csound_face(km1)/s% r(km1))
      !         Nred_km1 = sqrt(s% brunt_N2(km1)) / J
      ! ignoring factor sigma as it cancels later when doing powerlaw
      ! also using P and Q so we use smoothed values
      Sred_km1 = sqrt(2d0 * J / (2d0 * J - s% xtra4_array(km1)))
      Nred_km1 = sqrt(1d0 - s% xtra5_array(km1) / J)

      kp1 = k_mid + 1
      J = s% xtra1_array(kp1)
      !         Sred_k = J * (sqrt(2d0)*s% csound_face(kp1)/s% r(kp1))
      !         Nred_k = sqrt(s% brunt_N2(kp1)) / J
      Sred_k = sqrt(2d0 * J / (2d0 * J - s% xtra4_array(kp1)))
      Nred_k = sqrt(1d0 - s% xtra5_array(kp1) / J)

      beta_Sl = - (log(Sred_km1) - log(Sred_k)) / (s% lnr(km1) - s% lnr(kp1))
      beta_N = - (log(Nred_km1) - log(Nred_k)) / (s% lnr(km1) - s% lnr(kp1))
      alfa_Sl = (r_2 / r_1) ** beta_Sl
      alfa_N = (r_2 / r_1) ** beta_N

      dlnc_ds_s0_part_Sl = -beta_Sl * (alfa_Sl / (1d0 - alfa_Sl) + 1d0 / log(alfa_Sl))  ! first term in eq D.10 in Pincon 2019
      dlnc_ds_s0_part_N = -beta_N * (alfa_N / (1d0 - alfa_N) + 1d0 / log(alfa_N))  ! first term in eq D.10 in Pincon 2019
   end subroutine calc_dlnc_ds_s0_part_ap


   subroutine do_strong_coupling(id, nu_q, out, out_int)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s

      real(dp), intent(in) :: nu_q
      real(dp) :: sigma

      integer :: k_P, k_Q, k_u, k_l, k_u2, k_l2, k
      real(dp) :: r_1, r_2, m_1, m_2
      real(dp) :: r_1b, r_2b, m_1b, m_2b
      real(dp) :: PQ_integral, PQ_integral_b

      real(dp) :: G2_div_2kappa_s0, kappa_s0, NuAJ_s0, dlnPds_s0, dlnQds_s0
      real(dp) :: dlnc_ds_s0_part_grd, dlnc_ds_s0_part_ip, dlnc_ds_s0_part_Sl, dlnc_ds_s0_part_N
      real(dp) :: G2_div_2kappa_s0b, kappa_s0b, NuAJ_s0b, dlnPds_s0b, dlnQds_s0b
      real(dp) :: dlnc_ds_s0_part_grd_b, dlnc_ds_s0_part_ip_b, dlnc_ds_s0_part_Sl_b, dlnc_ds_s0_part_N_b

      real(dp) :: X, Xb, q, qb
      logical :: two_evn_zones, use_G2_v2, on_CHeB, save_mesh
      real(dp), intent(out) :: out(30, 2)
      integer, intent(out) :: out_int(2, 2)
      real(dp) :: a_grd(2, 4), a_int(2, 4)

      X = 0d0
      out = 0d0 / X
      out_int = 0
      dlnc_ds_s0_part_grd = 0d0 / X
      dlnc_ds_s0_part_grd_b = 0d0 / X

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      sigma = (2d0 * pi / 1d6) * nu_q ! nu_q in uHz
      save_mesh = s% x_logical_ctrl(i_save_mesh)

      call calc_JANu(id)
      call calc_PQ(id, sigma)
      call get_evanescent_zones(id, sigma, k_P, k_Q, k_u2, k_l2)

      call calc_evan_radii(id, k_P, k_Q, .true., r_1, r_2, m_1, m_2)
      call calc_evan_radii(id, k_u2, k_l2, .false., r_1b, r_2b, m_1b, m_2b)

      if ((k_P == 0) .or. (k_Q) == 0) then
         out(1, 1) = r_1 / rsun
         out(2, 1) = r_2 / rsun
         out(3, 1) = m_1
         out(4, 1) = m_2
         out_int(1, 1) = k_P
         out_int(2, 1) = k_Q
         return
      end if

      two_evn_zones = .false.
      if ((k_u2 /= 0) .and. (k_l2 /= 0)) then
         two_evn_zones = .true.
      end if

      PQ_integral = 0d0
      G2_div_2kappa_s0 = 0d0

      k_u = min(k_P, k_Q)
      k_l = max(k_P, k_Q)

      on_CHeB = ((s% mixing_type(s% nz) == convective_mixing) .and. (s% center_h1 < 1d-6)  .and. (s% center_he4 > 1d-6))
      use_G2_v2 = .true.
      do k = k_u, k_l
         if ((s% brunt_N2(k) <= 0d0)) then
            use_G2_v2 = .false.
            exit
         end if
      end do

      call calc_quadr_coeffs_grd(id, r_1, r_2, k_u, k_l, a_grd, save_mesh)
      out(17:20, 1) = a_grd(1, :)
      out(21:24, 1) = a_grd(2, :)
      call calc_quadr_coeffs_int(id, r_1, r_2, k_u, k_l, .true., a_int)
      out(25:27, 1) = a_int(1, 1:3)
      out(28:30, 1) = a_int(2, 1:3)
      call calc_PQ_integral_v2(id, r_1, r_2, k_u, k_l, PQ_integral, a_int)
      ! works better than v2 when close intermediate/large regime, but call anyway for checking later
      call calc_G2_div_2kappa_s0(id, r_1, r_2, k_u, k_l, G2_div_2kappa_s0, &
            kappa_s0, NuAJ_s0, dlnPds_s0, dlnQds_s0, &
            dlnc_ds_s0_part_ip, dlnc_ds_s0_part_Sl, dlnc_ds_s0_part_N, a_int)
      if (use_G2_v2) then
         call calc_G2_div_2kappa_s0_v2(id, r_1, r_2, k_u, k_l, G2_div_2kappa_s0, &
               kappa_s0, NuAJ_s0, dlnPds_s0, dlnQds_s0, &
               dlnc_ds_s0_part_grd, dlnc_ds_s0_part_Sl, dlnc_ds_s0_part_N, a_grd, a_int)
      end if

      X = PQ_integral / pi + G2_div_2kappa_s0
      q = (1d0 - sqrt(1d0 - exp(-2d0 * pi * X))) / (1d0 + sqrt(1d0 - exp(-2d0 * pi * X)))

      out_int(1, 1) = k_P
      out_int(2, 1) = k_Q

      out(1, 1) = r_1 / rsun
      out(2, 1) = r_2 / rsun
      out(3, 1) = m_1
      out(4, 1) = m_2

      out(5, 1) = q
      out(6, 1) = PQ_integral / pi
      out(7, 1) = G2_div_2kappa_s0
      out(8, 1) = dlnc_ds_s0_part_ip
      out(9, 1) = dlnc_ds_s0_part_Sl
      out(10, 1) = dlnc_ds_s0_part_N

      out(11, 1) = kappa_s0
      out(12, 1) = NuAJ_s0

      out(13, 1) = dlnPds_s0
      out(14, 1) = dlnQds_s0

      PQ_integral_b = 0d0
      G2_div_2kappa_s0b = 0d0

      if (two_evn_zones) then
         use_G2_v2 = .true.
         do k = k_u2, k_l2 - 1
            if ((s% brunt_N2(k) <= 0d0) .and. (s% brunt_N2(k + 1) <= 0d0)) then
               use_G2_v2 = .false.
               exit
            end if
         end do
         call calc_quadr_coeffs_grd(id, r_1b, r_2b, k_u2, k_l2, a_grd, .false.)  ! works better for gradients
         out(17:20, 2) = a_grd(1, :)
         out(21:24, 2) = a_grd(2, :)
         call calc_quadr_coeffs_int(id, r_1b, r_2b, k_u2, k_l2, .false., a_int)  ! works better for integral
         out(25:27, 2) = a_int(1, 1:3)
         out(28:30, 2) = a_int(2, 1:3)
         call calc_PQ_integral_v2(id, r_1b, r_2b, k_u2, k_l2, PQ_integral_b, a_int)
         call calc_G2_div_2kappa_s0(id, r_1b, r_2b, k_u2, k_l2, G2_div_2kappa_s0b, &
               kappa_s0b, NuAJ_s0b, dlnPds_s0b, dlnQds_s0b, &
               dlnc_ds_s0_part_ip_b, dlnc_ds_s0_part_Sl_b, dlnc_ds_s0_part_N_b, a_int)
         if (use_G2_v2) then  ! works better on small evan zones
            call calc_G2_div_2kappa_s0_v2(id, r_1b, r_2b, k_u2, k_l2, G2_div_2kappa_s0b, &
                  kappa_s0b, NuAJ_s0b, dlnPds_s0b, dlnQds_s0b, &
                  dlnc_ds_s0_part_grd_b, dlnc_ds_s0_part_Sl_b, dlnc_ds_s0_part_N_b, a_grd, &
                  a_int)
         end if

         Xb = PQ_integral_b / pi + G2_div_2kappa_s0b
         qb = (1d0 - sqrt(1d0 - exp(-2d0 * pi * Xb))) / (1d0 + sqrt(1d0 - exp(-2d0 * pi * Xb)))

         out_int(1, 2) = k_u2
         out_int(2, 2) = k_l2

         out(1, 2) = r_1b / rsun
         out(2, 2) = r_2b / rsun
         out(3, 2) = m_1b
         out(4, 2) = m_2b

         out(5, 2) = qb
         out(6, 2) = PQ_integral_b / pi
         out(7, 2) = G2_div_2kappa_s0b
         out(8, 2) = dlnc_ds_s0_part_ip_b
         out(9, 2) = dlnc_ds_s0_part_Sl_b
         out(10, 2) = dlnc_ds_s0_part_N_b

         out(11, 2) = kappa_s0b
         out(12, 2) = NuAJ_s0b

         out(13, 2) = dlnPds_s0b
         out(14, 2) = dlnQds_s0b
      else
         out(:, 2) = 0d0 / PQ_integral_b
         out_int(1, 2) = k_u2
         out_int(2, 2) = k_l2
      end if

      if (s% x_logical_ctrl(i_verbose_coupling)) then
         write(*, *) out_int(1, 1), out_int(2, 1), 'X parts', out(6, 1), out(7, 1), out(5, 1)
         write(*, *) out_int(1, 2), out_int(2, 2), 'X parts', out(6, 2), out(7, 2), out(5, 2)
      end if

   end subroutine do_strong_coupling

   subroutine pen_overshoot(id, ierr)
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type (star_info), pointer :: s

      logical :: core_convective
      integer :: k, k_conv, k_ovhe
      real(dp) :: aovhe, f_aovhe, ovhe_r, soft_r, a_soft
      real(dp) :: Ymin, Ymax

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      core_convective = (s%mixing_type(s%nz) == convective_mixing .and. s%mixing_type(s%nz - 1) == convective_mixing)
      aovhe = s%x_ctrl(i_aovhe_mod_pen_conv)
      a_soft = s%x_ctrl(i_aovhe_mod_pen_conv_exp)

      ! Grow convective core as He depletes
      if (s% x_logical_ctrl(i_grow_aovhe)) then
         Ymin = 0.2d0
         Ymax = 0.9d0
         f_aovhe = (Ymax - s% center_he4) / (Ymax - Ymin)
         f_aovhe = max(f_aovhe, 0d0)
         f_aovhe = min(f_aovhe, 1d0)
         aovhe = aovhe * f_aovhe
      end if
      a_soft = min(a_soft, aovhe)
      if (aovhe <= 0d0) then
         return
      end if

      if (core_convective .and. s% center_h1 < 1d-6 .and. s% center_he4 > 0) then
         k_conv = s% conv_bdy_loc(1) + 1  ! top of convective core, inside convective zone
         s% mixing_type(k_conv:s% nz) = convective_mixing
         s% D_mix(k_conv:s% nz) = s% D_mix(k_conv:s% nz)

         if (s% r(k_conv)<=s% scale_height(k_conv)) then
            ovhe_r = min(s% r(k_conv) * 2, s% r(k_conv) + aovhe * s% scale_height(k_conv))
            soft_r = min(s% r(k_conv) * 2, s% r(k_conv) + (aovhe - a_soft) * s% scale_height(k_conv))
         else
            ovhe_r = s% r(k_conv) + aovhe * s% scale_height(k_conv)
            soft_r = s% r(k_conv) + (aovhe - a_soft) * s% scale_height(k_conv)
         endif

         k_ovhe = k_conv
         do while (s% r(k_ovhe) <= ovhe_r)
            k_ovhe = k_ovhe - 1
         end do

         if (k_ovhe < k_conv) then
            s% mixing_type(k_ovhe:k_conv - 1) = overshoot_mixing
            s% D_mix(k_ovhe:k_conv + 1) = s% D_mix(k_conv + 2)
         end if

         ! Soften end of overshooting
         if (a_soft > 0d0) then
            do k = k_ovhe, k_conv-1
               if (s% r(k) < soft_r) then
                  cycle
               end if
               s% D_mix(k) = s% D_mix(k_conv) * 10**(log10(s% D_mix(k_conv) / s%overshoot_D_min) * (soft_r - s%r(k)) / a_soft)
            end do
         end if

         do k = k_ovhe, k_conv
            s% conv_vel(k) = 3 * s% D_mix(k) / (s% alpha_mlt(k) * (s% scale_height(k - 1) + s% scale_height(k + 1)) / 2)
         end do
!         call other_adjust_mlt_gradT_fraction_Peclet(id, ierr)
      end if

   end subroutine pen_overshoot

   subroutine turbulent_mixing(id, ierr)
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type (star_info), pointer :: s

      real(dp) :: f_turb, n_turb, T_turb, new_Dmix, DHe_0, Rho_0
      real(dp) :: alfa, beta
      integer :: k

      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      f_turb = s% x_ctrl(i_turb_constant)
      n_turb = s% x_ctrl(i_turb_exponent)
      T_turb = s% x_ctrl(i_turb_reference)

      do k = 1, s% nz
         if (s% lnT(k) > log(T_turb)) then
            exit
         end if
      end do
      alfa = (s% lnT(k) - log(T_turb)) / (s% lnT(k) - s% lnT(k - 1))
      beta = 1d0 - alfa

      Rho_0 = (beta * s% Rho(k) + alfa * s% Rho(k - 1))
      DHe_0 = (3.3d-15 * T_turb**2.5d0) / (4d0 * Rho_0 * log(1 + 1.125d-16 * T_turb**3 * Rho_0))

      do k = 1, s% nz
         new_Dmix = f_turb * DHe_0 * (s% Rho(k)/Rho_0)**n_turb
         new_Dmix = min(new_Dmix, 1d10)  ! Upper limit for numerical stability if using relax_tau_factor
         if ((new_Dmix > s% D_mix(k)) .and. (new_Dmix > 1d-3)) then
            s% D_mix(k) = new_Dmix
            s% mixing_type(k) = minimum_mixing
         end if
      end do

   end subroutine turbulent_mixing

   ! Michielsen 2023
   subroutine IGW_D_mix_rho(id, ierr)
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      real(dp) :: rho_switch, D_ext, new_Dmix, igw_exponent
      real(dp) :: f_IGW, f_flash, L_He_full_on, L_He_full_off, f_Mcc, M_cc_full_on, M_cc_full_off
      integer :: k, k_DmixMin
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      k_DmixMin = 0

      ! Turn off IGW mixing during He flash
      L_He_full_on = 7d0
      L_He_full_off = 8d0
      f_flash = min(1d0, max(0d0, (L_He_full_on - log10(s% power_he_burn))/(L_He_full_off - L_He_full_on)))

      ! Turn off IGW for very small convective cores
      M_cc_full_on = 0.05d0
      M_cc_full_off = 0.01d0
      f_Mcc = min(1d0, max(0d0, (M_cc_full_off - s% mass_conv_core)/(M_cc_full_off - M_cc_full_on)))

      f_IGW = min(f_flash, f_Mcc)
      if (f_IGW <= 0d0) then
         return
      end if

      ! IGW profile will not be implemented if MESA is doing the startup of a run
      if ((.not. s% doing_first_model_of_run) .and. (.not. s% doing_relax ) .and. (s% mixing_type(s%nz) .eq. 1)) then
         !------- Determining position where to switch to IGW profile -------
         ! Going from core towards the surface of the star/model, check when Dmix <= s% min_D_mix
         ! and store the index at this position to be used for rescaling of the y-axis below

         if (s% center_h1 < 1d-9) then  ! Lower IGW mixing intensity during CHeB, avoids breathing pulses
            D_ext = s% x_ctrl(i_IGW_D_ext_postMS)
         else
            D_ext = s% x_ctrl(i_IGW_D_ext)
         end if
         D_ext = D_ext * f_IGW

         do k = s% nz, 1, -1
             if (s% D_mix(k) <= D_ext) then
                 k_DmixMin = k
                 exit
             end if
         end do
         !------------------------------------------------------------------
         ! Rescaling along the y-axis of the Dmix profile according to the min_D_mix
         igw_exponent = s% x_ctrl(i_IGW_exponent)
         if (k_DmixMin > 0) then
             rho_switch = s%Rho(k_DmixMin)**igw_exponent
             do k = 1, k_DmixMin
                 new_Dmix = s%Rho(k)**igw_exponent/rho_switch * D_ext
                 if (new_Dmix > s% D_mix(k)) then
                     s% D_mix(k) = new_Dmix
                     s% mixing_type(k) = minimum_mixing
                 end if
             end do
         end if
      end if

   end subroutine IGW_D_mix_rho


   ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   ! Michielsen+ 2023
   ! Adjust temperature gradient in overshoot zone to be adiabatic (Pe>1d2) or radiative (Pe<1d-2) based upon the Peclet number,
   ! with a gradual transition between the two regimes.
   ! Only works if conv_premix = .true. since the last iteration in a time step has NaNs in D_mix if conv_premix = .false.
   ! The other hook on the next line needs to be included in run_star_extras to use this routine.
   ! s% other_adjust_mlt_gradT_fraction => other_adjust_mlt_gradT_fraction_Peclet
   subroutine other_adjust_mlt_gradT_fraction_Peclet(id, ierr)
   integer, intent(in) :: id
   integer, intent(out) :: ierr
   type(star_info), pointer :: s
   real(dp) :: fraction, Peclet_number, conductivity, Hp       ! f is fraction to compose grad_T = f*grad_ad + (1-f)*grad_rad
   integer :: k
   logical, parameter :: DEBUG = .FALSE.

   ierr = 0
   call star_ptr(id, s, ierr)
   if (ierr /= 0) return

   if (s%D_mix(1) .ne. s%D_mix(1)) return  ! To ignore iterations where Dmix and gradT are NaNs

   if (s%num_conv_boundaries .lt. 1) then  ! Is zero at initialisation of the run
   if (DEBUG) then
      write(*,*) ' skip since there are no convective boundaries'
   end if
   return
   endif

   do k= s%nz, 1, -1
      if (s%D_mix(k) <= s% min_D_mix) exit

      conductivity = 16.0_dp * boltz_sigma * pow3(s% T(k)) / ( 3.0_dp * s% opacity(k) * pow2(s% rho(k)) * s% cp(k) )
      Hp = s% Peos(k)/(s% rho(k)*s% grav(k)) ! Pressure scale height
      Peclet_number = s% conv_vel(k) * Hp * s% mixing_length_alpha / conductivity

      if (Peclet_number >= 100.0_dp) then
          fraction = 1.0_dp
      else if (Peclet_number .le. 0.01_dp) then
          fraction = 0.0_dp
      else
          fraction = (safe_log10(Peclet_number)+2.0_dp)/4.0_dp
      end if

      s% adjust_mlt_gradT_fraction(k) = fraction
   end do

   end subroutine other_adjust_mlt_gradT_fraction_Peclet

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   subroutine extended_convective_penetration(id, i, j, k_a, k_b, D, vc, ierr)
   integer, intent(in) :: id, i, j
   integer, intent(out) :: k_a, k_b
   real(dp), intent(out), dimension(:) :: D, vc
   integer, intent(out) :: ierr
   type (star_info), pointer :: s

   logical, parameter :: DEBUG = .FALSE.
   real(dp) :: f, f2, f0
   real(dp) :: D0, Delta0
   real(dp) :: w
   real(dp) :: factor
   real(dp) :: r_cb, Hp_cb
   real(dp) :: r_ob, D_ob, vc_ob
   logical  :: outward
   integer  :: dk, k, k_ob
   real(dp) :: r, dr, r_step
   real(dp) :: Ymin, Ymax, f_aovhe, L_He_full_on, L_He_full_off, f_flash

   ! Evaluate the overshoot diffusion coefficient D(k_a:k_b) and
   ! mixing velocity vc(k_a:k_b) at the i'th convective boundary,
   ! using the j'th set of overshoot parameters. The overshoot
   ! follows the extended convective penetration scheme description by Mathias
   ! Michielsen, "Probing the shape of the mixing profile and of the thermal
   ! structure at the convective core boundary through asteroseismology",
   ! A&A, 628, 76 (2019)

   ierr = 0
   call star_ptr(id, s, ierr)
   if (ierr /= 0) return

   ! Extract parameters
   f = s%overshoot_f(j)        ! extend of step function (a_ov)
   f0 = s%overshoot_f0(j)
   f2 = s%x_ctrl(j)            ! exponential decay (f_ov)

   D0 = s%overshoot_D0(j)
   Delta0 = s%overshoot_Delta0(j)

   if (f < 0.0_dp .OR. f0 <= 0.0_dp .OR. f2 < 0.0_dp) then
      write(*,*) 'ERROR: for extended convective penetration, must set f0 > 0, and f and f2 >= 0'
      write(*,*) 'see description of overshooting in star/defaults/control.defaults'
      ierr = -1
      return
   end if

   ! Grow convective core as He depletes and turn off if flash
   if (s% x_logical_ctrl(i_grow_aovhe)) then
      if (s% overshoot_zone_type(j) == 'burn_He') then
         if (s% overshoot_zone_loc(j) == 'core') then
            if (s% overshoot_bdy_loc(j) == 'top') then
               ! Turn off IGW mixing during He flash
               L_He_full_on = 7d0
               L_He_full_off = 8d0
               f_flash = min(1d0, max(0d0, (L_He_full_on - log10(s% power_he_burn))/(L_He_full_off - L_He_full_on)))

               Ymin = 0.2d0
               Ymax = 0.9d0
               f_aovhe = (Ymax - s% center_he4) / (Ymax - Ymin)
               f_aovhe = max(f_aovhe, 0d0)
               f_aovhe = min(f_aovhe, 1d0)
               f = f * f_aovhe * f_flash
               f0 = f0 * f_aovhe * f_flash
               f2 = f2 * f_aovhe * f_flash
            end if
         end if
      end if
   end if

   ! Apply mass limits
   if (s%star_mass < s%overshoot_mass_full_on(j)) then
      if (s%star_mass > s%overshoot_mass_full_off(j)) then
         w = (s%star_mass - s%overshoot_mass_full_off(j)) / &
         (s%overshoot_mass_full_on(j) - s%overshoot_mass_full_off(j))
         factor = 0.5_dp*(1.0_dp - cospi(w))
         f = f*factor
         f0 = f0*factor
         f2 = f2*factor
      else
         f = 0.0_dp
         f0 = 0.0_dp
         f2 = 0.0_dp
      endif
   endif

   ! Evaluate convective boundary (_cb) parameters
   call star_eval_conv_bdy_r(s, i, r_cb, ierr)
   if (ierr /= 0) return

   call star_eval_conv_bdy_Hp(s, i, Hp_cb, ierr)
   if (ierr /= 0) return

   ! Evaluate overshoot boundary (_ob) parameters
   call star_eval_over_bdy_params(s, i, f0, k_ob, r_ob, D_ob, vc_ob, ierr)
   if (ierr /= 0) return

   ! Loop over cell faces, adding overshoot until D <= overshoot_D_min
   outward = s%top_conv_bdy(i)

   if (outward) then
      k_a = k_ob
      k_b = 1
      dk = -1
   else
      k_a = k_ob+1
      k_b = s%nz
      dk = 1
   endif

   if (f > 0.0_dp) then
      r_step = f*Hp_cb
   else
      r_step = 0.0_dp
   endif

   face_loop : do k = k_a, k_b, dk
      ! Evaluate the extended convective penetration factor
      r = s%r(k)
      if (outward) then
          dr = r - r_ob
      else
          dr = r_ob - r
      endif

      if (dr < r_step .AND. f > 0.0_dp) then  ! step factor
          factor = 1.0_dp
      else
          if ( f2 > 0.0_dp) then                ! exponential factor
              factor = exp(-2.0_dp*(dr-r_step)/(f2*Hp_cb))
          else
              factor = 0.0_dp
          endif
      endif

      ! Store the diffusion coefficient and velocity
      D(k) = (D0 + Delta0*D_ob)*factor
      vc(k) = (D0/D_ob + Delta0)*vc_ob*factor

      ! Check for early overshoot completion
      if (D(k) < s%overshoot_D_min) then
          k_b = k
          exit face_loop
      endif

   end do face_loop

   if (DEBUG) then
      write(*,*) 'step exponential overshoot:'
      write(*,*) '  k_a, k_b   =', k_a, k_b
      write(*,*) '  r_a, r_b   =', s%r(k_a), s%r(k_b)
      write(*,*) '  r_ob, r_cb =', r_ob, r_cb
      write(*,*) '  Hp_cb      =', Hp_cb
   end if

   end subroutine extended_convective_penetration

   subroutine saving_routine(id, save_now, extras_finish_step)
      integer, intent(in) :: id
      logical, intent(out) :: save_now
      integer, intent(inout) :: extras_finish_step
      type(star_info), pointer :: s
      integer :: ierr

      real(dp) :: dT, dlogT, dlogL, dHc, dHec, dlogHc, dlogHec

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      save_now = .false.

      if (first_step) then
         prev_Teff = s% Teff
         prev_L = s% photosphere_L
         prev_Hc = s% center_h1
         prev_Hec = s% center_he4
         first_step = .false.
      end if

      ! Skip if in PMS
      if (s% power_h_burn / s% photosphere_L < 1 .and. s% center_h1 > 0.6) then
         return
      end if

      dT = abs(s% Teff - prev_Teff)
      dlogT = abs(log10(s% Teff / prev_Teff))
      dlogL = abs(log10(s% photosphere_L / prev_L))
      dHc = abs(s% center_h1 - prev_Hc)
      dHec = abs(s% center_he4 - prev_Hec)
      dlogHc = abs(log10(s% center_h1 / prev_Hc))
      dlogHec = abs(log10(s% center_he4 / prev_Hec))

      if (s% x_ctrl(i_save_dT) < dT .and. s% x_ctrl(i_save_dT) > 0d0) then
         save_now = .true.
      else if (s% x_ctrl(i_save_dlogL) < dlogL .and. s% x_ctrl(i_save_dlogL) > 0d0) then
         save_now = .true.
      else if (s% x_ctrl(i_save_dHc) < dHc .and. s% x_ctrl(i_save_dHc) > 0d0) then
         save_now = .true.
      else if (s% x_ctrl(i_save_dHec) < dHec .and. s% x_ctrl(i_save_dHec) > 0d0) then
         save_now  = .true.
      else if (s% center_h1 < dHc .and. s% x_ctrl(i_save_dlogHc) > 0d0 .and. &
         s% x_ctrl(i_save_min_logHc) < log10(s% center_h1)) then
         if (s% x_ctrl(i_save_dlogHc) < dlogHc) then
         save_now = .true.
         end if
      else if (s% center_he4 < dHec .and. s% x_ctrl(i_save_dlogHec) > 0d0 .and. &
         s% x_ctrl(i_save_min_logHec) < log10(s% center_he4)) then
         if (s% x_ctrl(i_save_dlogHec) < dlogHec) then
            save_now = .true.
         end if
      else
         save_now = .false.
      end if

      if (save_now) then
         prev_Teff = s% Teff
         prev_L = s% photosphere_L
         prev_Hc = s% center_h1
         prev_Hec = s% center_he4
      end if

   end subroutine saving_routine

end module run_star_extras
