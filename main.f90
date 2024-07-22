module global_variables
  implicit none
  complex(8),parameter :: zi = (0d0, 1d0)
  real(8),parameter :: pi = 4d0*atan(1d0)
  real(8),parameter :: a_B = 0.529177210903d-10
  real(8),parameter :: ev = 27.2114d0
  real(8),parameter :: fs = 0.024189d0
  real(8),parameter :: clight = 137d0


  integer,parameter :: num_drude = 1
  integer,parameter :: num_lorentz = 2

! spatial grid
  real(8) :: left_boundary, right_boundary, matter_thickness, dx
  integer :: nx_l, nx_r, mx
  real(8),allocatable :: xx_cor(:)

! time grid
  real(8) :: Tprop, dt
  integer :: nt

! electric fields
  real(8),allocatable :: Elec_x(:),Elec_x_old(:),Elec_x_new(:)
  real(8),allocatable :: Lap_Elec_x(:)
  real(8),allocatable :: acc_dns_x(:)


! base matter
  real(8) :: eps0, sigma_re, eps_D
! Drude parameters
  real(8) :: mass_drude(num_drude), gamma_drude(num_drude), density_drude(num_drude)
  real(8),allocatable :: vt_drude(:,:),vt_drude_old(:,:),vt_drude_new(:,:)

! Lorentz parameters
  real(8) :: mass_lorentz(num_lorentz), gamma_lorentz(num_lorentz), density_lorentz(num_lorentz)
  real(8) :: kconst_lorentz(num_lorentz)
  real(8),allocatable :: xt_lorentz(:,:),vt_lorentz(:,:)
  real(8),allocatable :: xt_lorentz_old(:,:),vt_lorentz_old(:,:)
  real(8),allocatable :: xt_lorentz_new(:,:),vt_lorentz_new(:,:)
  
! Averaged value in film
  real(8) :: current_ave, Efield_ave, dEdt_ave

! experimental data
  integer,parameter :: nt_exp =  11135 !111370
  real(8)  :: tt_exp(nt_exp)
  real(8) :: E_exp_in(5,nt_exp)
  real(8) :: E_exp(nt_exp)

! Fourier-filtering (experimental data)
  integer :: nt_exp_f
  real(8),allocatable  :: tt_exp_f(:)
  real(8),allocatable :: E_exp_f(:), E_exp_f_m_dt(:)
  real(8) :: tshift_exp_f


end module global_variables
!------------------------------------------------------------------------
program main
  use global_variables
  implicit none


  call set_model_parameters
  call check_dielectric_function
  call time_propergation


end program main
!------------------------------------------------------------------------
subroutine set_model_parameters
  use global_variables
  implicit none
  integer :: ix
  real(8) :: omega_p, factor_vac

  factor_vac = 1d0
! time propagation
  Tprop = 60d0/fs
!  dt = 0.05d0
  dt = 0.05d0
  nt = aint(Tprop/dt)+1
  write(*,*)"nt = ",nt

! base matter
  eps0 = 1d0 ! vacuum dielectric constant
  eps_D = 1d0  ! Matter dielectric constant
  sigma_re = 9.5d-2*factor_vac

! material parameters
  omega_p = 8.5d0/ev
  mass_drude(1) = 1d0
  gamma_drude(1) = 1d0/(14d0/0.024189d0)
  density_drude(1) = omega_p**2*0.080d0
  density_drude(1) = 0d0

  mass_lorentz(1) = 1d0
  gamma_lorentz(1) = (0.2d0/ev)
  density_lorentz(1) = 1.42d-3*factor_vac
!  density_lorentz(1) = 0d0
  kconst_lorentz(1) = (2.2d0/ev)**2


  mass_lorentz(2) = 1d0
  gamma_lorentz(2) = gamma_lorentz(1)
  density_lorentz(2) = density_lorentz(1)*0.20d0*factor_vac
!  density_lorentz(2) = 0d0
  kconst_lorentz(2) = (2.47d0/ev)**2



! spatial grid
  left_boundary = -20d-6/a_B
  right_boundary = 20d-6/a_B
  matter_thickness = 20d-9/a_B
!  dx = 1d-9/a_B
  dx = 5d-10/a_B

  mx = aint(matter_thickness/dx)+1
  write(*,*)'mx = ',mx
  dx = matter_thickness/mx

  nx_l = -aint( abs(left_boundary)/dx ) -1
  nx_r =  aint( abs(right_boundary)/dx ) +1

  allocate(xx_cor(nx_l:nx_r))

  do ix = nx_l, nx_r
    xx_cor(ix) = dx*(dble(ix)-0.5d0)
  end do


  allocate(Elec_x(nx_l:nx_r),Elec_x_old(nx_l:nx_r))
  allocate(Elec_x_new(nx_l:nx_r))
  Elec_x = 0d0
  Elec_x_old = 0d0

  allocate(Lap_Elec_x(nx_l:nx_r))
  
  allocate(acc_dns_x(mx))

  allocate(vt_drude(num_drude, mx),vt_drude_old(num_drude, mx))
  allocate(vt_drude_new(num_drude, mx))
  vt_drude = 0d0; vt_drude_old = 0d0; vt_drude_new = 0d0
  allocate(xt_lorentz(num_lorentz, mx),xt_lorentz_old(num_lorentz, mx))
  allocate(xt_lorentz_new(num_lorentz, mx))
  xt_lorentz = 0d0; xt_lorentz_old = 0d0; xt_lorentz_new = 0d0
  allocate(vt_lorentz(num_lorentz, mx),vt_lorentz_old(num_lorentz, mx))
  allocate(vt_lorentz_new(num_lorentz, mx))
  vt_lorentz = 0d0; vt_lorentz_old = 0d0; vt_lorentz_new = 0d0




end subroutine set_model_parameters
!------------------------------------------------------------------------
subroutine set_initial_laser
  use global_variables
  implicit none
  real(8) :: omega_ev, pulse_width_fs, E0
  real(8) :: omega, pulse_width
  real(8) :: xx, tt, velocity, ss
  integer :: ix


  E0 = 1d0
  velocity = clight/sqrt(eps0)

  omega_ev = 2.0d0
!  pulse_width_fs = 20d0
  pulse_width_fs = 10d0

  omega = omega_ev/ev
  pulse_width = pulse_width_fs/fs

  Elec_x = 0d0
  Elec_x_old = 0d0

  do ix = nx_l, nx_r
    xx = xx_cor(ix)
    tt = -xx/velocity
    
    ss = tt - 0.5d0*pulse_width
    if(abs(ss) < 0.5d0*pulse_width)then
      Elec_x(ix) = E0*cos(omega*ss)*cos(pi*ss/pulse_width)**4
    end if

    ss = tt - 0.5d0*pulse_width - dt
    if(abs(ss) < 0.5d0*pulse_width)then
      Elec_x_old(ix) = E0*cos(omega*ss)*cos(pi*ss/pulse_width)**4
    end if
    

  end do
  

  

end subroutine set_initial_laser
!------------------------------------------------------------------------
subroutine set_initial_laser_exp
  use global_variables
  implicit none
  real(8) :: xx, tt, velocity, ss, Et
  integer :: ix

  call read_exp_data

  velocity = clight/sqrt(eps0)

  Elec_x = 0d0
  Elec_x_old = 0d0

  do ix = nx_l, nx_r
    xx = xx_cor(ix)
    tt = -xx/velocity

    call calc_field_strength(Et,tt)
    Elec_x(ix) = Et

    xx = xx_cor(ix)
    tt = -xx/velocity -dt
    call calc_field_strength(Et,tt)
    Elec_x_old(ix) = Et
    

  end do
  
  


contains
  subroutine calc_field_strength(Et,tt)
    implicit none
    real(8) :: tt, Et
    integer :: it
    real(8) :: r1,r2

    if(tt<minval(tt_exp_f) .or. tt>= maxval(tt_exp_f))then
      Et = 0d0
      return
    end if

    
    do it = 0, nt_exp_f-1
      if(tt < tt_exp_f(it))then
        r1 = (tt - tt_exp_f(it-1))/(tt_exp_f(it) - tt_exp_f(it-1))
        r2 = 1d0 - r1
        Et = r1 *E_exp_f(it) + r2*E_exp_f(it-1)
        return
      end if
    end do
  end subroutine calc_field_strength

end subroutine set_initial_laser_exp
!------------------------------------------------------------------------
subroutine time_propergation
  use global_variables
  implicit none
  integer :: it
  real(8) :: tshift

  tshift = matter_thickness/(clight/sqrt(eps0))

  call set_initial_laser
!  call set_initial_laser_exp

  open(101,file="Et_vac.out")
  write(101,"(A)")"# t (a.u.), t-shifted (a.u.), E_front(t), E_rear(t), current_ave, Efield_ave, dEdt_ave"

  do it = 0, nt
    call calc_average_in_film
    write(101,"(999e26.16e3)")it*dt,it*dt-tshift, Elec_x(0), Elec_x(mx+1), &
        current_ave, Efield_ave, dEdt_ave
    call dt_propagation

!    if(mod(it, 200) == 0) call output_field(it)


  end do
  close(101)

end subroutine time_propergation
!------------------------------------------------------------------------
subroutine dt_propagation
  use global_variables
  implicit none

  call dt_newton
  call calc_acc
  call dt_maxwell

  Elec_x_old = Elec_x
  Elec_x = Elec_x_new

  vt_drude_old = vt_drude
  vt_drude = vt_drude_new

  xt_lorentz_old = xt_lorentz
  xt_lorentz = xt_lorentz_new
  
  vt_lorentz_old = vt_lorentz
  vt_lorentz = vt_lorentz_new

end subroutine dt_propagation
!------------------------------------------------------------------------
subroutine dt_newton
  use global_variables
  implicit none
  integer :: ix, imodel
  real(8) :: force, acc_t

! Drude model
  do ix = 1, mx
    do imodel = 1, num_drude
      acc_t = -gamma_drude(imodel)*vt_drude(imodel, ix) &
             + Elec_x(ix)/mass_drude(imodel)

      vt_drude_new(imodel, ix) = vt_drude_old(imodel, ix) + 2d0*dt*acc_t

    end do
  end do


! Lorentz model
  do ix = 1, mx
    do imodel = 1, num_lorentz
      acc_t = -gamma_lorentz(imodel)*vt_lorentz(imodel, ix) &
              -(kconst_lorentz(imodel)/mass_lorentz(imodel))*xt_lorentz(imodel, ix) &
             + Elec_x(ix)/mass_lorentz(imodel)

      vt_lorentz_new(imodel, ix) = vt_lorentz_old(imodel, ix) + 2d0*dt*acc_t
      xt_lorentz_new(imodel, ix) = 2d0*xt_lorentz(imodel, ix) - xt_lorentz_old(imodel, ix) &
                                 + acc_t*dt**2

    end do
  end do

!  write(*,*)vt_lorentz_new(1,1),vt_drude_new(1,1)


end subroutine dt_newton
!------------------------------------------------------------------------
subroutine calc_acc
  use global_variables
  implicit none
  integer :: ix, imodel
  real(8) :: acc_t

  acc_dns_x = 0d0

  do ix = 1, mx

! Drude model
    do imodel = 1, num_drude
      acc_t = 0.5d0*(vt_drude_new(imodel,ix)-vt_drude_old(imodel,ix))/dt
      acc_dns_x(ix) = acc_dns_x(ix) + density_drude(imodel)*acc_t
    end do

! Lorentz model
    do imodel = 1, num_lorentz
      acc_t = 0.5d0*(vt_lorentz_new(imodel,ix)-vt_lorentz_old(imodel,ix))/dt
      acc_dns_x(ix) = acc_dns_x(ix) + density_lorentz(imodel)*acc_t
    end do

  end do


end subroutine calc_acc
!------------------------------------------------------------------------
subroutine calc_average_in_film
  use global_variables
  implicit none
  integer :: ix, imodel


  current_ave = 0d0
  Efield_ave = 0d0
  dEdt_ave = 0d0


  do ix = 1, mx

    Efield_ave = Efield_ave + Elec_x(ix)
    dEdt_ave = dEdt_ave + (Elec_x(ix)- Elec_x_old(ix))/dt

! Drude model
    do imodel = 1, num_drude
      current_ave = current_ave + vt_drude(imodel, ix)*density_drude(imodel)
    end do

! Lorentz model
    do imodel = 1, num_lorentz
      current_ave = current_ave + vt_lorentz(imodel, ix)*density_lorentz(imodel)
    end do

! metalic current
      current_ave = current_ave + sigma_re*Elec_x(ix)

  end do

  current_ave = current_ave*dx
  Efield_ave = Efield_ave*dx
  dEdt_ave = dEdt_ave*dx

end subroutine calc_average_in_film
!------------------------------------------------------------------------
subroutine dt_maxwell
  use global_variables
  implicit none
  integer :: ix
  real(8) :: velocity_c, factor

! calc laplacian
  ix = nx_l
  Lap_Elec_x(ix) = (Elec_x(ix+1)-2d0*Elec_x(ix))/dx**2
  do ix = nx_l+1, nx_r-1
    Lap_Elec_x(ix) = (Elec_x(ix+1)-2d0*Elec_x(ix)+Elec_x(ix-1))/dx**2
  end do
  ix = nx_r
  Lap_Elec_x(ix) = (-2d0*Elec_x(ix)+Elec_x(ix-1))/dx**2

  velocity_c = clight/sqrt(eps0)


! vacuum
  Elec_x_new(nx_l:0) = 2d0*Elec_x(nx_l:0) - Elec_x_old(nx_l:0) &
                     + velocity_c**2*dt**2*Lap_Elec_x(nx_l:0)
  Elec_x_new(mx+1:nx_r) = 2d0*Elec_x(mx+1:nx_r) - Elec_x_old(mx+1:nx_r) &
                     + velocity_c**2*dt**2*Lap_Elec_x(mx+1:nx_r)

! matter
  factor=(eps_D/dt**2+2d0*pi*sigma_re/dt)

  Elec_x_new(1:mx) = velocity_c**2*Lap_Elec_x(1:mx) -4d0*pi*acc_dns_x(1:mx) &
      + 2d0*eps_D*Elec_x(1:mx)/(dt**2) &
      - (eps_D/dt**2-2d0*pi*sigma_re/dt)*Elec_x_old(1:mx)

  Elec_x_new(1:mx) = Elec_x_new(1:mx)/factor


!  Elec_x_new = 2d0*Elec_x -Elec_x_old +velocity_c**2*dt**2*Lap_Elec_x


!  write(*,*)nx_l, nx_r, mx
!  Elec_x_new(1:mx) = Elec_x_new(1:mx) &
!      -4d0*pi*(velocity_c/clight)**2*dt**2*acc_dns_x(1:mx)
  


end subroutine dt_maxwell
!------------------------------------------------------------------------
subroutine output_field(it)
  use global_variables
  implicit none
  integer,intent(in) :: it
  character(256) :: cit, cfilename
  integer :: ix

  write(cit, "(I9.9)")it
  cfilename = "Efields_x_"//trim(cit)//".out"

  open(20,file=cfilename)
  
  do ix = nx_l, nx_r
    write(20,"(99e26.16e3)")xx_cor(ix), Elec_x_old(ix)
  end do

  close(20)
  


end subroutine output_field
!------------------------------------------------------------------------
subroutine check_dielectric_function
  use global_variables
  implicit none
  integer,parameter :: nw = 600
  real(8),parameter :: wi = 0.1d0/ev, wf =30d0/ev, dw =(wf-wi)/nw
  integer :: iw, imodel
  real(8) :: ww, w0
  complex(8) :: zeps


  open(30,file="epsilon.out")
  do iw = 0, nw
    ww = wi + dw*iw

    zeps = eps0

! Drude
    do imodel = 1, num_drude
      zeps = zeps + (4d0*pi*zi/ww) &
                  * (density_drude(imodel)/mass_drude(imodel)) &
                  * 1d0/(gamma_drude(imodel)-zi*ww)
    end do

! Lorentz
    do imodel = 1, num_lorentz
      w0 = sqrt(kconst_lorentz(imodel)/mass_lorentz(imodel))
      zeps = zeps + (4d0*pi*density_lorentz(imodel)/mass_lorentz(imodel)) &
                  * 1d0/(w0**2-ww**2-zi*gamma_lorentz(imodel)*ww)
    end do

    write(30,"(99e26.16e3)")ww, zeps
    
  end do
  close(30)


end subroutine check_dielectric_function
!------------------------------------------------------------------------
subroutine read_exp_data
  use global_variables
  implicit none
  real(8),parameter :: fact_common = 1d35
  real(8) :: t0
  integer :: it

  tt_exp = 0d0
  open(20,file="time_axis.txt")
  read(20,*)tt_exp(1:nt_exp)
  close(20)

  open(20,file="no_np.txt")
  do it = 1, nt_exp
    read(20,*)E_exp_in(1:5,it)
  end do
  close(20)

  do it = 1, nt_exp
    E_exp(it) = fact_common*sum(E_exp_in(1:5,it))/5d0
  end do

  t0 = minval(tt_exp)
  tt_exp = tt_exp - t0
  tt_exp = tt_exp/fs
  

  open(30,file="efields.out")
  do it = 1, nt_exp
    write(30, "(999e26.16e3)")tt_exp(it)*fs, E_exp(it)
  end do
  close(30)

  call Fourier_filter_exp_data

end subroutine read_exp_data
!------------------------------------------------------------------------
subroutine Fourier_filter_exp_data
  use global_variables
  implicit none
  real(8),parameter :: ww_cut = 3d0/ev, sigma_ww = 0.25d0/ev
  integer :: nt_shift
  complex(8),allocatable :: zEw(:)
  complex(8) :: zs
  integer :: iw, it
  real(8) :: ww, ss, xx, mask, dt_tmp
  real(8),allocatable :: Et_org(:)

  tshift_exp_f = 5d0/fs
  dt_tmp = (tt_exp(2)-tt_exp(1))
  nt_shift = tshift_exp_f/dt_tmp


  
  nt_exp_f = nt_exp + 2*nt_shift
  allocate(tt_exp_f(0:nt_exp_f-1))
  allocate(E_exp_f(0:nt_exp_f-1))
  allocate(Et_org(0:nt_exp_f-1))
  allocate(E_exp_f_m_dt(0:nt_exp_f-1))

  E_exp_f = 0d0
  E_exp_f(nt_shift:nt_shift+nt_exp-1) = E_exp(1:nt_exp)



  Et_org = E_exp_f
  do it = 0, nt_exp_f-1
    tt_exp_f(it) = it*dt_tmp
  end do

  allocate(zEw(-nt_exp_f/2:nt_exp_f/2))

! Forward Fourier transform
  do iw = -nt_exp_f/2, nt_exp_f/2
    zs = 0d0
    do it = 0, nt_exp_f-1
      zs = zs + E_exp_f(it)*exp(zi*2d0*pi*dble(iw*it)/nt_exp_f)
    end do
    zEw(iw) = zs
  end do

! Filtering in Fourier space
    do iw = -nt_exp_f/2, nt_exp_f/2
      ww = 2d0*pi*dble(iw)/(nt_exp_f*dt_tmp)
      mask = 1d0/(exp((abs(ww)-ww_cut)/sigma_ww)+1d0)
      zEw(iw) = zEw(iw)*mask
    end do


! Backward Fourier transform
    do it = 0, nt_exp_f-1
      ss = 0d0
      do iw = -nt_exp_f/2, nt_exp_f/2
        ss = ss + zEw(iw)*exp(-zi*2d0*pi*dble(iw*it)/nt_exp_f)
      end do
      E_exp_f(it) = ss/nt_exp_f
    end do

! Real-time filter
    do it = 0, nt_shift
      xx = dble(it)/nt_shift
      ss = sin(0.5d0*pi*xx)**2
      E_exp_f(it) = E_exp_f(it)*ss
      E_exp_f(nt_exp_f-1-it) = E_exp_f(nt_exp_f-1-it)*ss
    end do


    open(101,file="Et_filtering_exp.out")
    do it = 0, nt_exp_f-1
      write(101, "(999e26.16e3)")tt_exp_f(it), E_exp_f(it), Et_org(it)
    end do
    close(101)


end subroutine Fourier_filter_exp_data
!------------------------------------------------------------------------
