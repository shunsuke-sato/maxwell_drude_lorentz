module global_variables
  implicit none
  real(8),parameter :: pi = 4d0*atan(1d0)
  real(8),parameter :: a_B = 0.529177210903d-10
  real(8),parameter :: ev = 27.2114d0
  real(8),parameter :: fs = 0.024189d0
  real(8),parameter :: clight = 137d0


  integer,parameter :: num_drude = 1
  integer,parameter :: num_lorentz = 1

! spatial grid
  real(8) :: left_boundary, right_boundary, matter_thickness, dx
  real(8) :: nx_l, nx_r, mx
  real(8),allocatable :: xx_cor(:)

! time grid
  real(8) :: Tprop, dt
  integer :: nt

! electric fields
  real(8),allocatable :: Elec_x(:),Elec_x_old(:),Elec_x_new(:)
  real(8),allocatable :: Lap_Elec_x(:)
  real(8),allocatable :: acc_dns_x(:)


! base matter
  real(8) :: eps0
! Drude parameters
  real(8) :: mass_drude(num_drude), gamma_drude(num_drude), density_drude(num_drude)
  real(8),allocatable :: vt_drude(:,:),vt_drude_old(:,:),vt_drude_new(:,:)

! Lorentz parameters
  real(8) :: mass_lorentz(num_drude), gamma_lorentz(num_drude), density_lorentz(num_drude)
  real(8) :: kconst_lorentz(num_drude)
  real(8),allocatable :: xt_lorentz(:,:),vt_lorentz(:,:)
  real(8),allocatable :: xt_lorentz_old(:,:),vt_lorentz_old(:,:)
  real(8),allocatable :: xt_lorentz_new(:,:),vt_lorentz_new(:,:)
  


end module global_variables
!------------------------------------------------------------------------
program main
  use global_variables
  implicit none


  call set_model_parameters

  call time_propergation


end program main
!------------------------------------------------------------------------
subroutine set_model_parameters
  use global_variables
  implicit none
  integer :: ix

! time propagation
  Tprop = 40d0/fs
  dt = 0.1d0
  nt = aint(Tprop/dt)+1

! base matter
  eps0 = 1d0

! material parameters
  mass_drude(1) = 1d0
  gamma_drude(1) = 1d0
  density_drude(1) = 0d0

  mass_lorentz(1) = 1d0
  gamma_lorentz(1) = 1d0
  density_lorentz(1) = 0d0
  kconst_lorentz(1) = 1d0


! spatial grid
  left_boundary = -10d-6/a_B
  right_boundary = 10d-6/a_B
  matter_thickness = 100d-9/a_B
  dx = 10d-9/a_B

  mx = aint(matter_thickness/dx)+1
  dx = matter_thickness/mx

  nx_l = -aint( abs(left_boundary)/dx ) -1
  nx_r =  aint( abs(right_boundary)/dx ) +1

  allocate(xx_cor(nx_l:nx_r))

  do ix = nx_l, nx_r
    xx_cor(ix) = dx*(dble(ix)-0.5d0)
  end do


  allocate(Elec_x(nx_l:nx_r),Elec_x_old(nx_l:nx_r))
  Elec_x = 0d0
  Elec_x_old = 0d0

  allocate(Lap_Elec_x(nx_l:nx_r))
  
  allocate(acc_dns_x(mx))

  allocate(vt_drude(num_drude, mx),vt_drude_old(num_drude, mx))
  allocate(vt_drude_new(num_drude, mx))
  vt_drude = 0d0; vt_drude_old = 0d0; vt_drude_new = 0d0
  allocate(xt_lorentz(num_drude, mx),xt_lorentz_old(num_drude, mx))
  allocate(xt_lorentz_new(num_drude, mx))
  xt_lorentz = 0d0; xt_lorentz_old = 0d0; xt_lorentz_new = 0d0
  allocate(vt_lorentz(num_drude, mx),vt_lorentz_old(num_drude, mx))
  allocate(vt_lorentz_new(num_drude, mx))
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


  velocity = clight/sqrt(eps0)

  omega_ev = 1.55d0
  pulse_width_fs = 20d0

  omega = omega_ev/ev
  pulse_width = pulse_width_fs/fs

  Elec_x = 0d0
  Elec_x_old = 0d0

  do ix = nx_l, nx_r
    xx = xx_cor(ix)
    tt = -xx/velocity
    
    ss = tt - 0.5d0*tpulse_width
    if(abs(ss) < 0.5d0*tpulse_width)then
      Elec_x(ix) = E0*cos(omega*ss)*cos(0.5d0*pi*ss/tpulse)**4
    end if

    ss = tt - 0.5d0*tpulse_width + dt
    if(abs(ss) < 0.5d0*tpulse_width)then
      Elec_x_old(ix) = E0*cos(omega*ss)*cos(0.5d0*pi*ss/tpulse)**4
    end if
    

  end do
  

  

end subroutine set_initial_laser
!------------------------------------------------------------------------
subroutine time_propergation
  use global_variables
  implicit none
  integer :: it


  call set_initial_laser
  do it = 0, nt

    call dt_propagation

    if(mod(it, nt/200) == 0) call output_field(it)

  end do

end subroutine time_propergation
!------------------------------------------------------------------------
subroutine dt_propagation
  use global_variables
  implicit none

  call dt_newton
  call calc_acc
  call dt_maxwell

  Elec_x_old = Eelec_x
  Eelec_x = Elec_x_new

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
      acc_t = -gamma_drude(imode)*vt_drude(imodel, ix) &
             + Eelec_x(ix)/mass_drude(imode)

      vt_drude_new(imode, ix) = vt_drude_old(imode, ix) + 2d0*dt*acc_t

    end do
  end do


! Lorentz model
  do ix = 1, mx
    do imodel = 1, num_drude
      acc_t = -gamma_lorentz(imode)*vt_lorentz(imodel, ix) &
              -(kconst_lorentz(imode)/mass_lorentz(imode)))*xt_lorentz(imodel, ix) &
             + Eelec_x(ix)/mass_lorentz(imode)

      vt_lorentz_new(imode, ix) = vt_lorentz_old(imode, ix) + 2d0*dt*acc_t
      xt_lorentz_new(imode, ix) = 2d0*xt_lorentz(imode, ix) - xt_lorentz_old(imode, ix) &
                                 + acc_t*dt**2

    end do
  end do

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
      acc_dns_x(ix) = acc_dns_x(ix) + density_drude(num_drude)*acc_t
    end do

! Lorentz model
    do imodel = 1, num_lorentz
      acc_t = 0.5d0*(vt_lorentz_new(imodel,ix)-vt_lorentz_old(imodel,ix))/dt
      acc_dns_x(ix) = acc_dns_x(ix) + density_lorentz(num_drude)*acc_t
    end do

  end do


end subroutine calc_acc
!------------------------------------------------------------------------
subroutine dt_maxwell
  use global_variables
  implicit none
  integer :: ix
  real(8) :: velocity_c

! calc laplacian
  ix = nx_l
  Lap_Elec_x(ix) = (Elec_x(ix+1)-2d0*Elec_x(ix+1))/dx**2
  do ix = nx_l+1, nx_r-1
    Lap_Elec_x(ix) = (Elec_x(ix+1)-2d0*Elec_x(ix+1)+Elec_x(ix-1))/dx**2
  end do
  ix = nx_r
  Lap_Elec_x(ix) = (-2d0*Elec_x(ix+1)+Elec_x(ix-1))/dx**2

  velocity_c = clight/sqrt(eps0)

  Eelec_x_new = 2d0*Elec_x -Elec_x_old +velocity_c**2*dt**2*Lap_Elec_x


  Elec_x_new(1:mx) = -4d0*pi/clight**2*dt**2*acc_dns_x(1:mx)
  


end subroutine dt_maxwell
!------------------------------------------------------------------------
subroutine output_field(it)
  use global_variables
  implicit none
  integer,intent(in) :: it
  character(256) :: cit, cfilename

  write(cit, "(I9.9)")it
  cfilename = "Efields_x_"//trim(cit)

  open(20,file=cfilename)
  
  do ix = nx_l, nx_r
    write(20,"(99e26.16e3)")xx_cor(ix), Elec_x(ix)
  end do

  close(20)
  


end subroutine output_field
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
