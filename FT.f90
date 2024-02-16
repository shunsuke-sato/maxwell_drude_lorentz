program main
  implicit none
  complex(8),parameter :: zi = (0d0, 1d0)
  integer,parameter :: nt = 16537, nw = 100
  real(8),parameter :: ev = 27.2114d0
  real(8),parameter :: wi = 1d0/ev, wf = 3d0/ev
  real(8),parameter :: dw = (wf-wi)/nw
  real(8) :: Et_vac_w_np(0:nt), Et_vac_wo_np(0:nt)
  complex(8) :: zEw_vac_w_np(0:nw), zEw_vac_wo_np(0:nw)
  real(8) :: tt(0:nt), dt, f1
  complex(8) :: zt1, zt2
  real(8) :: ww
  integer :: it, iw

  open(20,file="Et_vac_w_np.out")
  read(20,*)
  do it = 0, nt
    read(20,*)tt(it), f1, Et_vac_w_np(it)
  end do
  close(20)
  dt = tt(1)-tt(0)


  open(20,file="Et_vac_wo_np.out")
  read(20,*)
  do it = 0, nt
    read(20,*)tt(it), f1, Et_vac_wo_np(it)
  end do
  close(20)
  dt = tt(1)-tt(0)
  

  do iw = 0, nw
    ww = wi + dw*iw
    zt1 = 0d0; zt2 = 0d0
    do it = 0,nt
      zt1 = zt1 + exp(zi*ww*tt(it))*Et_vac_w_np(it)
      zt2 = zt2 + exp(zi*ww*tt(it))*Et_vac_wo_np(it)
    end do
    zEw_vac_w_np(iw)  = zt1*dw
    zEw_vac_wo_np(iw) = zt2*dw
  end do


  open(30,file="absorbance.out")
  do iw = 0, nw
    ww = wi + dw*iw
    write(30,"(999e26.16e3)")ww,abs(zEw_vac_wo_np(iw))**2/abs(zEw_vac_w_np(iw))**2
  end do
  close(30)
  
end program main
