program check_transport
  use patmo_volc, only: volc_transport, volc_form_sulfate
  implicit none
  integer,parameter::n=120
  real*8::q(n),initial(n),dz(n),air(n),k(n),v(0:n),z(n),p(n),a(n)
  real*8::loss,mean0,variance0,mean1,variance1,errors(3)
  integer::i,j
  dz=1d5
  z=[(dble(i)*1d5,i=1,n)]
  air=exp(-z/7d5)
  k=6d3
  v=0d0
  q=air*dz
  initial=q
  call volc_transport(q,dz,air,k,v,1d6,loss)
  call check(maxval(abs(q-initial))/maxval(initial)<1d-12,'constant mixing ratio stationary')
  dz=[(0.5d5+dble(i)*1d3,i=1,n)]
  q=air*dz
  initial=q
  call volc_transport(q,dz,air,k,v,1d6,loss)
  call check(maxval(abs(q-initial))/maxval(initial)<1d-12,'nonuniform layer volumes')
  dz=1d5
  air=1d0
  k=0d0
  q=0d0
  q(1)=1d0
  v=-0.2d0
  call volc_transport(q,dz,air,k,v,1d6,loss)
  call check(abs(sum(q)+loss-1d0)<1d-12.and.minval(q)>=0d0,'bottom loss budget and positivity')
  q=exp(-0.5d0*((z-60d5)/4d5)**2)
  q=q/sum(q)
  initial=q
  mean0=sum(q*z)
  variance0=sum(q*(z-mean0)**2)
  call volc_transport(q,dz,air,k,v,86400d0*10d0,loss)
  mean1=sum(q*z)/sum(q)
  variance1=sum(q*(z-mean1)**2)/sum(q)
  call check(abs(mean1-mean0+0.2d0*864000d0)<2d2,'settling centroid displacement')
  call check(abs(sum(q)+loss-1d0)<1d-12.and.minval(q)>=0d0,'settling conservation')
  print *, 'Numerical variance increase [km2]: ',(variance1-variance0)/1d10
  call check(abs(variance1-variance0)<0.2d10,'limited numerical spreading')
  q=initial
  v=0d0
  k=1d4
  call volc_transport(q,dz,air,k,v,864000d0,loss)
  mean1=sum(q*z)
  variance1=sum(q*(z-mean1)**2)
  call check(abs((variance1-variance0)/(2d0*1d4*864000d0)-1d0)<1d-6,'physical diffusion variance')
  p=0d0
  a=0d0
  p(60)=1d0
  do j=1,100
     call volc_form_sulfate(p,a,864d0,86400d0)
  end do
  call check(abs(sum(p)-exp(-1d0))<1d-12,'precursor e-folding')
  call check(abs(sum(p+a)-1d0)<1d-12,'formation budget conservation')
  do j=1,3
     call convergence_error(40*2**j,errors(j))
  end do
  print *, 'Advection L1 errors [1, 0.5, 0.25 km]: ',errors
  call check(errors(2)<0.65d0*errors(1).and.errors(3)<0.65d0*errors(2),'grid convergence')
  print *, 'All volcanic transport checks passed.'
contains
  subroutine check(ok,label)
    logical,intent(in)::ok
    character(len=*),intent(in)::label
    if(.not.ok) then
       print *, 'FAIL: ',label
       error stop 1
    end if
    print *, 'PASS: ',label
  end subroutine check
  subroutine convergence_error(nc,error)
    integer,intent(in)::nc
    real*8,intent(out)::error
    real*8::state(nc),width(nc),density(nc),diff(nc),speed(0:nc),center(nc),exact(nc),sigma,shift
    integer::m
    width=80d5/nc
    center=[((m-0.5d0)*width(1),m=1,nc)]
    density=1d0
    diff=0d0
    speed=-0.2d0
    sigma=2d5
    shift=-0.2d0*864000d0
    state=0.5d0*(erf((center+0.5d0*width-40d5)/(sqrt(2d0)*sigma)) &
         -erf((center-0.5d0*width-40d5)/(sqrt(2d0)*sigma)))
    exact=0.5d0*(erf((center+0.5d0*width-40d5-shift)/(sqrt(2d0)*sigma)) &
         -erf((center-0.5d0*width-40d5-shift)/(sqrt(2d0)*sigma)))
    call volc_transport(state,width,density,diff,speed,864000d0)
    error=sum(abs(state-exact))
  end subroutine convergence_error
end program check_transport
