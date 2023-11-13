*
       program init
      implicit real*8 (a-h,o-z)
      parameter (Imax=513, Jmax=65, Kmax=129)
      character*12 fncp,fndat
      dimension
     > u(0:Imax,0:Jmax,0:Kmax)
     >,v(0:Imax,0:Jmax,0:Kmax)
     >,w(0:Imax,0:Jmax,0:Kmax)
     >,p(0:Imax,0:Jmax,0:Kmax)
     >,u1(0:Imax,0:Jmax,0:Kmax)
     >,v1(0:Imax,0:Jmax,0:Kmax)
     >,w1(0:Imax,0:Jmax,0:Kmax)
     >,p1(0:Imax,0:Jmax,0:Kmax)
      common
     >/dim/Xmax,epsr
     >/dimx/hx,Im,lx
     >/dimr/rt(0:128),rt1(0:128),yt(129),yt1(129),hr,Jm
     >/dimt/ht,Km,lt
     >/Re/Re
     >/pi/pi
*
      open(5,file='init.car')
      read(5,*) Xmax
      read(5,*) lx
      read(5,*) Jm, epsr
      read(5,*) lt
      read(5,*) Re
      read(5,*) amp
      read(5,*) dt
      read(5,*) fncp
*
      call com
      open(9,file=fncp,form='unformatted',status='new',err=1)
      t=0.d0
      do k=1,Km
        tt=(1.d0*k)/Km
        do j=1,Jm
          r=yt(j)
          do i=1,Im
            x=(1.d0*i)/Im
            u(i,j,k)=1.-r**2
            v(i,j,k)=0.
            w(i,j,k)=0.
            u1(i,j,k)=amp*r*(sin(2.*pi*(2.*x-tt)
     >       +r*sin(2.*pi*(x+2.*tt))))
            v1(i,j,k)=0.
            w1(i,j,k)=0.
          end do
        end do
      end do
      dd=0.
      do k=1,Km
        do j=1,Jm
          do i=1,Im
            call div(i,j,k,u1,v1,w1,d,Imax,Jmax)
            dd=max(dd,abs(d))
          end do
        end do
      end do
      write(*,*)' Div0puls=',dd
      p1(0,0,1)=0.d0
      call pres(u1,v1,w1,p1,Imax,Jmax)
      dd=0.
      do k=1,Km
        do j=1,Jm
          do i=1,Im
            call div(i,j,k,u1,v1,w1,d,Imax,Jmax)
            dd=max(dd,abs(d))
          end do
        end do
      end do
      write(*,*)' Div1 puls =',dd
      amp_prev = 0.
      call calc_amp(u1,v1,w1,amp_prev, Imax, Jmax)
      write(*,*) "amp= ", amp_prev
      call correct_amp(u1,v1,w1,amp_prev, amp, Imax, Jmax)
      call calc_amp(u1,v1,w1,amp_prev, Imax, Jmax)
      write(*,*) "amp_new", amp_prev

      do k=1,Km
        do j=1,Jm
          do i=1,Im
            u(i,j,k) = u(i,j,k) + u1(i,j,k)
            v(i,j,k) = v(i,j,k) + v1(i,j,k)
            w(i,j,k) = w(i,j,k) + w1(i,j,k)
          end do
        end do
      end do

      p(0,0,1)=0.5d0
      call pres(u,v,w,p,Imax,Jmax)
      dd=0.
      do k=1,Km
        do j=1,Jm
          do i=1,Im
            call div(i,j,k,u,v,w,d,Imax,Jmax)
            dd=max(dd,abs(d))
          end do
        end do
      end do
      write(*,*)' Div flow =',dd


      call ubulk(u,Imax,Jmax,ub)
      write(*,*) "ubulk=", ub
      Dp=4.d0/Re
      write(9)t,dt,Dp,Re,Xmax,epsr,lx,Jm,lt
      do k=1,Km
        do j=1,Jm
          write(9)(u(i,j,k),i=1,Im)
          write(9)(v(i,j,k),i=1,Im)
          write(9)(w(i,j,k),i=1,Im)
        end do
      end do
      close(9)
      stop
1     write(*,*)'  File ',fncp,' already exists'
      stop
      end program init
*
      subroutine gradp(u,v,w,p,Imax,Jmax)
      implicit real*8 (a-h,o-z)
      dimension
     > u(0:Imax,0:Jmax,0:*)
     >,v(0:Imax,0:Jmax,0:*)
     >,w(0:Imax,0:Jmax,0:*)
     >,p(0:Imax,0:Jmax,0:*)
      common
     >/dimx/hx,Im,lx
     >/dimr/rt(0:128),rt1(0:128),yt(129),yt1(129),hr,Jm
     >/dimt/ht,Km,lt
      Dp=p(0,0,0)
      do k=1,Km
        do j=1,Jm
          do i=1,Im-1
            u(i,j,k)=u(i,j,k)-(p(i+1,j,k)-p(i,j,k))/hx+Dp
          end do
          i=Im
            u(i,j,k)=u(i,j,k)-(p(1,j,k)-p(i,j,k))/hx+Dp
        end do
      end do
*
      do k=1,Km
        do j=1,Jm-1
          do i=1,Im
            v(i,j,k)=v(i,j,k)-(p(i,j+1,k)-p(i,j,k))/rt1(j)
          end do
        end do
      end do
*
      do j=1,Jm
        do i=1,Im
          do k=1,Km-1
            w(i,j,k)=w(i,j,k)-(p(i,j,k+1)-p(i,j,k))/(yt(j)*ht)
          end do
          k=Km
            w(i,j,k)=w(i,j,k)-(p(i,j,1)-p(i,j,k))/(yt(j)*ht)
        end do
      end do
      return
      end
*
*
      real*8 function rrt(x,i)
      implicit real*8 (a-h,o-z)
      common/set/y0,y01,aset,bset,iset
      yp(a,b)=0.01*sqrt(0.3*(b+1))*(a*(58.*b-21.-21.*b*b)+100.)
      y1p(a,b)=1.+a*(b-0.45*(b+1.)**2)
      dyda(a,b)=0.01*sqrt(0.3*(b+1))*(58.*b-21.-21.*b*b)
      dydb(a,b)=0.01*sqrt(0.3*(b+1))*
     *         (0.5*(a*(58.*b-21.-21.*b*b)+100.)/(b+1)+a*(58.-42*b))
      dy1da(a,b)=b-0.45*(b+1.)**2
      dy1db(a,b)=a*(1.-0.9*(b+1.))
      if(iset.ne.0) goto 2
      a=0
      b=10./3.*y0*y0-1.
      nit=0
1     continue
      c1=dyda(a,b)
      c2=dydb(a,b)
      c3=y0-yp(a,b)
      c4=dy1da(a,b)
      c5=dy1db(a,b)
      c6=y01-y1p(a,b)
      d=c1*c5-c2*c4
      da=(c3*c5-c2*c6)/d
      db=(c1*c6-c3*c4)/d
      nit=nit+1
      a=a+da
      b=b+db
      ! write(*,*)' nit=',nit,' da=',da,' db=',db
      ! write(*,*)'     a =',a,'   b=',b
      if(abs(da).gt.1.e-5.or.abs(db).gt.1.e-5) goto 1
      aset=a
      bset=b
      iset=1
2     continue
      x2=x*x
      if(i.eq.0)rrt=x*(1.+aset*(1.-x2)*(bset-x2))
      if(i.eq.1)rrt=1.+aset*(bset-3.*x2*(1.+bset)+5.*x2*x2)
      return
      end
*
      subroutine calc_amp(u,v,w, amp, Imax, Jmax)
        implicit real*8 (a-h,o-z)
        dimension
     > u(0:Imax,0:Jmax,0:*)
     >,v(0:Imax,0:Jmax,0:*)
     >,w(0:Imax,0:Jmax,0:*)
       common
     >/dimx/hx,Im,lx
     >/dimr/rt(0:128),rt1(0:128),yt(129),yt1(129),hr,Jm
     >/dimt/ht,Km,lt
        amp = 0.d0
        Ss = 0.d0
        do j=1,Jm
          u0=0.
          Ss=Ss+yt(j)*yt1(j)
          do k=1,Km
            do i=1,Im
              u0 = u0 + u(i,j,k)
            end do
          end do
          u0 = u0/(Im*Km)
          do k = 1, Km
            do i = 1, Im
              amp = amp + 
     >          ((u(i,j,k)-u0)**2+w(i,j,k)**2+v(i,j,k)**2)*yt(j)*yt1(j)
            end do
          end do
        end do
        amp = sqrt(amp/(Im*Ss*Km))
        end
*
      subroutine correct_amp(u,v,w, amp_prev, amp_cur, Imax, Jmax)
      implicit real*8 (a-h,o-z)
      dimension
     > u(0:Imax,0:Jmax,0:*),
     > v(0:Imax,0:Jmax,0:*),
     > w(0:Imax,0:Jmax,0:*)
      common
     >/dimx/hx,Im,lx
     >/dimr/rt(0:128),rt1(0:128),yt(129),yt1(129),hr,Jm
     >/dimt/ht,Km,lt
        c = amp_cur / amp_prev
        do j=1,Jm
          u0=0.
          do k=1,Km
            do i=1,Im
              u0 = u0 + u(i,j,k)
            end do
          end do
          u0 = u0/(Im*Km)
          do k = 1, Km
            do i = 1, Im
              d = u(i,j,k) - u0
              u(i,j,k) = u0 + c * d
              v(i,j,k) = c * v(i,j,k)
              w(i,j,k) = c * w(i,j,k)
            end do
          end do
        end do
        end
*
      subroutine ubulk(u, Imax, Jmax, ub)
      implicit real*8 (a-h,o-z)
      dimension
     > u(0:Imax,0:Jmax,0:*)
       common
     >/dimx/hx,Im,lx
     >/dimr/rt(0:128),rt1(0:128),yt(129),yt1(129),hr,Jm
     >/dimt/ht,Km,lt
      ub=0.
      Ss = 0.d0
      do j=1,Jm
        Ss=Ss+yt(j)*yt1(j)
        do k=1,Km
          ub=ub+u(1,j,k)*yt(j)*yt1(j)
        end do
      end do
      ub = ub/(Ss*Km)
          end
*