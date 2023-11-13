*
      subroutine lin(tau,u,v,w,Imax,Jmax)
      implicit real*8 (a-h,o-z)
      dimension
     > u(0:Imax,0:Jmax,0:*)
     >,v(0:Imax,0:Jmax,0:*)
     >,w(0:Imax,0:Jmax,0:*)
     >,ap(256),bp(256),cp(256),dp(256),ep(256)
      common
     >/dimx/hx,Im,Imm,lx
     >/dimr/rt(0:128),rt1(0:128),yt(129),yt1(129),hr,Jm
     >/dimt/ht,Km,lt
     >/pry/apy(128),bpy(128),cpy(128)
     >/prv/apv(128),bpv(128),cpv(128)
     >/prw/apw(128),bpw(128),cpw(128)
     >/Re/Re
      ct=Re/tau
*

*
* At
*   U,W
      do i=1,Im
        do j=1,Jm
          do k=1,Km
            ap(k)=1.d0/ht**2
            cp(k)=ap(k)
            bp(k)=-ap(k)-cp(k)-ct*yt(j)**2
            dp(k)=-ct*u(i,j,k)*yt(j)**2
          end do
          call prg3cycl(ap,bp,cp,dp,ep,Km)
          do k=1,Km
            u(i,j,k)=ep(k)
          end do
          do k=1,Km
            ap(k)=1.d0/ht**2
            cp(k)=ap(k)
            bp(k)=-ap(k)-cp(k)-ct*yt(j)**2
            dp(k)=-ct*w(i,j,k)*yt(j)**2
          end do
          call prg3cycl(ap,bp,cp,dp,ep,Km)
          do k=1,Km
            w(i,j,k)=ep(k)
          end do
            w(i,j,0)=w(i,j,Km)
        end do
      end do
*   V
      do i=1,Im
        do j=1,Jm-1
          do k=1,Km
            ap(k)=1.d0/ht**2
            cp(k)=ap(k)
            bp(k)=-ap(k)-cp(k)-ct*rt(j)**2
            wt0=(w(i,j,k)-w(i,j,k-1))/(yt(j)*ht)
            wt1=(w(i,j+1,k)-w(i,j+1,k-1))/(yt(j+1)*ht)
            f1=(wt1-wt0)/rt1(j)
            wr0=(yt(j+1)*w(i,j+1,k-1)-yt(j)*w(i,j,k-1))/(rt(j)*rt1(j))
            wr1=(yt(j+1)*w(i,j+1,k)-yt(j)*w(i,j,k))/(rt(j)*rt1(j))
            f2=(wr1-wr0)/(rt(j)*ht)
            dp(k)=(-ct*v(i,j,k)-f1+f2)*rt(j)**2
          end do
          call prg3cycl(ap,bp,cp,dp,ep,Km)
          do k=1,Km
            v(i,j,k)=ep(k)
          end do
        end do
      end do
*
* Ar
*   U,V
      do i=1,Im
        do k=1,Km
          do j=1,Jm
            bp(j)=bpy(j)-ct
            dp(j)=-ct*u(i,j,k)
          end do
          bp(Jm)=bp(Jm)-cpy(Jm)
          call prog3(apy,bp,cpy,dp,ep,Jm)
          do j=1,Jm
            u(i,j,k)=ep(j)
          end do
          do j=1,Jm-1
            bp(j)=bpv(j)-ct
            dp(j)=-ct*rt(j)*v(i,j,k)
          end do
          call prog3(apv,bp,cpv,dp,ep,Jm-1)
          do j=1,Jm-1
            v(i,j,k)=ep(j)/rt(j)
          end do
        end do
        do j=1,Jm
          v(i,j,Km+1)=v(i,j,1)
        end do
      end do
*   W
      do i=1,Im
*    Axisymmetric part
        do j=1,Jm
          ap(j)=0.d0
          do k=1,Km
            ap(j)=ap(j)+w(i,j,k)
          end do
          ap(j)=ap(j)/Km
          do k=1,Km
            w(i,j,k)=w(i,j,k)-ap(j)
          end do
          bp(j)=bpw(j)-ct
          dp(j)=-ct*ap(j)*yt(j)
        end do
        bp(1)=bp(1)-2.d0/(yt(1)*yt1(1))
        bp(Jm)=bp(Jm)-cpw(Jm)
        call prog3(apw,bp,cpw,dp,cp,Jm)
*   Unsymmetric part
        do k=1,Km
          j=1
            bp(j)=bpw(j)-ct
            vr0=(rt(j)*v(i,j,k)-rt(j-1)*v(i,j-1,k))/(yt(j)*yt1(j))
            vr1=(rt(j)*v(i,j,k+1)-rt(j-1)*v(i,j-1,k+1))/(yt(j)*yt1(j))
            f1=(vr1-vr0)/(yt(j)*ht)
            vt0=0.d0
            vt1=(v(i,j,k+1)-v(i,j,k))/(rt(j)*ht)
            f2=(vt1-vt0)/yt1(j)
            dp(j)=(-ct*w(i,j,k)-f1+f2)*yt(j)
          do j=2,Jm
            bp(j)=bpw(j)-ct
            vr0=(rt(j)*v(i,j,k)-rt(j-1)*v(i,j-1,k))/(yt(j)*yt1(j))
            vr1=(rt(j)*v(i,j,k+1)-rt(j-1)*v(i,j-1,k+1))/(yt(j)*yt1(j))
            f1=(vr1-vr0)/(yt(j)*ht)
            vt0=(v(i,j-1,k+1)-v(i,j-1,k))/(rt(j-1)*ht)
            vt1=(v(i,j,k+1)-v(i,j,k))/(rt(j)*ht)
            f2=(vt1-vt0)/yt1(j)
            dp(j)=(-ct*w(i,j,k)-f1+f2)*yt(j)
          end do
          bp(Jm)=bp(Jm)-cpw(Jm)
          call prog3(apw,bp,cpw,dp,ep,Jm)
          do j=1,Jm
            w(i,j,k)=(ep(j)+cp(j))/yt(j)
          end do
        end do
      end do
      return
      end