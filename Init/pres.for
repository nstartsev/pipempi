*
      subroutine pres(u,v,w,p,Imax,Jmax)
        implicit real*8 (a-h,o-z)
        dimension
     > u(0:Imax,0:Jmax,0:*)
     >,v(0:Imax,0:Jmax,0:*)
     >,w(0:Imax,0:Jmax,0:*)
     >,p(0:Imax,0:Jmax,0:*)
        common
     > a1(256),a2(256),a3(256),a4(256)
     >,b1(256),b2(256),b3(256),b4(256)
     >,c1(256),c2(256),c3(256),c4(256)
     >,ap(256),bp(256),cp(256)
     >/dimx/hx,Im,lx
     >/dimr/rt(0:128),rt1(0:128),yt(129),yt1(129),hr,Jm
     >/dimt/ht,Km,lt
     >/rl/rlx(0:128),rlt(0:128)
     >/pry/apy(128),bpy(128),cpy(128)
*
      Im2=Im/2
      Km2=Km/2
      cik=4./(Im2*Km2)
*
      do k=1,Km
        do i=1,Im
          v(i,Jm,k)=0.d0
          do j=1,Jm
            call div(i,j,k,u,v,w,d,Imax,Jmax)
            p(i,j,k)=d
          end do
        end do
      end do
*
      do j=1,Jm
*   Symmetrization
        do k=1,Km2
          k1=Km+1-k
          do i=1,Im2
            i1=Im+1-i
            p1=p(i,j,k)
            p2=p(i1,j,k)
            p3=p(i,j,k1)
            p4=p(i1,j,k1)
            p(i,j,k)=0.25d0*(p1+p2+p3+p4)
            p(i1,j,k)=0.25d0*(p1-p2+p3-p4)
            p(i,j,k1)=0.25d0*(p1+p2-p3-p4)
            p(i1,j,k1)=0.25d0*(p1-p2-p3+p4)
          end do
        end do
*   Direct Fourier transformation
        do k=1,Km2
          k1=Km+1-k
          do i=1,Im2
            i1=Im+1-i
            a1(i)=p(i,j,k)
            a2(i)=p(i1,j,k)
            a3(i)=p(i,j,k1)
            a4(i)=p(i1,j,k1)
          end do
          call ftc05d(a1,b1,lx-1)
          call fts05d(a2,b2,lx-1)
          call ftc05d(a3,b3,lx-1)
          call fts05d(a4,b4,lx-1)
          do i=1,Im2
            i1=Im+1-i
            p(i,j,k)=b1(i)
            p(i1,j,k)=b2(i)
            p(i,j,k1)=b3(i)
            p(i1,j,k1)=b4(i)
          end do
        end do
        do i=1,Im2
          i1=Im+1-i
          do k=1,Km2
            k1=Km+1-k
            a1(k)=p(i,j,k)
            a2(k)=p(i1,j,k)
            a3(k)=p(i,j,k1)
            a4(k)=p(i1,j,k1)
          end do
          call ftc05d(a1,b1,lt-1)
          call ftc05d(a2,b2,lt-1)
          call fts05d(a3,b3,lt-1)
          call fts05d(a4,b4,lt-1)
          do k=1,Km2
            k1=Km+1-k
            p(i,j,k)=b1(k)*cik
            p(i1,j,k)=b2(k)*cik
            p(i,j,k1)=b3(k)*cik
            p(i1,j,k1)=b4(k)*cik
          end do
        end do
      end do
*   Solution in normal co-ordinate
      do k=1,Km2
        k1=Km+1-k
        do i=1,Im2
          i1=Im+1-i
          do j=1,Jm
            a1(j)=p(i,j,k)
            a2(j)=p(i1,j,k)
            a3(j)=p(i,j,k1)
            a4(j)=p(i1,j,k1)
            b1(j)=bpy(j)-rlx(i-1)-rlt(k-1)/yt(j)**2
            b2(j)=bpy(j)-rlx(i)-rlt(k-1)/yt(j)**2
            b3(j)=bpy(j)-rlx(i-1)-rlt(k)/yt(j)**2
            b4(j)=bpy(j)-rlx(i)-rlt(k)/yt(j)**2
            ap(j)=apy(j)
            cp(j)=cpy(j)
          end do
          b1(Jm)=b1(Jm)+cpy(Jm)
          b2(Jm)=b2(Jm)+cpy(Jm)
          b3(Jm)=b3(Jm)+cpy(Jm)
          b4(Jm)=b4(Jm)+cpy(Jm)
          Jm1=Jm
          if(i.eq.1.and.k.eq.1) then
            Jm1=Jm-1
            c1(Jm)=0.d0
          end if
          call prog3(ap,b1,cp,a1,c1,Jm1)
          call prog3(ap,b2,cp,a2,c2,Jm)
          call prog3(ap,b3,cp,a3,c3,Jm)
          call prog3(ap,b4,cp,a4,c4,Jm)
          do j=1,Jm
            p(i,j,k)=c1(j)
            p(i1,j,k)=c2(j)
            p(i,j,k1)=c3(j)
            p(i1,j,k1)=c4(j)
          end do
        end do
      end do
*
      do j=1,Jm
*   Inverse Fourier transformation
        do k=1,Km2
          k1=Km+1-k
          do i=1,Im2
            i1=Im+1-i
            a1(i)=p(i,j,k)
            a2(i)=p(i1,j,k)
            a3(i)=p(i,j,k1)
            a4(i)=p(i1,j,k1)
          end do
          call ftc05b(a1,b1,lx-1)
          call fts05b(a2,b2,lx-1)
          call ftc05b(a3,b3,lx-1)
          call fts05b(a4,b4,lx-1)
          do i=1,Im2
            i1=Im+1-i
            p(i,j,k)=b1(i)
            p(i1,j,k)=b2(i)
            p(i,j,k1)=b3(i)
            p(i1,j,k1)=b4(i)
          end do
        end do
        do i=1,Im2
          i1=Im+1-i
          do k=1,Km2
            k1=Km+1-k
            a1(k)=p(i,j,k)
            a2(k)=p(i1,j,k)
            a3(k)=p(i,j,k1)
            a4(k)=p(i1,j,k1)
          end do
          call ftc05b(a1,b1,lt-1)
          call ftc05b(a2,b2,lt-1)
          call fts05b(a3,b3,lt-1)
          call fts05b(a4,b4,lt-1)
          do k=1,Km2
            k1=Km+1-k
            p(i,j,k)=b1(k)
            p(i1,j,k)=b2(k)
            p(i,j,k1)=b3(k)
            p(i1,j,k1)=b4(k)
          end do
        end do
*   Desymmetrization
        do k=1,Km2
          k1=Km+1-k
          do i=1,Im2
            i1=Im+1-i
            p1=p(i,j,k)
            p2=p(i1,j,k)
            p3=p(i,j,k1)
            p4=p(i1,j,k1)
            p(i,j,k)=p1+p2+p3+p4
            p(i1,j,k)=p1-p2+p3-p4
            p(i,j,k1)=p1+p2-p3-p4
            p(i1,j,k1)=p1-p2-p3+p4
          end do
        end do
      end do
*
*  Mean pressure gradient
      Ub=p(0,0,1)
      ss=0.d0
      su=0.d0
      do j=1,Jm
        ss=ss+yt(j)*yt1(j)
        ssu=0.d0
        do k=1,Km
          do i=1,Im
            ssu=ssu+u(i,j,k)
          end do
        end do
        su=su+ssu*yt(j)*yt1(j)
      end do
      Dp=Ub-su/(Im*Km*ss)
      p(0,0,0)=Dp
*
      call gradp(u,v,w,p,Imax,Jmax)
      return
      end