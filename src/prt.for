*
      subroutine prt(t,dt,u,v,w,p,Imax,Jmax)
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      dimension
     > u(0:Imax,0:Jmax,0:*)
     >,v(0:Imax,0:Jmax,0:*)
     >,w(0:Imax,0:Jmax,0:*)
     >,p(0:Imax,0:Jmax,0:*)
     >,u0(0:Jmax),u0s(0:Jmax),u0mpi(0:Jmax)
      common
     >/dimx/hx,Im,Imm,lx
     >/dimr/rt(0:128),rt1(0:128),yt(129),yt1(129),hr,Jm
     >/dimt/ht,Km,lt
     >/Re/Re
     >/proc/Np,Npm
*

* compute u0n on each node
      do j=1,Jm
        u0(j) = 0.d0
        u0s(j) = 0.d0
      end do
      do j=1,Jm
        do k=1,Km
          do i=1,Im
            u0s(j)=u0s(j)+u(i,j,k)
          end do
        end do
        u0s(j)=u0s(j)/(Im*Npm*Km)
      end do
*
      call MPI_REDUCE(u0s(1),u0(1),Jm,MPI_DOUBLE_PRECISION,MPI_SUM,0
     >               ,MPI_COMM_WORLD,ier)
      Dp=p(0,0,0)
      ubulk=0.
      Ss=0.
      amp=0.
      dd=0.
*
      do j=1,Jm
        Ss=Ss+yt(j)*yt1(j)
        ! u0=0.
        do k=1,Km
          ubulk=ubulk+u(1,j,k)*yt(j)*yt1(j)
          ! do i=1,Im
          !   u0=u0+u(i,j,k)
          ! end do
        end do
        ! u0=u0/(Im*Km)
        do k=1,Km
          do i=1,Im
            amp=amp
     >       +((u(i,j,k)-u0(j))**2+w(i,j,k)**2+v(i,j,k)**2)*yt(j)*yt1(j)
            call div(i,j,k,u,v,w,d,Imax,Jmax)
            dd=max(dd,abs(d))
          end do
        end do
      end do
      ubulk=ubulk/(Ss*Km)
      amp=amp/(Ss*Im*Km)
      ucl=0.
      do k=1,Km
        do i=1,Im
          ucl=ucl+u(i,1,k)
        end do
      end do
      ucl=ucl/(Im*Km)
*
      call MPI_REDUCE(Dp,Dps,1,MPI_DOUBLE_PRECISION,MPI_SUM,0
     >               ,MPI_COMM_WORLD,ier)
      call MPI_REDUCE(amp,amps,1,MPI_DOUBLE_PRECISION,MPI_SUM,0
     >               ,MPI_COMM_WORLD,ier)
      call MPI_REDUCE(ucl,ucls,1,MPI_DOUBLE_PRECISION,MPI_SUM,0
     >               ,MPI_COMM_WORLD,ier)
      call MPI_REDUCE(dd,dds,1,MPI_DOUBLE_PRECISION,MPI_MAX,0
     >               ,MPI_COMM_WORLD,ier)
*
      if(Np.eq.0)then
        Dp=Dps/Npm
        amp=sqrt(amps/Npm)
        ucl=ucls/Npm
        i=1
        k=1
        j1=1
        j2=Jm/4
        j3=3*Jm/4
        u1=u(i,j1,k)
        u2=u(i,j2,k)
        u3=u(i,j3,k)
        v1=v(i,j1,k)
        v2=v(i,j2,k)
        v3=v(i,j3,k)
        w1=w(i,j1,k)
        w2=w(i,j2,k)
        w3=w(i,j3,k)
        write(8,120)t,dt,Dp,amp,ucl,dds,u1,u2,u3,v1,v2,v3,w1,w2,w3
*        write(8,120)t,dt,Dp,amp,ucl,dds
        write(*,110)Dp,amp,ucl,dds
      end if 
120   format(1pe14.6,15e12.4)
110   format('    Dp =',1pe12.4,'  amp =',e12.4,'  Ucl =',e12.4,
     >       '  Div =',e12.4)
      return
      end