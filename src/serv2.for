*
      subroutine servis(t,u,v,w,ox,oy,oz,p,ii,Imax,Jmax)
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      dimension
     > u(0:Imax,0:Jmax,0:*)
     >,v(0:Imax,0:Jmax,0:*)
     >,w(0:Imax,0:Jmax,0:*)
     >,p(0:Imax,0:Jmax,0:*)
     >,ox(0:Imax,0:Jmax,0:*)
     >,oy(0:Imax,0:Jmax,0:*)
     >,oz(0:Imax,0:Jmax,0:*)
     >,ums(128),u2s(128),v2s(128),w2s(128),uvs(128)
     >,ozms(128),ox2s(128),oy2s(128),oz2s(128)
      common
     >/dim/Xmax,epsy
     >/dimx/hx,Im,Imm,lx
     >/dimr/rt(0:128),rt1(0:128),yt(129),yt1(129),hr,Jm
     >/dimt/hz,Km,lz
     >/Re/Re
     >/proc/Np,Npm
*
     >/servst/iserv
     >/serv/tend,dt
     >,um(128),u2(128),v2(128),w2(128),uv(128)
     >,ozm(128),ox2(128),oy2(128),oz2(128)
      cxz=1./(Im*Km)
*
      if(iserv.eq.0) then
        iserv=1
        tend=t
        dt=0.
        do j=1,Jm
          um(j)=0.
          u2(j)=0.
          v2(j)=0.
          w2(j)=0.
          uv(j)=0.
          ozm(j)=0.
          ox2(j)=0.
          oy2(j)=0.
          oz2(j)=0.
        end do
        return
      end if
*
      if(ii.eq.0) then
        tau=t-tend
        dt0=dt
        dt=dt0+tau
        tend=t
        if(dt.eq.0.)return
        do j=1,Jm
          uu=0.
          uu2=0.
          vv2=0.
          ww2=0.
          uuv=0.
          ooz=0.
          oox2=0.
          ooy2=0.
          ooz2=0.
          do k=1,Km
            do i=1,Im
              uu=uu+u(i,j,k)
              uu2=uu2+u(i,j,k)**2
              vv2=vv2+v(i,j,k)**2
              ww2=ww2+w(i,j,k)**2
              uuv=uuv+0.25*(u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))
              ooz=ooz+oz(i,j,k)
              oox2=oox2+ox(i,j,k)**2
              ooy2=ooy2+oy(i,j,k)**2
              ooz2=ooz2+oz(i,j,k)**2
            end do
          end do
          um(j)=(dt0*um(j)+tau*uu*cxz)/dt
          u2(j)=(dt0*u2(j)+tau*uu2*cxz)/dt
          v2(j)=(dt0*v2(j)+tau*vv2*cxz)/dt
          w2(j)=(dt0*w2(j)+tau*ww2*cxz)/dt
          uv(j)=(dt0*uv(j)+tau*uuv*cxz)/dt
          ozm(j)=(dt0*ozm(j)+tau*ooz*cxz)/dt
          ox2(j)=(dt0*ox2(j)+tau*oox2*cxz)/dt
          oy2(j)=(dt0*oy2(j)+tau*ooy2*cxz)/dt
          oz2(j)=(dt0*oz2(j)+tau*ooz2*cxz)/dt
        end do
        return
      end if
*
      if(ii.eq.1) return
*
      if(ii.eq.2) then
        call MPI_REDUCE(um(0),ums(0),Jm+1,MPI_DOUBLE_PRECISION,MPI_SUM
     >                 ,0,MPI_COMM_WORLD,ier)
        call MPI_REDUCE(u2(0),u2s(0),Jm+1,MPI_DOUBLE_PRECISION,MPI_SUM
     >                 ,0,MPI_COMM_WORLD,ier)
        call MPI_REDUCE(v2(0),v2s(0),Jm+1,MPI_DOUBLE_PRECISION,MPI_SUM
     >                 ,0,MPI_COMM_WORLD,ier)
        call MPI_REDUCE(w2(0),w2s(0),Jm+1,MPI_DOUBLE_PRECISION,MPI_SUM
     >                 ,0,MPI_COMM_WORLD,ier)
        call MPI_REDUCE(uv(0),uvs(0),Jm+1,MPI_DOUBLE_PRECISION,MPI_SUM
     >                 ,0,MPI_COMM_WORLD,ier)
        call MPI_REDUCE(ozm(0),ozms(0),Jm+1,MPI_DOUBLE_PRECISION,MPI_SUM
     >                 ,0,MPI_COMM_WORLD,ier)
        call MPI_REDUCE(ox2(0),ox2s(0),Jm+1,MPI_DOUBLE_PRECISION,MPI_SUM
     >                 ,0,MPI_COMM_WORLD,ier)
        call MPI_REDUCE(oy2(0),oy2s(0),Jm+1,MPI_DOUBLE_PRECISION,MPI_SUM
     >                 ,0,MPI_COMM_WORLD,ier)
        call MPI_REDUCE(oz2(0),oz2s(0),Jm+1,MPI_DOUBLE_PRECISION,MPI_SUM
     >                 ,0,MPI_COMM_WORLD,ier)
        call MPI_BARRIER(MPI_COMM_WORLD,ier)
        if(Np.eq.0) then
	    open(10,file='mean.dat',form='unformatted',status='old',err=1)
          read(10)ttend,dtt,Re0,Xmax0,epsy0,lx0,Jm0,lz0
          read(10)(um(j),j=1,Jm0)
          read(10)(u2(j),j=1,Jm0)
          read(10)(v2(j),j=1,Jm0)
          read(10)(w2(j),j=1,Jm0)
          read(10)(uv(j),j=1,Jm0)
          read(10)(ozm(j),j=1,Jm0)
          read(10)(ox2(j),j=1,Jm0)
          read(10)(oy2(j),j=1,Jm0)
          read(10)(oz2(j),j=1,Jm0)
          close(10)
	    goto 2
1         continue
          dtt=0.
          do j=1,Jm
            um(j)=0.
            u2(j)=0.
            v2(j)=0.
            w2(j)=0.
            uv(j)=0.
            ozm(j)=0.
            ox2(j)=0.
            oy2(j)=0.
            oz2(j)=0.
          end do
2         continue
          dtt0=dtt
	    dtt=dtt+dt
	    do j=1,Jm
            um(j)=(dtt0*um(j)+dt*ums(j)/Npm)/dtt
            u2(j)=(dtt0*u2(j)+dt*u2s(j)/Npm)/dtt
            v2(j)=(dtt0*v2(j)+dt*v2s(j)/Npm)/dtt
            w2(j)=(dtt0*w2(j)+dt*w2s(j)/Npm)/dtt
            uv(j)=(dtt0*uv(j)+dt*uvs(j)/Npm)/dtt
            ozm(j)=(dtt0*ozm(j)+dt*ozms(j)/Npm)/dtt
            ox2(j)=(dtt0*ox2(j)+dt*ox2s(j)/Npm)/dtt
            oy2(j)=(dtt0*oy2(j)+dt*oy2s(j)/Npm)/dtt
            oz2(j)=(dtt0*oz2(j)+dt*oz2s(j)/Npm)/dtt
	    end do 
	    open(10,file='mean.dat',form='unformatted')
          write(10)tend,dtt,Re,Xmax,epsy,lx,Jm,lz
          write(10)(um(j),j=1,Jm)
          write(10)(u2(j),j=1,Jm)
          write(10)(v2(j),j=1,Jm)
          write(10)(w2(j),j=1,Jm)
          write(10)(uv(j),j=1,Jm)
          write(10)(ozm(j),j=1,Jm)
          write(10)(ox2(j),j=1,Jm)
          write(10)(oy2(j),j=1,Jm)
          write(10)(oz2(j),j=1,Jm)
          close(10)
        end if
        dt=0.
        do j=1,Jm
          um(j)=0.
          u2(j)=0.
          v2(j)=0.
          w2(j)=0.
          uv(j)=0.
          ozm(j)=0.
          ox2(j)=0.
          oy2(j)=0.
          oz2(j)=0.
        end do
        return
      end if
*
      end