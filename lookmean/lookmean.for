      program lookmean
      implicit real*8 (a-h,o-z)
      parameter (Imax=2049, Jmax=129, Kmax=129)
      dimension
     > u(0:Imax,0:Jmax,0:Kmax)
     >,v(0:Imax,0:Jmax,0:Kmax)
     >,w(0:Imax,0:Jmax,0:Kmax)
     >,p(0:Imax,0:Jmax,0:Kmax)
     >,ox(0:Imax,0:Jmax,0:Kmax)
     >,oy(0:Imax,0:Jmax,0:Kmax)
     >,oz(0:Imax,0:Jmax,0:Kmax)
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

      open(10,file='mean.dat',form='unformatted',status='old',err=1)
      read(10)tend,dtt,Re,Xmax,epsy,lx,Jm,lz   
      call com

      Im = 2**(lx)
      Km = 2**(lz)
      write(*,200) tend,dtt,Re,Xmax,epsy,Im,Jm,Km
      read(10)(um(j),j=1,Jm)
      read(10)(u2(j),j=1,Jm)
      read(10)(v2(j),j=1,Jm)
      read(10)(w2(j),j=1,Jm)
      read(10)(uv(j),j=1,Jm)
      read(10)(ozm(j),j=1,Jm)
      read(10)(ox2(j),j=1,Jm)
      read(10)(oy2(j),j=1,Jm)
      read(10)(oz2(j),j=1,Jm)
      close(10)
      
      open(11, file='um.dat')
      do j=1, Jm
            write(11,100) rt(j), um(j)
      end do
      close(11)
      open(11, file='ozm.dat')
      do j=1, Jm
            write(11,100) rt(j), ozm(j)
      end do
      close(11)
      
*
      goto 1
200   format('    t=',1pe10.3,' dt=',e9.2,/,
     >'    Re=',e9.2,/,
     >'    Xmax=',e9.2,/,
     >'    epsr=',e9.2,' Im=',i4,' Jm=',i4,' Km=',i4)
101   format(10(i4))
100   format(10(1pe12.4))
1     end