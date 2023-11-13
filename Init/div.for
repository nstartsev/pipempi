      subroutine div(i,j,k,u,v,w,d,Imax,Jmax)
          implicit real*8 (a-h,o-z)
          dimension
     > u(0:Imax,0:Jmax,0:*)
     >,v(0:Imax,0:Jmax,0:*)
     >,w(0:Imax,0:Jmax,0:*)
          common
     >/dimx/hx,Im,lx
     >/dimr/rt(0:128),rt1(0:128),yt(129),yt1(129),hr,Jm
     >/dimt/ht,Km,lt
          i1=mod(Im+i-2,Im)+1
          k1=mod(Km+k-2,Km)+1
          d=(u(i,j,k)-u(i1,j,k))/hx
     > +(rt(j)*v(i,j,k)-rt(j-1)*v(i,j-1,k))/(yt(j)*yt1(j))
     > +(w(i,j,k)-w(i,j,k1))/(yt(j)*ht)
          return
          end