*
      block data bdftc
      common
     >/ftcst/ist
     >/ftnl/nl(13)
      data ist/0/
      data nl/1,2,4,8,16,32,64,128,256,512,1024,2048,4096/
      end
*
      subroutine stftc(l)
      implicit real*8 (a-h,o-z)
      common
     >/ftcst/ist
     >/ftco4/e4(15),e3(7),e2(3)
     >/ftco12/e12(4095)
     >/ftnl/nl(13)
      ist=l
      n=nl(l+1)
      pi=atan(1.d0)/2048.
      jd=nl(13-l)
      do j=jd,4095,jd
        e12(j)=cos(j*pi)
      end do
      do j=1,15
        e4(j)=e12(256*j)
      end do
      do j=1,7
        e3(j)=e12(512*j)
      end do
      do j=1,3
        e2(j)=e12(1024*j)
      end do
      return
      end
*
      subroutine ftc05d(x,y,l)
      implicit real*8 (a-h,o-z)
      dimension x(*),y(*)
      common
     >/ftcst/ist
     >/ftco4/e4(15),e3(7),e2(3)
     >/ftco12/e12(4095)
      if(l.lt.1.or.l.gt.12) then
        write(*,*)'  ***** fftcos:  l is out of available range *****'
        stop
      end if
      if(ist.lt.l) call stftc(l)
      goto(10,20,30,40,50,60,70,80,90,100,110,120) l
10    call ftcod1(x,y)
      goto 130
20    call ftcod2(x,y)
      goto 130
30    call ftcod3(x,y)
      goto 130
40    call ftcod4(x,y)
      goto 130
50    call ftcod5(x,y)
      goto 130
60    call ftcod6(x,y)
      goto 130
70    call ftcod7(x,y)
      goto 130
80    call ftcod8(x,y)
      goto 130
90    call ftcod9(x,y)
      goto 130
100   call ftcod10(x,y)
      goto 130
110   call ftcod11(x,y)
      goto 130
120   call ftcod12(x,y)
      goto 130
130   continue
      return
      end
*
      subroutine ftc05b(x,y,l)
      implicit real*8 (a-h,o-z)
      dimension x(*),y(*)
      common
     >/ftcst/ist
     >/ftco4/e4(15),e3(7),e2(3)
     >/ftco12/e12(4095)
      if(l.lt.1.or.l.gt.12) then
        write(*,*)'  ***** fftcos:  l is out of available range *****'
        stop
      end if
      if(ist.lt.l) call stftc(l)
      goto(10,20,30,40,50,60,70,80,90,100,110,120) l
10    call ftcob1(x,y)
      goto 130
20    call ftcob2(x,y)
      goto 130
30    call ftcob3(x,y)
      goto 130
40    call ftcob4(x,y)
      goto 130
50    call ftcob5(x,y)
      goto 130
60    call ftcob6(x,y)
      goto 130
70    call ftcob7(x,y)
      goto 130
80    call ftcob8(x,y)
      goto 130
90    call ftcob9(x,y)
      goto 130
100   call ftcob10(x,y)
      goto 130
110   call ftcob11(x,y)
      goto 130
120   call ftcob12(x,y)
      goto 130
130   continue
      return
      end
*
      subroutine fts05d(x,y,l)
      implicit real*8 (a-h,o-z)
      dimension x(*),y(*)
      common
     >/ftnl/nl(13)
      N=nl(l+1)
      do j=1,N,2
        y(j)=x(j)
        y(j+1)=-x(j+1)
      end do
      call ftc05d(y,x,l)
      do k=1,N
        y(k)=x(N+1-k)
      end do
      return
      end
*
      subroutine fts05b(x,y,l)
      implicit real*8 (a-h,o-z)
      dimension x(*),y(*)
      common
     >/ftnl/nl(13)
      N=nl(l+1)
      do k=1,N
        y(k)=x(N+1-k)
      end do
      call ftc05b(y,x,l)
      do j=1,N,2
        y(j)=x(j)
        y(j+1)=-x(j+1)
      end do
      return
      end
*
      subroutine ftc(x,y,l)
      implicit real*8 (a-h,o-z)
      dimension x(*),y(*)
      common
     >/ftcst/ist
      if(l.lt.1.or.l.gt.12) then
        write(*,*)'  ***** ftc:  l is out of available range *****'
        stop
      end if
      if(ist.lt.l) call stftc(l)
      goto(10,20,30,40,50,60,70,80,90,100,110,120) l
10    call ftc1(x,y)
      goto 130
20    call ftc2(x,y)
      goto 130
30    call ftc3(x,y)
      goto 130
40    call ftc4(x,y)
      goto 130
50    call ftc5(x,y)
      goto 130
60    call ftc6(x,y)
      goto 130
70    call ftc7(x,y)
      goto 130
80    call ftc8(x,y)
      goto 130
90    call ftc9(x,y)
      goto 130
100   call ftc10(x,y)
      goto 130
110   call ftc11(x,y)
      goto 130
120   call ftc12(x,y)
      goto 130
130   continue
      return
      end
*
      subroutine fts(x,y,l)
      implicit real*8 (a-h,o-z)
      dimension x(*),y(*)
      common
     >/ftcst/ist
      if(l.lt.1.or.l.gt.12) then
        write(*,*)'  ***** fts:  l is out of available range *****'
        stop
      end if
      if(ist.lt.l) call stftc(l)
      goto(10,20,30,40,50,60,70,80,90,100,110,120) l
10    call fts1(x,y)
      goto 130
20    call fts2(x,y)
      goto 130
30    call fts3(x,y)
      goto 130
40    call fts4(x,y)
      goto 130
50    call fts5(x,y)
      goto 130
60    call fts6(x,y)
      goto 130
70    call fts7(x,y)
      goto 130
80    call fts8(x,y)
      goto 130
90    call fts9(x,y)
      goto 130
100   call fts10(x,y)
      goto 130
110   call fts11(x,y)
      goto 130
120   call fts12(x,y)
      goto 130
130   continue
      return
      end
*
      subroutine ftc025(x,y,l)
      implicit real*8 (a-h,o-z)
      dimension x(*),y(*)
      common
     >/ftcst/ist
     >/ftco12/e12(4095)
     >/ftnl/nl(13)
      if(l.lt.1.or.l.gt.11) then
        write(*,*)'  ***** ftc025:  l is out of available range *****'
        stop
      end if
      if(ist.lt.l+1) call stftc(l+1)
      n=nl(l+1)
      n0=nl(12-l)
      nd=n0+n0
      do j=1,n
        y(j)=x(j)*e12(j*nd-n0)
      end do
      call ftc05d(y,x,l)
      y(1)=x(1)
      do k=2,n
        y(k)=x(k)+x(k)-y(k-1)
      end do
      return
      end
*
      subroutine ftcod1(x,y)
      implicit real*8 (a-h,o-z)
      dimension x(2),y(2)
      common/ftco4/e4(15),e3(7),e2(3)
*
      c=e2(2)
      y(1)= x(1)+x(2)
      y(2)=(x(1)-x(2))*c
*
      return
      end
*
      subroutine ftcod2(x,y)
      implicit real*8 (a-h,o-z)
      dimension x(4),y(4)
      common/ftco4/e4(15),e3(7),e2(3)
*  Step N1
      y(1) = x(1)+x(4)
      y(2) =(x(1)-x(4))*e2(1)
      y(3) = x(2)+x(3)
      y(4) =(x(2)-x(3))*e2(3)
*  Step N2
      c=e2(2)
      x(3)=(y(1)-y(3))*c
      x(4)=(y(2)-y(4))*c
      y(1)= y(1)+y(3)
      y(2)= y(2)+y(4)
*  Step N5
      y(4)=x(4)+x(4)-y(2)
      y(3)=x(3)
*
      return
      end
*
      subroutine ftcod3(x,y)
      implicit real*8 (a-h,o-z)
      dimension x(8),y(8)
      common/ftco4/e4(15),e3(7),e2(3)
*  Step N1
      y(1) = x(1)+x(8)
      y(2) =(x(1)-x(8))*e3(1)
      y(3) = x(2)+x(7)
      y(4) =(x(2)-x(7))*e3(3)
      y(5) = x(3)+x(6)
      y(6) =(x(3)-x(6))*e3(5)
      y(7) = x(4)+x(5)
      y(8) =(x(4)-x(5))*e3(7)
*  Step N2
      x(1) = y(1)+y(7)
      x(2) = y(2)+y(8)
      x(3) =(y(1)-y(7))*e3(2)
      x(4) =(y(2)-y(8))*e3(2)
      x(5) = y(3)+y(5)
      x(6) = y(4)+y(6)
      x(7) =(y(3)-y(5))*e3(6)
      x(8) =(y(4)-y(6))*e3(6)
*  Step N3
      c=e3(4)
      y(1)= x(1)+x(5)
      y(2)= x(2)+x(6)
      y(3)= x(3)+x(7)
      y(4)= x(4)+x(8)
      y(5)=(x(1)-x(5))*c
      y(6)=(x(2)-x(6))*c
      y(7)=(x(3)-x(7))*c
      y(8)=(x(4)-x(8))*c
*  Step N4
      y(7)=y(7)+y(7)-y(3)
      y(8)=y(8)+y(8)-y(4)
*  Step N5
      y(4)=y(4)+y(4)-y(2)
      y(6)=y(6)+y(6)-y(4)
      y(8)=y(8)+y(8)-y(6)
*
      return
      end
*
      subroutine ftcod4(x,y)
      implicit real*8 (a-h,o-z)
      dimension x(16),y(16)
      common/ftco4/e4(15),e3(7),e2(3)
*  Step N1
      y( 1) = x( 1)+x(16)
      y( 2) =(x( 1)-x(16))*e4( 1)
      y( 3) = x( 2)+x(15)
      y( 4) =(x( 2)-x(15))*e4( 3)
      y( 5) = x( 3)+x(14)
      y( 6) =(x( 3)-x(14))*e4( 5)
      y( 7) = x( 4)+x(13)
      y( 8) =(x( 4)-x(13))*e4( 7)
      y( 9) = x( 5)+x(12)
      y(10) =(x( 5)-x(12))*e4( 9)
      y(11) = x( 6)+x(11)
      y(12) =(x( 6)-x(11))*e4(11)
      y(13) = x( 7)+x(10)
      y(14) =(x( 7)-x(10))*e4(13)
      y(15) = x( 8)+x( 9)
      y(16) =(x( 8)-x( 9))*e4(15)
*  Step N2
      x( 1) = y( 1)+y(15)
      x( 2) = y( 2)+y(16)
      x( 3) =(y( 1)-y(15))*e4( 2)
      x( 4) =(y( 2)-y(16))*e4( 2)
      x( 5) = y( 3)+y(13)
      x( 6) = y( 4)+y(14)
      x( 7) =(y( 3)-y(13))*e4( 6)
      x( 8) =(y( 4)-y(14))*e4( 6)
      x( 9) = y( 5)+y(11)
      x(10) = y( 6)+y(12)
      x(11) =(y( 5)-y(11))*e4(10)
      x(12) =(y( 6)-y(12))*e4(10)
      x(13) = y( 7)+y( 9)
      x(14) = y( 8)+y(10)
      x(15) =(y( 7)-y( 9))*e4(14)
      x(16) =(y( 8)-y(10))*e4(14)
*  Step N3
      y( 1) = x( 1)+x(13)
      y( 2) = x( 2)+x(14)
      y( 3) = x( 3)+x(15)
      y( 4) = x( 4)+x(16)
      y( 5) =(x( 1)-x(13))*e4( 4)
      y( 6) =(x( 2)-x(14))*e4( 4)
      y( 7) =(x( 3)-x(15))*e4( 4)
      y( 8) =(x( 4)-x(16))*e4( 4)
      y( 9) = x( 5)+x( 9)
      y(10) = x( 6)+x(10)
      y(11) = x( 7)+x(11)
      y(12) = x( 8)+x(12)
      y(13) =(x( 5)-x( 9))*e4(12)
      y(14) =(x( 6)-x(10))*e4(12)
      y(15) =(x( 7)-x(11))*e4(12)
      y(16) =(x( 8)-x(12))*e4(12)
*  Step N4
      c=e4( 8)
      x( 1)= y( 1)+y( 9)
      x( 2)= y( 2)+y(10)
      x( 3)= y( 3)+y(11)
      x( 4)= y( 4)+y(12)
      x( 5)= y( 5)+y(13)
      x( 6)= y( 6)+y(14)
      x( 7)= y( 7)+y(15)
      x( 8)= y( 8)+y(16)
      x( 9)=(y( 1)-y( 9))*c
      x(10)=(y( 2)-y(10))*c
      x(11)=(y( 3)-y(11))*c
      x(12)=(y( 4)-y(12))*c
      x(13)=(y( 5)-y(13))*c
      x(14)=(y( 6)-y(14))*c
      x(15)=(y( 7)-y(15))*c
      x(16)=(y( 8)-y(16))*c
*  Step N5
      y(13)=x(13)+x(13)-x( 5)
      y(14)=x(14)+x(14)-x( 6)
      y(15)=x(15)+x(15)-x( 7)
      y(16)=x(16)+x(16)-x( 8)
*  Step N6
      y( 7)=x( 7)+x( 7)-x( 3)
      y( 8)=x( 8)+x( 8)-x( 4)
      y(11)=x(11)+x(11)-y( 7)
      y(12)=x(12)+x(12)-y( 8)
      y(15)=y(15)+y(15)-y(11)
      y(16)=y(16)+y(16)-y(12)
*  Step N7
      y( 4)=x( 4)+x( 4)-x( 2)
      y( 6)=x( 6)+x( 6)-y( 4)
      y( 8)=y( 8)+y( 8)-y( 6)
      y(10)=x(10)+x(10)-y( 8)
      y(12)=y(12)+y(12)-y(10)
      y(14)=y(14)+y(14)-y(12)
      y(16)=y(16)+y(16)-y(14)
*
      y( 1)=x( 1)
      y( 2)=x( 2)
      y( 3)=x( 3)
      y( 5)=x( 5)
      y( 9)=x( 9)
      return
      end
*
      subroutine ftcod5(x,y)
      implicit real*8 (a-h,o-z)
      dimension x(*),y(*)
      common
     >/ftco12/e12(4095)
      do j=1,16
        y(j)=x(j)+x(33-j)
        y(16+j)=(x(j)-x(33-j))*e12(128*(2*j-1))
      end do
      call ftcod4(y,x)
      call ftcod4(y(17),x(17))
      do j=18,32
        x(j)=x(j)+x(j)-x(j-1)
      end do
      do j=1,16
        y(j+j-1)=x(j)
        y(j+j)=x(16+j)
      end do
      return
      end
*
      subroutine ftcod6(x,y)
      implicit real*8 (a-h,o-z)
      dimension x(*),y(*)
      common
     >/ftco12/e12(4095)
      do j=1,32
        y(j)=x(j)+x(65-j)
        y(32+j)=(x(j)-x(65-j))*e12(64*(2*j-1))
      end do
      call ftcod5(y,x)
      call ftcod5(y(33),x(33))
      do j=34,64
        x(j)=x(j)+x(j)-x(j-1)
      end do
      do j=1,32
        y(j+j-1)=x(j)
        y(j+j)=x(32+j)
      end do
      return
      end
*
      subroutine ftcod7(x,y)
      implicit real*8 (a-h,o-z)
      dimension x(*),y(*)
      common
     >/ftco12/e12(4095)
      do j=1,64
        y(j)=x(j)+x(129-j)
        y(64+j)=(x(j)-x(129-j))*e12(32*(2*j-1))
      end do
      call ftcod6(y,x)
      call ftcod6(y(65),x(65))
      do j=66,128
        x(j)=x(j)+x(j)-x(j-1)
      end do
      do j=1,64
        y(j+j-1)=x(j)
        y(j+j)=x(64+j)
      end do
      return
      end
*
      subroutine ftcod8(x,y)
      implicit real*8 (a-h,o-z)
      dimension x(*),y(*)
      common
     >/ftco12/e12(4095)
      do j=1,128
        y(j)=x(j)+x(257-j)
        y(128+j)=(x(j)-x(257-j))*e12(16*(2*j-1))
      end do
      call ftcod7(y,x)
      call ftcod7(y(129),x(129))
      do j=130,256
        x(j)=x(j)+x(j)-x(j-1)
      end do
      do j=1,128
        y(j+j-1)=x(j)
        y(j+j)=x(128+j)
      end do
      return
      end
*
      subroutine ftcod9(x,y)
      implicit real*8 (a-h,o-z)
      dimension x(*),y(*)
      common
     >/ftco12/e12(4095)
      do j=1,256
        y(j)=x(j)+x(513-j)
        y(256+j)=(x(j)-x(513-j))*e12(8*(2*j-1))
      end do
      call ftcod8(y,x)
      call ftcod8(y(257),x(257))
      do j=258,512
        x(j)=x(j)+x(j)-x(j-1)
      end do
      do j=1,256
        y(j+j-1)=x(j)
        y(j+j)=x(256+j)
      end do
      return
      end
*
      subroutine ftcod10(x,y)
      implicit real*8 (a-h,o-z)
      dimension x(*),y(*)
      common
     >/ftco12/e12(4095)
      do j=1,512
        y(j)=x(j)+x(1025-j)
        y(512+j)=(x(j)-x(1025-j))*e12(4*(2*j-1))
      end do
      call ftcod9(y,x)
      call ftcod9(y(513),x(513))
      do j=514,1024
        x(j)=x(j)+x(j)-x(j-1)
      end do
      do j=1,512
        y(j+j-1)=x(j)
        y(j+j)=x(512+j)
      end do
      return
      end
*
      subroutine ftcod11(x,y)
      implicit real*8 (a-h,o-z)
      dimension x(*),y(*)
      common
     >/ftco12/e12(4095)
      do j=1,1024
        y(j)=x(j)+x(2049-j)
        y(1024+j)=(x(j)-x(2049-j))*e12(2*(2*j-1))
      end do
      call ftcod10(y,x)
      call ftcod10(y(1025),x(1025))
      do j=1026,2048
        x(j)=x(j)+x(j)-x(j-1)
      end do
      do j=1,1024
        y(j+j-1)=x(j)
        y(j+j)=x(1024+j)
      end do
      return
      end
*
      subroutine ftcod12(x,y)
      implicit real*8 (a-h,o-z)
      dimension x(*),y(*)
      common
     >/ftco12/e12(4095)
      do j=1,2048
        y(j)=x(j)+x(4097-j)
        y(2048+j)=(x(j)-x(4097-j))*e12(2*j-1)
      end do
      call ftcod11(y,x)
      call ftcod11(y(2049),x(2049))
      do j=2050,4096
        x(j)=x(j)+x(j)-x(j-1)
      end do
      do j=1,2048
        y(j+j-1)=x(j)
        y(j+j)=x(2048+j)
      end do
      return
      end
*
      subroutine ftcob1(x,y)
      implicit real*8 (a-h,o-z)
      dimension x(16),y(16)
      common/ftco4/e4(15),e3(7),e2(3)
      x(1)=0.5*x(1)
      x(2)=e2(2)*x(2)
      y(1)=x(1)+x(2)
      y(2)=x(1)-x(2)
      return
      end
*
      subroutine ftcob2(x,y)
      implicit real*8 (a-h,o-z)
      dimension x(4),y(4)
      common/ftco4/e4(15),e3(7),e2(3)
*  Step N1
      x(4)=x(4)+x(4)
      x(2)=x(2)+x(2)-x(4)
*
      y(1)=0.5*x(1)
      y(2)=0.5*x(2)
      y(3)=e2(2)*x(3)
      y(4)=e2(2)*x(4)
*  Step N2
      x(1)= y(1)+y(3)
      x(3)= y(1)-y(3)
      x(2)=(y(2)+y(4))*e2(1)
      x(4)=(y(2)-y(4))*e2(3)
*  Step N3
      y(1)=x(1)+x(2)
      y(4)=x(1)-x(2)
      y(2)=x(3)+x(4)
      y(3)=x(3)-x(4)
*
      return
      end
*
      subroutine ftcob3(x,y)
      implicit real*8 (a-h,o-z)
      dimension x(8),y(8)
      common/ftco4/e4(15),e3(7),e2(3)
*  Step N1
      x(8)=x(8)+x(8)
      x(6)=x(6)+x(6)-x(8)
      x(4)=x(4)+x(4)-x(6)
      x(2)=x(2)+x(2)-x(4)
*  Step N2
      x(7)=x(7)+x(7)
      x(3)=x(3)+x(3)-x(7)
      x(8)=x(8)+x(8)
      x(4)=x(4)+x(4)-x(8)
*
      x(1)=0.5*x(1)
      x(2)=0.5*x(2)
      x(3)=0.5*x(3)
      x(4)=0.5*x(4)
      x(5)=e3(4)*x(5)
      x(6)=e3(4)*x(6)
      x(7)=e3(4)*x(7)
      x(8)=e3(4)*x(8)
*  Step N3
      y(1)= x(1)+x(5)
      y(5)= x(1)-x(5)
      y(3)=(x(3)+x(7))*e3(2)
      y(7)=(x(3)-x(7))*e3(6)
      y(2)= x(2)+x(6)
      y(6)= x(2)-x(6)
      y(4)=(x(4)+x(8))*e3(2)
      y(8)=(x(4)-x(8))*e3(6)
*  Step N4
      x(1)= y(1)+y(3)
      x(7)= y(1)-y(3)
      x(3)= y(5)+y(7)
      x(5)= y(5)-y(7)
      x(2)=(y(2)+y(4))*e3(1)
      x(8)=(y(2)-y(4))*e3(7)
      x(4)=(y(6)+y(8))*e3(3)
      x(6)=(y(6)-y(8))*e3(5)
*  Step N5
      y(1)=x(1)+x(2)
      y(8)=x(1)-x(2)
      y(2)=x(3)+x(4)
      y(7)=x(3)-x(4)
      y(3)=x(5)+x(6)
      y(6)=x(5)-x(6)
      y(4)=x(7)+x(8)
      y(5)=x(7)-x(8)
      return
      end
*
      subroutine ftcob4(x,y)
      implicit real*8 (a-h,o-z)
      dimension x(16),y(16)
      common/ftco4/e4(15),e3(7),e2(3)
*  Step N1
      x(16)=x(16)+x(16)
      x(14)=x(14)+x(14)-x(16)
      x(12)=x(12)+x(12)-x(14)
      x(10)=x(10)+x(10)-x(12)
      x( 8)=x( 8)+x( 8)-x(10)
      x( 6)=x( 6)+x( 6)-x( 8)
      x( 4)=x( 4)+x( 4)-x( 6)
      x( 2)=x( 2)+x( 2)-x( 4)
*  Step N2
      x(15)=x(15)+x(15)
      x(11)=x(11)+x(11)-x(15)
      x( 7)=x( 7)+x( 7)-x(11)
      x( 3)=x( 3)+x( 3)-x( 7)
      x(16)=x(16)+x(16)
      x(12)=x(12)+x(12)-x(16)
      x( 8)=x( 8)+x( 8)-x(12)
      x( 4)=x( 4)+x( 4)-x( 8)
*  Step N3
      x(13)=x(13)+x(13)
      x( 5)=x( 5)+x( 5)-x(13)
      x(15)=x(15)+x(15)
      x( 7)=x( 7)+x( 7)-x(15)
      x(14)=x(14)+x(14)
      x( 6)=x( 6)+x( 6)-x(14)
      x(16)=x(16)+x(16)
      x( 8)=x( 8)+x( 8)-x(16)
*
      y( 1)=0.5*x( 1)
      y( 2)=0.5*x( 2)
      y( 3)=0.5*x( 3)
      y( 4)=0.5*x( 4)
      y( 5)=0.5*x( 5)
      y( 6)=0.5*x( 6)
      y( 7)=0.5*x( 7)
      y( 8)=0.5*x( 8)
      y( 9)=e4(8)*x( 9)
      y(10)=e4(8)*x(10)
      y(11)=e4(8)*x(11)
      y(12)=e4(8)*x(12)
      y(13)=e4(8)*x(13)
      y(14)=e4(8)*x(14)
      y(15)=e4(8)*x(15)
      y(16)=e4(8)*x(16)
*  Step N4
      x( 1)= y( 1)+y( 9)
      x( 9)= y( 1)-y( 9)
      x( 5)=(y( 5)+y(13))*e4( 4)
      x(13)=(y( 5)-y(13))*e4(12)
      x( 3)= y( 3)+y(11)
      x(11)= y( 3)-y(11)
      x( 7)=(y( 7)+y(15))*e4( 4)
      x(15)=(y( 7)-y(15))*e4(12)
      x( 2)= y( 2)+y(10)
      x(10)= y( 2)-y(10)
      x( 6)=(y( 6)+y(14))*e4( 4)
      x(14)=(y( 6)-y(14))*e4(12)
      x( 4)= y( 4)+y(12)
      x(12)= y( 4)-y(12)
      x( 8)=(y( 8)+y(16))*e4( 4)
      x(16)=(y( 8)-y(16))*e4(12)
*  Step N5
      y( 1)= x( 1)+x( 5)
      y(13)= x( 1)-x( 5)
      y( 3)=(x( 3)+x( 7))*e4( 2)
      y(15)=(x( 3)-x( 7))*e4(14)
      y( 2)= x( 2)+x( 6)
      y(14)= x( 2)-x( 6)
      y( 4)=(x( 4)+x( 8))*e4( 2)
      y(16)=(x( 4)-x( 8))*e4(14)
      y( 5)= x( 9)+x(13)
      y( 9)= x( 9)-x(13)
      y( 7)=(x(11)+x(15))*e4( 6)
      y(11)=(x(11)-x(15))*e4(10)
      y( 6)= x(10)+x(14)
      y(10)= x(10)-x(14)
      y( 8)=(x(12)+x(16))*e4( 6)
      y(12)=(x(12)-x(16))*e4(10)
*  Step N6
      x( 1)= y( 1)+y( 3)
      x(15)= y( 1)-y( 3)
      x( 2)=(y( 2)+y( 4))*e4( 1)
      x(16)=(y( 2)-y( 4))*e4(15)
      x( 3)= y( 5)+y( 7)
      x(13)= y( 5)-y( 7)
      x( 4)=(y( 6)+y( 8))*e4( 3)
      x(14)=(y( 6)-y( 8))*e4(13)
      x( 5)= y( 9)+y(11)
      x(11)= y( 9)-y(11)
      x( 6)=(y(10)+y(12))*e4( 5)
      x(12)=(y(10)-y(12))*e4(11)
      x( 7)= y(13)+y(15)
      x( 9)= y(13)-y(15)
      x( 8)=(y(14)+y(16))*e4( 7)
      x(10)=(y(14)-y(16))*e4( 9)
*  Step N7
      y( 1)=x( 1)+x( 2)
      y(16)=x( 1)-x( 2)
      y( 2)=x( 3)+x( 4)
      y(15)=x( 3)-x( 4)
      y( 3)=x( 5)+x( 6)
      y(14)=x( 5)-x( 6)
      y( 4)=x( 7)+x( 8)
      y(13)=x( 7)-x( 8)
      y( 5)=x( 9)+x(10)
      y(12)=x( 9)-x(10)
      y( 6)=x(11)+x(12)
      y(11)=x(11)-x(12)
      y( 7)=x(13)+x(14)
      y(10)=x(13)-x(14)
      y( 8)=x(15)+x(16)
      y( 9)=x(15)-x(16)
*
      return
      end
*
      subroutine ftcob5(x,y)
      implicit real*8 (a-h,o-z)
      dimension x(*),y(*)
      common
     >/ftco12/e12(4095)
      do j=1,16
        y(j)=x(j+j-1)
        y(16+j)=x(j+j)
      end do
      y(32)=y(32)+y(32)
      do j=31,17,-1
        y(j)=y(j)+y(j)-y(j+1)
      end do
      call ftcob4(y,x)
      call ftcob4(y(17),x(17))
      do j=1,16
        x(16+j)=x(16+j)*e12(128*(j+j-1))
      end do
      do j=1,16
        y(j)=x(j)+x(16+j)
        y(33-j)=x(j)-x(16+j)
      end do
      return
      end
*
      subroutine ftcob6(x,y)
      implicit real*8 (a-h,o-z)
      dimension x(*),y(*)
      common
     >/ftco12/e12(4095)
      do j=1,32
        y(j)=x(j+j-1)
        y(32+j)=x(j+j)
      end do
      y(64)=y(64)+y(64)
      do j=63,33,-1
        y(j)=y(j)+y(j)-y(j+1)
      end do
      call ftcob5(y,x)
      call ftcob5(y(33),x(33))
      do j=1,32
        x(32+j)=x(32+j)*e12(64*(j+j-1))
      end do
      do j=1,32
        y(j)=x(j)+x(32+j)
        y(65-j)=x(j)-x(32+j)
      end do
      return
      end
*
      subroutine ftcob7(x,y)
      implicit real*8 (a-h,o-z)
      dimension x(*),y(*)
      common
     >/ftco12/e12(4095)
      do j=1,64
        y(j)=x(j+j-1)
        y(64+j)=x(j+j)
      end do
      y(128)=y(128)+y(128)
      do j=127,65,-1
        y(j)=y(j)+y(j)-y(j+1)
      end do
      call ftcob6(y,x)
      call ftcob6(y(65),x(65))
      do j=1,64
        x(64+j)=x(64+j)*e12(32*(j+j-1))
      end do
      do j=1,64
        y(j)=x(j)+x(64+j)
        y(129-j)=x(j)-x(64+j)
      end do
      return
      end
*
      subroutine ftcob8(x,y)
      implicit real*8 (a-h,o-z)
      dimension x(*),y(*)
      common
     >/ftco12/e12(4095)
      do j=1,128
        y(j)=x(j+j-1)
        y(128+j)=x(j+j)
      end do
      y(256)=y(256)+y(256)
      do j=255,129,-1
        y(j)=y(j)+y(j)-y(j+1)
      end do
      call ftcob7(y,x)
      call ftcob7(y(129),x(129))
      do j=1,128
        x(128+j)=x(128+j)*e12(16*(j+j-1))
      end do
      do j=1,128
        y(j)=x(j)+x(128+j)
        y(257-j)=x(j)-x(128+j)
      end do
      return
      end
*
      subroutine ftcob9(x,y)
      implicit real*8 (a-h,o-z)
      dimension x(*),y(*)
      common
     >/ftco12/e12(4095)
      do j=1,256
        y(j)=x(j+j-1)
        y(256+j)=x(j+j)
      end do
      y(512)=y(512)+y(512)
      do j=511,257,-1
        y(j)=y(j)+y(j)-y(j+1)
      end do
      call ftcob8(y,x)
      call ftcob8(y(257),x(257))
      do j=1,256
        x(256+j)=x(256+j)*e12(8*(j+j-1))
      end do
      do j=1,256
        y(j)=x(j)+x(256+j)
        y(513-j)=x(j)-x(256+j)
      end do
      return
      end
*
      subroutine ftcob10(x,y)
      implicit real*8 (a-h,o-z)
      dimension x(*),y(*)
      common
     >/ftco12/e12(4095)
      do j=1,512
        y(j)=x(j+j-1)
        y(512+j)=x(j+j)
      end do
      y(1024)=y(1024)+y(1024)
      do j=1023,513,-1
        y(j)=y(j)+y(j)-y(j+1)
      end do
      call ftcob9(y,x)
      call ftcob9(y(513),x(513))
      do j=1,512
        x(512+j)=x(512+j)*e12(4*(j+j-1))
      end do
      do j=1,512
        y(j)=x(j)+x(512+j)
        y(1025-j)=x(j)-x(512+j)
      end do
      return
      end
*
      subroutine ftcob11(x,y)
      implicit real*8 (a-h,o-z)
      dimension x(*),y(*)
      common
     >/ftco12/e12(4095)
      do j=1,1024
        y(j)=x(j+j-1)
        y(1024+j)=x(j+j)
      end do
      y(2048)=y(2048)+y(2048)
      do j=2047,1025,-1
        y(j)=y(j)+y(j)-y(j+1)
      end do
      call ftcob10(y,x)
      call ftcob10(y(1025),x(1025))
      do j=1,1024
        x(1024+j)=x(1024+j)*e12(2*(j+j-1))
      end do
      do j=1,1024
        y(j)=x(j)+x(1024+j)
        y(2049-j)=x(j)-x(1024+j)
      end do
      return
      end
*
      subroutine ftcob12(x,y)
      implicit real*8 (a-h,o-z)
      dimension x(*),y(*)
      common
     >/ftco12/e12(4095)
      do j=1,2048
        y(j)=x(j+j-1)
        y(2048+j)=x(j+j)
      end do
      y(4096)=y(4096)+y(4096)
      do j=4095,2049,-1
        y(j)=y(j)+y(j)-y(j+1)
      end do
      call ftcob11(y,x)
      call ftcob11(y(2049),x(2049))
      do j=1,2048
        x(2048+j)=x(2048+j)*e12(j+j-1)
      end do
      do j=1,2048
        y(j)=x(j)+x(2048+j)
        y(4097-j)=x(j)-x(2048+j)
      end do
      return
      end
*
      subroutine ftc1(x,y)
      implicit real*8 (a-h,o-z)
      dimension x(3),y(3)
      y(2)=0.5*(x(1)+x(3))
      y(1)=y(2)+x(2)
      y(3)=y(2)-x(2)
      y(2)=0.5*(x(1)-x(3))
      return
      end
*
      subroutine ftc2(x,y)
      implicit real*8 (a-h,o-z)
      dimension x(5),y(5)
      do k=1,2
        y(k)=x(k)-x(6-k)
        y(2+k)=x(k)+x(6-k)
      end do
      y(5)=x(3)+x(3)
      call ftcob1(y,x)
      call ftc1(y(3),x(3))
      do j=1,2
        y(j+j)=x(j)
        y(j+j-1)=x(2+j)
      end do
      y(5)=x(5)
      return
      end
*
      subroutine ftc3(x,y)
      implicit real*8 (a-h,o-z)
      dimension x(9),y(9)
      do k=1,4
        y(k)=x(k)-x(10-k)
        y(4+k)=x(k)+x(10-k)
      end do
      y(9)=x(5)+x(5)
      call ftcob2(y,x)
      call ftc2(y(5),x(5))
      do j=1,4
        y(j+j)=x(j)
        y(j+j-1)=x(4+j)
      end do
      y(9)=x(9)
      return
      end
*
      subroutine ftc4(x,y)
      implicit real*8 (a-h,o-z)
      dimension x(17),y(17)
      do k=1,8
        y(k)=x(k)-x(18-k)
        y(8+k)=x(k)+x(18-k)
      end do
      y(17)=x(9)+x(9)
      call ftcob3(y,x)
      call ftc3(y(9),x(9))
      do j=1,8
        y(j+j)=x(j)
        y(j+j-1)=x(8+j)
      end do
      y(17)=x(17)
      return
      end
*
      subroutine ftc5(x,y)
      implicit real*8 (a-h,o-z)
      dimension x(33),y(33)
      do k=1,16
        y(k)=x(k)-x(34-k)
        y(16+k)=x(k)+x(34-k)
      end do
      y(33)=x(17)+x(17)
      call ftcob4(y,x)
      call ftc4(y(17),x(17))
      do j=1,16
        y(j+j)=x(j)
        y(j+j-1)=x(16+j)
      end do
      y(33)=x(33)
      return
      end
*
      subroutine ftc6(x,y)
      implicit real*8 (a-h,o-z)
      dimension x(65),y(65)
      do k=1,32
        y(k)=x(k)-x(66-k)
        y(32+k)=x(k)+x(66-k)
      end do
      y(65)=x(33)+x(33)
      call ftcob5(y,x)
      call ftc5(y(33),x(33))
      do j=1,32
        y(j+j)=x(j)
        y(j+j-1)=x(32+j)
      end do
      y(65)=x(65)
      return
      end
*
      subroutine ftc7(x,y)
      implicit real*8 (a-h,o-z)
      dimension x(129),y(129)
      do k=1,64
        y(k)=x(k)-x(130-k)
        y(64+k)=x(k)+x(130-k)
      end do
      y(129)=x(65)+x(65)
      call ftcob6(y,x)
      call ftc6(y(65),x(65))
      do j=1,64
        y(j+j)=x(j)
        y(j+j-1)=x(64+j)
      end do
      y(129)=x(129)
      return
      end
*
      subroutine ftc8(x,y)
      implicit real*8 (a-h,o-z)
      dimension x(257),y(257)
      do k=1,128
        y(k)=x(k)-x(258-k)
        y(128+k)=x(k)+x(258-k)
      end do
      y(257)=x(129)+x(129)
      call ftcob7(y,x)
      call ftc7(y(129),x(129))
      do j=1,128
        y(j+j)=x(j)
        y(j+j-1)=x(128+j)
      end do
      y(257)=x(257)
      return
      end
*
      subroutine ftc9(x,y)
      implicit real*8 (a-h,o-z)
      dimension x(513),y(513)
      do k=1,256
        y(k)=x(k)-x(514-k)
        y(256+k)=x(k)+x(514-k)
      end do
      y(513)=x(257)+x(257)
      call ftcob8(y,x)
      call ftc8(y(257),x(257))
      do j=1,256
        y(j+j)=x(j)
        y(j+j-1)=x(256+j)
      end do
      y(513)=x(513)
      return
      end
*
      subroutine ftc10(x,y)
      implicit real*8 (a-h,o-z)
      dimension x(1025),y(1025)
      do k=1,512
        y(k)=x(k)-x(1026-k)
        y(512+k)=x(k)+x(1026-k)
      end do
      y(1025)=x(513)+x(513)
      call ftcob9(y,x)
      call ftc9(y(513),x(513))
      do j=1,512
        y(j+j)=x(j)
        y(j+j-1)=x(512+j)
      end do
      y(1025)=x(1025)
      return
      end
*
      subroutine ftc11(x,y)
      implicit real*8 (a-h,o-z)
      dimension x(2049),y(2049)
      do k=1,1024
        y(k)=x(k)-x(2050-k)
        y(1024+k)=x(k)+x(2050-k)
      end do
      y(2049)=x(1025)+x(1025)
      call ftcob10(y,x)
      call ftc10(y(1025),x(1025))
      do j=1,1024
        y(j+j)=x(j)
        y(j+j-1)=x(1024+j)
      end do
      y(2049)=x(2049)
      return
      end
*
      subroutine ftc12(x,y)
      implicit real*8 (a-h,o-z)
      dimension x(4097),y(4097)
      do k=1,2048
        y(k)=x(k)-x(4098-k)
        y(2048+k)=x(k)+x(4098-k)
      end do
      y(4097)=x(2049)+x(2049)
      call ftcob11(y,x)
      call ftc11(y(2049),x(2049))
      do j=1,2048
        y(j+j)=x(j)
        y(j+j-1)=x(2048+j)
      end do
      y(4097)=x(4097)
      return
      end
*
      subroutine fts1(x,y)
      implicit real*8 (a-h,o-z)
      y=x
      return
      end
*
      subroutine fts2(x,y)
      implicit real*8 (a-h,o-z)
      dimension x(3),y(3)
      common
     >/ftco4/e4(15),e3(7),e2(3)
      y(2)=e2(2)*(x(1)+x(3))
      y(1)=y(2)+x(2)
      y(3)=y(2)-x(2)
      y(2)=x(1)-x(3)
      return
      end
*
      subroutine fts3(x,y)
      implicit real*8 (a-h,o-z)
      dimension x(7),y(7)
      y(1)=x(4)+x(4)
      do k=1,3
        y(k+1)=x(4+k)+x(4-k)
        y(4+k)=x(4+k)-x(4-k)
      end do
      call ftcob2(y,x)
      call fts2(y(5),x(5))
      do j=1,4,2
        y(j+j-1)=x(j)
        y(j+j+1)=-x(j+1)
      end do
      y(2)=-x(5)
      y(4)=x(6)
      y(6)=-x(7)
      return
      end
*
      subroutine fts4(x,y)
      implicit real*8 (a-h,o-z)
      dimension x(15),y(15)
      y(1)=x(8)+x(8)
      do k=1,7
        y(k+1)=x(8+k)+x(8-k)
        y(8+k)=x(8+k)-x(8-k)
      end do
      call ftcob3(y,x)
      call fts3(y(9),x(9))
      do j=1,8,2
        y(j+j-1)=x(j)
        y(j+j+1)=-x(j+1)
      end do
      do j=1,6,2
        y(j+j)=-x(8+j)
        y(j+j+2)=x(9+j)
      end do
      y(14)=-x(15)
      return
      end
*
      subroutine fts5(x,y)
      implicit real*8 (a-h,o-z)
      dimension x(31),y(31)
      y(1)=x(16)+x(16)
      do k=1,15
        y(k+1)=x(16+k)+x(16-k)
        y(16+k)=x(16+k)-x(16-k)
      end do
      call ftcob4(y,x)
      call fts4(y(17),x(17))
      do j=1,16,2
        y(j+j-1)=x(j)
        y(j+j+1)=-x(j+1)
      end do
      do j=1,14,2
        y(j+j)=-x(16+j)
        y(j+j+2)=x(17+j)
      end do
      y(30)=-x(31)
      return
      end
*
      subroutine fts6(x,y)
      implicit real*8 (a-h,o-z)
      dimension x(63),y(63)
      y(1)=x(32)+x(32)
      do k=1,31
        y(k+1)=x(32+k)+x(32-k)
        y(32+k)=x(32+k)-x(32-k)
      end do
      call ftcob5(y,x)
      call fts5(y(33),x(33))
      do j=1,32,2
        y(j+j-1)=x(j)
        y(j+j+1)=-x(j+1)
      end do
      do j=1,30,2
        y(j+j)=-x(32+j)
        y(j+j+2)=x(33+j)
      end do
      y(62)=-x(63)
      return
      end
*
      subroutine fts7(x,y)
      implicit real*8 (a-h,o-z)
      dimension x(127),y(127)
      y(1)=x(64)+x(64)
      do k=1,63
        y(k+1)=x(64+k)+x(64-k)
        y(64+k)=x(64+k)-x(64-k)
      end do
      call ftcob6(y,x)
      call fts6(y(65),x(65))
      do j=1,64,2
        y(j+j-1)=x(j)
        y(j+j+1)=-x(j+1)
      end do
      do j=1,62,2
        y(j+j)=-x(64+j)
        y(j+j+2)=x(65+j)
      end do
      y(126)=-x(127)
      return
      end
*
      subroutine fts8(x,y)
      implicit real*8 (a-h,o-z)
      dimension x(255),y(255)
      y(1)=x(128)+x(128)
      do k=1,127
        y(k+1)=x(128+k)+x(128-k)
        y(128+k)=x(128+k)-x(128-k)
      end do
      call ftcob7(y,x)
      call fts7(y(129),x(129))
      do j=1,128,2
        y(j+j-1)=x(j)
        y(j+j+1)=-x(j+1)
      end do
      do j=1,126,2
        y(j+j)=-x(128+j)
        y(j+j+2)=x(129+j)
      end do
      y(254)=-x(255)
      return
      end
*
      subroutine fts9(x,y)
      implicit real*8 (a-h,o-z)
      dimension x(511),y(511)
      y(1)=x(256)+x(256)
      do k=1,255
        y(k+1)=x(256+k)+x(256-k)
        y(256+k)=x(256+k)-x(256-k)
      end do
      call ftcob8(y,x)
      call fts8(y(257),x(257))
      do j=1,256,2
        y(j+j-1)=x(j)
        y(j+j+1)=-x(j+1)
      end do
      do j=1,254,2
        y(j+j)=-x(256+j)
        y(j+j+2)=x(257+j)
      end do
      y(510)=-x(511)
      return
      end
*
      subroutine fts10(x,y)
      implicit real*8 (a-h,o-z)
      dimension x(1023),y(1023)
      y(1)=x(512)+x(512)
      do k=1,511
        y(k+1)=x(512+k)+x(512-k)
        y(512+k)=x(512+k)-x(512-k)
      end do
      call ftcob9(y,x)
      call fts9(y(513),x(513))
      do j=1,512,2
        y(j+j-1)=x(j)
        y(j+j+1)=-x(j+1)
      end do
      do j=1,510,2
        y(j+j)=-x(512+j)
        y(j+j+2)=x(513+j)
      end do
      y(1022)=-x(1023)
      return
      end
*
      subroutine fts11(x,y)
      implicit real*8 (a-h,o-z)
      dimension x(2047),y(2047)
      y(1)=x(1024)+x(1024)
      do k=1,1023
        y(k+1)=x(1024+k)+x(1024-k)
        y(1024+k)=x(1024+k)-x(1024-k)
      end do
      call ftcob10(y,x)
      call fts10(y(1025),x(1025))
      do j=1,1024,2
        y(j+j-1)=x(j)
        y(j+j+1)=-x(j+1)
      end do
      do j=1,1022,2
        y(j+j)=-x(1024+j)
        y(j+j+2)=x(1025+j)
      end do
      y(2046)=-x(2047)
      return
      end
*
      subroutine fts12(x,y)
      implicit real*8 (a-h,o-z)
      dimension x(4095),y(4095)
      y(1)=x(2048)+x(2048)
      do k=1,2047
        y(k+1)=x(2048+k)+x(2048-k)
        y(2048+k)=x(2048+k)-x(2048-k)
      end do
      call ftcob11(y,x)
      call fts11(y(2049),x(2049))
      do j=1,2048,2
        y(j+j-1)=x(j)
        y(j+j+1)=-x(j+1)
      end do
      do j=1,2046,2
        y(j+j)=-x(2048+j)
        y(j+j+2)=x(2049+j)
      end do
      y(4094)=-x(4095)
      return
      end