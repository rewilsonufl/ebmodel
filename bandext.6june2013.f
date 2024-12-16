      subroutine bandext (extin,iband,xx,ff,extout)
c  Version of June 7, 2013
      implicit real*8(a-h,o-z)
      dimension xw(15),xx(*),ff(*),extout(*)
      common /effwave/ effwvl(100)
      data xw(5),xw(6),xw(7),xw(8),xw(9),xw(10),xw(11),xw(12),xw(13)/
     $2.78d0,2.27d0,1.82d0,1.43d0,1.11d0,0.80d0,0.46d0,0.29d0,0.63d0/
c
c  The approximation functions applied here are from Cardelli, Clayton, & Mathis (1989)
c  BANDEXT computes ratios of extinctions in each of 93 photometric
c    bands to extinction in the Johnson V  band.  
c    It outputs extinction in stellar magnitudes for all the bands in the 
c    EXTOUT array. Adopted effective wavelengths were computed by W. Van Hamme 
c    from response curves in the literature except for bands 5 thtough 12
c    which are from Cardelli, Clayton, & Mathis. Those bands are: 
c 
c       5        U      Buser, R. 1978, Ang, 62, 411
c       6        B      Azusienis and Straizys 1969, Sov. Astron., 13, 316
c       7        V          "             "                "
c       8        R      Johnson, H.L. 1965, ApJ, 141, 923
c       9        I         "            "    "
c      10        J         "            "    "
c      11        K         "            "    "
c      12        L         "            "    "

c  Input 'extin' is extinction in band 'iband' and dimensioned variable 'ff' contains
c     the extinction ratios. 
c
c
      do 23 iyn=1,93
      xx(iyn)=1.d3/effwvl(iyn)
      if(iyn.ge.5.and.iyn.le.12) xx(iyn)=xw(iyn)
   23 continue
      rrrv=3.1d0
      do 18 ix=1,93
      xm=xx(ix)
      if(xm.ge.3.30d0) goto 14
      if(xm.le.1.10d0) goto 11
      yyy=xm-1.82d0
      ysq=yyy*yyy
      ycb=ysq*yyy
      y4=ycb*yyy
      y5=y4*yyy
      y6=y5*yyy
      y7=y6*yyy
      aofx=1.d0+.17699d0*yyy-.50447d0*ysq-.02427d0*ycb+.72085d0*y4
     $+.01979d0*y5-.77530d0*y6+.32999d0*y7
      bofx=1.41338d0*yyy+2.28305d0*ysq+1.07233d0*ycb-5.38434d0*y4
     $-.62251d0*y5+5.30260d0*y6-2.09002d0*y7
      goto 12
   11 continue
      aofx=.574d0*xm**1.61d0
      bofx=-.527d0*xm**1.61d0
      goto 12
   14 fa=0.d0
      fb=0.d0
      if(xm.lt.5.9d0) goto 15
      xmf2=(xm-5.9d0)**2
      xmf3=(xm-5.9d0)**3
      fa=-0.04473d0*xmf2-0.009779d0*xmf3
      fb=0.2130d0*xmf2+0.1207d0*xmf3
   15 continue
      aofx=1.752d0-0.316d0*xm-0.104d0/((xm-4.67d0)**2+0.341d0)+fa
      bofx=-3.090d0+1.825d0*xm+1.206d0/((xm-4.62d0)**2+0.263d0)+fb
   12 ff(ix)=aofx+bofx/rrrv
   18 continue
c
c  Calculate V-band extinction (as an intermediary scaling quantity)
c
      ev=extin/ff(iband)
c
c  Calculate extinction in all bands
c
      do 25 iyn=1,93
   25 extout(iyn)=ff(iyn)*ev
      return
      end
