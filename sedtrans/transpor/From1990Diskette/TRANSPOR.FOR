$debug
C     PROGRAM TRANSPOR
C
C     COMPUTATION OF CONCENTRATION VERTICALS AND SUSPENDED LOAD TRANSPORT
C     THROUGH NUMERICAL INTEGRATION AND FORMULAE APPROXIMATION
C
      LOGICAL BO
      DIMENSION RMAT(25,6)
      CHARACTER KAR*1,line*80
      common/niorez/tp,vr,phi,g,pi,hd
      external fwl,tolx
      logical log
C
C  READING INPUT FROM SCREEN AND WRITING TO FILE TRANSPOR.OUT
C
      OPEN(2,FILE='TRANSPOR.HLP',STATUS='UNKNOWN')
      itel=0
   10 read(2,'(a)',end=20)line
      itel=itel+1
      if(mod(itel,18).eq.0)then
       write(*,'(a)')' Type C to continue '
       read(*,'(a)')kar
       itel=0
      endif
      write(*,'(1x,a)')line(1:79)
      goto 10
   20 continue
      close(2)
      OPEN(2,FILE='TRANSPOR.OUT',STATUS='UNKNOWN')
      WRITE(2,'(A)')CHAR(12)
      WRITE(2,'('' INPUT :'',/)')
      WRITE(*,'('' HD   = WATER DEPTH [ M ] =                      '')')
      WRITE(2,'('' HD   = WATER DEPTH [ M ] =                      '')')
      READ (*,*)HD
      WRITE(2,*)HD
      WRITE(*,'('' VR   = MEAN VEL. IN CURRENT DIR. [ M/S ] =      '')')
      WRITE(2,'('' VR   = MEAN VEL. IN CURRENT DIR. [ M/S ] =      '')')
      READ (*,*)VR
      WRITE(2,*)VR
      WRITE(*,'('' UR   = MEAN VEL. IN WAVE DIR(-BACK) [ M/S ] =   '')')
      WRITE(2,'('' UR   = MEAN VEL. IN WAVE DIR(-BACK) [ M/S ] =   '')')
      READ (*,*)UR
      WRITE(2,*)UR
      WRITE(*,'('' UB   = NEAR-BED VEL IN WAVE DIR(-BACK) [ M/S ]= '')')
      WRITE(2,'('' UB   = NEAR-BED VEL IN WAVE DIR(-BACK) [ M/S ]= '')')
      READ (*,*)UB
      WRITE(2,*)UB
      WRITE(*,'('' HS   = SIGNIFICANT WAVE HEIGHT [ M ] =          '')')
      WRITE(2,'('' HS   = SIGNIFICANT WAVE HEIGHT [ M ] =          '')')
      READ (*,*)HS
      WRITE(2,*)HS
      WRITE(*,'('' TP   = PEAK WAVE PERIOD        [ S ] =          '')')
      WRITE(2,'('' TP   = PEAK WAVE PERIOD        [ S ] =          '')')
      READ (*,*)TP
      WRITE(2,*)TP
      WRITE(*,'('' PHI  = ANGLE CURRENT AND WAVES 0-360 [ DEG ]    '')')
      WRITE(2,'('' PHI  = ANGLE CURRENT AND WAVES 0-360 [ DEG ]    '')')
      READ (*,*)PHI
      WRITE(2,*)PHI
      WRITE(*,'('' D50  = MEDIAN PARTICLE SIZE OF BED [ M ] =      '')')
      WRITE(2,'('' D50  = MEDIAN PARTICLE SIZE OF BED [ M ] =      '')')
      READ (*,*)D50
      WRITE(2,*)D50
      WRITE(*,'('' D90  = 90 0/0 PARTICLE SIZE OF BED [ M ] =      '')')
      WRITE(2,'('' D90  = 90 0/0 PARTICLE SIZE OF BED [ M ] =      '')')
      READ (*,*)D90
      WRITE(2,*)D90
      WRITE(*,'('' DSS  = SUSPENDED SEDIMENT SIZE   [ M ] =        '')')
      WRITE(2,'('' DSS  = SUSPENDED SEDIMENT SIZE   [ M ] =        '')')
      READ (*,*)DSS
      WRITE(2,*)DSS
      WRITE(*,'('' RC   = CURRENT-RELATED ROUGHNESS [ M ] =        '')')
      WRITE(2,'('' RC   = CURRENT-RELATED ROUGHNESS [ M ] =        '')')
      READ (*,*)RC
      WRITE(2,*)RC
      WRITE(*,'('' RW   = WAVE-RELATED ROUGHNESS [ M ] =           '')')
      WRITE(2,'('' RW   = WAVE-RELATED ROUGHNESS [ M ] =           '')')
      READ (*,*)RW
      WRITE(2,*)RW
      WRITE(*,'('' TE   = WATER TEMPERATURE [ CELSIUS ] =         '')')
      WRITE(2,'('' TE   = WATER TEMPERATURE [ CELSIUS ] =         '')')
      READ (*,*)TE
      WRITE(2,*)TE
      WRITE(*,'('' SA   = SALINITY OF FLUID [ PROMILLE ] =         '')')
      WRITE(2,'('' SA   = SALINITY OF FLUID [ PROMILLE ] =         '')')
      READ (*,*)SA
      WRITE(2,*)SA
C
C  INITIALIZING SOME NUMERICAL AND PHYSICAL PARAMETERS
C
      if(hs.eq.0.)then
        ur=0.
        ub=0.
      endif
      if(tp.gt.0..and.tp.lt.1.)tp=1.
      if(phi.eq.0.)phi=0.000001
      if(phi.eq.180.)phi=179.99999
      if(phi.eq.360.)phi=359.99999
      NN     = 12
      G      = 9.81
      PI     = 4.*ATAN(1.)
      phi    = phi/180.*pi
      cl=(sa-0.03)/1.805
      rhow=1000.+1.455*cl-0.0065*(te-4.+0.4*cl)**2.
      rhos=2650.
      rnu=(4.e-5)/(20.+te)
      DEL=(RHOS-RHOW)/RHOW
      dsh=0.01*g*del*dss**3./rnu/rnu
      d50h=0.01*g*del*d50**3./rnu/rnu
      if(dss.lt.0.0001)ws=(del*g*dss*dss)/(18.*rnu)
      if(d50.lt.0.0001)wsb=(del*g*d50*d50)/(18.*rnu)
      if(dss.ge.0.0001.and.dss.lt.0.001)
     * ws=(10.*rnu/dss)*((1.+dsh)**0.5-1.)
      if(d50.ge.0.0001.and.d50.lt.0.001)
     * wsb=(10.*rnu/d50)*((1.+d50h)**0.5-1.)
      if(dss.ge.0.001)ws=1.1*(del*g*dss)**0.5
      if(d50.ge.0.001)wsb=1.1*(del*g*d50)**0.5
      RKAP=.4
      if(rc.lt.0.01)rc=0.01
      if(rw.lt.0.01)rw=0.01
      DSTER=D50*(DEL*G/RNU**2)**(1./3.)
C
C   COMPUTATION OF WAVE PARAMETERS
C
      IF(HD.LE.0.)STOP 'WATER DEPTH .LE. 0'
      IF (TP.GE.1.) THEN
         xx=.1
         yy=1.8*tp*tp
         uugg=VR
         VR=0.
         call zeroin(xx,yy,log,fwl,tolx)
         if(.not.log)then
           stop 'zeroin no wave length'
         else
           yy=2.0*xx
           xx=.2*xx
           VR=uugg
           call zeroin(xx,yy,log,fwl,tolx)
           if(.not.log)then
              stop 'zeroin no wave length'
           else
              rls=xx
           endif
         endif
         arg=2.*pi*hd/rls
         tp1=tp/(1.-(VR*tp*cos(phi))/rls)
       IF (ARG.GT.50.) THEN
          ABW=0.
          FW=0.
          FW1=0.
          UBW=0.
        ELSE
          ABW=HS/(2.*SINH(ARG))
          FW=0.
          FW1=0.
          IF(ABW.GT.0.)FW=EXP(-6.+5.2*(ABW/RW)**(-0.19))
          IF(ABW.GT.0.)FW1=EXP(-6.+5.2*(ABW/(3.*D90))**(-0.19))
          IF(FW.GT.0.3 )FW=0.3
          IF(FW1.GT.0.3)FW1=0.3
          UBW  = 2.*PI/TP1*ABW
        ENDIF
      ELSE
        ABW=0.
        RLS=0.
        FW=0.
        FW1=0.
        UBW=0.
        tp1=tp
      ENDIF
        ubw1=0.
        if(rls.gt.0.)
     *  ubw1=(3.*pi*pi*hs*hs)/(4.*tp1*rls*(sinh(arg)**4.))
        rhs=hs/hd
      h1=1.+0.3*rhs
      hlimit=0.01*g*tp1*tp1
      if(hd.ge.hlimit)then
        ubwfor=ubw+ubw1
        ubwback=ubw-ubw1
      endif
      if(hd.lt.hlimit)then
        ubwfor=h1*ubw
        ubwback=(2.-h1)*ubw
      endif
      ufmax=1.2*ubw
      ubmin=0.8*ubw
      if(ubwfor.ge.ufmax)ubwfor=ufmax
      if(ubwback.le.ubmin)ubwback=ubmin
      asym=0.
      if(abw.gt.0.)asym=ubwfor/(ubwfor+ubwback)
      htrough=hd*(0.95-0.35*(hs/hd))
      if(ur.eq.9.)ur=-0.125*g**0.5*hs**2./(hd**0.5*htrough)
      if(ub.eq.9.)ub=(0.05-(asym-0.5))*ubw
C
C  CRITICAL SHEAR STRESS SHIELDS
C
      IF(DSTER.LE.4.)THETCR=.24/DSTER
      IF(4. .LT.DSTER.AND.DSTER.LE.10.)THETCR=.14*DSTER**(-.64)
      IF(10..LT.DSTER.AND.DSTER.LE.20.)THETCR=.04*DSTER**(-.1 )
      IF(20..LT.DSTER.AND.DSTER.LE.150.)THETCR=.013*DSTER**(.29 )
      IF(DSTER.GT.150.)THETCR=.055
      TAUCR=(RHOS-RHOW)*G*D50*THETCR
      umcr=5.75*(del*g*d50)**0.5*(thetcr)**0.5*alog10(4.*hd/d90)
      ubwcr=0.
      if(tp.gt.0.)then
       if(d50.le.0.0005)then
             ubwcr=(0.12*del*g*d50**0.5*tp**0.5)**0.667
       else
             ubwcr=(1.09*del*g*d50**0.75*tp**0.25)**0.57
       endif
      endif
C
C  COMPUTATION OF REFERENCE CONCENTRATION CA
C
      a=max(rc,rw)
      CC=18.*ALOG10(12.*HD/RC)
      CC1=18.*ALOG10(12.*HD/3./D90)
      FC=0.24*ALOG10(12.*HD/RC)**(-2)
      FC1=0.24*ALOG10(12.*HD/3./D90)**(-2)
      vrr=(VR**2+ur**2)**0.5
      RMUC=FC1/FC
      RMUWA=0.6/DSTER
      if(dster.gt.10.)RMUWA=0.06
      rmuw=0.
      if(fw.ne.0.)rmuw=fw1/fw
      UST=G**0.5*ABS(vrr)/CC
      DELm=0.
      delw=0.
      IF(ABW.GT.0.)THEN
         DELW=0.072*ABW*(ABW/RW)**(-0.25)
         delm=3.*delw
      ENDIF
      delm=max(delm,rc)
      GAMMA=0.
      h2=phi
      if(h2.gt.pi)h2=2.*pi-phi
      gamma=0.8+h2-0.3*h2*h2
      uratio=ubw/vrr
      if(uratio.gt.5.)uratio=5.
      RA=EXP(GAMMA*uratio)*RC
      RCC=10.*RC
      IF(RA.GE.RCC)RA=RCC
      FCA=0.24*ALOG10(12.*HD/RA)**(-2)
      IF(DELm.LE.Ra/29.9)delm=ra/29.9
      h3=(ALOG(30.*DELM/RA)/ALOG(30.*DELM/RC))**2
      h4=((-1.+alog(30.*hd/rc))/(-1.+alog(30.*hd/ra)))**2
      alfaw=h3*h4
      TAUC=0.125*RHOW*FC*vrr*vrr
      TAUW=0.25*RHOW*FW*UBW*UBW
      tauce=rmuc*tauc
      TAUCEF=RMUC*ALFAW*TAUC
      TAUWEF=RMUWA*TAUW
      TAUCWE=TAUCEF+TAUWEF
      if(taucwe.le.tauce)taucwe=tauce
      uster1=sqrt(taucef/rhow)
      uster2=sqrt(tauce/rhow)
      IF(USTER1.LE.USTER2)USTER1=USTER2
      TA=(TAUCWE-TAUCR)/TAUCR
      TA=MAX(.0001,TA)
      CA=0.015*D50/A*DSTER**(-.3)*TA**1.5
      if(ca.ge.0.60)ca=0.6
      taub=rmuc*alfaw*tauc+rmuw*tauw
      if(taub.le.tauce)taub=tauce
      t=(taub-taucr)/taucr
      t=max(0.0001,t)
C
C SUSPENDED LOAD TRANSPORT IN CURRENT DIRECTION
C
C COMPUTATION OF NUMERICAL PARAMETERS AT Z=A
C
      JTAL= 8
      NN  = JTAL*NN
      DYM = CA/NN
      DXM = HD/NN
      DYX = DYM/DXM
      BO  = .FALSE.
      HULP30=-1.+ALOG(30.*HD/RA)
      IF(DELM.GT.0.)THEN
         UDEL=VR*ALOG(30.*DELM/RA)/HULP30
      ENDIF
C
C  COMPUTATION OF DERIVATIVE DC/DY OR DC/DX
C
      beta=1.+2.*(ws/ust)**2.
      if(beta.ge.1.5)beta=1.5
      ds=0.3*hd*(hs/hd)**0.5
      if(ds.gt.0.2)ds=0.2
      if(ds.le.0.05)ds=0.05
      EBW=.004*DSTER*DS*UBW
      IF (TP.GE.1.) THEN
        EMAXW=0.035*HD*HS/TP
      ELSE
        EMAXW=0.
      ENDIF
      IF(EMAXW.LE.EBW)EMAXW=EBW
      EMAXC=0.25*RKAP*UST*HD*BETA
      C=CA
      Z=A
      IF(Z.LE.DS)ESW=EBW
      IF(Z.GT.DS.AND.Z.LE.0.5*HD)ESW=EBW+(EMAXW-EBW)*((Z-DS)/
     *(0.5*HD-DS))
      IF(Z.GE.0.5*HD)ESW=EMAXW
      IF(Z.GE.0.5*HD)ESC=EMAXC
      IF(Z.LT.0.5*HD)ESC=EMAXC-EMAXC*(1.-2.*Z/HD)**2
      ES=(ESW**2.+ESC**2.)**0.5
      cm=0.65
      fcc=0.
      IF (C.GT.1.E-8) THEN
        IF(Z.GE.A)then
            fi =1.+(c/cm)**.8-2.*(c/cm)**.4
            fcc=-WS/ES/fi*C*(1.-C)**5.
        endif
      ENDIF
      YPRIME=fcc
      FF = 1./CA*YPRIME
      IF(DELM.GT.0.)THEN
         UC=UDEL*ALOG(30.*A/RC)/ALOG(30.*DELM/RC)
      ENDIF
      IF(A.GE.DELM)UC=VR*ALOG(30.*A/RA)/HULP30
      if(a.le.rc/30.)uc=0.
C
C  STORAGE OF RESULTS IN MATRIX RMAT
C
      RMAT(1,1)=Z
      RMAT(1,2)=CA*RHOS
      RMAT(1,3)=UC
      RMAT(1,4)=0.
C
C INTEGRATION FROM Z=A TO SURFACE
C
      Y = CA
      TERM1=UC*Y
      XEND=A
      SSC=0.
      NTEL = 0
      IT   = 2
C
  100 CONTINUE
      NTEL = NTEL+1
      XOLD = XEND
      YOLD = Y
      IF (-YPRIME .GT. DYX) THEN
         Y = YOLD-DYM
         IF (Y .LT. 2./3.*YOLD) Y = 2./3.*YOLD
         XEND = XOLD+ALOG(Y/YOLD)/FF
      ELSE
         XEND = XOLD+DXM
         IF (XEND .GE. HD) THEN
            XEND = HD
            BO = .TRUE.
         ENDIF
         Y = EXP(ALOG(YOLD)+(XEND-XOLD)*FF)
      ENDIF
      C=Y
      Z=XEND
      IF(Z.LE.DS)ESW=EBW
      IF(Z.GT.DS.AND.Z.LE.0.5*HD)ESW=EBW+(EMAXW-EBW)*((Z-DS)/
     *(0.5*HD-DS))
      IF(Z.GE.0.5*HD)ESW=EMAXW
      IF(Z.GE.0.5*HD)ESC=EMAXC
      IF(Z.LT.0.5*HD)ESC=EMAXC-EMAXC*(1.-2.*Z/HD)**2
      ES=(ESW**2.+ESC**2.)**0.5
      fcc=0.
      IF (C.GT.1.E-8) THEN
        IF(Z.GE.A)then
            fi =1.+(c/cm)**.8-2.*(c/cm)**.4
            fcc=-WS/ES/fi*C*(1.-C)**5.
        endif
      ENDIF
      YPRIME=fcc
      FF = 1./Y*YPRIME
      IF(DELM.GT.0.)THEN
         UC=UDEL*ALOG(30.*XEND/RC)/ALOG(30.*DELM/RC)
      ENDIF
      IF(XEND.GE.DELM)UC=VR*ALOG(30.*XEND/RA)/HULP30
      if(xend.le.rc/30.)uc=0.
      TERM2=UC*Y
      TERM1=TERM2
C
C  STORAGE OF RESULTS IN MATRIX RMAT
C
      IF (NTEL .EQ. NN/JTAL.OR.BO) THEN
         RMAT(IT,1)=Z
         dz=rmat(it,1)-rmat(it-1,1)
         RMAT(IT,2)=Y*RHOS
         ssc=ssc+(rmat(it,2)*uc+rmat(it-1,2)*rmat(it-1,3))/2.*dz
         RMAT(IT,3)=UC
         RMAT(IT,4)=SSC
         IT = IT+1
         NTEL      = 0
      ENDIF
      IF (.NOT. BO) GOTO 100
c
c  SUSPENDED LOAD TRANSPORT IN WAVE DIRECTION
c
       ssw=0.
       const=(ur/VR)
       uz1=const*rmat(1,3)
       rmat(1,5)=uz1
       rmat(1,6)=0.
       do 250 i=2,it-1
           uz2=const*rmat(i,3)
           rmat(i,5)=uz2
           dz=rmat(i,1)-rmat(i-1,1)
           torm=(uz1*rmat(i-1,2)+uz2*rmat(i,2))/2.*dz/rhos
           ssw=ssw+torm*rhos
           rmat(i,6)=ssw
           uz1=uz2
  250 continue
c
c  INSTANTANEOUS BED LOAD TRANSPORT 
c
      vrdelm=vr*alog(30.*delm/ra)/(-1.+alog(30.*hd/ra))
      urdelm=(ur/vr)*vrdelm
      IF(DELM.LE.RA/30.)THEN
        vrdelm=0.
        urdelm=0.
      endif
      ubtot=ub+urdelm
      rmfor=(ubwfor+ubtot)**2./(del*g*d50)
      rmback=(-ubwback+ubtot)**2./(del*g*d50)
      rmcr=ubwcr**2./(del*g*d50)
      ntime=201
      dtt=tp1/(ntime-1)
      tfor=0.
      if(abw.gt.0.)tfor=ubwback/(ubwfor+ubwback)*tp1
      tback=tp1-tfor
      profact=((-1.+alog(30.*hd/rc))/alog(30.*delm/rc))**2.
      fc11=0.25*fc1*profact
      acw=abs(vrdelm)/(abs(ubw)+abs(vrdelm))
      fcw1=acw*fc11+(1.-acw)*fw1
      sbx=0.
      sby=0.
      tau1som=0.
      tau1x=0.
      tau1y=0.
      nnn=0
  150 continue
      time=nnn*dtt
      udt=0.
      if(abw.gt.0.)then
      if(hd.ge.hlimit)then
         udt=ubw*sin(2.*pi*time/tp1)+ubw1*sin(4.*pi*time/tp1-0.5*pi)
      else
         if(time.lt.tfor)then
         udt=ubwfor*sin(pi*time/tfor)
      else
         udt=-ubwback*sin(pi*(time-tfor)/tback)
      endif
      endif
      endif
      udtx=udt*cos(phi)
      udty=udt*sin(phi)
      ubtotx=ubtot*cos(phi)
      ubtoty=ubtot*sin(phi)
      uxt=vrdelm+udtx+ubtotx
      uyt=udty+ubtoty
      utvec=(uxt**2.+uyt**2.)**0.5
      tau1som=tau1som+0.5*rhow*fcw1*utvec*utvec
      tau1t=0.5*rhow*fcw1*utvec*utvec
      tt=(tau1t-taucr)/taucr
      tt=max(0.0001,tt)
      uster1t=(tau1t/rhow)**0.5
      coeff=0.25*(1.-(hs/hd)**0.5)
      if(coeff.le.0.075)coeff=0.075
      sbt=coeff*d50*rhos*uster1t*tt**1.5/(dster**0.3)
      sbx=sbx+uxt/utvec*sbt
      sby=sby+uyt/utvec*sbt
      tau1x=tau1x+uxt/utvec*tau1t
      tau1y=tau1y+uyt/utvec*tau1t
      nnn=nnn+1
      if(nnn.le.ntime)goto 150
      tau1m=tau1som/(ntime-1)
      tau1x=tau1x/(ntime-1)
      tau1y=tau1y/(ntime-1)
      tau1net=(tau1x**2.+tau1y**2.)**0.5
      phi1t=(180./pi)*atan(tau1y/tau1x)
      sbx=sbx/(ntime-1)
      sby=sby/(ntime-1)
c
c   TRANSPORT VECTORS
c
c   SUSPENDED LOAD TRANSPORT
c
      ssx=ssw*cos(phi)+ssc
      ssy=ssw*sin(phi)
      ssvec=(ssx**2.+ssy**2.)**0.5
      phis1=(180./pi)*atan(ssy/ssx)
      phis2=(180./pi)*phi-phis1
c
c   BED LOAD TRANSPORT
c
      sbvec=(sbx**2.+sby**2.)**0.5
      sbc=sbx-sby/tan(phi)
      sbw=sby/sin(phi)
      phib1=(180./pi)*atan(sby/sbx)
      phib2=(180./pi)*phi-phib1
c
c   TOTAL LOAD TRANSPORT
c
      stx=ssx+sbx
      sty=ssy+sby
      stvec=(stx**2.+sty**2.)**0.5
      phit1=(180./pi)*atan(sty/stx)
      phit2=(180./pi)*phi-phit1
      stc=stx-sty/tan(phi)
      stw=sty/sin(phi)
c
c   BED FORM COMPUTATION
c
      fr=vr/sqrt(g*hd)
      rm=ubw*ubw/(del*g*d50)
      if(fr.gt.1.)fr=1.
c
c   RIVERS AND ESTUARIES
c
      d1=0.
      d2=0.
      d3=0.
      d4=0.
      d5=0.
      d6=0.
      d7=0.
      d8=0.
      d9=0.
      d10=0.
      d11=0.
      d12=0.
      rl1=0.
      rl2=0.
      rl3=0.
      rl4=0.
      rl5=0.
      rl6=0.
      rl7=0.
      rl8=0.
      rl9=0.
      rl10=0.
      rl11=0.
      rl12=0.
      if(hs.le.0.01.and.dster.le.10..and. t.le.3.)then
        d1=100.*d50
        rl1=750.*d50
        d5=d1
        rl5=rl1
      endif
      if(hs.le.0.01.and.dster.le.10..and.t.gt.3..and.t.le.10.)then
        d2=0.02*hd*(10.-t)*(1.-exp(-0.1*t))
        rl2=0.5*hd
        d3=0.11*hd*(d50/hd)**0.3*(25.-t)*(1.-exp(-0.5*t))
        rl3=7.3*hd
        d6=d2
        rl6=rl2
        d7=0.05*hd
        rl7=2.*hd
      endif
      if(hs.le.0.01.and.t.gt.10..and.t.le.15.)then
        d3=0.11*hd*(d50/hd)**0.3*(25.-t)*(1.-exp(-0.5*t))
        rl3=7.3*hd
      endif
      if(hs.le.0.01.and.t.gt.15..and.t.le.25.)then
        d3=0.11*hd*(d50/hd)**0.3*(25.-t)*(1.-exp(-0.5*t))
        rl3=7.3*hd
        d4=0.15*hd*(1.-fr**2.)*(1.-exp(-0.5*(t-15.)))
        rl4=10.*hd
        d7=0.05*hd
        rl7=2.*hd
        d8=0.1*hd
        rl8=3.*hd
      endif
      if(hs.le.0.01.and.t.gt.25..and.fr.lt.0.8)then
        d4=0.15*hd*(1.-fr**2.)*(1.-exp(-0.5*(t-15.)))
        rl4=10.*hd
        d8=0.1*hd
        rl8=3.*hd
      endif
      if(hs.le.0.01.and.dster.gt.10..and.t.le.15.)then
        d3=0.11*hd*(d50/hd)**0.3*(25.-t)*(1.-exp(-0.5*t))
        rl3=7.3*hd
        d7=0.05*hd
        rl7=2.*hd
      endif
C
C    COASTAL SEA
C
      if(hs.gt.0.01.and.rm.le.10.)then
        d9=0.22*abw
        rl9=5.6*d9
      endif
      if(hs.gt.0.01.and.rm.gt.10..and.rm.le.250.)then
        d9=2.8e-13*(250.-rm)**5.*abw
        rl9=(1./(2.e-7*(250.-rm)**2.5))*d9
      endif
      if(hs.gt.0.01.and.rm.gt.250.)then
        d9=0.
        rl9=0.
      endif
      if(hs.gt.0.01.and.rm.le.10..and.vr.ge.0.3)then
        d9=0.22*abw
        rl9=5.6*d9
        d10=0.02*hd
        rl10=0.5*hd
      endif
      if(hs.gt.0.01.and.rm.gt.10..and.rm.le.250..and.vr.ge.0.3)then
        d9=2.8e-13*(250.-rm)**5.*abw
        rl9=(1./(2.e-7*(250.-rm)**2.5))*d9
        d10=0.02*hd
        rl10=0.5*hd
      endif
      if(hs.gt.0.01.and.rm.le.10..and.vr.ge.0.5.and.hd.gt.10.)then
        d9=0.22*abw
        rl9=5.6*d9
        d10=0.02*hd
        rl10=0.5*hd
        d12=0.15*hd
        rl12=10.*hd
      endif
      if(hs.gt.0.01.and.rm.gt.10..and.rm.le.250..and.vr.ge.0.5.and.
     *hd.gt.10.)then
        d9=2.8e-13*(250.-rm)**5.*abw
        rl9=(1./(2.e-7*(250.-rm)**2.5))*d9
        d10=0.02*hd
        rl10=0.5*hd
        d12=0.15*hd
        rl12=10.*hd
      endif
      if(hs.gt.0.01.and.rm.ge.250..and.vr.ge.0.5.and.hd.gt.10.)then
        d9=0.
        rl9=0.
        d10=0.02*hd
        rl10=0.5*hd
        d12=0.15*hd
        rl12=10.*hd
      endif
c   RESULTS OUTPUT TO FILE TRANSPOR.OUT
C
      WRITE(2,'(A)')CHAR(12)
      WRITE(2,'('' PHYSICAL PARAMETERS :'',/)')
      WRITE(2,'('' DSTER = PARTICLE PARAMETER               [   -  ] '',
     *E10.4)')DSTER
      WRITE(2,'('' L     = WAVE LENGTH                      [  M   ] '',
     *E10.4)')RLS
      WRITE(2,'('' TR    = RELATIVE WAVE PERIOD             [  S   ] '',
     *E10.4)')TP1
      WRITE(2,'('' UBW   = PEAK ORBITAL VELOCITY            [ M/S  ] '',
     *E10.4)')UBW
      WRITE(2,'('' UBWF  = PEAK ORBITAL VELOCITY FORWARD    [ M/S  ] '',
     *E10.4)')UBWFOR
      WRITE(2,'('' UBWB  = PEAK ORBITAL VELOCITY BACKWARD   [ M/S  ] '',
     *E10.4)')UBWBACK
      WRITE(2,'('' ASYM  = ASYMMETRY FACTOR PEAK ORB. VEL.  [  -   ] '',
     *E10.4)')ASYM
      WRITE(2,'('' ABW   = PEAK ORBITAL EXCURSION AT BED    [  M   ] '',
     *E10.4)')ABW
      WRITE(2,'('' DELW  = THICKNESS WAVE BOUNDARY LAYER    [  M   ] '',
     *E10.4)')DELW
      WRITE(2,'('' DELM  = THICKNESS WAVE MIXING LAYER      [  M   ] '',
     *E10.4)')DELM
      WRITE(2,'('' DS    = THICKNESS SEDIMENT MIXING LAYER  [  M   ] '',
     *E10.4)')DS
      WRITE(2,'('' TAUW  = WAVE-RELATED BED-SHEAR STRESS    [ N/M2 ] '',
     *E10.4)')TAUW
      WRITE(2,'('' TAUC  = CURRENT-RELATED BED-SHEAR STRESS [ N/M2 ] '',
     *E10.4)')TAUC
      WRITE(2,'('' FW    = WAVE-RELATED FRICTION COEFFICIENT[  -   ] '',
     *E10.4)')FW
      WRITE(2,'('' FW1   = WAVE-RELATED GRAIN FRIC. COEFF.  [  -   ] '',
     *E10.4)')FW1 
      WRITE(2,'('' FC    = CURR-RELATED FRICTION COEFFICIENT[  -   ] '',
     *E10.4)')FC
      WRITE(2,'('' FC1   = CURR-RELATED GRAIN FRIC.COEFF.   [  -   ] '',
     *E10.4)')FC1
      WRITE(2,'('' C     = CHEZY COEFFICIENT                [M0.5/S] '',
     *E10.4)')CC
      WRITE(2,'('' C1    = CHEZY GRAIN COEFFICIENT          [M0.5/S] '',
     *E10.4)')CC1
      WRITE(2,'('' RA    = APPARENT ROUGHNESS               [  M   ] '',
     *E10.4)')RA
      WRITE(2,'('' ALFAW = WAVE-CURRENT COEFFICIENT         [  -   ] '',
     *E10.4)')ALFAW
      WRITE(2,'('' RHOW  = FLUID DENSITY                    [KG/M3 ] '',
     *E10.4)')RHOW
      WRITE(2,'('' RHOS  = DENSITY SEDIMENT MATERIAL        [KG/M3 ] '',
     *E10.4)')RHOS
      WRITE(2,'('' WS    = FALL VELOCITY SUSP. MATERIAL     [ M/S  ] '',
     *E10.4)')WS
      WRITE(2,'('' WSB   = FALL VELOCITY BED MATERIAL       [ M/S  ] '',
     *E10.4)')WSB
      WRITE(2,'('' TAUCR = CRITICAL BED-SHEAR STRESS        [ N/M2 ] '',
     *E10.4)')TAUCR
      WRITE(2,'('' UCR   = CRITICAL DEPTH-AVERAGED VELOCITY [ M/S  ] '',
     *E10.4)')UMCR
      WRITE(2,'('' UBWCR = CRITICAL PEAK ORBITAL VELOCITY   [ M/S  ] '',
     *E10.4)')UBWCR
      WRITE(2,'('' UC    = C-RELATED EFFICIENCY FACTOR      [  -   ] '',
     *E10.4)')RMUC
      WRITE(2,'('' UWA   = W-RELATED EFF.FAC. CONCENTRATION [  -   ] '',
     *E10.4)')RMUWA
      WRITE(2,'('' UW    = W-RELATED EFF.FAC. BED LOAD TR.  [  -   ] '',
     *E10.4)')RMUW
      WRITE(2,'('' TA    = BED-SHEAR STRESS PAR.CONCENTR.   [  -   ] '',
     *E10.4)')TA
      WRITE(2,'('' T     = BED-SHEAR STRESS PAR.BED LOAD TR [  -   ] '',
     *E10.4)')T
      WRITE(2,'('' MFOR  = MOBILITY NUMBER FORWARD          [  -   ] '',
     *E10.4)')RMFOR
      WRITE(2,'('' MBACK = MOBILITY NUMBER BACKWARD         [  -   ] '',
     *E10.4)')RMBACK
      WRITE(2,'('' UR    = RETURN VELOCITY IN WAVE DIR.      [ M/S ] '',
     *E10.4)')UR
      WRITE(2,'('' UB    = NEAR-BED VELOCITY IN WAVE DIR.    [ M/S ] '',
     *E10.4)')UB
      WRITE(2,'('' VRDELM= VEL.AT EDGE W.MIX.LAYER, CURR DIR [ M/S ] '',
     *E10.4)')VRDELM
      WRITE(2,'('' URDELM= VEL.AT EDGE W.MIX.LAYER, WAVE DIR [ M/S ] '',
     *E10.4)')URDELM
      WRITE(2,'('' UBTOT = UB+URDELM=TOTAL VEL. AT EDGE W.M.L[ M/S ] '',
     *E10.4)')UBTOT
      WRITE(2,'(///)')
      WRITE(2,'(''1. BED LOAD TRANSPORT '')')
      WRITE(2,'(/)')
      WRITE(2,'(''   NET GRAIN BED SHEAR STRESS VECTOR    (N/M2)   ='',
     *E10.4)')TAU1NET
      WRITE(2,'(''   ANGLE BETWEEN CURRENT DIR. AND STRESS VECTOR  ='',
     *E10.4)')PHI1T
      WRITE(2,'(''   BED LOAD TRANSPORT VECTOR            (KG/SM)  ='',
     *E10.4)')SBVEC
      WRITE(2,'(''   ANGLE BETWEEN CURRENT DIR. AND B.L.TR. VECTOR ='',
     *E10.4)')PHIB1
      WRITE(2,'(''   ANGLE BETWEEN B.L.TR.VECTOR AND WAVE DIRECTION='',
     *E10.4)')PHIB2
      WRITE(2,'(''   BED L. TRANSP.COMPONENT IN CURR. DIR.(KG/SM)  ='',
     *E10.4)')SBC
      WRITE(2,'(''   BED L. TRANSP.COMPONENT IN WAVE DIR. (KG/SM)  ='',
     *E10.4)')SBW
      WRITE(2,'('' ( POS. ANGLE = ANTI-CLOCKWISE,NEG.=CLOCKWISE ) '')')
      WRITE(2,'(//)')
      WRITE(2,'(''2. SUSPENDED LOAD TRANSPORT '')')
      WRITE(2,'(/)')
      WRITE(2,'(''   SUSPENDED LOAD TRANSPORT VECTOR      (KG/SM)  ='',
     *E10.4)')SSVEC
      WRITE(2,'(''   ANGLE BETWEEN CURRENT DIR. AND S.L.TR. VECTOR ='',
     *E10.4)')PHIS1
      WRITE(2,'(''   ANGLE BETWEEN S.L.TR.VECTOR AND WAVE DIRECTION='',
     *E10.4)')PHIS2
      WRITE(2,'(''   SUSP. L. TR. COMPONENT IN CURR. DIR. (KG/SM)  ='',
     *E10.4)')SSC
      WRITE(2,'(''   SUSP. L. TR. COMPONENT IN WAVE DIR.  (KG/SM)  ='',
     *E10.4)')SSW
      WRITE(2,'(//)')
      WRITE(2,'(''3. TOTAL LOAD TRANSPORT '')')
      WRITE(2,'(/)')
      WRITE(2,'(''   TOTAL LOAD TRANSPORT VECTOR           (KG/SM) ='',
     *E10.4)')STVEC
      WRITE(2,'(''   ANGLE BETWEEN CURRENT DIR. AND T.L.TR.VECTOR  ='',
     *E10.4)')PHIT1
      WRITE(2,'(''   ANGLE BETWEEN T.L.TR.VECTOR AND WAVE DIRECTION='',
     *E10.4)')PHIT2
      WRITE(2,'(''   TOTAL L. TR. COMPONENT IN CURRENT DIR.(KG/SM) ='',
     *E10.4)')STC
      WRITE(2,'(''   TOTAL L. TR. COMPONENT IN WAVE DIR.   (KG/SM) ='',
     *E10.4)')STW
      WRITE(2,'(////)')   
      WRITE(2,'(''BED FORM HEIGHT AND LENGTH ( IN M.)'')')
      WRITE(2,'(1x)')
      line=' '
      line(1:10)='RIVER'
      line(15:28)='ESTUARY'
      line(29:40)='COASTAL SEA'
      write(2,'(22x,a)')line
      line=' '
      if(d1.ne.0.)write(line(1:10),'(e10.4)')d1
      if(d5.ne.0.)write(line(15:24),'(e10.4)')d5
      if(d9 .ne.0.)write(line(29:38),'(e10.4)')d9
      write(2,'('' Ripple height      = '',a50)')line
      line=' '
      if(rl1.ne.0.)write(line(1:10),'(e10.4)')rl1
      if(rl5.ne.0.)write(line(15:24),'(e10.4)')rl5
      if(rl9 .ne.0.)write(line(29:38),'(e10.4)')rl9
      write(2,'('' Ripple length      = '',a50)')line
      line=' '
      if(d2.ne.0.)write(line(1:10),'(e10.4)')d2
      if(d6.ne.0.)write(line(15:24),'(e10.4)')d6
      if(d10.ne.0.)write(line(29:38),'(e10.4)')d10
      write(2,'('' Mega ripple height = '',a50)')line
      line=' '
      if(rl2.ne.0.)write(line(1:10),'(e10.4)')rl2
      if(rl6.ne.0.)write(line(15:24),'(e10.4)')rl6
      if(rl10.ne.0.)write(line(29:38),'(e10.4)')rl10
      write(2,'('' Mega ripple length = '',a50)')line
      line=' '
      if(d3.ne.0.)write(line(1:10),'(e10.4)')d3
      if(d7.ne.0.)write(line(15:24),'(e10.4)')d7
      if(d11.ne.0.)write(line(29:38),'(e10.4)')d11
      write(2,'('' Dune height        = '',a50)')line
      line=' '
      if(rl3.ne.0.)write(line(1:10),'(e10.4)')rl3
      if(rl7.ne.0.)write(line(15:24),'(e10.4)')rl7
      if(rl11.ne.0.)write(line(29:38),'(e10.4)')rl11
      write(2,'('' Dune length        = '',a50)')line
      line=' '
      if(d4.ne.0.)write(line(1:10),'(e10.4)')d4
      if(d4.ne.0.)write(line(15:24),'(e10.4)')d8
      if(d12.ne.0.)write(line(29:38),'(e10.4)')d12
      write(2,'('' Sand wave height   = '',a50)')line
      line=' '
      if(rl4.ne.0.)write(line(1:10),'(e10.4)')rl4
      if(rl8.ne.0.)write(line(15:24),'(e10.4)')rl8
      if(rl12.ne.0.)write(line(29:38),'(e10.4)')rl12
      write(2,'('' Sand wave length   = '',a50)')line
c
      WRITE(2,'(A)')CHAR(12)
      WRITE(2,'('' PROFILES :'',//)')
      WRITE(2,'('' HEIGHT Z    CONCENTRATION    VELOCITY     TRANSPORT
     * VELOCITY     TRANSPORT'')')
      WRITE(2,'(''                                ( CURRENT DIRECTION)
     *  ( WAVE DIRECTION )'')')
      WRITE(2,'(''  [ M ]        [KG/M3]         [M/S]        [KG/SM2]
     *  [M/S]       [KG/SM2]'')')
      DO 200 I=1,IT-1
        WRITE(2,'(6(1X,E10.4,2x))')RMAT(I ,1),RMAT(I ,2),RMAT(I ,3),
     *  RMAT(I ,4),RMAT(I,5),RMAT(I,6)
  200 CONTINUE
C
C  RESULTS OUTPUT TO SCREEN
C
      rewind 2
      itel=0
  300 read(2,'(a)',end=350)line
      itel=itel+1
      if(mod(itel,18).eq.0)then
      write(*,'(a)')' Type C to continue.'
      read(*,'(a)')kar
      itel=0
      endif
      write(*,'(1x,a)')line(1:79)
      goto 300
  350 continue
      close(2)
      STOP 'NORMAL END TRANSPOR'
      END
      SUBROUTINE ZEROIN(X,Y,LOG,F,TOLX)
C***********************************************************************
C
C     ZEROIN LEVERT IN X DE BESTE BENADER ING VAN DE WORTEL VAN
C     DE VERGELIJKING F(X)=0. OP HET INTERVAL (X,Y) GEGEVEN IN DE BEGIN
C     WAARDEN VAN X ENY.
C     INDIEN ER GEEN WORTEL IS LOG=.FALSE.
C     ZEROIN ROEPT FLAMBDA EN TOLLD
C
C***********************************************************************
      LOGICAL LOG
      A=X
      FA=F(X)
      B=Y
      X=Y
      FB=F(X)
C INTERPOLATE
    1 C=A
      FC=FA
C EXTRAPOLATE:
    2 IF(ABS(FC).GE.ABS(FB))GOTO 19
      A=B
      FA=FB
      X=C
      B=C
      FB=FC
      C=A
      FC=FA
C END INTERCHANGE
   19 TOL=TOLX(X)
      RM=(C+B)*.5
      IF(ABS(RM-B).LE.TOL)GOTO 200
      P=(B-A)*FB
      IF(P.LT.0.)GOTO 5
      Q=FA-FB
      GOTO 6
    5 Q=FB-FA
      P=-P
    6 A=B
      FA=FB
      IF(P.GT.ABS(Q)*TOL)GOTO 8
      IF(C.LE.B)GOTO 10
      B=B+TOL
      X=B
      GOTO 11
   10 B=B-TOL
      X=B
   11 GOTO 12
    8 IF(P.GE.(RM-B)*Q)GOTO 14
      B=P/Q+B
      X=B
      GOTO 12
   14 B=RM
      X=B
   12 FB=F(X)
      IF(FC.LT.0.)GOTO 16
      IF(FB.LT.0.)GOTO 18
      GOTO 1
   18 GOTO 2
   16 IF(FB.GT.0.)GOTO 21
      GOTO 1
   21 GOTO 2
  200 Y=C
      IF(FC.LT.0.)GOTO 31
      IF(FB.GT.0.)GOTO 33
      LOG=.TRUE.
      GOTO 300
   33 LOG=.FALSE.
      GOTO 300
   31 IF(FB.LT.0.)GOTO 35
      LOG=.TRUE.
      GOTO 300
   35 LOG=.FALSE.
  300 RETURN
      END
      function fwl(rls)
      common/niorez/tp,vr,phi,g,pi,hd
      fwl=(rls/tp-vr*cos(phi))**2-g*rls/2./pi*tanh(2.*pi*hd/rls)
      return
      end
      function tolx(x)
      tolx=.001*abs(x)+1.e-5
      return
      end
