
        bits of java from              FoilSim III  - Airfoil  mode

                           A Java Applet
               to perform Kutta-Joukowski Airfoil analysis
                including drag from wind tunnel tests

                              Written by

                               Tom Benson
                       NASA Glenn Research Center

                                 and

                               Anthony Vila
                          Vanderbilt University

>                                NOTICE
>This software is in the Public Domain.  It may be freely copied and used in
>non-commercial products, assuming proper credit to the author is given.  IT
>MAY NOT BE RESOLD.


import java.awt.*
import java.lang.Math

public class Foil extends java.applet.Applet {

  dim convdr = 3.1415926/180. as double
   dim pid2 = 3.1415926/2.0 as double
   dim  rval,ycval,xcval,gamval,alfval,thkval,camval,chrd,clift as double
  dim dragCoeff,drag,liftOverDrag,reynolds,viscos as double
  dim alfd,thkd,camd,dragco as double
   dim thkinpt,caminpt as double            'MODS 10 Sep 99
   dim leg,teg,lem,tem as double
   dim usq,vsq,alt,altmax,area,armax,armin as double
  dim chord,span,aspr,arold,chrdold,spnold  as double ' Mod 13 Jan 00
   dim g0,q0,ps0,pt0,ts0,rho,rlhum,temf,presm as double
  dim lyg,lrg,lthg,lxgt,lygt,lrgt,lthgt as double ' MOD 20 Jul
   dim  lxm,lym,lxmt,lymt,vxdir as double ' MOD 20 Jul
  dim deltb,xflow as double             ' MODS  20 Jul 99
  dim delx,delt,vfsd,spin,spindr,yoff,radius as double
   dim vel,pres,lift,side,omega,radcrv,relsy,angr as double

   dim rg(20,40) as double     'do same for others
   dim thg(20,40) as double
   dim xg(20,40) as double
   dim yg(20,40) as double
   dim xm(20,40) as double
   dim ym(20,40) as double
   dim xpl(20,40) as double
   dim ypl(20,40) as double
   dim plp(40) as double
   dim plv(40) as double

   dim inptopt,outopt as integer   ' do same for others
   dim nptc,npt2,nlnc,nln2,rdflag,browflag,probflag,anflag as integer
   dim foil,flflag,lunits,lftout,planet,dragOut as integer
   dim displ,viewflg,dispp,dout,doutb,antim,ancol,sldloc as integer
   dim calcrange,arcor,indrag,recor,bdragflag as integer


       ' units data
  dim vmn,almn,angmn,vmx,almx,angmx as double
   dim camn,thkmn,camx,thkmx as double
   dim chrdmn,spanmn,armn,chrdmx,spanmx,armx as double
   dim radmn,spinmn,radmx,spinmx as  double

  dim vconv,vmax as double
   dim pconv,pmax,pmin,lconv,rconv,fconv,fmax,fmaxb as double
   dim lflag,gflag,plscale,nond as double
       '  plot & probe data
   dim fact,xpval,ypval,pbval,factp as double
   dim prg,pthg,pxg,pyg,pxm,pym,pxpl,pypl as double
   dim pboflag,xt,yt,ntikx,ntiky,npt,xtp,ytp as integer
   dim xt1,yt1,xt2,yt2,spanfac as integer
   dim lines,nord,nabs,ntr as integer
   dim begx,endx,begy,endy as double
   static String labx,labxu,laby,labyu
   dim pltx(3,40) as double
   dim plty(3,40) as double
   dim plthg(2) as double

   Solver solve
   Viewer view
   Con con
     Out out
   CardLayout layin,layout,layplt
   Image offImg1
   Graphics off1Gg
   Image offImg2
   Graphics off2Gg
   Image offImg3
   Graphics off3Gg

   private sub init()
     dim i as integer
     solve = new Solver()

     offImg1 = createImage(this.size().width,
                      this.size().height)
     off1Gg = offImg1.getGraphics()
     offImg2 = createImage(this.size().width,
                      this.size().height)
     off2Gg = offImg2.getGraphics()
     offImg3 = createImage(this.size().width,
                      this.size().height)
     off3Gg = offImg3.getGraphics()

     setLayout(new GridLayout(2,2,5,5))

     solve.setDefaults ()

     view  = new Viewer(this)
     con = new Con(this)
     out = new Out(this)

     add(view)

     solve.getFreeStream ()
     computeFlow ()
     view.start()
     out.plt.start()
  end sub

  private sub insets()
     return new Insets(10,10,10,10)
  }

  private sub computeFlow()

     if (flflag == 1) then
         solve.getFreeStream ()
         solve.getCirc ()                   ' get circulation
         solve.genFlow ()
     end if

     if (foil <= 3) then

         reynolds = vfsd/vconv * chord/lconv * rho / viscos

     else

         reynolds = vfsd/vconv * 2 * radius/lconv * rho / viscos
      end if

       thkd = thkinpt
     camd = caminpt
     alfd = alfval

     out.plt.loadPlot()
  end sub

  private sub filter0(double inumbr)
        '  output only to .
       dim number as integer
       dim intermed as integer

       number = (int) (inumbr)
       return number
  end sub

  private sub filter1(double inumbr)
     '  output only to .1
       float number
       int intermed

       intermed = (int) (inumbr * 10.)
       number = (float) (intermed / 10. )
       return number
  end sub

  private sub filter3(double inumbr)
     '  output only to .001
       float number
       dim intermed as integer

       intermed = (int) (inumbr * 1000.)
       number = (float) (intermed / 1000. )
       return number
  end sub


  private sub setUnits()    ' Switching Units
       dim ovs,chords,spans,aros,chos,spos,rads as double
       dim alts,ares as double

       alts = alt / lconv
       chords = chord / lconv
       spans = span / lconv
       ares = area /lconv/lconv
       aros = arold /lconv/lconv
       chos = chrdold / lconv
       spos = spnold / lconv
       ovs = vfsd / vconv
       rads = radius / lconv

       switch (lunits) {
          case 0: {                             ' English
            lconv = 1.                      '  feet
            vconv = .6818 vmax = 250.   '  mph
            fconv = 1.0 fmax = 100000. fmaxb = .5  ' pounds
            pconv = 14.7                     ' lb/sq in
            break
          }
       }

       alt = alts * lconv
       chord = chords * lconv
       span = spans * lconv
       area = ares * lconv * lconv
       arold = aros * lconv * lconv
       chrdold = chos * lconv
       spnold = spos * lconv
       vfsd  = ovs * vconv
       radius  = rads * lconv

       return
  end sub

  private sub loadInput()    ' load the input panels

       computeFlow()
       return
  end sub


  class Solver {    '????????????

     Solver () {
     }

     private sub setDefaults()

        dragOut = 0
        arcor = 1
        indrag = 1
        recor = 1
        bdragflag = 1  ' smooth ball
        planet = 0
        lunits = 0
        lftout = 0
        inptopt = 0
        outopt = 0
        nlnc = 15
        nln2 = nlnc/2 + 1
        nptc = 37
        npt2 = nptc/2 + 1
        deltb = .5
        foil = 1
        flflag = 1
        thkval = .5
        thkinpt = 12.5                    ' MODS 10 SEP 99
        camval = 0.0
        caminpt = 0.0
        alfval = 5.0
        gamval = 0.0
        radius = 1.0
        spin = 0.0
        spindr = 1.0
        rval = 1.0
        ycval = 0.0
        xcval = 0.0
        displ   = 1
        viewflg = 0
        dispp = 20
        calcrange = 0
        dout = 0
        doutb = 0

        dragCoeff = 0

        xpval = 2.1
        ypval = -.5
        pboflag = 0
        xflow = -10.0                             ' MODS  20 Jul 99

        pconv = 14.7
        pmin = .5
        pmax = 1.0
        fconv = 1.0
        fmax = 100000.
        fmaxb = .50
        vconv = .6818
        vfsd = 100.
        vmax = 250.
        lconv = 1.0

        alt = 0.0
        altmax = 50000.
        chrdold = chord = 5.0
        spnold = span = 20.0
        aspr = 4.0
        arold = area = 100.0
        armax = 2500.01
        armin = .01                  ' MODS 9 SEP 99

        xt = 170  yt = 105 fact = 40.0
        sldloc = 50
        xtp = 95 ytp = 165 factp = 30.0
        spanfac = (int)(2.0*fact*aspr*.3535)
        xt1 = xt + spanfac
        yt1 = yt - spanfac
        xt2 = xt - spanfac
        yt2 = yt + spanfac
        plthg(1) = 0.0

        probflag = 0
        anflag = 1
        vmn = 0.0     vmx = 250.0
        almn = 0.0    almx = 50000.0
        angmn = -20.0 angmx = 20.0
        camn = -20.0  camx = 20.0
        thkmn = 1.0 thkmx = 20.0
        chrdmn = .1   chrdmx = 20.1
        spanmn = .1   spanmx = 125.1
        armn = .1   armx = 2500.1
        spinmn = -1500.0   spinmx = 1500.0
        radmn = .05   radmx = 5.0

        return
     end sub

     private sub getFreeStream()     '  free stream conditions
       dim hite,pvap,rgas,gama,mu0 as double        ' MODS  19 Jan 00  whole routine

       g0 = 32.2
       rgas = 1718.                 ' ft2/sec2 R
       gama = 1.4
       hite = alt/lconv
       mu0 = .000000362
            ts0 = 518.6 - 3.56 * hite/1000.
            ps0 = 2116. * Math.pow(ts0/518.6,5.256)


       rho = ps0/(rgas * ts0)
       viscos = mu0 * 717.408/(ts0 + 198.72)*Math.pow(ts0/518.688,1.5)

       q0  = .5 * rho * vfsd * vfsd / (vconv * vconv)
       pt0 = ps0 + q0

       return
     end sub

     private sub getCirc()    ' circulation from Kutta condition
       double thet,rdm,thtm
       double beta
       int index

       xcval = 0.0
               ' Juokowski geometry
              ycval = camval / 2.0
              rval = thkval/4.0 +Math.sqrt(thkval*thkval/16.0+ycval*ycval +1.0)
              xcval = 1.0 - Math.sqrt(rval*rval - ycval*ycval)
              beta = Math.asin(ycval/rval)/convdr      ' Kutta condition
              gamval = 2.0*rval*Math.sin((alfval+beta)*convdr)

                             ' geometry
       for (index =1 index <= nptc ++index) {
           thet = (index -1)*360./(nptc-1)
           xg(0,index) = rval * Math.cos(convdr * thet) + xcval
           yg(0,index) = rval * Math.sin(convdr * thet) + ycval
           rg(0,index) = Math.sqrt(xg(0,index)*xg(0,index) +
                                yg(0,index)*yg(0,index))
           thg(0,index) = Math.atan2(yg(0,index),xg(0,index))/convdr
           xm(0,index) = (rg(0,index) + 1.0/rg(0,index))*
                    Math.cos(convdr*thg(0,index))
           ym(0,index) = (rg(0,index) - 1.0/rg(0,index))*
                    Math.sin(convdr*thg(0,index))
           rdm = Math.sqrt(xm(0,index)*xm(0,index) +
                           ym(0,index)*ym(0,index))
           thtm = Math.atan2(ym(0,index),xm(0,index))/convdr
           xm(0,index) = rdm * Math.cos((thtm - alfval)*convdr)
           ym(0,index) = rdm * Math.sin((thtm - alfval)*convdr)
           getVel(rval,thet)
           plp(index) = ((ps0 + pres * q0)/2116.) * pconv
           plv(index) = vel * vfsd
       }

       xt1 = xt + spanfac
       yt1 = yt - spanfac
       xt2 = xt - spanfac
       yt2 = yt + spanfac

       return
     end sub

     private sub genFlow()    ' generate flowfield
       double rnew,thet,psv,fxg
       int k,index
                              ' all lines of flow  except stagnation line
       for (k=1 k<=nlnc ++k) {
         psv = -.5*(nln2-1) + .5*(k-1)
         fxg = xflow
         for (index =1 index <=nptc ++ index) {
           solve.getPoints (fxg,psv)
           xg(k,index)  = lxgt
           yg(k,index)  = lygt
           rg(k,index)  = lrgt
           thg(k,index) = lthgt
           xm(k,index)  = lxmt
           ym(k,index)  = lymt

           solve.getVel(lrg,lthg)
           fxg = fxg + vxdir*deltb
         }
       }
                                            '  stagnation line
       k = nln2
       psv = 0.0
                                              '  incoming flow
       for (index =1 index <= npt2 ++ index) {
           rnew = 10.0 - (10.0 - rval)*Math.sin(pid2*(index-1)/(npt2-1))
           thet = Math.asin(.999*(psv - gamval*Math.log(rnew/rval))/
                                   (rnew - rval*rval/rnew))
           fxg =  - rnew * Math.cos(thet)
           solve.getPoints (fxg,psv)
           xg(k,index)  = lxgt
           yg(k,index)  = lygt
           rg(k,index)  = lrgt
           thg(k,index) = lthgt
           xm(k,index)  = lxmt
           ym(k,index)  = lymt
       }
                                              '  downstream flow
       for (index = 1 index <= npt2 ++ index) {
           rnew = 10.0 + .01 - (10.0 - rval)*Math.cos(pid2*(index-1)/(npt2-1))
           thet = Math.asin(.999*(psv - gamval*Math.log(rnew/rval))/
                                      (rnew - rval*rval/rnew))
           fxg =   rnew * Math.cos(thet)
           solve.getPoints (fxg,psv)
           xg(k,npt2+index)  = lxgt
           yg(k,npt2+index)  = lygt
           rg(k,npt2+index)  = lrgt
           thg(k,npt2+index) = lthgt
           xm(k,npt2+index)  = lxmt
           ym(k,npt2+index)  = lymt
       }
                                              '  stagnation point
       xg(k,npt2)  = xcval
       yg(k,npt2)  = ycval
       rg(k,npt2)  = Math.sqrt(xcval*xcval+ycval*ycval)
       thg(k,npt2) = Math.atan2(ycval,xcval)/convdr
       xm(k,npt2)  = (xm(k,npt2+1) + xm(k,npt2-1))/2.0
       ym(k,npt2)  = (ym(0,nptc/4+1) + ym(0,nptc/4*3+1))/2.0
                                '  compute lift coefficient
       leg = xcval - Math.sqrt(rval*rval - ycval*ycval)
       teg = xcval + Math.sqrt(rval*rval - ycval*ycval)
       lem = leg + 1.0/leg
       tem = teg + 1.0/teg
       chrd = tem - lem
       clift = gamval*4.0*3.1415926/chrd

       return
     end sub

     private sub getPoints(double fxg, double psv)    ' flow in x-psi
       double radm,thetm
       double fnew,ynew,yold,rfac,deriv
       double xold,xnew,thet
       double rmin,rmax
       int iter,isign

       ynew = 10.0
       yold = 10.0
       if (psv < 0.0) ynew = -10.0
       if (Math.abs(psv) < .001 && alfval < 0.0) ynew = rval
       if (Math.abs(psv) < .001 && alfval >= 0.0) ynew = -rval
       fnew = 0.1
       iter = 1
       while (Math.abs(fnew) >= .00001 && iter < 25) {
           ++iter
           rfac = fxg*fxg + ynew*ynew
           if (rfac < rval*rval) rfac = rval*rval + .01
           fnew = psv - ynew*(1.0 - rval*rval/rfac)
                  - gamval*Math.log(Math.sqrt(rfac)/rval)
           deriv = - (1.0 - rval*rval/rfac)
               - 2.0 * ynew*ynew*rval*rval/(rfac*rfac)
               - gamval * ynew / rfac
           yold = ynew
           ynew = yold  - .5*fnew/deriv
       }
       lyg = yold
                                     ' rotate for angle of attack
       lrg = Math.sqrt(fxg*fxg + lyg*lyg)
       lthg = Math.atan2(lyg,fxg)/convdr
       lxgt = lrg * Math.cos(convdr*(lthg + alfval))
       lygt = lrg * Math.sin(convdr*(lthg + alfval))
                              ' translate cylinder to generate airfoil
       lxgt = lxgt + xcval
       lygt = lygt + ycval
       lrgt = Math.sqrt(lxgt*lxgt + lygt*lygt)
       lthgt = Math.atan2(lygt,lxgt)/convdr
                               '  Kutta-Joukowski mapping
       lxm = (lrgt + 1.0/lrgt)*Math.cos(convdr*lthgt)
       lym = (lrgt - 1.0/lrgt)*Math.sin(convdr*lthgt)
                              ' tranforms for view fixed with free stream
                ' take out rotation for angle of attack mapped and cylinder
       radm = Math.sqrt(lxm*lxm+lym*lym)
       thetm = Math.atan2(lym,lxm)/convdr
       lxmt = radm*Math.cos(convdr*(thetm-alfval))
       lymt = radm*Math.sin(convdr*(thetm-alfval))

       lxgt = lxgt - xcval
       lygt = lygt - ycval
       lrgt = Math.sqrt(lxgt*lxgt + lygt*lygt)
       lthgt = Math.atan2(lygt,lxgt)/convdr
       lxgt = lrgt * Math.cos((lthgt - alfval)*convdr)
       lygt = lrgt * Math.sin((lthgt - alfval)*convdr)

       return
     end sub

     private sub getVel(double rad, double theta)   'velocity and pressure
      double ur,uth,jake1,jake2,jakesq
      double xloc,yloc,thrad,alfrad

      thrad = convdr * theta
      alfrad = convdr * alfval
                                ' get x, y location in cylinder plane
      xloc = rad * Math.cos(thrad)
      yloc = rad * Math.sin(thrad)
                                ' velocity in cylinder plane
      ur  = Math.cos(thrad-alfrad)*(1.0-(rval*rval)/(rad*rad))
      uth = -Math.sin(thrad-alfrad)*(1.0+(rval*rval)/(rad*rad))
                            - gamval/rad
      usq = ur*ur + uth*uth
      vxdir = ur * Math.cos(thrad) - uth * Math.sin(thrad)
                                ' translate to generate airfoil
      xloc = xloc + xcval
      yloc = yloc + ycval

      rad = Math.sqrt(xloc*xloc + yloc*yloc)
      thrad  = Math.atan2(yloc,xloc)
                                   ' compute Joukowski Jacobian
      jake1 = 1.0 - Math.cos(2.0*thrad)/(rad*rad)
      jake2 = Math.sin(2.0*thrad)/(rad*rad)
      jakesq = jake1*jake1 + jake2*jake2
      if (Math.abs(jakesq) <= .01) jakesq = .01   ' protection
      vsq = usq / jakesq
          ' vel is velocity ratio - pres is coefficient  (p-p0)/q0

           vel = Math.sqrt(vsq)
           pres = 1.0 - vsq

          return
    end sub



  } ' end Solver class

  class Con extends Panel {
     Foil outerparent

     Con (Foil target) {

     }


  } ' Con


  class Out extends Panel {
     Foil outerparent
     Plt plt

     Out (Foil target) {
        outerparent = target
        layout = new CardLayout()
        setLayout(layout)

        plt = new Plt(outerparent)

        add ("first", plt)
       }

     class Plt extends Canvas
         implements Runnable{
        Foil outerparent
        Thread run2
        Point locp,ancp

        Plt (Foil target) {
           setBackground(Color.blue)
           run2 = null
        }

        private sub start()
           if (run2 == null) then
              run2 = new Thread(this)
              run2.start()
           end if
        end sub

        private sub run()
          int timer

          timer = 100
          while (true) {
             try { Thread.sleep(timer) }
             catch (InterruptedException e) {}
             out.plt.repaint()
          }
        end sub

        private sub loadPlot()
          double rad,ang,xc,yc,lftref,clref,drgref,cdref
          double del,spd,awng,ppl,tpl,hpl,angl,thkpl,campl,clpl,cdpl
          int index,ic

          lines = 1
          clref =  getClplot(camval,thkval,alfval)
          if (Math.abs(clref) <= .001) clref = .001     ' protection
          lftref = clref * q0 * area/lconv/lconv
          alfd = alfval
          thkd = thkinpt
          camd = caminpt
          cdref = dragco
          drgref = cdref * q0 * area/lconv/lconv

' load up the view image
          for (ic = 0 ic <= nlnc ++ ic) {
             for (index = 0 index <= nptc ++ index) {
                   xpl(ic,index) = xm(ic,index)
                   ypl(ic,index) = ym(ic,index)

             }
          }

         end sub

        public double getClplot (double camb, double thic, double angl) {
           double beta,xc,yc,rc,gamc,lec,tec,lecm,tecm,crdc
           double stfact,number

           xc = 0.0
           yc = camb / 2.0
           rc = thic/4.0 + Math.sqrt( thic*thic/16.0 + yc*yc + 1.0)
           xc = 1.0 - Math.sqrt(rc*rc - yc*yc)
           beta = Math.asin(yc/rc)/convdr        ' Kutta condition
           gamc = 2.0*rc*Math.sin((angl+beta)*convdr)
           lec = xc - Math.sqrt(rc*rc - yc*yc)
           tec = xc + Math.sqrt(rc*rc - yc*yc)
           lecm = lec + 1.0/lec
           tecm = tec + 1.0/tec
           crdc = tecm - lecm
                                      ' stall model 1
           stfact = 1.0

           number = stfact*gamc*4.0*3.1415926/crdc

           if (arcor == 1) then  ' correction for low aspect ratio
               number = number /(1.0 + number/(3.14159*aspr))
           end if

           return (number)
        }

        private sub update(Graphics g)
           out.plt.paint(g)
        end sub

        private sub paint(Graphics g)
           dim i,j,k,n,index as integer
           dim xlabel,ylabel,ind,inmax,inmin as integer
           dim  exes(8)  as integer
           dim whys(8) as integer
           double offx,scalex,offy,scaley,waste,incy,incx
           double xl,yl
           double liftab,dragab
           dim camx(19) as integer
           dim camy(19) as integer
           Color col

           if (ntikx < 2) ntikx = 2      ' protection 13June96
           if (ntiky < 2) ntiky = 2
           offx = 0.0 - begx
           scalex = 6.0/(endx-begx)
           incx = (endx-begx)/(ntikx-1)
           offy = 0.0 - begy
           scaley = 4.5/(endy-begy)
           incy = (endy-begy)/(ntiky-1)


                g.drawImage(offImg2,0,0,this)
       end sub
     }     ' Plt

      } ' Out

  class Viewer extends Canvas
         implements Runnable{
     Foil outerparent
     Thread runner
     Point locate,anchor

     Viewer (Foil target) {
         setBackground(Color.black)
         runner = null
     }

     public Insets insets() {
        return new Insets(0,10,0,10)
     }

     public boolean mouseDown(Event evt, int x, int y) {
        anchor = new Point(x,y)
        return true
     }

     private sub start()
        if (runner == null) {
           runner = new Thread(this)
           runner.start()
        }
        antim = 0
        ancol = 1
     end sub

     private sub run()
       int timer

       timer = 100
       while (true) {
          ++ antim
          try { Thread.sleep(timer) }
          catch (InterruptedException e) {}
          view.repaint()
          if (antim == 3) {
             antim = 0
             ancol = - ancol
          }
          timer = 135 - (int) (.227 *vfsd/vconv)
        }
     end sub

     private sub update(Graphics g)
        view.paint(g)
     end sub

     private sub paint(Graphics g)
        dim i,j,k,n as integer
        dim xlabel,ylabel,ind,inmax,inmin as integer
        dim exes(8) as integer
        dim whys(8) as integer
        dim offx,scalex,offy,scaley,waste,incy,incx as double
        dim xl,yl,slope,radvec,xvec,yvec as double
        dim camx(19) as integer
        dim camy(19) as integer
        Color col

        col = new Color(0,0,0)
        if(planet == 0) col = Color.cyan
        off1Gg.setColor(Color.black)
        off1Gg.fillRect(0,0,500,500)

     '   do NOT need viewflg 1 cases             ' Top View
     'lots deleted


        if (viewflg == 0 || viewflg == 2) then   ' edge View
         if (vfsd > .01) then
                                            ' plot airfoil flowfield
          radvec = .5
          for (j=1 j<=nln2-1 ++j) {           ' lower half (AIRFLOW)
             for (i=1  i<= nptc-1 ++i) {
                exes(0) = (int) (fact*xpl(j)(i)) + xt
                whys(0) = (int) (fact*(-ypl(j)(i))) + yt
                slope = (ypl(j)(i+1)-ypl(j)(i))/(xpl(j)(i+1)-xpl(j)(i))
                xvec = xpl(j)(i) + radvec / Math.sqrt(1.0 + slope*slope)
                yvec = ypl(j)(i) + slope * (xvec - xpl(j)(i))
                exes(1) = (int) (fact*xvec) + xt
                whys(1) = (int) (fact*(-yvec)) + yt

                if (displ == 2 && (i/3*3 == i) ) then
                  off1Gg.setColor(col)
                  for (n=1  n <= 4  ++n) {
                     if(i == 6 + (n-1)*9) off1Gg.setColor(Color.yellow)
                  }
                  if(i/9*9 == i) off1Gg.setColor(Color.white)
                  off1Gg.drawLine(exes(0),whys(0),exes(1),whys(1))
                end if
                if (displ == 1 && ((i-antim)/3*3 == (i-antim)) ) {
                  if (ancol == -1) {          ' MODS  27 JUL 99
                    if((i-antim)/6*6 == (i-antim))off1Gg.setColor(col)
                    if((i-antim)/6*6 != (i-antim))off1Gg.setColor(Color.white)
                  }
                  if (ancol == 1) {          ' MODS  27 JUL 99
                    if((i-antim)/6*6 == (i-antim))off1Gg.setColor(Color.white)
                    if((i-antim)/6*6 != (i-antim))off1Gg.setColor(col)
                  }
                  off1Gg.drawLine(exes(0),whys(0),exes(1),whys(1))
                }
             }
          }

          off1Gg.setColor(Color.white)  ' stagnation
          exes(1) = (int) (fact*xpl(nln2)(1)) + xt
          whys(1) = (int) (fact*(-ypl(nln2)(1))) + yt
          for (i=2  i<= npt2-1 ++i) {
                exes(0) = exes(1)
                whys(0) = whys(1)
                exes(1) = (int) (fact*xpl(nln2)(i)) + xt
                whys(1) = (int) (fact*(-ypl(nln2)(i))) + yt
                if (displ <= 2) {             ' MODS  21 JUL 99
                  off1Gg.drawLine(exes(0),whys(0),exes(1),whys(1))
                }
          }
          exes(1) = (int) (fact*xpl(nln2)(npt2+1)) + xt
          whys(1) = (int) (fact*(-ypl(nln2)(npt2+1))) + yt
          for (i=npt2+2 to nptc) 
                exes(0) = exes(1)
                whys(0) = whys(1)
                exes(1) = (int) (fact*xpl(nln2)(i)) + xt
                whys(1) = (int) (fact*(-ypl(nln2)(i))) + yt
                if (displ <= 2) then                         ' MODS  21 JUL 99
                  off1Gg.drawLine(exes(0),whys(0),exes(1),whys(1))
                end if
          next i


            for (j=nln2+1 j<=nlnc ++j)           ' upper half (AIRFLOW)
             for (i=1  i<= nptc-1 ++i) 
                exes(0) = (int) (fact*xpl(j)(i)) + xt
                whys(0) = (int) (fact*(-ypl(j)(i))) + yt
                slope = (ypl(j)(i+1)-ypl(j)(i))/(xpl(j)(i+1)-xpl(j)(i))
                xvec = xpl(j,i) + radvec / Math.sqrt(1.0 + slope*slope)
                yvec = ypl(j)(i) + slope * (xvec - xpl(j)(i))
                exes(1) = (int) (fact*xvec) + xt
                whys(1) = (int) (fact*(-yvec)) + yt
                if (displ == 0) {
                  off1Gg.setColor(col)
                  exes(1) = (int) (fact*xpl(j)(i+1)) + xt
                  whys(1) = (int) (fact*(-ypl(j)(i+1))) + yt
                  off1Gg.drawLine(exes(0),whys(0),exes(1),whys(1))
                }
                if (displ == 2 && (i/3*3 == i) ) {
                  off1Gg.setColor(col)
                  for (n=1  n <= 4  ++n) {
                     if(i == 6 + (n-1)*9) off1Gg.setColor(Color.yellow)
                  }
                  if(i/9*9 == i) off1Gg.setColor(Color.white)
                  off1Gg.drawLine(exes(0),whys(0),exes(1),whys(1))
                }
                if (displ == 1 && ((i-antim)/3*3 == (i-antim)) ) then
                  if (ancol == -1) {
                    if((i-antim)/6*6 == (i-antim))off1Gg.setColor(col)
                    if((i-antim)/6*6 != (i-antim))off1Gg.setColor(Color.white)
                  }
                  if (ancol == 1) then
                    if((i-antim)/6*6 == (i-antim)) then 
						  off1Gg.setColor(Color.white)
						  end if
                    if((i-antim)/6*6 != (i-antim))then 
						  off1Gg.setColor(col)
						  end if
                  end if
                  off1Gg.drawLine(exes(0),whys(0),exes(1),whys(1))
                end if
             next i
          next j
        end if

         if (viewflg == 0) then
  ' draw the airfoil geometry
             off1Gg.setColor(Color.white)
             exes(1) = (int) (fact*(xpl(0,npt2))) + xt
             whys(1) = (int) (fact*(-ypl(0,npt2))) + yt
             exes(2) = (int) (fact*(xpl(0,npt2))) + xt
             whys(2) = (int) (fact*(-ypl(0,npt2))) + yt
             for (i=1  i<= npt2-1 ++i) 
                exes(0) = exes(1)
                whys(0) = whys(1)
                exes(1) = (int) (fact*(xpl(0,npt2-i))) + xt
                whys(1) = (int) (fact*(-ypl(0,npt2-i))) + yt
                exes(3) = exes(2)
                whys(3) = whys(2)
                exes(2) = (int) (fact*(xpl(0,npt2+i))) + xt
                whys(2) = (int) (fact*(-ypl(0,npt2+i))) + yt
                camx(i) = (exes(1) + exes(2)) / 2
                camy(i) = (whys(1) + whys(2)) / 2
                off1Gg.setColor(Color.white)
                off1Gg.fillPolygon(exes,whys,4)    ' THIS DRAWS AIRFOIL
                next i
             end if
          end if

        g.drawImage(offImg1,0,0,this)      ' draws airflow!!!!

     end sub ' end paint

  } ' end Viewer

} ' end all
