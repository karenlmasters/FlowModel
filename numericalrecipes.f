CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C NUMERICAL RECIPES ROUTINES
      SUBROUTINE lfit(x,y,sig,ndat,a,ia,ma,covar,npc,chisq,funcs)
      INTEGER ma,ia(ma),npc,ndat,MMAX
      REAL chisq,a(ma),covar(npc,npc),sig(ndat),x(ndat),y(ndat)
      EXTERNAL funcs
      PARAMETER (MMAX=50)
CU    USES covsrt,gaussj
      INTEGER i,j,k,l,m,mfit
      REAL sig2i,sum,wt,ym,afunc(MMAX),beta(MMAX)
      mfit=0
      do 11 j=1,ma
        if(ia(j).ne.0) mfit=mfit+1
11    continue
      if(mfit.eq.0) pause 'lfit: no parameters to be fitted'
      do 13 j=1,mfit
        do 12 k=1,mfit
          covar(j,k)=0.
12      continue
        beta(j)=0.
13    continue
      do 17 i=1,ndat
        call funcs(x(i),afunc,ma)
        ym=y(i)
        if(mfit.lt.ma) then
          do 14 j=1,ma
            if(ia(j).eq.0) ym=ym-a(j)*afunc(j)
14        continue
        endif
        sig2i=1./sig(i)**2
        j=0
        do 16 l=1,ma
          if (ia(l).ne.0) then
            j=j+1
            wt=afunc(l)*sig2i
            k=0
            do 15 m=1,l
              if (ia(m).ne.0) then
                k=k+1
                covar(j,k)=covar(j,k)+wt*afunc(m)
              endif
15          continue
            beta(j)=beta(j)+ym*wt
          endif
16      continue
17    continue
      do 19 j=2,mfit
        do 18 k=1,j-1
          covar(k,j)=covar(j,k)
18      continue
19    continue
      call gaussj(covar,mfit,npc,beta,1,1)
      j=0
      do 21 l=1,ma
        if(ia(l).ne.0) then
          j=j+1
          a(l)=beta(j)
        endif
21    continue
      chisq=0.
      do 23 i=1,ndat
        call funcs(x(i),afunc,ma)
        sum=0.
        do 22 j=1,ma
          sum=sum+a(j)*afunc(j)
22      continue
        chisq=chisq+((y(i)-sum)/sig(i))**2
23    continue
      call covsrt(covar,npc,ma,ia,mfit)
      return
      END

      SUBROUTINE fit(x,y,ndata,sig,mwt,a,b,siga,sigb,chi2,q)
      INTEGER mwt,ndata
      REAL a,b,chi2,q,siga,sigb,sig(ndata),x(ndata),y(ndata)
CU    USES gammq
      INTEGER i
      REAL sigdat,ss,st2,sx,sxoss,sy,t,wt,gammq
      sx=0.
      sy=0.
      st2=0.
      b=0.
      if(mwt.ne.0) then
        ss=0.
        do 11 i=1,ndata
          wt=1./(sig(i)**2)
          ss=ss+wt
          sx=sx+x(i)*wt
          sy=sy+y(i)*wt
11      continue
      else
        do 12 i=1,ndata
          sx=sx+x(i)
          sy=sy+y(i)
12      continue
        ss=float(ndata)
      endif
      sxoss=sx/ss
      if(mwt.ne.0) then
        do 13 i=1,ndata
          t=(x(i)-sxoss)/sig(i)
          st2=st2+t*t
          b=b+t*y(i)/sig(i)
13      continue
      else
        do 14 i=1,ndata
          t=x(i)-sxoss
          st2=st2+t*t
          b=b+t*y(i)
14      continue
      endif
      b=b/st2
      a=(sy-sx*b)/ss
      siga=sqrt((1.+sx*sx/(ss*st2))/ss)
      sigb=sqrt(1./st2)
      chi2=0.
      q=1.
      if(mwt.eq.0) then
        do 15 i=1,ndata
          chi2=chi2+(y(i)-a-b*x(i))**2
15      continue
        sigdat=sqrt(chi2/(ndata-2))
        siga=siga*sigdat
        sigb=sigb*sigdat
      else
        do 16 i=1,ndata
          chi2=chi2+((y(i)-a-b*x(i))/sig(i))**2
16      continue
        if(ndata.gt.2) q=gammq(0.5*(ndata-2),0.5*chi2)
      endif
      return
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE fitexy(x,y,ndat,sigx,sigy,a,b,siga,sigb,chi2,q)
      INTEGER ndat,NMAX
      REAL x(ndat),y(ndat),sigx(ndat),sigy(ndat),a,b,siga,sigb,chi2,q,
     *POTN,PI,BIG,ACC
      PARAMETER (NMAX=1000,POTN=1.571000,BIG=1.e30,PI=3.14159265,
     *ACC=1.e-3)
CU    USES avevar,brent,chixy,fit,gammq,mnbrak,zbrent
      INTEGER j,nn
      REAL xx(NMAX),yy(NMAX),sx(NMAX),sy(NMAX),ww(NMAX),swap,amx,amn,
     *varx,vary,aa,offs,ang(6),ch(6),scale,bmn,bmx,d1,d2,r2,dum1,dum2,
     *dum3,dum4,dum5,brent,chixy,gammq,zbrent
      COMMON /fitxyc/ xx,yy,sx,sy,ww,aa,offs,nn
      EXTERNAL chixy
      if (ndat.gt.NMAX) pause 'NMAX too small in fitexy'
      call avevar(x,ndat,dum1,varx)
      call avevar(y,ndat,dum1,vary)
      scale=sqrt(varx/vary)
      nn=ndat
      do 11 j=1,ndat
        xx(j)=x(j)
        yy(j)=y(j)*scale
        sx(j)=sigx(j)
        sy(j)=sigy(j)*scale
        ww(j)=sqrt(sx(j)**2+sy(j)**2)
11    continue
      call fit(xx,yy,nn,ww,1,dum1,b,dum2,dum3,dum4,dum5)
      offs=0.
      ang(1)=0.
      ang(2)=atan(b)
      ang(4)=0.
      ang(5)=ang(2)
      ang(6)=POTN
      do 12 j=4,6
        ch(j)=chixy(ang(j))
12    continue
      call mnbrak(ang(1),ang(2),ang(3),ch(1),ch(2),ch(3),chixy)
      chi2=brent(ang(1),ang(2),ang(3),chixy,ACC,b)
      chi2=chixy(b)
      a=aa
      q=gammq(0.5*(nn-2),0.5*chi2)
      r2=0.
      do 13 j=1,nn
        r2=r2+ww(j)
13    continue
      r2=1./r2
      bmx=BIG
      bmn=BIG
      offs=chi2+1.
      do 14 j=1,6
        if (ch(j).gt.offs) then
          d1=mod(abs(ang(j)-b),PI)
          d2=PI-d1
          if(ang(j).lt.b)then
            swap=d1
            d1=d2
            d2=swap
          endif
          if (d1.lt.bmx) bmx=d1
          if (d2.lt.bmn) bmn=d2
        endif
14    continue
      if (bmx.lt. BIG) then
        bmx=zbrent(chixy,b,b+bmx,ACC)-b
        amx=aa-a
        bmn=zbrent(chixy,b,b-bmn,ACC)-b
        amn=aa-a
        sigb=sqrt(0.5*(bmx**2+bmn**2))/(scale*cos(b)**2)
        siga=sqrt(0.5*(amx**2+amn**2)+r2)/scale
      else
        sigb=BIG
        siga=BIG
      endif
      a=a/scale
      b=tan(b)/scale
      return
      END

      FUNCTION chixy(bang)
      REAL chixy,bang,BIG
      INTEGER NMAX
      PARAMETER (NMAX=1000,BIG=1.E30)
      INTEGER nn,j
      REAL xx(NMAX),yy(NMAX),sx(NMAX),sy(NMAX),ww(NMAX),aa,offs,avex,
     *avey,sumw,b
      COMMON /fitxyc/ xx,yy,sx,sy,ww,aa,offs,nn
      b=tan(bang)
      avex=0.
      avey=0.
      sumw=0.
      do 11 j=1,nn
        ww(j)=(b*sx(j))**2+sy(j)**2
        if(ww(j).lt.1./BIG) then
          ww(j)=BIG
        else
          ww(j)=1./ww(j)
        endif
        sumw=sumw+ww(j)
        avex=avex+ww(j)*xx(j)
        avey=avey+ww(j)*yy(j)
11    continue
      avex=avex/sumw
      avey=avey/sumw
      aa=avey-b*avex
      chixy=-offs
      do 12 j=1,nn
        chixy=chixy+ww(j)*(yy(j)-aa-b*xx(j))**2
12    continue
      return
      END

      SUBROUTINE covsrt(covar,npc,ma,ia,mfit)
      INTEGER ma,mfit,npc,ia(ma)
      REAL covar(npc,npc)
      INTEGER i,j,k
      REAL swap
      do 12 i=mfit+1,ma
        do 11 j=1,i
          covar(i,j)=0.
          covar(j,i)=0.
11      continue
12    continue
      k=mfit
      do 15 j=ma,1,-1
        if(ia(j).ne.0)then
          do 13 i=1,ma
            swap=covar(i,k)
            covar(i,k)=covar(i,j)
            covar(i,j)=swap
13        continue
          do 14 i=1,ma
            swap=covar(k,i)
            covar(k,i)=covar(j,i)
            covar(j,i)=swap
14        continue
          k=k-1
        endif
15    continue
      return
      END

      FUNCTION gammq(a,x)
      REAL a,gammq,x
CU    USES gcf,gser
      REAL gammcf,gamser,gln
      if(x.lt.0..or.a.le.0.)pause 'bad arguments in gammq'
      if(x.lt.a+1.)then
        call gser(gamser,a,x,gln)
        gammq=1.-gamser
      else
        call gcf(gammcf,a,x,gln)
        gammq=gammcf
      endif
      return
      END

      SUBROUTINE gcf(gammcf,a,x,gln)
      INTEGER ITMAX
      REAL a,gammcf,gln,x,EPS,FPMIN
      PARAMETER (ITMAX=100,EPS=3.e-7,FPMIN=1.e-30)
CU    USES gammln
      INTEGER i
      REAL an,b,c,d,del,h,gammln
      gln=gammln(a)
      b=x+1.-a
      c=1./FPMIN
      d=1./b
      h=d
      do 11 i=1,ITMAX
        an=-i*(i-a)
        b=b+2.
        d=an*d+b
        if(abs(d).lt.FPMIN)d=FPMIN
        c=b+an/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1./d
        del=d*c
        h=h*del
        if(abs(del-1.).lt.EPS)goto 1
11    continue
      pause 'a too large, ITMAX too small in gcf'
1     gammcf=exp(-x+a*log(x)-gln)*h
      return
      END

      SUBROUTINE gser(gamser,a,x,gln)
      INTEGER ITMAX
      REAL a,gamser,gln,x,EPS
      PARAMETER (ITMAX=1000,EPS=3.e-7)
CU    USES gammln
      INTEGER n
      REAL ap,del,sum,gammln
      gln=gammln(a)
      if(x.le.0.)then
        if(x.lt.0.)pause 'x < 0 in gser'
        gamser=0.
        return
      endif
      ap=a
      sum=1./a
      del=sum
      do 11 n=1,ITMAX
        ap=ap+1.
        del=del*x/ap
        sum=sum+del
        if(abs(del).lt.abs(sum)*EPS)goto 1
11    continue
      pause 'a too large, ITMAX too small in gser'
1     gamser=sum*exp(-x+a*log(x)-gln)
      return
      END

      FUNCTION gammln(xx)
      REAL gammln,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
      END

      SUBROUTINE gaussj(a,n,np,b,m,mp)
      INTEGER m,mp,n,np,NMAX
      REAL a(np,np),b(np,mp)
      PARAMETER (NMAX=50)
      INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
      REAL big,dum,pivinv
      do 11 j=1,n
        ipiv(j)=0
11    continue
      do 22 i=1,n
        big=0.
        do 13 j=1,n
          if(ipiv(j).ne.1)then
            do 12 k=1,n
              if (ipiv(k).eq.0) then
                if (abs(a(j,k)).ge.big)then
                  big=abs(a(j,k))
                  irow=j
                  icol=k
                endif
              endif
12          continue
          endif
13      continue
        ipiv(icol)=ipiv(icol)+1
        if (irow.ne.icol) then
          do 14 l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
14        continue
          do 15 l=1,m
            dum=b(irow,l)
            b(irow,l)=b(icol,l)
            b(icol,l)=dum
15        continue
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol).eq.0.) pause 'singular matrix in gaussj'
        pivinv=1./a(icol,icol)
        a(icol,icol)=1.
        do 16 l=1,n
          a(icol,l)=a(icol,l)*pivinv
16      continue
        do 17 l=1,m
          b(icol,l)=b(icol,l)*pivinv
17      continue
        do 21 ll=1,n
          if(ll.ne.icol)then
            dum=a(ll,icol)
            a(ll,icol)=0.
            do 18 l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
18          continue
            do 19 l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
19          continue
          endif
21      continue
22    continue
      do 24 l=n,1,-1
        if(indxr(l).ne.indxc(l))then
          do 23 k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
23        continue
        endif
24    continue
      return
      END

      SUBROUTINE mnbrak(ax,bx,cx,fa,fb,fc,func)
      REAL ax,bx,cx,fa,fb,fc,func,GOLD,GLIMIT,TINY
      EXTERNAL func
      PARAMETER (GOLD=1.618034, GLIMIT=100., TINY=1.e-20)
      REAL dum,fu,q,r,u,ulim
      fa=func(ax)
      fb=func(bx)
      if(fb.gt.fa)then
        dum=ax
        ax=bx
        bx=dum
        dum=fb
        fb=fa
        fa=dum
      endif
      cx=bx+GOLD*(bx-ax)
      fc=func(cx)
1     if(fb.ge.fc)then
        r=(bx-ax)*(fb-fc)
        q=(bx-cx)*(fb-fa)
        u=bx-((bx-cx)*q-(bx-ax)*r)/(2.*sign(max(abs(q-r),TINY),q-r))
        ulim=bx+GLIMIT*(cx-bx)
        if((bx-u)*(u-cx).gt.0.)then
          fu=func(u)
          if(fu.lt.fc)then
            ax=bx
            fa=fb
            bx=u
            fb=fu
            return
          else if(fu.gt.fb)then
            cx=u
            fc=fu
            return
          endif
          u=cx+GOLD*(cx-bx)
          fu=func(u)
        else if((cx-u)*(u-ulim).gt.0.)then
          fu=func(u)
          if(fu.lt.fc)then
            bx=cx
            cx=u
            u=cx+GOLD*(cx-bx)
            fb=fc
            fc=fu
            fu=func(u)
          endif
        else if((u-ulim)*(ulim-cx).ge.0.)then
          u=ulim
          fu=func(u)
        else
          u=cx+GOLD*(cx-bx)
          fu=func(u)
        endif
        ax=bx
        bx=cx
        cx=u
        fa=fb
        fb=fc
        fc=fu
        goto 1
      endif
      return
      END

      FUNCTION brent(ax,bx,cx,f,tol,xmin)
      INTEGER ITMAX
      REAL brent,ax,bx,cx,tol,xmin,f,CGOLD,ZEPS
      EXTERNAL f
      PARAMETER (ITMAX=100,CGOLD=.3819660,ZEPS=1.0e-10)
      INTEGER iter
      REAL a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
      a=min(ax,cx)
      b=max(ax,cx)
      v=bx
      w=v
      x=v
      e=0.
      fx=f(x)
      fv=fx
      fw=fx
      do 11 iter=1,ITMAX
        xm=0.5*(a+b)
        tol1=tol*abs(x)+ZEPS
        tol2=2.*tol1
        if(abs(x-xm).le.(tol2-.5*(b-a))) goto 3
        if(abs(e).gt.tol1) then
          r=(x-w)*(fx-fv)
          q=(x-v)*(fx-fw)
          p=(x-v)*q-(x-w)*r
          q=2.*(q-r)
          if(q.gt.0.) p=-p
          q=abs(q)
          etemp=e
          e=d
          if(abs(p).ge.abs(.5*q*etemp).or.p.le.q*(a-x).or.p.ge.q*(b-x)) 
     *goto 1
          d=p/q
          u=x+d
          if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
          goto 2
        endif
1       if(x.ge.xm) then
          e=a-x
        else
          e=b-x
        endif
        d=CGOLD*e
2       if(abs(d).ge.tol1) then
          u=x+d
        else
          u=x+sign(tol1,d)
        endif
        fu=f(u)
        if(fu.le.fx) then
          if(u.ge.x) then
            a=x
          else
            b=x
          endif
          v=w
          fv=fw
          w=x
          fw=fx
          x=u
          fx=fu
        else

          if(u.lt.x) then
            a=u
          else
            b=u
          endif
          if(fu.le.fw .or. w.eq.x) then
            v=w
            fv=fw
            w=u
            fw=fu
          else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
            v=u
            fv=fu
          endif
        endif
11    continue
      pause 'brent exceed maximum iterations'
3     xmin=x
      brent=fx
      return
      END

      FUNCTION zbrent(func,x1,x2,tol)
      INTEGER ITMAX
      REAL zbrent,tol,x1,x2,func,EPS
      EXTERNAL func
      PARAMETER (ITMAX=100,EPS=3.e-8)
      INTEGER iter
      REAL a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
      a=x1
      b=x2
      fa=func(a)
      fb=func(b)
      if((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.))pause
     *'root must be bracketed for zbrent'
      c=b
      fc=fb
      do 11 iter=1,ITMAX
        if((fb.gt.0..and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))then
          c=a
          fc=fa
          d=b-a
          e=d
        endif
        if(abs(fc).lt.abs(fb)) then
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
        endif
        tol1=2.*EPS*abs(b)+0.5*tol
        xm=.5*(c-b)
        if(abs(xm).le.tol1 .or. fb.eq.0.)then
          zbrent=b
          return
        endif
        if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
          s=fb/fa
          if(a.eq.c) then
            p=2.*xm*s
            q=1.-s
          else
            q=fa/fc
            r=fb/fc
            p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
            q=(q-1.)*(r-1.)*(s-1.)
          endif
          if(p.gt.0.) q=-q
          p=abs(p)
          if(2.*p .lt. min(3.*xm*q-abs(tol1*q),abs(e*q))) then
            e=d
            d=p/q
          else
            d=xm
            e=d
          endif
        else
          d=xm
          e=d
        endif
        a=b
        fa=fb
        if(abs(d) .gt. tol1) then
          b=b+d
        else
          b=b+sign(tol1,xm)
        endif
        fb=func(b)
11    continue
      pause 'zbrent exceeding maximum iterations'
      zbrent=b
      return
      END

      SUBROUTINE avevar(data,n,ave,var)
      INTEGER n
      REAL ave,var,data(n)
      INTEGER j
      REAL s,ep
      ave=0.0
      do 11 j=1,n
        ave=ave+data(j)
11    continue
      ave=ave/n
      var=0.0
      ep=0.0
      do 12 j=1,n
        s=data(j)-ave
        ep=ep+s
        var=var+s*s
12    continue
      var=(var-ep**2/n)/(n-1)
      return
      END
