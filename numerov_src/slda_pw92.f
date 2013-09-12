cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c correlation energy and its derivative w.r.t. rs and z at mu=infinity
c Perdew & Wang PRB 45, 13244 (1992)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ecPW(x,y,ec,ecd,ecz)
c in Hartree; ec=ec(rs,zeta)
c x -> rs; y -> zeta
ccc ecd is d/drs ec
ccc ecz is d/dz ec
      implicit none
      double precision pi,f02,ff,x,y,ec,ecd,ec0,ec0d,ec1,ec1d,
     $     aaa,G,Gd,alfac,alfacd,ecz
      pi=dacos(-1.d0)
      
      f02=4.d0/(9.d0*(2.d0**(1.d0/3.d0)-1.d0))

      ff=((1.d0+y)**(4.d0/3.d0)+(1.d0-y)**(4.d0/3.d0)-
     $     2.d0)/(2.d0**(4.d0/3.d0)-2.d0)

      aaa=(1.d0-log(2.d0))/pi**2
      call  GPW(x,aaa,0.21370d0,7.5957d0,3.5876d0,
     $     1.6382d0,0.49294d0,G,Gd)
      ec0=G
      ec0d=Gd

      aaa=aaa/2.d0
      call GPW(x,aaa,0.20548d0,14.1189d0,6.1977d0,
     $     3.3662d0,0.62517d0,G,Gd)
      ec1=G
      ec1d=Gd
      call GPW(x,0.016887d0,0.11125d0,10.357d0,3.6231d0,
     $     0.88026d0,0.49671d0,G,Gd)
      alfac=-G
      alfacd=-Gd

      ec=ec0+alfac*ff/f02*(1.d0-y**4)+(ec1-ec0)*ff*y**4
      ecd=ec0d+alfacd*ff/f02*(1.d0-y**4)+(ec1d-ec0d)*
     $     ff*y**4
      ecz=alfac*(-4.d0*y**3)*ff/f02+alfac*(1.d0-y**4)/f02*
     $     4.d0/3.d0*((1.d0+y)**(1.d0/3.d0)-(1.d0-y)**(1.d0/3.d0))/
     $     (2.d0**(4.d0/3.d0)-2.d0)+(ec1-ec0)*(4.d0*y**3*ff+
     $     4.d0/3.d0*((1.d0+y)**(1.d0/3.d0)-(1.d0-y)**(1.d0/3.d0))/
     $     (2.d0**(4.d0/3.d0)-2.d0)*y**4)

      return
      end

      subroutine GPW(x,Ac,alfa1,beta1,beta2,beta3,beta4,G,Gd)
ccc Gd is d/drs G
      implicit none
      double precision G,Gd,Ac,alfa1,beta1,beta2,beta3,beta4,x
      G=-2.d0*Ac*(1.d0+alfa1*x)*dlog(1.d0+1.d0/(2.d0*
     $     Ac*(beta1*x**0.5d0+
     $     beta2*x+beta3*x**1.5d0+beta4*x**2)))
      Gd=(1.d0+alfa1*x)*(beta2+beta1/(2.d0*sqrt(x))+3.d0*beta3*
     $     sqrt(x)/2.d0+2.d0*beta4*x)/((beta1*sqrt(x)+beta2*x+
     $     beta3*x**(3.d0/2.d0)+beta4*x**2)**2*(1.d0+1.d0/
     $     (2.d0*Ac*(beta1*sqrt(x)+beta2*x+beta3*x**(3.d0/2.d0)+
     $     beta4*x**2))))-2.d0*Ac*alfa1*dlog(1.d0+1.d0/(2.d0*Ac*
     $     (beta1*sqrt(x)+beta2*x+beta3*x**(3.d0/2.d0)+
     $     beta4*x**2)))
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
