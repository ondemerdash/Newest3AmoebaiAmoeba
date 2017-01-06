      tauxx1d=tauxx1d+ (0.5d0*rr5*dsc5*dixr(1)*xr 
     &        +0.5d0*2.0d0*rr5*dsc5*(yr*qi(3)-zr*qi(2)) 
     &         - rr7*dsc7*rxqir(1)*xr)*(-1.0d0)
      tauxy1d=tauxy1d+(-0.5d0*rr3*dsc3*(-di(3))
     &             +0.5d0*rr5*dsc5*dixr(1)*yr
     &          + 0.5d0*2.0d0*rr5*dsc5*(qir(3) + yr*qi(6)-zr*qi(5))
     &          -rr7*dsc7*rxqir(1)*yr)*(-1.0d0)
      tauxz1d=tauxz1d+ (-0.5d0*rr3*dsc3*di(2) 
     &           + 0.5d0*rr5*dsc5*dixr(1)*zr
     &          + 0.5d0*2.0d0*rr5*dsc5*(-qir(2) + yr*qi(9)-zr*qi(8))
     &          -rr7*dsc7*rxqir(1)*zr)*(-1.0d0)


      tauyx1d= tauyx1d+ (-0.5d0*rr3*dsc3*di(3) 
     &             + 0.5d0*rr5*dsc5*dixr(2)*xr
     &          +0.5d0*2.0d0*rr5*dsc5*( -qir(3) + zr*qi(1) - xr*qi(3))
     &          -rr7*dsc7*rxqir(2)*xr)*(-1.0d0)
      tauyy1d= tauyy1d+ (0.5d0*rr5*dsc5*dixr(2)*yr
     &         +0.5d0*2.0d0*rr5*dsc5*( zr*qi(4)-xr*qi(6))
     &         -rr7*dsc7*rxqir(2)*yr)*(-1.0d0)
      tauyz1d=tauyz1d+ (-0.5d0*rr3*dsc3*(-di(1))
     &               +0.5d0*rr5*dsc5*dixr(2)*zr
     &         +0.5d0*2.0d0*rr5*dsc5*( qir(1) + zr*qi(7)-xr*qi(9) )
     &         -rr7*dsc7*rxqir(2)*zr)*(-1.0d0)

      tauzx1d=tauzx1d+ (-0.5d0*rr3*dsc3*(-di(2))
     &         + 0.5d0*rr5*dsc5*dixr(3)*xr
     &         +0.5d0*2.0d0*rr5*dsc5*(qir(2) + xr*qi(2)-yr*qi(1) )
     &         -rr7*dsc7*rxqir(3)*xr)*(-1.0d0)
      tauzy1d=tauzy1d+ (-0.5d0*rr3*dsc3*di(1) 
     &         + 0.5d0*rr5*dsc5*dixr(3)*yr
     &        +0.5d0*2.0d0*rr5*dsc5*(-qir(1) + xr*qi(5)-yr*qi(4) )
     &        -rr7*dsc7*rxqir(3)*yr)*(-1.0d0)
      tauzz1d=tauzz1d+ (0.5d0*rr5*dsc5*dixr(3)*zr
     &        +0.5d0*2.0d0*rr5*dsc5*(xr*qi(8)-yr*qi(7))
     &        -rr7*dsc7*rxqir(3)*zr)*(-1.0d0)


      tauxx2d=tauxx2d + (0.5d0*rr5*dsc5*dkxr(1)*xr
     &        -0.5d0*2.0d0*rr5*dsc5*(yr*qk(3)-zr*qk(2))
     &        -(-rr7*dsc7)*rxqkr(1)*xr)*(-1.0d0)
      tauxy2d=tauxy2d+ (-0.5d0*rr3*dsc3*(-dk(3))
     &           +0.5d0*rr5*dsc5*dkxr(1)*yr
     &       -0.5d0*2.0d0*rr5*dsc5*(qkr(3) + yr*qk(6)-zr*qk(5))
     &       -(-rr7*dsc7)*rxqkr(1)*yr)*(-1.0d0)
      tauxz2d=tauxz2d+ (-0.5d0*rr3*dsc3*dk(2) 
     &          + 0.5d0*rr5*dsc5*dkxr(1)*zr
     &        -0.5d0*2.0d0*rr5*dsc5*(-qkr(2) +yr*qk(9)-zr*qk(8))
     &        -(-rr7*dsc7)*rxqkr(1)*zr)*(-1.0d0)


      tauyx2d= tauyx2d+ (-0.5d0*rr3*dsc3*dk(3) 
     &                 + 0.5d0*rr5*dsc5*dkxr(2)*xr
     &          -0.5d0*2.0d0*rr5*dsc5*( -qkr(3) + zr*qk(1) - xr*qk(3))
     &          -(-rr7*dsc7)*rxqkr(2)*xr)*(-1.0d0)
      tauyy2d= tauyy2d+ (0.5d0*rr5*dsc5*dkxr(2)*yr
     &         -0.5d0*2.0d0*rr5*dsc5*( zr*qk(4)-xr*qk(6))
     &         -(-rr7*dsc7)*rxqkr(2)*yr)*(-1.0d0)
      tauyz2d=tauyz2d+ (-0.5d0*rr3*dsc3*(-dk(1))
     &                  +0.5d0*rr5*dsc5*dkxr(2)*zr
     &         -0.5d0*2.0d0*rr5*dsc5*( qkr(1) + zr*qk(7)-xr*qk(9) )
     &         -(-rr7*dsc7)*rxqkr(2)*zr)*(-1.0d0)

      tauzx2d=tauzx2d+ (-0.5d0*rr3*dsc3*(-dk(2))
     &                     +0.5d0*rr5*dsc5*dkxr(3)*xr
     &         -0.5d0*2.0d0*rr5*dsc5*(qkr(2) + xr*qk(2)-yr*qk(1) )
     &         -(-rr7*dsc7)*rxqkr(3)*xr)*(-1.0d0)
      tauzy2d=tauzy2d+ (-0.5d0*rr3*dsc3*dk(1) 
     &                     + 0.5d0*rr5*dsc5*dkxr(3)*yr
     &        -0.5d0*2.0d0*rr5*dsc5*(-qkr(1) + xr*qk(5)-yr*qk(4) )
     &        -(-rr7*dsc7)*rxqkr(3)*yr)*(-1.0d0)
      tauzz2d=tauzz2d+(0.5d0*rr5*dsc5*dkxr(3)*zr
     &        -0.5d0*2.0d0*rr5*dsc5*(xr*qk(8)-yr*qk(7))
     &        -(-rr7*dsc7)*rxqkr(3)*zr)*(-1.0d0)

