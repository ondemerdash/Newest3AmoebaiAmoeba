      tauxx1p=tauxx1p+ (0.5d0*rr5*psc5*dixr(1)*xr 
     &        +0.5d0*2.0d0*rr5*psc5*(yr*qi(3)-zr*qi(2)) 
     &         - rr7*psc7*rxqir(1)*xr)*(-1.0d0)
      tauxy1p=tauxy1p+(-0.5d0*rr3*psc3*(-di(3))
     &             +0.5d0*rr5*psc5*dixr(1)*yr
     &          + 0.5d0*2.0d0*rr5*psc5*(qir(3) + yr*qi(6)-zr*qi(5))
     &          -rr7*psc7*rxqir(1)*yr)*(-1.0d0)
      tauxz1p=tauxz1p+ (-0.5d0*rr3*psc3*di(2) 
     &           + 0.5d0*rr5*psc5*dixr(1)*zr
     &          + 0.5d0*2.0d0*rr5*psc5*(-qir(2) + yr*qi(9)-zr*qi(8))
     &          -rr7*psc7*rxqir(1)*zr)*(-1.0d0)


      tauyx1p= tauyx1p+ (-0.5d0*rr3*psc3*di(3) 
     &             + 0.5d0*rr5*psc5*dixr(2)*xr
     &          +0.5d0*2.0d0*rr5*psc5*( -qir(3) + zr*qi(1) - xr*qi(3))
     &          -rr7*psc7*rxqir(2)*xr)*(-1.0d0)
      tauyy1p= tauyy1p+ (0.5d0*rr5*psc5*dixr(2)*yr
     &         +0.5d0*2.0d0*rr5*psc5*( zr*qi(4)-xr*qi(6))
     &         -rr7*psc7*rxqir(2)*yr)*(-1.0d0)
      tauyz1p=tauyz1p+ (-0.5d0*rr3*psc3*(-di(1))
     &               +0.5d0*rr5*psc5*dixr(2)*zr
     &         +0.5d0*2.0d0*rr5*psc5*( qir(1) + zr*qi(7)-xr*qi(9) )
     &         -rr7*psc7*rxqir(2)*zr)*(-1.0d0)

      tauzx1p=tauzx1p+ (-0.5d0*rr3*psc3*(-di(2))
     &         + 0.5d0*rr5*psc5*dixr(3)*xr
     &         +0.5d0*2.0d0*rr5*psc5*(qir(2) + xr*qi(2)-yr*qi(1) )
     &         -rr7*psc7*rxqir(3)*xr)*(-1.0d0)
      tauzy1p=tauzy1p+ (-0.5d0*rr3*psc3*di(1) 
     &         + 0.5d0*rr5*psc5*dixr(3)*yr
     &        +0.5d0*2.0d0*rr5*psc5*(-qir(1) + xr*qi(5)-yr*qi(4) )
     &        -rr7*psc7*rxqir(3)*yr)*(-1.0d0)
      tauzz1p=tauzz1p+ (0.5d0*rr5*psc5*dixr(3)*zr
     &        +0.5d0*2.0d0*rr5*psc5*(xr*qi(8)-yr*qi(7))
     &        -rr7*psc7*rxqir(3)*zr)*(-1.0d0)


      tauxx2p=tauxx2p + (0.5d0*rr5*psc5*dkxr(1)*xr
     &        -0.5d0*2.0d0*rr5*psc5*(yr*qk(3)-zr*qk(2))
     &        -(-rr7*psc7)*rxqkr(1)*xr)*(-1.0d0)
      tauxy2p=tauxy2p+ (-0.5d0*rr3*psc3*(-dk(3))
     &           +0.5d0*rr5*psc5*dkxr(1)*yr
     &       -0.5d0*2.0d0*rr5*psc5*(qkr(3) + yr*qk(6)-zr*qk(5))
     &       -(-rr7*psc7)*rxqkr(1)*yr)*(-1.0d0)
      tauxz2p=tauxz2p+ (-0.5d0*rr3*psc3*dk(2) 
     &          + 0.5d0*rr5*psc5*dkxr(1)*zr
     &        -0.5d0*2.0d0*rr5*psc5*(-qkr(2) +yr*qk(9)-zr*qk(8))
     &        -(-rr7*psc7)*rxqkr(1)*zr)*(-1.0d0)


      tauyx2p= tauyx2p+ (-0.5d0*rr3*psc3*dk(3) 
     &                 + 0.5d0*rr5*psc5*dkxr(2)*xr
     &          -0.5d0*2.0d0*rr5*psc5*( -qkr(3) + zr*qk(1) - xr*qk(3))
     &          -(-rr7*psc7)*rxqkr(2)*xr)*(-1.0d0)
      tauyy2p= tauyy2p+ (0.5d0*rr5*psc5*dkxr(2)*yr
     &         -0.5d0*2.0d0*rr5*psc5*( zr*qk(4)-xr*qk(6))
     &         -(-rr7*psc7)*rxqkr(2)*yr)*(-1.0d0)
      tauyz2p=tauyz2p+ (-0.5d0*rr3*psc3*(-dk(1))
     &                  +0.5d0*rr5*psc5*dkxr(2)*zr
     &         -0.5d0*2.0d0*rr5*psc5*( qkr(1) + zr*qk(7)-xr*qk(9) )
     &         -(-rr7*psc7)*rxqkr(2)*zr)*(-1.0d0)

      tauzx2p=tauzx2p+ (-0.5d0*rr3*psc3*(-dk(2))
     &                     +0.5d0*rr5*psc5*dkxr(3)*xr
     &         -0.5d0*2.0d0*rr5*psc5*(qkr(2) + xr*qk(2)-yr*qk(1) )
     &         -(-rr7*psc7)*rxqkr(3)*xr)*(-1.0d0)
      tauzy2p=tauzy2p+ (-0.5d0*rr3*psc3*dk(1) 
     &                     + 0.5d0*rr5*psc5*dkxr(3)*yr
     &        -0.5d0*2.0d0*rr5*psc5*(-qkr(1) + xr*qk(5)-yr*qk(4) )
     &        -(-rr7*psc7)*rxqkr(3)*yr)*(-1.0d0)
      tauzz2p=tauzz2p+(0.5d0*rr5*psc5*dkxr(3)*zr
     &        -0.5d0*2.0d0*rr5*psc5*(xr*qk(8)-yr*qk(7))
     &        -(-rr7*psc7)*rxqkr(3)*zr)*(-1.0d0)

