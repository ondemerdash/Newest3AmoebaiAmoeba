!$OMP PARALLEL default(private) shared(demi,demk,fieldt,fieldtp,
!$OMP& viremrealt,emrealt,npole,electric,dielec,off2,toffset,
!$OMP& toffset0,dEindex,ntpair_dE,dEd1,dEd2,dEp1,dEp2,maxlocal,
!$OMP& nthread,maxelst,n,x,y,z,rpole,ipole,pdamp,thole,n12,i12,
!$OMP& m2scale,p2scale,n13,m3scale,p3scale,i13,n14,i14,p4scale,
!$OMP& m4scale,np11,ip11,p41scale,n15,m5scale,p5scale,i15,
!$OMP& np12,ip12,d1scale,d2scale,np13,ip13,d3scale,np14,ip14,d4scale,
!$OMP& nelst,elst,xaxis,yaxis,zaxis)
!$OMP& firstprivate(pscale,dscale,mscale,nlocal)

