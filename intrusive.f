c intrusive light routine for a single lamp
c
c 
       subroutine intrusive(intlu,d_o,h_o,hwindow,h_l,pvalno,srei,stype
     + ,Irad)
       real pvalno(181,120),d_o,h_o,hwindow,h_l,Irad,z,pi,srei            ! srei = suface reflectance
       real inteo,integ,intlu
       integer stype,i,iw
       pi=3.14159
       z_o=pi/2.-atan((h_o-h_l)/d_o)
       z_g=atan(d_o/h_l)
       z_w=pi/2.+atan((h_l-h_w)/d_o)

c integrate the LOP from obstacle base to obstacle top (inteo)
c and LOP from nadir to obstacle base (integ)
       inteo=0.
       integ=0.
       dz=pi/180.
       for i=1,181
         z=real(i-1)*pi/180.
         if ((z.ge.z_o).and.(z.lt.z_g)) then
            inteo=inteo+pvalno(i,ntype)*sin(z)*dz
         endif
         if ((z.ge.z_g).and.(z.lt.pi)) then
            integ=integ+pvalno(i,ntype)*sin(z)*dz
         endif
         if (z.eq.z_w) iw=i
       enddo
       Irad=intlu*(pvalno(iw,ntype)+2.*srei*integ+srei/4.*inteo
       return
       end
         

