c produce phase functions
        real ang,valeur(13,95),norm(43)
        real data(43,181),down,up,scat(43)
        integer i,j,k,wl,ii,jj,iid,iiu,wavscat(43)
        character*12 outfile
        character*3 wave
        open(unit=1,file='opacfog.txt',status='old')
          read(1,*) ((valeur(i,j),i=1,13),j=1,95)
        close(unit=1)
        open(unit=1,file='scatfog.txt',status='old')
          read(1,*)
          do i=1,43
            read(1,*) wavscat(i),scat(i)
          enddo
        close(unit=1)
        do jj=1,181
        do ii=1,43
           data(ii,jj)=-1.
        enddo
        enddo
        ang=-1.
        do j=2,95
           jj=nint(valeur(1,j))+1
           do i=2,13
              if ((valeur(i,1)*1000..ge.335.).and.
     +        (valeur(i,1)*1000..le.760.)) then
                 ii=nint((valeur(i,1)*1000.-340.)/10.)+1
                 data(ii,jj)=valeur(i,j)
              endif
           enddo
        enddo
c        print*,data
        do jj=1,181
           up=-1.
           down=-1.
           do k=1,11
              if ((data(1+k,jj).ge.0.).and.(down.lt.0.)) then
                 down=data(1+k,jj)
                 iid=1+k
              endif
              if ((data(1+k,jj).ge.0.).and.(down.ge.0.)) then
                 up=data(1+k,jj)
                 iiu=1+k
              endif
           enddo             
           data(1,jj)=down-real(iid-1)*(up-down)/(real(iiu-iid))
           up=-1.
           down=-1.
           do k=1,11
              if ((data(43-k,jj).ge.0.).and.(up.lt.0.)) then
                 up=data(43-k,jj)
                 iiu=43-k
              endif
              if ((data(43-k,jj).ge.0.).and.(up.ge.0.)) then
                 down=data(43-k,jj)
                 iid=43-k
              endif
           enddo             
           data(43,jj)=up+real(43-iiu)*(up-down)/(real(iiu-iid))


           
           do ii=2,42
              if (data(ii,jj).lt.0.) then
                 down=-1.
                 up=-1.
                 do k=1,4
                    if (data(ii+k,jj).ge.0.) then
                       up=data(ii+k,jj)
                       iiu=ii+k
                    endif
                    if (data(ii-k,jj).ge.0.) then
                       down=data(ii-k,jj)
                       iid=ii-k
                    endif
                 enddo
                 if ((down.lt.0.).and.(up.lt.0.)) then
                    data(ii,jj)=-1.
                 else
                  data(ii,jj)=down+(up-down)*real(ii-iid)/real(iiu-iid)
                 endif
              endif
           enddo     
        enddo       
c         print*,data
c         stop
  
            up=-1.
           down=-1.       
        do ii=1,43
           do jj=2,180
              if (data(ii,jj).lt.0.) then
                 data(ii,jj)=(data(ii,jj+1)+data(ii,jj-1))/2.
              endif
           enddo
        enddo
c normalisation de la fonction de phase integrale = 4pi
        pi=3.141592654
        do ii=1,43
           norm(ii)=0.
           do jj=1,181
              ang=real(jj-1)
              norm(ii)=norm(ii)+
     +        data(ii,jj)*2.*pi*sin(ang*pi/180.)*pi/180.
           enddo
           do jj=1,181
              ang=real(jj-1)
              data(ii,jj)=4.*pi*data(ii,jj)/norm(ii)
           enddo          
        enddo        
        wl=330       
        do ii=1,43
           wl=wl+10        
           write(wave, '(I3.3)' ) wl
           outfile='fog_'//wave//'.txt'
           open(unit=2,file=outfile,status='unknown')
           write(2,*) scat(ii), '   #    scattering/extinction'
           write(2,*) '  ScatAngle   PhaseFct'  
           do jj=1,181
              ang=real(jj-1)
              write(2,*) int(ang),data(ii,jj)
           enddo
           close(unit=2) 
        enddo
        stop
        end
