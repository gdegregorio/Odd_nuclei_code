       Program density
 
       use technical
       use geom

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc' 
     
  
      call input_den

      eps =1.d9 
     
      write(*,*) 'calling densities'

      call rohp(-1,jmax,levp,eps)
      call rohn(1,jmax,levn,eps)
      call ropp(-1,jmax,levp,eps)
      call ropn(1,jmax,levn,eps)
!************************************************
      deallocate(levn,levp)
      deallocate(dime,ph)
!************************************************      
      end
!subroutine functions
!************************************************
!
      subroutine ropp(ipart,isimax,lev,en_cut)
!     density for proton particle       
!
!*************************************************
      use technical
      use geom
 
      include 'commons.inc'
       
      integer :: lam1,lam2               !index for phonons do cycle
      integer :: p1,p2                   !indeces for particles do cycles
      integer :: j_sig                   !index for coupling momentum 
      integer :: ih                      !index for hole do cycle  
      integer :: Jp1,Jp2,ja,jb,jh        !angular mumenta 
      integer :: ilamp,ilam              !indeces of phonons
      integer :: ipmin,ipmax,ihmin,ihmax !space dimension indices
      integer :: irho                    !demension of density array
      integer :: ipp1,ipp2,inn1,inn2     

      double precision :: f1,rac           !service variable
      double precision :: rp

        character   :: ch1
        character*2 :: ch2
        character*5 :: ch3

      type(level_type),dimension(*)   :: lev
     
          ndrho=1000000000    
 
          en_cut=1.d9
          isimax=jmax

          ipmin=ippmin 
          ipmax=ippmax
          ihmin=ihpmin    
          ihmax=ihpmax

      allocate(den(ndrho))

      do lam1=1,num_ph
       ip1=ph(lam1)%par
       Jp1=ph(lam1)%Jph
       ilamp=ph(lam1)%lam
       if(allocated(C1p)) deallocate(C1p)
       if(allocated(C1n)) deallocate(C1n)
       call cph(ilamp,ip1,Jp1,ipp1,inn1,imax)
       allocate(C1p(ipp1),C1n(inn1))
       C1p=cp
       C1n=cn
        irho=0
       deallocate(cp,cn)

        if(ph(lam1)%el.lt.en_cut) then

      write(*,*) 'Generating dens. for J=',Jp1,'pi=',ip1,'lam=',ilamp

       if(ip1.eq.-1) ch1='0'
       if(ip1.eq.1)  ch1='1'

612    format(i2.2) 

       write(ch2,612) Jp1

613    format(i5.5)

       write(ch3,613) ilamp


      write(*,*)'Calculation of 1-phonon particle proton densities'

      open(1,file='dens/Prot_den_'//ch1//'_'//ch2//'_'//ch3//'.out',
     &                     status='unknown',form='unformatted')


      do lam2=1,num_ph

       ip2=ph(lam2)%par
       Jp2=ph(lam2)%Jph
       ilam=ph(lam2)%lam
       if(allocated(C2p)) deallocate(C2p)
       if(allocated(C2n)) deallocate(C2n)
       call cph(ilam,ip2,Jp2,ipp2,inn2,imax)
       Allocate(C2p(ipp2),C2n(inn2)) 
       C2p=cp
       C2n=cn
       deallocate(cp,cn)

         do p1=ipmin,ipmax
          ja=lev(p1)%j2
           do p2=ipmin,ipmax
            jb=lev(p2)%j2

           do j_sig=0,jmax
           rp=0.d0

            do ii=1,ipp2
              do jj=1,ipp1
          if(C1p(jj)%par.eq.p1.and.C2p(ii)%par.eq.p2)then
           
          do ih=ihmin,ihmax
                      
          if(C1p(jj)%hol.eq.ih.and.C2p(ii)%hol.eq.ih)then
            jh=lev(ih)%j2
            f1=(dble(2*Jp1+1)*dble(2*Jp2+1)*dble(2*j_sig+1))**0.5d0
            rac=racah(dble(Jp2),dble(j_sig),0.5d0*dble(jh),0.5d0*(ja),
     &                 dble(Jp1),0.5d0*dble(jb))
            rp=rp+rac*C1p(jj)%cph*C2p(ii)%cph*f1 
           end if !ih
          end do !ih
          end if!p1,p2=par1,par2
         end do !ii
         end do

         if(dabs(rp).gt.1.d-10)then
         irho=irho+1
         if(irho.gt.n_max_dim) then
         write(*,*) 'Dimension of den_p overflew n_max_dim !!!'
            stop
           endif
           den(irho)%lam=lam2
           den(irho)%aa=p1
           den(irho)%bb=p2
           den(irho)%sig=j_sig
           den(irho)%val=rp 
          endif

        end do !p2
       end do !p1
      end do !jsig
      end do !lam2
      end if!en
       write(*,*)'dimension of density',irho
       write(1)irho

       if(irho.gt.0) then 

       write(1) (den(k)%lam,k=1,irho)
       write(1) (den(k)%aa,k=1,irho)
       write(1) (den(k)%bb,k=1,irho)
       write(1) (den(k)%sig,k=1,irho)
       write(1) (den(k)%val,k=1,irho)
       end if
       close(1)
      end do !lam1

      deallocate(den)

      return
      end subroutine
!************************************************
      subroutine ropn(ipart,isimax,lev,en_cut)
!     density for neutron particle       
!*************************************************
      use technical
      use geom
 
      include 'commons.inc'
       
      integer :: lam1,lam2               !index for phonons do cycle
      integer :: p1,p2                   !indeces for particles do cycles
      integer :: j_sig                   !index for coupling momentum 
      integer :: ih                      !index for hole do cycle  
      integer :: Jp1,Jp2,ja,jb,jh        !angular mumenta 
      integer :: ilamp,ilam              !indeces of phonons
      integer :: ipmin,ipmax,ihmin,ihmax !space dimension indices
      integer :: irho                    !demension of density array
      integer :: ipp1,ipp2,inn1,inn2     

      double precision :: f1,rac           !service variable
      double precision :: rp
        character   :: ch1
        character*2 :: ch2
        character*5 :: ch3



      type(level_type),dimension(*)   :: lev
     
          ndrho=1000000000
 
          en_cut=1.d9
          isimax=jmax

          ipmin=ipnmin 
          ipmax=ipnmax
          ihmin=ihnmin    
          ihmax=ihnmax

      allocate(den(ndrho))

      do lam1=1,num_ph
       ip1=ph(lam1)%par
       Jp1=ph(lam1)%Jph
       ilamp=ph(lam1)%lam
       if(allocated(C1p)) deallocate(C1p)
       if(allocated(C1n)) deallocate(C1n)
       call cph(ilamp,ip1,Jp1,ipp1,inn1,imax)
       allocate(C1p(ipp1),C1n(inn1))
       C1p=cp
       C1n=cn
        irho=0
       deallocate(cp,cn)

       if(ph(lam1)%el.lt.en_cut) then

      write(*,*) 'Generating dens. for J=',Jp1,'pi=',ip1,'lam=',ilamp

       if(ip1.eq.-1) ch1='0'
       if(ip1.eq.1)  ch1='1'

612    format(i2.2) 

       write(ch2,612) Jp1

613    format(i5.5)

       write(ch3,613) ilamp


      write(*,*)'Calculation of 1-phonon particle neutron densities'

       open(1,file='dens/Neut_den_'//ch1//'_'//ch2//'_'//ch3//'.out',
     &                     status='unknown',form='unformatted')


      do lam2=1,num_ph
       ip2=ph(lam2)%par
       Jp2=ph(lam2)%Jph
       ilam=ph(lam2)%lam
       if(allocated(C2p)) deallocate(C2p)
       if(allocated(C2n)) deallocate(C2n)
       call cph(ilam,ip2,Jp2,ipp2,inn2,imax)
       Allocate(C2p(ipp2),C2n(inn2)) 
       C2p=cp
       C2n=cn
       deallocate(cp,cn)
         
         do p1=ipmin,ipmax
          ja=lev(p1)%j2
            do p2=ipmin,ipmax
           jb=lev(p2)%j2
         do j_sig=0,jmax!abs(Jp1-Jp2),(Jp1+Jp2)

           rp=0.d0
            do ii=1,inn2
              do jj=1,inn1
          if(C1n(jj)%par.eq.p1.and.C2n(ii)%par.eq.p2)then
           
          do ih=ihmin,ihmax
                      
          if(C1n(jj)%hol.eq.ih.and.C2n(ii)%hol.eq.ih)then
            jh=lev(ih)%j2
            f1=(dble(2*Jp1+1)*dble(2*Jp2+1)*dble(2*j_sig+1))**0.5d0
            rac=racah(dble(Jp2),dble(j_sig),0.5d0*dble(jh),0.5d0*(ja),
     &                 dble(Jp1),0.5d0*dble(jb))
            rp=rp+rac*C1n(jj)%cph*C2n(ii)%cph*f1 
           end if !ih
          end do !ih
          end if!p1,p2=par1,par2
         end do !ii
         end do
         if(dabs(rp).gt.1.d-10)then
         irho=irho+1
         if(irho.gt.n_max_dim) then
         write(*,*) 'Dimension of den_p overflew n_max_dim !!!'
            stop
           endif
           den(irho)%lam=lam2
           den(irho)%aa=p1
           den(irho)%bb=p2
           den(irho)%sig=j_sig
           den(irho)%val=rp
          endif
        end do !p2
       end do !p1
      end do !jsig
      end do !lam2
      end if!en
       write(*,*)'dimension of density',irho
       write(1)irho
       if(irho.gt.0) then 
       !write(1) irho
       write(1) (den(k)%lam,k=1,irho)
       write(1) (den(k)%aa,k=1,irho)
       write(1) (den(k)%bb,k=1,irho)
       write(1) (den(k)%sig,k=1,irho)
       write(1) (den(k)%val,k=1,irho)
       end if
       close(1)
      
      end do !lam1
      deallocate(den)
      return
      end subroutine

!************************************************
      subroutine rohp(ipart,isimax,lev,en_cut)
!     density for proton hole       
!*************************************************
      use technical
      use geom
 
      include 'commons.inc'
       
      integer :: lam1,lam2               !index for phonons do cycle
      integer :: h1,h2                   !indeces for particles do cycles
      integer :: j_sig                   !index for coupling momentum 
      integer :: ip                      !index for hole do cycle  
      integer :: Jp1,Jp2,ja,jb,jp        !angular momenta 
      integer :: ilamp,ilam              !indeces of phonons
      integer :: ipmin,ipmax,ihmin,ihmax !space dimension indices
      integer :: irho                    !demension of density array
      integer :: ipp1,ipp2,inn1,inn2     

      double precision :: f1,rac           !service variable
      double precision :: rp
      double precision :: phase

        character   :: ch1
        character*2 :: ch2
        character*5 :: ch3

      type(level_type),dimension(*)   :: lev
     
          ndrho=1000000000
 
          en_cut=1.d9
          isimax=jmax

          ipmin=ippmin 
          ipmax=ippmax
          ihmin=ihpmin    
          ihmax=ihpmax

      allocate(den(ndrho))


      do lam1=1,num_ph
       ip1=ph(lam1)%par
       Jp1=ph(lam1)%Jph
       ilamp=ph(lam1)%lam
       if(allocated(C1p)) deallocate(C1p)
       if(allocated(C1n)) deallocate(C1n)
       call cph(ilamp,ip1,Jp1,ipp1,inn1,imax)
       allocate(C1p(ipp1),C1n(inn1))
       C1p=cp
       C1n=cn
       deallocate(cp,cn)
        irho=0
        if(ph(lam1)%el.lt.en_cut) then

      write(*,*) 'Generating dens. for J=',Jp1,'pi=',ip1,'lam=',ilamp

       if(ip1.eq.-1) ch1='0'
       if(ip1.eq.1)  ch1='1'

612    format(i2.2) 

       write(ch2,612) Jp1

613    format(i5.5)

       write(ch3,613) ilamp


      write(*,*)'Calculation of 1-phonon hole proton densities'

      open(1,file='dens/Prot_hol_den_'//ch1//'_'//ch2//'_'//ch3//'.out',
     &                     status='unknown',form='unformatted')


      do lam2=1,num_ph
       ip2=ph(lam2)%par
       Jp2=ph(lam2)%Jph
       ilam=ph(lam2)%lam
       if(allocated(C2p)) deallocate(C2p)
       if(allocated(C2n)) deallocate(C2n)
       call cph(ilam,ip2,Jp2,ipp2,inn2,imax)
       Allocate(C2p(ipp2),C2n(inn2)) 
       C2p=cp
       C2n=cn      
       deallocate(cp,cn)
         
         do h1=ihmin,ihmax
          ja=lev(h1)%j2
            do h2=ihmin,ihmax
           jb=lev(h2)%j2
         do j_sig=0,jmax!abs(Jp1-Jp2),(Jp1+Jp2)

           rp=0.d0
            do ii=1,ipp2
              do jj=1,ipp1
          if(C1p(jj)%hol.eq.h2.and.C2p(ii)%hol.eq.h1)then
                  
          do ip=ipmin,ipmax
                      
          if(C1p(jj)%par.eq.ip.and.C2p(ii)%par.eq.ip)then
            jp=lev(ip)%j2
            phase=dble((-1)**((ja-jb)/2-j_sig))
            f1=(dble(2*Jp1+1)*dble(2*Jp2+1)*dble(2*j_sig+1))**0.5d0
            rac=racah(dble(Jp1),dble(j_sig),0.5d0*dble(jp),0.5d0*(ja),
     &                 dble(Jp2),0.5d0*dble(jb))
            rp=rp-phase*rac*C1p(jj)%cph*C2p(ii)%cph*f1 

           end if !ip
          end do !ip
          end if!p1,p2=par1,par2
         end do !ii
         end do
         if(dabs(rp).gt.1.d-10)then
         irho=irho+1
         if(irho.gt.n_max_dim) then
         write(*,*) 'Dimension of den_p overflew n_max_dim !!!'
            stop
           endif
           den(irho)%lam=lam2
           den(irho)%aa=h1
           den(irho)%bb=h2
           den(irho)%sig=j_sig
           den(irho)%val=rp
          endif
        end do !p2
       end do !p1
      end do !jsig
      end do !lam2
      end if!en
       write(*,*)'dimension of density',irho
       write(1)irho
       if(irho.gt.0) then 

       write(1) (den(k)%lam,k=1,irho)
       write(1) (den(k)%aa,k=1,irho)
       write(1) (den(k)%bb,k=1,irho)
       write(1) (den(k)%sig,k=1,irho)
       write(1) (den(k)%val,k=1,irho)
       end if
       close(1)
      
      end do !lam1
      deallocate(den)
      return
      end subroutine
!************************************************
      subroutine rohn(ipart,isimax,lev,en_cut)
!     density for neutron hole      
!*************************************************
      use technical
      use geom
 
      include 'commons.inc'
       
      integer :: lam1,lam2               !index for phonons do cycle
      integer :: h1,h2                   !indeces for particles do cycles
      integer :: j_sig                   !index for coupling momentum 
      integer :: ip                      !index for hole do cycle  
      integer :: Jp1,Jp2,ja,jb,jp        !angular mumenta 
      integer :: ilamp,ilam              !indeces of phonons
      integer :: ipmin,ipmax,ihmin,ihmax !space dimension indices
      integer :: irho                    !demension of density array
      integer :: ipp1,ipp2,inn1,inn2     

      double precision :: f1,rac           !service variable
      double precision :: rp
      double precision :: phase
        character   :: ch1
        character*2 :: ch2
        character*5 :: ch3

      type(level_type),dimension(*)   :: lev
     
          ndrho=1000000000
 
          en_cut=1.d9
          isimax=jmax

          ipmin=ipnmin 
          ipmax=ipnmax
          ihmin=ihnmin    
          ihmax=ihnmax

      allocate(den(ndrho))


      do lam1=1,num_ph
       ip1=ph(lam1)%par
       Jp1=ph(lam1)%Jph
       ilamp=ph(lam1)%lam
       if(allocated(C1p)) deallocate(C1p)
       if(allocated(C1n)) deallocate(C1n)
       call cph(ilamp,ip1,Jp1,ipp1,inn1,imax)
       allocate(C1p(ipp1),C1n(inn1))

       C1p=cp
       C1n=cn

       deallocate(cp,cn)

       irho=0

        if(ph(lam1)%el.lt.en_cut) then

      write(*,*) 'Generating dens. for J=',Jp1,'pi=',ip1,'lam=',ilamp

       if(ip1.eq.-1) ch1='0'
       if(ip1.eq.1)  ch1='1'

612    format(i2.2) 

       write(ch2,612) Jp1

613    format(i5.5)

       write(ch3,613) ilamp


      write(*,*)'Calculation of 1-phonon neutron hole densities' 
      open(1,file='dens/Neut_hol_den_'//ch1//'_'//ch2//'_'//ch3//'.out',
     &                     status='unknown',form='unformatted')


      do lam2=1,num_ph
       ip2=ph(lam2)%par
       Jp2=ph(lam2)%Jph
       ilam=ph(lam2)%lam
       if(allocated(C2p)) deallocate(C2p)
       if(allocated(C2n)) deallocate(C2n)
       call cph(ilam,ip2,Jp2,ipp2,inn2,imax)
       Allocate(C2p(ipp2),C2n(inn2)) 
       C2p=cp
       C2n=cn      
       deallocate(cp,cn)

         
         do h1=ihmin,ihmax
          ja=lev(h1)%j2
            do h2=ihmin,ihmax
           jb=lev(h2)%j2
         do j_sig=0,jmax!abs(Jp1-Jp2),(Jp1+Jp2)

           rp=0.d0
            do ii=1,inn2
              do jj=1,inn1
          if(C1n(jj)%hol.eq.h2.and.C2n(ii)%hol.eq.h1)then
           
          do ip=ipmin,ipmax
                      
          if(C1n(jj)%par.eq.ip.and.C2n(ii)%par.eq.ip)then
            jp=lev(ip)%j2
            phase=dble((-1)**((ja-jb)/2-j_sig))
            f1=(dble(2*Jp1+1)*dble(2*Jp2+1)*dble(2*j_sig+1))**0.5d0
            rac=racah(dble(Jp1),dble(j_sig),0.5d0*dble(jp),0.5d0*(ja),
     &                 dble(Jp2),0.5d0*dble(jb))
            rp=rp-phase*rac*C1n(jj)%cph*C2n(ii)%cph*f1 
           end if !ip
          end do !ip
          end if!p1,p2=par1,par2
         end do !ii
         end do
         if(dabs(rp).gt.1.d-10)then
         irho=irho+1
         if(irho.gt.n_max_dim) then
         write(*,*) 'Dimension of den_p overflew n_max_dim !!!'
            stop
           endif
           den(irho)%lam=lam2
           den(irho)%aa=h1
           den(irho)%bb=h2
           den(irho)%sig=j_sig
           den(irho)%val=rp

         endif

        end do !p2
       end do !p1
      end do !jsig
      end do !lam2
      end if!en

       write(*,*)'dimension of density',irho
       write(1)irho
       if(irho.gt.0) then 

       write(1) (den(k)%lam,k=1,irho)
       write(1) (den(k)%aa,k=1,irho)
       write(1) (den(k)%bb,k=1,irho)
       write(1) (den(k)%sig,k=1,irho)
       write(1) (den(k)%val,k=1,irho)
       end if
       close(1)
      

      end do !lam1
      deallocate(den)
      return
      end subroutine
