      subroutine read_density(par,inf,ip,Jp,lam)
     

       USE technical

       include 'commons.inc'

       character :: ch1
       character*2 :: ch2
       character*5 :: ch3
    
       integer :: mm,par,inf,ip,Jp,lam

!ipart=1 for neutrons
!ipart=-1 for protons
!info=1 for particle
!info=-1 for hole

       if(allocated(den)) deallocate(den)
       if(allocated(dens)) deallocate(dens)

       if(ip.eq.-1) ch1='0'
       if(ip.eq.1) ch1='1'
612    format(i2.2)
       write(ch2,612) Jp
613    format(i5.5)
       write(ch3,613) lam
       !write(*,*) ip,Jp,lam,ch1,ch2,ch3
      if (par.eq.1.and.inf.eq.1) then
       open(1,file='../density/dens/Neut_den_'//ch1//'_'//ch2//'_'//
     & ch3//'.out',status='old',form='unformatted')
       end if

      if (par.eq.1.and.inf.eq.-1) then
       open(1,file='../density/dens/Neut_hol_den_'//ch1//'_'//ch2//'_'//
     & ch3//'.out',status='old',form='unformatted')
       end if

       if (par.eq.-1.and.inf.eq.1) then
       open(1,file='../density/dens/Prot_den_'//ch1//'_'//ch2//'_'//
     & ch3//'.out',status='old',form='unformatted')
       end if

       if (par.eq.-1.and.inf.eq.-1) then
       open(1,file='../density/dens/Prot_hol_den_'//ch1//'_'//ch2//'_'//
     & ch3//'.out',status='old',form='unformatted')
       end if
         allocate(den(num_ph,ipnmin:ipnmax,ipnmin:ipnmax,0:2*jjmax))
         den=0.d0
         read(1) mm
         if(mm.gt.0) then
         allocate(dens(mm))
         read(1) (dens(k)%lam,k=1,mm)
         read(1) (dens(k)%aa,k=1,mm)
         read(1) (dens(k)%bb,k=1,mm)
         read(1) (dens(k)%sig,k=1,mm)
         read(1) (dens(k)%val,k=1,mm)
         end if
        close(1)
         do i=1,mm
         den(dens(i)%lam,dens(i)%aa,dens(i)%bb,dens(i)%sig)=dens(i)%val
         end do

        deallocate(dens)
      return
      end subroutine
