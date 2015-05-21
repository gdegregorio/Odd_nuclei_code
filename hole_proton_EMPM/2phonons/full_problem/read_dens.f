      subroutine read_density(par,inf,ip,Jp,lam,mm)
     

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

       if(ip.eq.-1) ch1='0'
       if(ip.eq.1) ch1='1'
612    format(i2.2)
       write(ch2,612) Jp
613    format(i5.5)
       write(ch3,613) lam
         !write(*,*) ip,Jp,lam,par,inf
      if (par.eq.1.and.inf.eq.1) then
       open(1,file='Neut_den_'//ch1//'_'//ch2//'_'//ch3//'.out',
     &                     status='old',form='unformatted')
       end if

      if (par.eq.1.and.inf.eq.-1) then
       open(1,file='Neut_hol_den_'//ch1//'_'//ch2//'_'//ch3//'.out',
     &                     status='old',form='unformatted')
       end if

       if (par.eq.-1.and.inf.eq.1) then
       open(1,file='Prot_den_'//ch1//'_'//ch2//'_'//ch3//'.out',
     &                     status='old',form='unformatted')
       end if

       if (par.eq.-1.and.inf.eq.-1) then
       open(1,file='Prot_hole_den_'//ch1//'_'//ch2//'_'//ch3//'.out',
     &                     status='old',form='unformatted')
       end if
  
         read(1) mm
         allocate(den(mm))
         read(1) (den(k)%lam,k=1,mm)
         read(1) (den(k)%aa,k=1,mm)
         read(1) (den(k)%bb,k=1,mm)
         read(1) (den(k)%sig,k=1,mm)
         read(1) (den(k)%val,k=1,mm)
        close(1)

      return
      end subroutine
