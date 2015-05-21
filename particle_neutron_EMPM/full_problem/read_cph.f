      subroutine cph(nn,ipar,ijj,idp,idn,idmax)
     

       USE technical

       include 'commons.inc'

       character :: ch1
       character*2 :: ch2
       character*3 :: ch3

    
       integer :: idp,idn        !dimension of proton and neutron subspace
       integer :: idmax          !dimension of TDA matrix
       integer :: a,b,c,e,f,g
       double precision :: dd,hh
       integer :: i,j,k          !service variables
       double precision :: val,w !service variable       
       
       
       !type(cph_type),allocatable,save :: ccp(:),ccn(:)

!ipar=1 for parity=1
!ipart=0 for parity=-1
           if(allocated(cp)) deallocate(cp)
           if(allocated(cn)) deallocate(cn)
       !write(*,*)'par',ipar
       if(ipar.eq.-1) ch1='0'
       if(ipar.eq.1) ch1='1'

613    format(i2.2)
614    format(i2.3) 

       write(ch2,613) ijj
       write(ch3,614) nn

       !write(*,*)'c3=',ch3,'c1=',ch1,'c2=',ch2
         open(12,file='phon/phon_'//ch3//'_'//ch1//'_'//ch2//'.out',
     &                     status='old',form='unformatted')
       read(12) idmax,idp,idn
         !write(*,*) idmax,idp,idn 
           allocate(cn(idn),cp(idp)) 
           read(12) w
           !write(*,*)w
           do ii=1,idp
          read(12)a
          read(12)b
          read(12)c
          read(12)dd
          !write(*,*) a,b,c,dd
          cp(ii)%en=w
          cp(ii)%par=a
          cp(ii)%hol=b
          cp(ii)%tz=c
          cp(ii)%cph=dd
          end do
          do ii=1,idn
          read(12)e
          read(12)f
          read(12)g
          read(12)hh 
          cn(ii)%en=w
          cn(ii)%par=e
          cn(ii)%hol=f
          cn(ii)%tz=g
          cn(ii)%cph=hh
          end do
      
         ! if(ipar.eq.-1.and.ijj.eq.0)then
         ! write(145,*)ccp(1)%en,ccp(1)%par,ccp(1)%hol,ccp(1)%cph
         ! end if
         
        close(12)
       return
      end subroutine
