      subroutine cTDA(nn,ipar,ijj)!,idp,idn,idmax)
     

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
       

           if(allocated(cph)) deallocate(cph)
           !if(allocated(cn)) deallocate(cn)

       if(ipar.eq.-1) ch1='0'
       if(ipar.eq.1) ch1='1'

613    format(i2.2)
614    format(i2.3) 

       write(ch2,613) ijj
       write(ch3,614) nn
!        cph=0.d0

      open(12,file='../TDA/phon/phon_'
     &              //ch3//'_'//ch1//'_'//ch2//'.out',
     &                     status='old',form='unformatted')
       read(12) idmax,idp,idn
 
           allocate(cph(1:ipnmax,1:ipnmax,-1:1))!,cp(idp)) 
           cph=0.d0
           read(12) w

           do ii=1,idp
          read(12)a
          read(12)b
          read(12)c
          read(12)dd

          !cp(ii)%en=w
          !cp(ii)%par=a
          !cp(ii)%hol=b
          !cp(ii)%tz=c
          cph(a,b,c)=dd
          end do

          do ii=1,idn
          read(12)e
          read(12)f
          read(12)g
          read(12)hh 
          !cn(ii)%en=w
          !cn(ii)%par=e
          !cn(ii)%hol=f
          !cn(ii)%tz=g
          cph(e,f,g)=hh
          end do
         
        close(12)
       return
      end subroutine
