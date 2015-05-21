    
      subroutine read_density(ig,ijj,ipp,ll)


      use technical
      include 'commons.inc'

      integer :: lam,aa,bb,sig,ijj,ipp
      integer :: ig,ndgg,igg,ndro,ll
      double precision :: val
        character :: ch1
        character*2 :: ch2

604    format(i2.2)

       if(ipp.eq.-1) ch1='0'
       if(ipp.eq.1)  ch1='1'

       write(ch2,604) ijj
       open(15,file='2ph/2f_rph_'//ch1//'_'//ch2//'.out',
     &            status='unknown',form='unformatted')

      ndro=5000000
      ndgg=0


      if(allocated(den)) deallocate(den)
      allocate(den(num_2ph,ippmin:ippmax,ippmin:ippmax,0:2*jjmax2))

      if(allocated(dens)) deallocate(dens)

      do while (.not.eof(15))
       read(15)igg,ndgg

       if (ndgg.gt.ndro) then
                write(*,*)'WARNING: Increase dimension in readro'
                stop
       endif

      
       if(ndgg.gt.0) then 
        if(igg.eq.ig) then
       ll=ndgg
       allocate(dens(ndgg))
       read(15)(dens(ii)%lam,ii=1,ndgg)
       read(15)(dens(ii)%sig,ii=1,ndgg)
       read(15)(dens(ii)%aa,ii=1,ndgg)
       read(15)(dens(ii)%bb,ii=1,ndgg)
       read(15)(dens(ii)%val,ii=1,ndgg)

         else if(igg.ne.ig) then
       read(15)(lam,ii=1,ndgg)
       read(15)(sig,ii=1,ndgg)
       read(15)(aa,ii=1,ndgg)
       read(15)(bb,ii=1,ndgg)
       read(15)(val,ii=1,ndgg)
        end if
  
        end if
        
        end do
        close(15)
        den=0.d0
        do i=1,ll
          den(dens(i)%lam,dens(i)%aa,dens(i)%bb,dens(i)%sig)=dens(i)%val
        end do
      deallocate(dens)

      return
      end subroutine read_density

