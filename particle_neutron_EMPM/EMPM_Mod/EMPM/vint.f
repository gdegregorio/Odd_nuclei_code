       subroutine Vinteraction

!********************************************************************
!this program generate the interaction part of EMPM for odd nuclei
!********************************************************************

       use technical
       use geom

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'
 
       integer :: ii,jj,j_sig,sigma,sig,si,p1,p2,lam1,lam2
       integer :: pp,pl,pl1
       double precision :: fatt1,fatt2,fatt3,vint,kk

       double precision,allocatable ::Fcall(:,:,:,:)

       Allocate(Fcall(ipnmin:ipnmax,num_ph,ipnmin:ipnmax,num_ph))
       write(*,*) ipnmin,ipnmax,num_ph
       Allocate(Vcall(dimbas,dimbas))

       Vcall=0.d0

       Fcall=0.d0

       open(1,file='phon_int.dat',status='unknown',form='unformatted')
         do while(.not.eof(1))
         read(1) p1,p2,lam1,lam2,kk
         Fcall(p1,lam1,p2,lam2)=kk
         end do
        close(1)


        write(*,*) 'Fcall'
       fatt1=0.d0
       fatt2=0.d0
       fatt3=0.d0
    

       ii=0
       jj=0
       sig=0
       
       do ii=1,dimbas
         pp=bas(ii)%ir
         pl= bas(ii)%ilam
        do jj=1,dimbas
         p1=bas(jj)%ir
         pl1=bas(jj)%ilam
          write(119,*) pp,pl,p1,pl1
          Vcall(ii,jj)=Fcall(pp,pl,p1,pl1)
          end do
        end do
       deallocate(Fcall)
       return

        end subroutine
