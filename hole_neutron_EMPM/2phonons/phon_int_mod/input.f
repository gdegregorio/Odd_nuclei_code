       subroutine input
        
       use technical

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       integer :: i,il,ij             !single particle input
       integer :: a,b,cc,c,ff,ss              !phonons input
       integer :: ind,par,jj,inum,dimea
       integer :: int1,int2,int3,int4
       integer, allocatable :: ll(:),part(:,:)
       double precision :: en
       double precision :: e,qe,v     !single particle input

!**********************************************************************
!******************input file******************************************
!**********************************************************************

      open(1,file='input.dat',status='old',form='formatted')

      read(1,*)AN,AZ
      read(1,*)ihnmin,ihnmax
      read(1,*)ihpmin,ihpmax
      read(1,*)ipnmin,ipnmax
      read(1,*)ippmin,ippmax   
      read(1,*)iparmin,iparmax      !TDA
      read(1,*)jjmin,jjmax          !TDA
      read(1,*)iparmin2,iparmax2    !2ph
      read(1,*)jjmin2,jjmax2        !2ph
      read(1,*)iv
      close(1)

       write(*,*) 'H02 calculation for hole proton phonon bas.'
       write(*,*) 'TDA phonons from j=',jjmin,'to',jjmax
       write(*,*) '2 phonons from j=',jjmin2,'to',jjmax2
       write(*,*) 'phonon parity from',iparmin,iparmax
       write(*,*) '2 phonons from j=',iparmin2,iparmax2
       write(*,*) 'iv=',iv
!**********************************************************************
!********************single particle levels****************************
!**********************************************************************

      open(12,file='../TDA/HF_n.out',
     &status='old',form='formatted')
          read(12,*)
 
         do while (.not.eof(12))
        read(12,*)i,il,ij,e,qe,v
       enddo
        imax=i
        allocate(levn(imax))
        close(12)
        open(12,file='../TDA/HF_n.out',
     &    status='old',form='formatted')
           read(12,*)
         do while (.not.eof(12))
          read(12,*)i,il,ij,e,qe,v
          levn(i)%index=i
          levn(i)%l=il
          levn(i)%j2=ij
          levn(i)%ei=e
          levn(i)%ipar=(-1)**il
         enddo

        close(12)
        id=imax

       open(12,file='../TDA/HF_p.out',
     &          status='old',form='formatted')
      read(12,*)

       do while (.not.eof(12))
       read(12,*)i,il,ij,e,qe,v

       enddo
        allocate(levp(i))
       close(12)
        open(12,file='../TDA/HF_p.out',
     &                status='old',form='formatted')
        read(12,*)
        do while (.not.eof(12))
         read(12,*)i,il,ij,e,qe,v
         levp(i)%index=i
         levp(i)%l=il
         levp(i)%j2=ij
         levp(i)%ei=e
         levp(i)%ipar=(-1)**il

        enddo
        if (i.gt.imax) imax=i
          close(12)
       jmax=0

       do i=1,id
        if(levp(i)%j2.gt.jmax) jmax=levp(i)%j2
       enddo
       do i=1,id
        if(levn(i)%j2.gt.jmax) jmax=levn(i)%j2
       enddo

       Jtot=levp(iv)%j2
       ipartot=(-1)**levp(iv)%l

       write(*,*) 'J e parità dello stato da generare'
       write(*,*)  Jtot,ipartot

!*************************************************************************
!********************phonons**************** ******************************
!*************************************************************************
       open(1,file='../TDA/phon_inf.out',
     &                 status='old',form='formatted')
       num_ph=0 
       read(1,*) c
       allocate(dime(c))
       numb=c
       do is=1,c
       read(1,*) a,b,cc
       dime(is)%parity=a
       dime(is)%jph=b
       dime(is)%dimen=cc
       num_ph=num_ph+cc 
        end do

       close(1)
         write(*,*)'num_ph=',num_ph
       allocate(ph(num_ph))
!*******************************************************************************
        open(3,file='../TDA/list.out',
     &            status='old',form='unformatted')
          inum=0
          do while(.not.eof(3))
          inum=inum+1

         read(3) ind,par,jj,en

         ph(inum)%par=par
         ph(inum)%Jph=jj
         ph(inum)%lam=ind
         ph(inum)%el=en 
          end do

         close(3)
!***********************************************************************************
!**********************************************************************
!***********************************************************************************

       return
        end subroutine
