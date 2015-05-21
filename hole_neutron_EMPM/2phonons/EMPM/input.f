       subroutine input
        
       use technical

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       integer :: i,il,ij,lll,ak,row,col             !single particle input
       integer :: a,b,cc,c,ff,ss,dimensi               !"iphonons input
       integer :: ind,par,jj,inum,dimea
       integer :: int1,int2,int3,int4
       integer, allocatable :: ll(:),part(:,:)

       double precision :: en,dd
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

       write(*,*) 'couple hole protons-phonons'
       write(*,*) 'TDA phonons from j=',jjmin,'to',jjmax
       write(*,*) '2 phonons from j=',jjmin2,'to',jjmax2
       write(*,*) 'phonon parity from',iparmin,iparmax
       write(*,*) '2 phonons from j=',iparmin2,iparmax2
       write(*,*) 'iv=',iv
!**********************************************************************
!********************single particle levels****************************
!**********************************************************************

          open(12,file='../TDA/HF_n.out',status='old',form='formatted')
          read(12,*)
 
         do while (.not.eof(12))
        read(12,*)i,il,ij,e,qe,v
       enddo
        imax=i
        allocate(levn(imax))
        close(12)
          open(12,file='../TDA/HF_n.out',status='old',form='formatted')
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

       open(12,file='../TDA/HF_p.out',status='old',form='formatted')
      read(12,*)

       do while (.not.eof(12))
       read(12,*)i,il,ij,e,qe,v

       enddo
        allocate(levp(i))
       close(12)
        open(12,file='../TDA/HF_p.out',status='old',form='formatted')
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
!********************phonons**********************************************
!*************************************************************************
       open(4,file='infobase.dat',status='old',form='formatted')
      
       dimensi=0
       do while(.not.eof(4))
        read(4,*) as,ae,ar
         dimensi=dimensi+1
        end do
        close(4)
        dimbas=dimensi
         write(*,*) 'dimension of coup. basis=',dimbas
       open(5,file='infobase.dat',status='old',form='formatted')
          allocate(bas2ph(dimbas))
          do while(.not.eof(5))
         read(5,*) ii,ir,ilam
         bas2ph(ii)%ir=ir
         bas2ph(ii)%ilam=ilam
          end do
!**********************************************************************
!***********************************************************************************
      open(12,file='2ph/2f_x.dat',status='old',form='unformatted')
      allocate(dime2ph(-1:1,0:jjmax2))
      ss=0
      dime2ph=0

      do while(.not.eof(12))
       read(12) int1,int2,int3,int4
        dime2ph(int1,int2)=int3
         do j=1,int3
         read(12)(ila,ilap,tt,i=1,int4)
          end do
         ss=ss+int3
         end do
       close(12)
       write(*,*) ss
       allocate(ll(0:100))
       allocate(part(-1:1,0:jjmax2))
       ll=0
       part=0
       tt=0
       do ii=-1,1
        do jj=0,jjmax2
          if(ii.ne.0) then
          tt=tt+1
         ll(0)=0
         ll(tt)=ll(tt-1)+dime2ph(ii,jj)
         part(ii,jj)=ll(tt-1)
         end if
        end do
       end do

      open(1,file='2ph/2f_states.dat',status='old',form='unformatted')
         dimea=0
         do while(.not.eof(1))
          read(1) int1,int2,int3,en
          dimea=dimea+1
          end do
        close(1)
       num_2ph=dimea
          ll=0
      open(1,file='2ph/2f_states.dat',status='old',form='unformatted')
          allocate(ph2(dimea))
         do while(.not.eof(1))
          read(1) int1,int2,int3,en

          ph2(int1)%par=int2
          ph2(int1)%J=int3
          ph2(int1)%el=en
          ph2(int1)%alfa=int1-part(int2,int3)
          end do
        close(1)
        write(*,*)'num_2ph=',num_2ph


         write(*,*) '2 phon energy therdshold?'
         read(*,*) en_cut
         !en_cut=30.d0
         write(*,*) '2 phon energy therdshold=',en_cut
         lll=1
        Allocate(good(num_2ph))
        do i=1,num_2ph
          if(ph2(i)%el.le.en_cut) then
           good(lll)=i
           lll=lll+1
           write(12,*) i,ph2(i)%alfa
          end if
         end do
          dim2b=lll-1

           deallocate(part,ll)
!**********************************************************
!**********************************************************

       open(6,file='dmat.dat',status='old',form='formatted')
        allocate(dmat(dimbas,dimbas))
       dmat=0.d0
       do while(.not.eof(6))
        read(6,*) row,col,dd
        dmat(row,col)=dd
        
       end do
       close(6)
!*******************************************************************************

      open(7,file='d1matr.dat',status='old',form='formatted')
      read(7,*) ak,no
      write(*,*) ak,no
      allocate(d1(no,ak))
      do while(.not.eof(7))
      read(7,*) i,j,d1(i,j)
      if(dabs(d1(i,j)).gt.1.d0) then
      write(134,*) i,j,d1(i,j)
      end if
      end do
  
      close(7)





         return
        end subroutine
