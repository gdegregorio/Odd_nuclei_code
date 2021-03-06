       program Hamiltonian
  

       use technical
       use geom

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'
    
      Integer :: ii,jj,i,j,nor
      Integer :: jr,pr,pl,jl,ilam,ir

      Integer, allocatable :: mxt(:),mxt1(:),nxt(:)
      double precision ::er,el
     
      Double precision :: en,xval
      Double precision, allocatable :: ser1(:,:),ser2(:,:), ser3(:,:)
      Double precision, allocatable :: ser4(:,:), ser5(:),ser6(:)
      Double precision, allocatable :: Amat(:,:),ADM(:,:),ham(:,:)
      Double precision, allocatable :: AD1(:,:),ADD(:,:),ide(:,:)
      Double precision, allocatable :: wr(:),wi(:),ww(:),vr(:,:)

      call input
!***************************************************************
!****************A matrix***************************************
!***************************************************************
        open(11,file='Amat.dat',status='unknown',form='formatted')
       ALLOCATE(Amat(dimbas,dimbas),ser5(dimbas),ser6(dimbas)) 
       ser5=0.d0
       ser6=0.d0
       Amat=0.d0
       do ii=1,dimbas
        do jj=1,dimbas
       if (ii.eq.jj)then
       Amat(ii,ii)=levn(bas2ph(ii)%ir)%ei+ph2(bas2ph(ii)%ilam)%el
       ser5(ii)=Amat(ii,ii)
       end if
       if(Amat(ii,jj).gt.0.9) then
       write(11,*) ii,jj,Amat(ii,jj)
       end if!
       end do
        end do
        close(11)

        call order(ser5,dimbas,ser6)
       open(33,file='un_energies.dat',status='unknown',form='formatted')
       do i=1,dimbas
        write(33,*)ser6(i)
       end do
        close(33)
       deallocate(ser5,ser6)
       
!****************************************************************
  
       allocate(AD(dimbas,dimbas),ADM(dimbas,dimbas))
       allocate(ser1(dimbas,dimbas),ser2(dimbas,dimbas))
      
       call DGEMM('N','N',dimbas,dimbas,dimbas,1.d0,Amat,max(1,dimbas),
     &        dmat,max(1,dimbas),0.d0,AD,max(1,dimbas))

        open(12,file='AD_notord.dat',status='unknown',form='formatted')
       do i=1,dimbas
        do j=1,dimbas
       ADM(i,j)=AD(i,j)
       ser1(i,j)=AD(i,j)
       ser2(i,j)=dmat(i,j)
       if(AD(i,j).gt.0.9) then
       write(12,*)i,j,AD(i,j)
       end if

        end do
        end do
        close(12)


          deallocate(AD)
!**************************************************************************
!**************************************************************************
       allocate(mxt(dimbas),mxt1(dimbas))   
       open(6,file='mixt.dat',status='old',form='formatted')

      do while(.not.eof(6))
      read(6,*)i,mm
       mxt(i)=mm
       mxt1(mm)=i
      enddo!

      close(6)
************************************************************
!*********** permutations ******************************************
!*******************************************************************

       call permut(dimbas,dmat,mxt)
  
       call permut(dimbas,ADM,mxt)


!******************************************************************
!******************************************************************
!******************************************************************

       write(*,*) 'perm done'  
       open(114,file='AD_ord.dat',status='unknown',form='formatted')
      allocate(ADD(no,no))!dimbas,dimbas))
      DD=0.d0
      do i=1,no!dimbas
       do j=1,no!dimbas
        ADD(i,j)=ADM(i,j)
        if(ADM(i,j).gt.0.9) then
        write(114,*) i,j,ADM(i,j)
        end if
      
        end do
      end do
      close(114)
      allocate(d1m(no,no))
       d1m=0.d0
      allocate(ser3(no,no))
       do i=1,no
        do j=i,no
         d1m(i,j)=dmat(i,j)
         ser3(i,j)=dmat(i,j)
         if(dabs(d1m(i,j)).gt.1.d0) then
         write(243,*) i,j,d1m(i,j)
         end if
        end do
       end do
!*****!**********************************************************
!*******************************check permutation***************
!***************************************************************


!****************************************************************
!****************************************************************
!****************************************************************
       nor=no

       call dpotrf('U',nor,d1m,nor,info)
        write(*,*)' Factorization info ',info
       
       call dpotri('U',nor,d1m,nor,info)
          write(*,*)' inverse info ',info   
  
      do i=1,nor
       do j=i+1,nor
         d1m(j,i)=d1m(i,j)
         ser3(j,i)=ser3(i,j)
       enddo
      enddo
       open(238,file='d1.dat',status='unknown',form='formatted')
       do i=1,nor
        do j=1,nor
         write(238,*) i,j,d1m(i,j)
         if(dabs(d1m(i,j)).gt.1.d0) then
         write(135,*) i,j,d1m(i,j)
         end if
        end do
       end do
        close(238)

      deallocate(ser3)

      allocate(ham(nor,nor))     
      allocate(ser4(nor,nor))!,ser5(nor,nor))
       ham=0.d0
      do i=1,nor
       do j=1,nor
         ham(i,j)=ADM(i,j)
         ser4(i,j)=ADM(i,j)
       enddo
      enddo
      


      allocate(DAD(nor,nor))


      call dgemm('N','N',nor,nor,nor,1.d0,d1m,nor,ham,nor,0.d0,
     &      DAD,nor)

       open(241,file='symDAD.dat',status='unknown',form='formatted')
      open(240,file='DAD.dat',status='unknown',form='formatted')
      open(242,file='DAD_big.dat',status='unknown',form='formatted')
      do i=1,no!dimbas
       do j=1,no!dimbas
        write(240,*) i,j,DAD(i,j)
       if(abs(DAD(i,j)-DAD(j,i)).gt.0.0001) then
       write(241,*)i,j,DAD(i,j)-DAD(j,i)
       end if
         if(abs(DAD(i,j)).gt.0.1) then
       write(242,*)i,j,DAD(i,j)
       end if


       end do
      end do 

      deallocate(d1m,ham)
!
      write(*,*) 'diagonalization' 

       allocate(wr(no),wi(no),ww(no))
       allocate(vr(no,no))

       wr=0.d0
       wi=0.d0
       vr=0.d0
       ww=0.d0

        call diagonalization(DAD,wr,wi,vr,no)
    
         open(43,file='eigenDAD.dat',status='unknown',form='formatted')
         open(448,file='IM_DAD.dat',status='unknown',form='formatted')
           do i=1,no
           if(wi(i).gt.0.0000001) then
          write(448,*) wi(i),i
           end if
          write(43,*) i,wr(i)
           end do
         close(43)
         close(448)
  
         open(44,file='eigvDAD.dat',status='unknown',form='formatted') 
         write(44,*) '|v>   |p alpha> '

         do i=1,no!dimbas
          do j=1,no!dimbas
           
           if(dabs(vr(j,i)).gt.0.5) then
           write(44,*) i,j,vr(j,i)
           end if
          end do
         end do
          close(44)

           open(22,file='struct.dat',status='unknown',form='formatted')
          do i=1,no!dimbas
           write(22,*)'***************************'
           write(22,*) i,wr(i)
           write(22,*)'***************************'
           write(22,*)'jr,pr,er,jl,pl,el,c,c^2'
            do j=1,nor
               if(dabs(vr(j,i)).gt.0.01) then
               ir=bas2ph(mxt(j))%ir
               ilam=bas2ph(mxt(j))%ilam
               jr=levp(ir)%j2
               pr=(-1)**levp(ir)%l
               er=levp(ir)%ei
               jl=ph2(ilam)%J
               pl=ph2(ilam)%par
               el=ph2(ilam)%el
          write(22,'(i2,i2,1x,f10.5,2x,i2,i2,1x,f10.5,4x,f10.5,f10.5)')
     &                  jr,pr,er,jl,pl,el,vr(j,i),vr(j,i)**2
              end if
             end do
           end do 
          


!*********************************************************
!*********************************************************
!*********************************************************       
          
         call order(wr,no,ww)
 
        deallocate(wr,wi,vr)

        
       open(486,file='ordeig_DAD.dat',status='unknown',form='formatted')
          do i=1,no!dimbas
          write(486,'(i3,2x,f10.5)') 
     &      i,ww(i)
           end do

        write(*,*) 'end diagonalization'

       contains

!***********************************************************************
!***********************************************************************
!***********************************************************************
       subroutine permut(ndim,matr,mxt)

      implicit double precision (a-h,o-z)

      double precision, dimension(:,:), allocatable :: matr
      double precision, dimension(:), allocatable :: work
      integer, dimension(:), allocatable :: mxt,nxt
      write(*,*) matr(1,1)
      !write(*,*) 'ndim',ndim,mxt(1)
      allocate (work(ndim))
      work=0.d0
      allocate(nxt(ndim))

      do i=1,ndim
       nxt(i)=i
      enddo

      do i=1,ndim
       work(:)=matr(i,:)
       do j=i,ndim
        if (nxt(j).eq.mxt(i)) then
        matr(i,:)=matr(j,:)
        matr(j,:)=work(:)
        nxt(j)=nxt(i)
        nxt(i)=mxt(i)
       endif
       enddo
      enddo

      do i=1,ndim
       nxt(i)=i
      enddo

      do i=1,ndim
       work(:)=matr(:,i)
       do j=i,ndim
        if (nxt(j).eq.mxt(i)) then
        matr(:,i)=matr(:,j)
        matr(:,j)=work(:)
        nxt(j)=nxt(i)
        nxt(i)=mxt(i)
        endif
       enddo
      enddo



      deallocate(work)

      return
      end subroutine
!***************************************************************************
!***************************************************************************
!***************************************************************************

      subroutine diagonalization(A,er,ei,v,i)

       include 'define.inc'

       double precision :: A(i,i),er(i),ei(i),v(i,i)
       double precision :: wor(max(1,6*i))

       icou1=max(1,6*i)
      er=0.d0
      ei=0.d0
      v=0.d0
      call DGEEV('N','V',i,A,max(1,i),er,ei,vl,
     &       max(1,i),v,max(1,i),wor,icou1,info)

       if(info.ne.0) then
        write(*,*) 'The diagonalization failed!!!'
        stop
       endif

       return
      end SUBROUTINE
!********************************************************************************
!********************************************************************************


        subroutine order(vec,mm,wro)

        integer ::mm,ipoz(mm) 
        double precision ::vec(mm)
        double precision, allocatable :: wro(:)
         !allocate(wro(mm))
        
        
      do i=1,mm

        ipoz(i)=i

        wro(i)=vec(i)

      enddo



      do i=1,mm

        do j=1,mm-i

          if (wro(j).gt.wro(j+1)) then

            xe=wro(j)

            wro(j)=wro(j+1)

            wro(j+1)=xe

            ipozz=ipoz(j)

            ipoz(j)=ipoz(j+1)

            ipoz(j+1)=ipozz

      endif

      enddo

      enddo
       return
      end subroutine  
!**************************************************************
!**************************************************************     
      end program
