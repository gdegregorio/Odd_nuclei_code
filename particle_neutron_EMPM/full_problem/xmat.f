       subroutine xmatr

!********************************************************************
!this program generate the interaction part of EMPM for odd nuclei
!********************************************************************

       use technical
       use geom

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'
 
       integer :: iii,jj
       double precision :: xx,kk

       double precision,allocatable ::eigv(:,:)
    
       allocate(eigv(no,no))
       allocate(xm(dimbas,no))    
       open(44,file='eigvDAD.dat',status='unknown',form='formatted')
        do while(.not.eof(44))
       read(44,*) ii,jj,kk
        eigv(ii,jj)=kk
        end do
        close(44)
 
        do i=1,dimbas
         do k=1,no
            xx=0.d0 
           do j=1,no
            xx=xx+dmat(i,j)*eigv(j,k)
            !write(333,*) i,k,j
            !write(333,*) dmat(i,j),eigv(j,k),xx
           end do 
           xm(i,k)=xx!dsqrt(dble(Jtot+1))*xx
          !write(333,*) i,k,xm(i,k)
         end do
        end do
   
       return

        end subroutine
