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
       double precision :: xxv,kk

       allocate(xm(dimbas,no))    

       open(4,file='xmat.dat',status='unknown',form='formatted')

        do i=1,dimbas
         do k=1,no
            xxv=0.d0 
           do j=1,no
            xxv=xxv+dmat(i,j)*eigv(j,k)
           end do 
           xm(i,k)=dsqrt(dble(Jtot+1))*xxv
          write(4,*) i,k,xm(i,k)
         end do
        end do
         close(4)
       
       return

        end subroutine
