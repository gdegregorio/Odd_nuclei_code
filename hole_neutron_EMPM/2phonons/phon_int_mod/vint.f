       subroutine Vinteraction

!********************************************************************
!this program generate the interaction part of EMPM for odd nuclei
!********************************************************************

       use technical
       use geom

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'
 
       integer :: i,ipl,ih,Jpl,jl,jh1,jp1,ilamp
       double precision ::somm
    
       allocate(Fcoup(ippmin:ippmax,ihpmin:ihpmax,num_ph))
       call Finteraction    
        do ilam=1,num_ph
         ipl=ph(ilam)%par
         Jpl=ph(ilam)%Jph
         ilamp=ph(ilam)%lam   
         write(*,*) ilam
       if(allocated(C1)) deallocate(C1)

       call cTDA(ilamp,ipl,Jpl)

       allocate(C1(1:ipnmax,1:ipnmax,-1:1)) 

       C1=0.d0
       C1=cph

       deallocate(cph)

       do ip=ippmin,ippmax
         do ih=ihpmin,ihpmax

          somm=0.d0  
!************************************************************     
        do ip1=ippmin,ippmax
          do ih1=ihpmin,ihpmax
        somm=somm+0.5d0*Fpp(ip,ih,ip1,ih1,Jpl)*C1(ip1,ih1,-1)
         end do
        end do

        do ip1=ipnmin,ipnmax
          do ih1=ihnmin,ihnmax      
         somm=somm+Fpn(ip,ih,ip1,ih1,Jpl)*C1(ip1,ih1,1)
         end do
        end do
!************************************************************
          Fcoup(ip,ih,ilam)=somm
         end do
        end do
       end do
  
       return

        end subroutine
