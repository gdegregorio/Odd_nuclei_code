       subroutine Vinteraction

!********************************************************************
!this program generate the interaction part of EMPM for odd nuclei
!********************************************************************

       use technical
       use geom

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'
 
       integer :: i,ip,ih,jp,jl,jh1,jp1
       double precision :: phase,somm

       double precision,allocatable ::Fcoup_prot(:),Fcoup_neut(:)
    
       allocate(Fcoup_prot(dimbas),Fcoup_neut(dimbas),Fcoup(dimbas))
      !need iv
       do i=1,dimbas
         ip=bas(i)%ir
         ilam=bas(i)%ilam       !1,numPh
         
         ipl=ph(ilam)%par
         Jpl=ph(ilam)%Jph
         ilamp=ph(ilam)%lam    !1,dimSubspace of given Jp and par
         jp=levp(ip)%j2


       if(allocated(C1p)) deallocate(C1p)
       if(allocated(C1n)) deallocate(C1n)

       call cph(ilamp,ipl,Jpl,ipp1,inn1,imax)

       allocate(C1p(ipp1),C1n(inn1))

       C1p=cp
       C1n=cn
         somm=0.d0
        do ip1=ipnmin,ipnmax
         jp1=levn(ip1)%j2
          do ih1=ihnmin,ihnmax
           jh1=levn(ih1)%j2
          do jj=1,inn1
          if(C1n(jj)%par.eq.ip1.and.C1n(jj)%hol.eq.ih1)then
            somm=somm+Fnn(ip1,ih1,iv,ip,Jpl)*C1n(jj)%cph
           end if
          end do
         end do
        end do
        Fcoup_neut(i)=0.5d0*somm
       end do

       do i=1,dimbas
         ip=bas(i)%ir
         ilam=bas(i)%ilam       !1,numPh

         ipl=ph(ilam)%par
         Jpl=ph(ilam)%Jph
         ilamp=ph(ilam)%lam    !1,dimSubspaceof given Jp and par
         jp=levp(ip)%j2
 

       if(allocated(C1p)) deallocate(C1p)
       if(allocated(C1n)) deallocate(C1n)

       call cph(ilamp,ipl,Jpl,ipp1,inn1,imax)

       allocate(C1p(ipp1),C1n(inn1))

       C1p=cp
       C1n=cn

        somm=0.d0
        do ip1=ippmin,ippmax
         jp1=levp(ip1)%j2
          do ih1=ihpmin,ihpmax
           jh1=levp(ih1)%j2
          do jj=1,inn1
          if(C1n(jj)%par.eq.ip1.and.C1n(jj)%hol.eq.ih1)then
            somm=somm+Fpn(ip1,ih1,iv,ip,Jpl)*C1p(jj)%cph
           end if
          end do
         end do
        end do
        Fcoup_prot(i)=somm
       end do

         do kk=1,dimbas
          Fcoup(kk)=Fcoup_prot(kk)+Fcoup_neut(kk)
          write(444,*) kk,Fcoup(kk)
         end do

        deallocate(Fcoup_prot,Fcoup_neut)
  
       return

        end subroutine
