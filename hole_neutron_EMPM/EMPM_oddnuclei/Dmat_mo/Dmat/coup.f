       Subroutine coup_basis
 
       use technical
       use geom

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc' 
  
       integer :: jlam,jp,il,ip
       integer :: ii,jj,kk
       double precision ::en_th

      write(*,*)'gen. of coup. basis for neutrons holes id=',id,
     &    'and num p=',num_ph
      open(1,file='infobase.dat',status='unknown',form='formatted')
      open(2,file='information.dat',status='unknown',form='formatted')

      write(*,*) 'energy threshold for 1 phonon states?'
      read(*,*) en_th
      write(*,*) 'energy=',en_th
      write(2,*)'neutron hole coupled to phonons'
      write(2,*)'energy threshold for phonons=',en_th

      ind=0
       do kk=1,num_ph
        if(ph(kk)%el.le.en_th) then
         ind=ind+1
        end if
       end do

       write(*,*) 'number of phonons taken=',ind,'of',num_ph
       write(2,*) 'number of phonons taken=',ind,'of',num_ph

       do ii=ihnmin,ihnmax
        do jj=1,num_ph        
         jlam=ph(jj)%Jph
         jp=levn(ii)%j2
         il=ph(jj)%par
         ip=levn(ii)%ipar
         if(ip*il.eq.ipartot.and.ph(jj)%el.le.en_th) then
          if(abs(2*jlam-jp).le.Jtot.and.Jtot.le.(2*jlam+jp)) then
           dimbas=dimbas+1
          endif
         end if
        end do
       end do

       allocate(indbase(dimbas))
       kk=1
       do ii=ihnmin,ihnmax
       do jj=1,num_ph
       jlam=ph(jj)%Jph
       jp=levn(ii)%j2
       il=ph(jj)%par
       ip=levn(ii)%ipar
       if(ip*il.eq.ipartot.and.ph(jj)%el.le.en_th) then
       if(abs(2*jlam-jp).le.Jtot.and.Jtot.le.(2*jlam+jp)) then
         indbase(kk)%ir=ii
         indbase(kk)%ilam=jj
       write(1,*)kk,ii,jj
       kk=kk+1
          end if
         end if
        end do
       end do
       close(1)
       write(*,*) 'basis of dimension =',dimbas,'generated'
       deallocate(dime)
        return
       end subroutine
  
