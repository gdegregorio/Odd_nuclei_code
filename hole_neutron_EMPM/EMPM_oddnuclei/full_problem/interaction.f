      subroutine Finteraction

       USE technical
       USE geom

       !include 'define.inc'
       include 'commons.inc'
      
       character :: ch
       double precision :: val
       integer :: a,b,c,d,J
       
       allocate(Fpp(ihpmin:ippmax,ihpmin:ippmax,ihpmin:ippmax,
     & ihpmin:ippmax,0:jmax))
       allocate(Fnn(ihnmin:ipnmax,ihnmin:ipnmax,ihnmin:ipnmax,
     & ihnmin:ipnmax,0:jmax))
       allocate(Fpn(ihpmin:ippmax,ihpmin:ippmax,ihnmin:ipnmax,
     & ihnmin:ipnmax,0:jmax))
       Fpp=0.d0
       Fnn=0.d0
       Fpn=0.d0
       a=0
       b=0
       c=0
       d=0
       J=0
       val=0.d0
       open(1,file='../TDA/fmat_pp.out',status='old',form='formatted')
       do while(.not.eof(1))
        read(1,*) a,b,c,d,J,val
        Fpp(a,b,c,d,J)=val
       end do
       close(1)
       a=0
       b=0
       c=0
       d=0
       J=0
       val=0.d0
       open(2,file='../TDA/fmat_nn.out',status='old',form='formatted')
       do while(.not.eof(2))
        read(2,*) a,b,c,d,J,val
        Fnn(a,b,c,d,J)=val
       end do
       close(2)
       a=0
       b=0
       c=0
       d=0
       J=0
       val=0.d0
       open(3,file='../TDA/fmat_pn.out',status='old',form='formatted')
       do while(.not.eof(3))
        read(3,*) a,b,c,d,J,val
        Fpn(a,b,c,d,J)=val
       end do
       close(3)

      
        return
      end
