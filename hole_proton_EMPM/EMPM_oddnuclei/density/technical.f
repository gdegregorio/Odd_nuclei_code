      MODULE technical

       include 'typedef.inc'

        PUBLIC :: levp,levn,cp,cn,eps

        type(level_type),dimension(:),allocatable, save :: levp,levn
        type(ph_type)   ,dimension(:),allocatable, save :: ph
        type(dens_type) ,dimension(:), allocatable, save :: den
        type(dimen_type),dimension(:),allocatable,save :: dime
        type(cph_type)  ,dimension(:),allocatable,save :: cp,cn,
     &                                                 C1p,C1n,C2p,C2n
         
        Integer, Save :: id,iv

        Integer, save :: Jtot,ipartot 
        Integer, save :: partc,state      
        Double precision,save ::eps

       END MODULE technical
