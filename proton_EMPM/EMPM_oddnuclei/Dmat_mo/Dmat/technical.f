      MODULE technical

       include 'typedef.inc'

        PUBLIC :: levp,levn,den 


        type(level_type),dimension(:),allocatable, save :: levp,levn
        type(ph_type),   dimension(:),allocatable, save :: ph
        type(dens_type) , dimension(:), allocatable, save :: dens,h_den
        type(dimen_type),dimension(:),allocatable,save :: dime
        type(coupbas_type),dimension(:),allocatable,save ::indbase 
        Integer, Save :: id,dimbas,iv

        Integer, save :: Jtot,ipartot 
        Integer, save :: partc,state      
        Double precision,save ::eps
      
        
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: Dmat(:,:)
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: den(:,:,:,:)
    
        
       END MODULE technical
