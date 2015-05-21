      MODULE technical

       include 'typedef.inc'

        PUBLIC :: levp,levn,dime,dime2ph,
     &            eps,bas,dmat,en_thr,dim2b

        type(level_type),dimension(:),allocatable, save :: levp,levn
        type(ph2_type)  ,dimension(:),allocatable, save :: ph2
        type(dens_type) ,dimension(:), allocatable, save :: dens
        !type(dens2_type) ,dimension(:), allocatable, save :: den2
        type(indbas_type), dimension(:),allocatable,save :: bas
        type(dimen_type),dimension(:),allocatable,save ::dime

        Integer, Save :: id,dimbas,no,iv,dim1ph,dimbs,dim2b
         
        Integer, allocatable, save :: dime2ph(:,:)
        Integer, save :: Jtot,ipartot,en_cut
        Integer, save :: partc,state  
    

        Double precision,save ::eps,en_thr
       Double precision,allocatable,save :: den(:,:,:,:)
       DOUBLE PRECISION, ALLOCATABLE, SAVE :: somma(:,:),ro(:,:),good(:)

        double precision,allocatable, save ::eigv(:,:) !2 ph eigv
       
        DOUBLE PRECISION, ALLOCATABLE, SAVE :: Dmat(:,:)
 
 
 
                
       END MODULE technical
