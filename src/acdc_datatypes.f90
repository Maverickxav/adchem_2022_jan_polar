module acdc_datatypes

    USE second_Precision, ONLY : dp    
    Use constants
    
    implicit NONE
    
    INTEGER                                :: n_clustering_vapors=2         ! Number of clustering vapors
    INTEGER                                :: n_clustering_vapors_3comp=3         ! Number of clustering vapors
    logical                                :: use_clustering_plugin        !! true if clusterin plugin is used 
    logical                                :: clust_evap, clust_firstcall


    type :: clustering_mod
    CHARACTER(LEN=11),allocatable          :: names_vapor(:)   ! Vapor names for clustering species
    
    REAL(dp) ,allocatable                  :: conc_vapor(:) , nmols_evap(:)
    real(dp)                               :: conc_out                  !! total outgoing conc from ACDC
    
    real(dp), allocatable                  :: conc_coag_molec(:,:)                            ! Concentrations of coagulated molecules on aerosol bins (size: bins, molec) (1/m^3)
    real(dp), allocatable                  :: conc_coag_clust(:,:)                            ! Concentrations of coagulated clusters on aerosol bins (size: bins, clusters) (1/m^3)
    integer, allocatable                   :: clust_molec(:,:)                                      ! Cluster composition as numbers of vapor molecules (not specifying their possible charge; size: clusters, molec)
    real(dp), allocatable                  :: conc_out_bin(:)                                ! Concentrations of outgrown clusters to different aerosol bins (size: bins) (1/m^3)
    real(dp), allocatable                  :: comp_out_bin(:,:)                            ! Composition of outgrown clusters to different aerosol bins (size: bins, molec) (molec; note: real, not integer)
    
    REAL(DP)                               :: mass_par(nr_bins)
    integer                                :: nclust_syst, neq_syst, nclust_out
    REAL(dp), ALLOCATABLE                  :: c_p_clust(:,:), v_clust(:)
    
    REAL(DP),allocatable                   :: Nconc_evap(:)
    REAL(DP),allocatable                   :: conc_out_all(:), comp_out_all(:,:) ! Concentrations and compositions of individual outgrown clusters (optional) clust_out_molec=nclust,number of vapors 
                                                                                          !  conc_out_all =nclust_molec
    INTEGER, ALLOCATABLE                   :: clust_out_molec(:,:)
    INTEGER, ALLOCATABLE                   :: ind_out_bin(:)

  
    end type
    

end module acdc_datatypes