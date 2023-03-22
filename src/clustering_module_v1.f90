module ClusterIn_plugin

USE second_Precision,  ONLY : dp    ! KPP Numerical type
USE acdc_datatypes
include 'cluster_chem_use.inc'

implicit none
PRIVATE
public:: clustering_subroutine, allocate_chem_dimensions

contains

subroutine allocate_chem_dimensions(chem_v, MX_v,qX_v, c_clusters,c_clusters_v, nz, n_clustering_vapors, kk)

    implicit NONE
    type(clustering_mod), intent(inout) ::chem_v
    real(dp), allocatable, intent(inout) :: Mx_v(:), qX_v(:),c_clusters(:,:),c_clusters_v(:)
    INTEGER, INTENT(IN):: nz, n_clustering_vapors,kk

    
    allocate(chem_v%names_vapor(n_clustering_vapors)) 
    allocate(chem_v%conc_vapor(n_clustering_vapors))   
    allocate(chem_v%nconc_evap(nz))
    allocate(Mx_v(n_clustering_vapors))
    allocate(qX_v(n_clustering_vapors))
    allocate(chem_v%nmols_evap(n_clustering_vapors))
   
    if (kk .eq. 1) then
        call get_system_size_1(neq_syst=chem_v%neq_syst)
    elseif (kk.eq.2) then
        call get_system_size_2(neq_syst=chem_v%neq_syst)
    elseif (kk.eq.3) then
        call get_system_size_3(neq_syst=chem_v%neq_syst)
    elseif (kk.eq.4) then
        call get_system_size_4(neq_syst=chem_v%neq_syst)
    elseif (kk.eq.5) then
        call get_system_size_4(neq_syst=chem_v%neq_syst)
    end if

    allocate(c_clusters(nz,chem_v%neq_syst))
    allocate(c_clusters_v(chem_v%neq_syst))
   
    c_clusters=0d0; c_clusters_v=0d0

end Subroutine allocate_chem_dimensions

subroutine initialize_chem(chem_v, kk, n_clustering_vapors)
    implicit none
    type(clustering_mod), intent(inout) ::chem_v
    INTEGER, INTENT(IN):: kk, n_clustering_vapors
    
    chem_v%Nconc_evap=0D0
    
    
    ALLOCATE(chem_v%conc_coag_molec(nr_bins,n_clustering_vapors))

    
    if (kk .eq. 1) then
        write(*,*) 'calling system 1'
        call get_system_size_1(neq_syst=chem_v%neq_syst,nclust_syst=chem_v%nclust_syst,nclust_out=chem_v%nclust_out)
    elseif (kk .eq. 2) then
        write(*,*) 'calling system 2'
        call get_system_size_2(neq_syst=chem_v%neq_syst,nclust_syst=chem_v%nclust_syst,nclust_out=chem_v%nclust_out)
    elseif (kk .eq. 3) then
        write(*,*) 'calling system 3'
        call get_system_size_3(neq_syst=chem_v%neq_syst,nclust_syst=chem_v%nclust_syst,nclust_out=chem_v%nclust_out)
    elseif (kk .eq. 4) then
        write(*,*) 'calling system 4'
        call get_system_size_4(neq_syst=chem_v%neq_syst,nclust_syst=chem_v%nclust_syst,nclust_out=chem_v%nclust_out)
    end if

    allocate(chem_v%clust_molec(chem_v%nclust_syst, n_clustering_vapors))
    ALLOCATE(chem_v%conc_coag_clust(nr_bins, chem_v%nclust_syst))
    allocate(chem_v%c_p_clust(NSPEC_P, chem_v%nclust_syst))
    allocate(chem_v%v_clust(chem_v%nclust_syst))      
    ALLOCATE(chem_v%conc_out_all(chem_v%nclust_out))
    ALLOCATE(chem_v%clust_out_molec(chem_v%nclust_out,n_clustering_vapors))
    ALLOCATE(chem_v%comp_out_all(chem_v%nclust_out,n_clustering_vapors))
    ALLOCATE(chem_v%ind_out_bin(chem_v%nclust_out))

    chem_v%clust_molec=0D0
    chem_v%c_p_clust = 0.D0
    chem_v%v_clust = 0.D0
    chem_v%conc_coag_molec=0d0
    chem_v%clust_out_molec=0.
    chem_v%conc_out_all=0d0
    chem_v%comp_out_all=0d0
    chem_v%ind_out_bin=0d0

    allocate(chem_v%conc_out_bin(nr_bins))
    allocate(chem_v%comp_out_bin(nr_bins,n_clustering_vapors))
    

end subroutine initialize_chem

Subroutine clustering_subroutine(chem_1, chem_2,chem_3, chem_4, clust_firstcall, n_clustering_vapors,n_clustering_vapors_3comp, nr_bins,&
     cH2SO4, cDMA, cNH3,cHIO2, cHIO3, cMSA,&
     Jnucl_N_out, Jnucl_D_out, Jnucl_I_out,Jnucl_Ms_out, &
     dia_n, dia_dma,dia_hio,dia_msa, CS_H2SO4,T,press, ipr,dt, d, d_p, m_p, N_bins_in, Mx, qX,n_evap, comp_evap, &
     clustering_systems, c_clusters1, c_clusters2,c_clusters3, c_clusters4, layer,&
      l_cond_evap, l_coag_loss, three_comp_sys, two_comp, AN, AD, IiIo, AMsD )
   
    IMPLICIT NONE

    type(clustering_mod), intent(inout),optional ::chem_1, chem_2,chem_3, chem_4
    REAL(DP) :: c_dma , C_ACID, C_BASE, c_Ii, c_IO, sa_old, dm_old, am_old, hio2old, hio3old, nconc_evap1,&
                 c_msa,    nconc_evap2, nconc_evap3, nconc_evap4
    logical, intent(inout) :: clust_firstcall
    integer, intent(in) :: n_clustering_vapors, nr_bins, clustering_systems, layer, n_clustering_vapors_3comp
    REAL(DP), intent(inout), optional :: cH2SO4,cDMA, cNH3,cHIO2, cHIO3 ,cMSA
    REAL(DP), intent(out), optional :: Jnucl_N_out, Jnucl_D_out, Jnucl_I_out, Jnucl_Ms_out, dia_n, dia_dma, dia_hio, dia_msa
    REAL(DP), intent(in) :: CS_H2SO4,T,press,ipr,dt, d_p(:),d(nr_bins+1), N_bins_in(:), m_p(nr_bins)
    REAL(dp), INTENT(in) :: MX(NSPEC_P),qX(NSPEC_P)
    REAL(DP), intent(in):: n_evap, comp_evap(:)
    REAL(dp), intent(inout), optional ::  c_clusters1(:),c_clusters2(:),c_clusters3(:),c_clusters4(:)
    
    ! real(dp), optional :: 
    
    integer:: i, kk    
    logical, intent(in):: l_cond_evap,l_coag_loss, three_comp_sys, two_comp,AN, AD, IiIo, AMsD
    

    if (clust_firstcall) then

        ! if (two_comp) then    
        if (AN) then    
            CALL initialize_chem(chem_1,1, n_clustering_vapors)
            write(*,*) 'nceq_syst', chem_1%neq_syst, '', chem_1%nclust_syst,'',chem_1%nclust_out
        end if    
        if (AD) then     
            CALL initialize_chem(chem_2,2, n_clustering_vapors)
            write(*,*) 'nceq_syst', chem_2%neq_syst, '', chem_2%nclust_syst,'',chem_2%nclust_out
        end if
        if (IiIo) then     
            CALL initialize_chem(chem_3,3, n_clustering_vapors)
            write(*,*) 'nceq_syst', chem_3%neq_syst, '', chem_3%nclust_syst,'',chem_3%nclust_out
        end if
        if (AMsD) then     
            CALL initialize_chem(chem_4,4, n_clustering_vapors_3comp)
            write(*,*) 'nceq_syst', chem_4%neq_syst, '', chem_4%nclust_syst,'',chem_4%nclust_out
        end if    
            
        clust_firstcall=.false.      
      

    end if

    ! nconc_evap1=chem_1%Nconc_evap(layer)
    ! nconc_evap2=chem_2%Nconc_evap(layer)
    ! if (two_comp) then
        c_acid = cH2SO4*1D6 ! from molec/cm3  -> molec/m3 
        c_base = cNH3*1D6 
        c_DMA  = cDMA*1D6
        c_Ii   = cHIO3*1D6 
        c_IO   = cHIO2*1D6
        c_MSA  = cMSA*1d6
        
    if (AN) then     
        !cH2SO4 = cH2SO4 + comp_evap(1)
        cNH3   = cNH3   + comp_evap(4)
    end if    
    if (AD) then 
        !cH2SO4 = cH2SO4 + comp_evap(1)
        cDMA   = cDMA   + comp_evap(11)
    end if
    if (IiIo) then 
        cHIO3  = cHIO3  + comp_evap(10)
        cHIO2  = cHIO2  + comp_evap(12)
    end if
    if (AMsD) THEN
        !cH2SO4 = cH2SO4 + comp_evap(1)
        cMSA   = cMSA   + comp_evap(9)
        ! cDMA   = cDMA   + comp_evap(11)
    end if
    if (AN .or. AD .or. AMsD )then 
        ! write(*,*) 'here'
        cH2SO4 = cH2SO4 + comp_evap(1)
    end if
    if ( AD .or. AMsD )then 
        ! write(*,*) 'here'
        cDMA = cDMA + comp_evap(1)
    end if
    ! ! write(*,*) 'in clusterin module nconc_evap1 and 2' , nconc_evap1, nconc_evap2,layer
 
    if (n_evap > 0D0) then 
        ! write(*,*) 'goes in here'
        if (AN) then
            ! write(*,*) 'here'
            chem_1%nmols_evap(1) = comp_evap(1)*1.D6/n_evap ! H2SO4
            chem_1%nmols_evap(2) = comp_evap(4)*1.D6/n_evap ! NH3
        end if
        if(AD) then    
            chem_2%nmols_evap(1) = comp_evap(1)*1.D6/n_evap ! NH3
            chem_2%nmols_evap(2) = comp_evap(11)*1.D6/n_evap ! NH3
        end if
        if(IiIo) then 
            ! chem_3%nmols_evap(1) = comp_evap(10)*1.D6/n_evap ! HIO3
            ! chem_3%nmols_evap(2) = comp_evap(12)*1.D6/n_evap ! HIO2
            chem_3%nmols_evap(1) = comp_evap(10)*1.D6/n_evap ! HIO3
            chem_3%nmols_evap(2) = comp_evap(12)*1.D6/n_evap ! HIO2
        end if
        if(AMsD) then    
            chem_4%nmols_evap(1) = comp_evap(1)*1.D6/n_evap
            chem_4%nmols_evap(2) = comp_evap(9)*1.D6/n_evap ! MSA
            chem_4%nmols_evap(3) = comp_evap(11)*1.D6/n_evap ! DMA
        end if
        !end if
        ! if (three_comp_sys) then
        !     chem_4%nmols_evap(1) = comp_evap(1)*1.D6/n_evap ! H2SO4
        ! end if
    end if    



  !!!!!!!!! AN system !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
    if (AN) then
        ! write(*,*) 'here AN'
        ! c_acid = cH2SO4*1D6 ! from molec/cm3  -> molec/m3 
        ! c_base = cNH3*1D6 
        chem_1%conc_vapor=(/c_acid ,c_base/)
        ! write(*,*) chem_1%conc_vapor
        
        if (l_cond_evap) then
            
            CALL cluster_dynamics_1(chem_1%names_vapor,chem_1%conc_vapor,CS_H2SO4,T,ipr,dt,0.d0,&
            &    Jnucl_N_out,dia_n,c_inout=c_clusters1,naero=nr_bins,dp_aero_lim=d ,mp_aero=m_p,c_aero=N_bins_in,dp_aero=d_p,pres=Press, &
            c_evap=n_evap,nmols_evap=chem_1%nmols_evap,&
            c_coag_molec=chem_1%conc_coag_molec,c_coag_clust=chem_1%conc_coag_clust,clust_molec=chem_1%clust_molec, &
            c_out_bin=chem_1%conc_out_bin,comp_out_bin=chem_1%comp_out_bin,&
            c_out_all=chem_1%conc_out_all,clust_out_molec=chem_1%clust_out_molec)
            
            ! write(*,*) , sum(REAL(chem_1%clust_molec(:,1),KIND=dp)),  sum(REAL(chem_1%clust_molec(:,2),KIND=dp))
            elseif (l_coag_loss) then
                
                CALL cluster_dynamics_1(chem_1%names_vapor,chem_1%conc_vapor,CS_H2SO4,T,ipr,dt,0.d0,&
                &    Jnucl_N_out,dia_n,c_inout=c_clusters1, naero=nr_bins,dp_aero_lim=d ,dp_aero=d_p,mp_aero=m_p,c_aero=N_bins_in,pres=Press,&
                c_evap=n_evap,nmols_evap=chem_1%nmols_evap,&
                c_coag_clust=chem_1%conc_coag_clust,clust_molec=chem_1%clust_molec,&
                c_out_bin=chem_1%conc_out_bin,comp_out_bin=chem_1%comp_out_bin,&
                c_out_all=chem_1%conc_out_all,clust_out_molec=chem_1%clust_out_molec) 
                ! write(*,*) , sum(REAL(chem_1%clust_molec(:,1),KIND=dp)),  sum(REAL(chem_1%clust_molec(:,2),KIND=dp))
            else
                write(*,*) 'No aerosol-cluster feedback selected'
            end if    
            
            
            cH2SO4=chem_1%conc_vapor(1)*1d-6
            cNH3=chem_1%conc_vapor(2)*1d-6
    end if
        !!!!!!!!!!!!!!!!!!!! end AN system!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (AD) then
            !! update c_acid for AD system
            ! write(*,*) 'here AD'
            c_acid=cH2SO4*1D6 ! from molec/cm3  -> molec/m3 
            
            !!!!!!!!!!!!!! AD system !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            chem_2%conc_vapor=(/c_acid ,c_dma/)
            ! write(*,*) chem_2%conc_vapor
        
        
        if (l_cond_evap) then
        
            CALL cluster_dynamics_2(chem_2%names_vapor,chem_2%conc_vapor,CS_H2SO4,T,ipr,dt,0.d0,&
            &    Jnucl_D_out,dia_dma,c_inout=c_clusters2,naero=nr_bins,mp_aero=m_p,c_aero=N_bins_in,dp_aero_lim=d ,dp_aero=d_p,pres=Press, &
            c_evap=n_evap,nmols_evap=chem_2%nmols_evap, &
            c_coag_molec=chem_2%conc_coag_molec,c_coag_clust=chem_2%conc_coag_clust,clust_molec=chem_2%clust_molec,&
            c_out_bin=chem_2%conc_out_bin,comp_out_bin=chem_2%comp_out_bin,&
            c_out_all=chem_2%conc_out_all,clust_out_molec=chem_2%clust_out_molec)
            ! write(*,*) , sum(REAL(chem_2%clust_molec(:,1),KIND=dp)),  sum(REAL(chem_2%clust_molec(:,2),KIND=dp))
        elseif (l_coag_loss) then

            CALL cluster_dynamics_2(chem_2%names_vapor,chem_2%conc_vapor,CS_H2SO4,T,ipr,dt,0.d0,&
            &    Jnucl_D_out,dia_n,c_inout=c_clusters2, naero=nr_bins,dp_aero_lim=d ,dp_aero=d_p,mp_aero=m_p,c_aero=N_bins_in,pres=Press,&
            c_evap=n_evap,nmols_evap=chem_2%nmols_evap,&
            c_coag_clust=chem_2%conc_coag_clust,clust_molec=chem_2%clust_molec,&
            c_out_bin=chem_2%conc_out_bin,comp_out_bin=chem_2%comp_out_bin,&
            c_out_all=chem_2%conc_out_all,clust_out_molec=chem_2%clust_out_molec) 
            
            ! write(*,*) , sum(REAL(chem_2%clust_molec(:,1),KIND=dp)),  sum(REAL(chem_2%clust_molec(:,2),KIND=dp))
        else
            write(*,*) 'No aerosol-cluster feedback selected'
        end if
        
        !!! update acid and DMA!!!!
        cH2SO4=chem_2%conc_vapor(1)*1d-6
        cDMA=chem_2%conc_vapor(2)*1d-6
    
    end if    
    
    If (IiIo) then    
        !!!!!! HIO3-HIO2 system
        c_Ii   = cHIO3*1D6 
        c_IO   = cHIO2*1D6
        chem_3%conc_vapor=(/c_Ii ,c_Io/)
        ! write(*,*) 'Here IiIo'

        if (l_cond_evap) then
            CALL cluster_dynamics_3(chem_3%names_vapor,chem_3%conc_vapor,CS_H2SO4,T,ipr,dt,0.d0,&
            &    Jnucl_I_out,dia_hio,c_inout=c_clusters3,naero=nr_bins,mp_aero=m_p,c_aero=N_bins_in,dp_aero_lim=d ,dp_aero=d_p,pres=Press, &
            c_evap=n_evap,nmols_evap=chem_3%nmols_evap, &
            c_coag_molec=chem_3%conc_coag_molec,c_coag_clust=chem_3%conc_coag_clust,clust_molec=chem_3%clust_molec,&
            c_out_bin=chem_3%conc_out_bin,comp_out_bin=chem_3%comp_out_bin,&
            c_out_all=chem_3%conc_out_all,clust_out_molec=chem_3%clust_out_molec)
            ! write(*,*) 'here in here'
            
        elseif (l_coag_loss) then
        
            
            CALL cluster_dynamics_3(chem_3%names_vapor,chem_3%conc_vapor,CS_H2SO4,T,ipr,dt,0.d0,&
            &    Jnucl_I_out,dia_n,c_inout=c_clusters3, naero=nr_bins,dp_aero_lim=d ,dp_aero=d_p,mp_aero=m_p,c_aero=N_bins_in,pres=Press,&
            c_evap=n_evap,nmols_evap=chem_3%nmols_evap,&
            c_coag_clust=chem_3%conc_coag_clust,clust_molec=chem_3%clust_molec,&
            c_out_bin=chem_3%conc_out_bin,comp_out_bin=chem_3%comp_out_bin,&
            c_out_all=chem_3%conc_out_all,clust_out_molec=chem_3%clust_out_molec) 
            

        else
            write(*,*) 'No aerosol-cluster feedback selected'
        end if

    !!!!! update HIO3, HIO2 conc
        cHIO3=chem_3%conc_vapor(1)*1d-6
        cHIO2=chem_3%conc_vapor(2)*1d-6
    end if
    
    if (AMsD) then

        c_acid=cH2SO4*1D6
        c_dma=cDMA*1D6
        c_msa=cMSA*1D6
        
        chem_4%conc_vapor=(/C_ACID ,c_msa, c_DMA/)
        ! chem_4%conc_vapor=(/cH2SO4*1D6 ,cMSA*1D6, cDMA*1D6/)

        if (l_cond_evap) then
        
            CALL cluster_dynamics_4(chem_4%names_vapor,chem_4%conc_vapor,CS_H2SO4,T,ipr,dt,0.d0,&
            &    Jnucl_Ms_out,dia_msa,c_inout=c_clusters4,naero=nr_bins,dp_aero_lim=d ,mp_aero=m_p,c_aero=N_bins_in,dp_aero=d_p,pres=Press, &
            c_evap=n_evap,nmols_evap=chem_4%nmols_evap,&
            c_coag_molec=chem_4%conc_coag_molec,c_coag_clust=chem_4%conc_coag_clust,clust_molec=chem_4%clust_molec, &
            c_out_bin=chem_4%conc_out_bin,comp_out_bin=chem_4%comp_out_bin,&
            c_out_all=chem_4%conc_out_all,clust_out_molec=chem_4%clust_out_molec)

        elseif (l_coag_loss) then
            
            CALL cluster_dynamics_4(chem_4%names_vapor,chem_4%conc_vapor,CS_H2SO4,T,ipr,dt,0.d0,&
            &    Jnucl_Ms_out,dia_msa,c_inout=c_clusters4, naero=nr_bins,dp_aero_lim=d ,dp_aero=d_p,mp_aero=m_p,c_aero=N_bins_in,pres=Press,&
            c_evap=n_evap,nmols_evap=chem_4%nmols_evap,&
            c_coag_clust=chem_4%conc_coag_clust,clust_molec=chem_4%clust_molec,&
            c_out_bin=chem_4%conc_out_bin,comp_out_bin=chem_4%comp_out_bin,&
            c_out_all=chem_4%conc_out_all,clust_out_molec=chem_4%clust_out_molec) 
            
        else
            write(*,*) 'No aerosol-cluster feedback selected'
        end if    
        
        
        cH2SO4=chem_4%conc_vapor(1)*1d-6
        cMSA=chem_4%conc_vapor(2)*1d-6
        cDMA=chem_4%conc_vapor(3)*1d-6
            
            
    end if    

 

end subroutine clustering_subroutine



end module ClusterIn_plugin