MODULE output

   USE second_Precision,  ONLY : dp    ! KPP Numerical type
   USE constants ! Set case specifics (module also includes physical constants)
   USE second_Parameters, ONLY : NSPEC, ind_H2SO4, ind_HNO3
   use netcdf 
   USE second_Monitor

   IMPLICIT NONE
   public :: open_write_nc_files


CONTAINS
    
    SUBROUTINE open_write_nc_files(openfname,writefname,closefile, ncids,general_ids,increment,RH,TEMP,pressure, species_names, time, conc,psat,&
                                    N_bins, pH, conc_pp, Jnucl_N, Jnucl_D, Jnucl_I, Jnucl_M, Jnucl_AND, CS_H2SO4)
        implicit NONE

        integer :: status
        INTEGER, PARAMETER :: N_FILES = 3
        INTEGER, PARAMETER :: shuff=1, compress=1, compression=9

        
        INTEGER, intent(out), optional      :: ncids(N_FILES)
        CHARACTER(220) :: ncfile_names(N_FILES) = (['General.nc  ', 'Chemistry.nc', 'Aerosols.nc '])
        integer :: i
        integer :: dt_id, nz_id, nbins_id, nspec_id, ncond_id, nspecp_id
        integer ::jnh3_id,jdma_id, jmsa_id, jhio3_id 
        integer ::rh_id,temp_id, cond_index_id, molar_mass_id, ph_id, psat_id, spec_name_id, dens_id, pres_id, time_id, conc_id
        integer ::nconc_id,jnucl_n_id,jnucl_d_id,jnucl_i_id,jnucl_m_id,jnucl_and_id, cs_id, concpp_id
        logical, intent(in):: openfname, writefname, closefile
        integer, intent(in), optional:: increment
        REAL(dp), DIMENSION(Nz), INTENT(in), optional :: RH,TEMP, pressure
        integer, intent(out), optional :: general_ids(nr_outputs)
        REAL(dp), INTENT(in), optional :: time
        character(len=15),INTENT(IN), optional:: species_names(nspec)
        REAL(dp), DIMENSION(Nz,NSPEC), INTENT(in), optional :: conc
        REAL(dp), DIMENSION(NCOND), INTENT(in), optional :: psat
        REAL(dp), DIMENSION(Nz,nr_bins), INTENT(in), optional :: N_bins,pH
        REAL(dp), DIMENSION(Nz,NSPEC_P,nr_bins), INTENT(in), optional :: conc_pp
        REAL(dp), DIMENSION(Nz), INTENT(in), optional :: Jnucl_N,Jnucl_D, Jnucl_I,Jnucl_M, Jnucl_AND, CS_H2SO4


        
        do i=1, N_files
          if (i==1) then
                IF (openfname) THEN    
                
                    status = nf90_create('output/netcdf/General.nc',IOR(NF90_NETCDF4, NF90_CLASSIC_MODEL), ncids(i))
                    call check(status, 'open')
                    
                    call check( nf90_def_dim(ncids(i), 'time', NF90_UNLIMITED, dt_id), 'def_dim' )
                    call check( nf90_def_dim(ncids(i), 'nz', nz, nz_id), 'def_dim nz' )
                    call check( nf90_def_dim(ncids(i), 'NPSEC',NSPEC, nspec_id), 'def_dim nspec' )
                    call check( nf90_def_dim(ncids(i), 'index_cond',NCOND, ncond_id), 'def_dim ncond' )

                    call check( nf90_def_var(ncids(i), 'RH', NF90_FLOAT, [nz_id,dt_id], rh_id),'def_var' )
                    call check( nf90_def_var(ncids(i), 'Pr', NF90_FLOAT, [nz_id,dt_id], pres_id),'def_var' )
                    call check( nf90_def_var(ncids(i), 'Temp', NF90_FLOAT, [nz_id,dt_id], temp_id),'def_var' )
                    call check( nf90_def_var(ncids(i), 'Species', NF90_CHAR, [nspec_id], spec_name_id),'def_var' )
                    call check( nf90_def_var(ncids(i), 'Time', NF90_FLOAT, [dt_id], time_id),'def_var' )
               
                    call check( nf90_put_att(ncids(i), NF90_GLOBAL, 'Information', 'This netcdf file contains general Meteorology variables saved along the trajectory'), 'att1')
                    call check( nf90_put_att(ncids(i), rh_id, 'units', '%'), 'rh')
                    call check( nf90_put_att(ncids(i), pres_id, 'units', 'Pa'), 'pressure')
                    call check( nf90_put_att(ncids(i), temp_id, 'units', 'K'), 'temp')
                    call check( nf90_put_att(ncids(i), spec_name_id, 'units', 'Character'), 'species_name')
                    call check( nf90_put_att(ncids(i), time_id, 'units', 's'), 'time')
                    call check( nf90_enddef(ncids(i)) )  
                    general_ids(1:4)=(/rh_id,pres_id,temp_id, time_id/)
                    
                    !! for variables saved only once
                    call check( nf90_put_var(ncids(i), spec_name_id,species_names, start=(/1/), count=(/NSPEC/) ),'write-spec' )
                    ! write(*,*) 'L52:', ncids(i),general_ids(1), openfname, writefname,i
                elseif (writefname) THEN
                  
                        call check( nf90_put_var(ncids(i), general_ids(1),RH,start=(/1, increment/), count=(/nz/) ),'write-rh' )
                        call check( nf90_put_var(ncids(i), general_ids(2),pressure,start=(/1, increment/), count=(/nz/) ),'write-pressure' )
                        call check( nf90_put_var(ncids(i), general_ids(3),temp,start=(/1, increment/), count=(/nz/) ),'write-temperature' )
                        call check( nf90_put_var(ncids(i), general_ids(4),time,start=(/increment/) ),'write-time' )
                end if  
            
            elseif (i == 2) then 
                IF (openfname) THEN    
                
                    status = nf90_create('output/netcdf/Chemistry.nc',IOR(NF90_NETCDF4, NF90_CLASSIC_MODEL), ncids(i))
                    call check(status, 'open')

                    call check( nf90_def_dim(ncids(i), 'time', NF90_UNLIMITED, dt_id), 'def_dim' )
                    call check( nf90_def_dim(ncids(i), 'nz', nz, nz_id), 'def_dim nz' )
                    call check( nf90_def_dim(ncids(i), 'NPSEC',NSPEC, nspec_id), 'def_dim nspec' )
                    call check( nf90_def_dim(ncids(i), 'ncond',NCOND, ncond_id), 'def_dim ncond' )
                    
                    call check( nf90_def_var(ncids(i), 'conc', NF90_FLOAT, [nz_id,nspec_id,dt_id], conc_id),'def_var' )
                    call check( nf90_def_var(ncids(i), 'psat', NF90_FLOAT, [ncond_id], psat_id),'def_var' )
                    
                    call check( nf90_def_var_deflate(ncids(i), conc_id, shuffle=1, deflate=1,deflate_level=9),'compression' )

                    call check( nf90_put_att(ncids(i), NF90_GLOBAL, 'Information', 'This netcdf file contains gas-phase concentrations of all species saved along the trajectory'), 'att1')
                    call check( nf90_put_att(ncids(i), conc_id, 'units', '[#/cm^3]'), 'conc')
                    call check( nf90_put_att(ncids(i), psat_id, 'units', '[atm]'), 'psat')
                    call check( nf90_enddef(ncids(i)) )  
                     !! for variables saved only once
                    
                    general_ids(5:6)=(/conc_id, psat_id/)
                    
                    elseif (writefname) THEN
                        ! write(*,*) 'L54:', ncids(i), general_ids(1), openfname, writefname,i
                        call check( nf90_put_var(ncids(i), general_ids(5),conc,start=(/1,1, increment/), count=(/nz, nspec,1/) ),'write-conc' )
                        call check( nf90_put_var(ncids(i), general_ids(6),psat, start=(/1/), count=(/Ncond/) ),'write-spec' )
                end if    
                
                elseif (i == 3) then 
                    IF (openfname) THEN    
                    
                        status = nf90_create('output/netcdf/Aerosols.nc',IOR(NF90_NETCDF4, NF90_CLASSIC_MODEL), ncids(i))
                        call check(status, 'open')
    
                        call check( nf90_def_dim(ncids(i), 'time', NF90_UNLIMITED, dt_id), 'def_dim' )
                        call check( nf90_def_dim(ncids(i), 'nz', nz, nz_id), 'def_dim nz' )
                        call check( nf90_def_dim(ncids(i), 'NPSEC',NSPEC, nspec_id), 'def_dim nspec' )
                        call check( nf90_def_dim(ncids(i), 'NPSEC_P',NSPEC_P, nspecp_id), 'def_dim nspec' )
                        call check( nf90_def_dim(ncids(i), 'nbins',nr_bins, nbins_id), 'def_dim nbins' )
                        
                        call check( nf90_def_var(ncids(i), 'nconc', NF90_FLOAT, [nz_id,nbins_id,dt_id], nconc_id),'def_var' )
                        call check( nf90_def_var(ncids(i), 'Jnucl_AN', NF90_FLOAT, [nz_id,dt_id], jnucl_n_id),'def_var' )
                        call check( nf90_def_var(ncids(i), 'Jnucl_AD', NF90_FLOAT, [nz_id,dt_id], jnucl_d_id),'def_var' )
                        call check( nf90_def_var(ncids(i), 'Jnucl_IiIo', NF90_FLOAT, [nz_id,dt_id], jnucl_i_id),'def_var' )
                        call check( nf90_def_var(ncids(i), 'Jnucl_AMsD', NF90_FLOAT, [nz_id,dt_id], jnucl_m_id),'def_var' )
                        call check( nf90_def_var(ncids(i), 'Jnucl_AND', NF90_FLOAT, [nz_id,dt_id], jnucl_and_id),'def_var' )
                        call check( nf90_def_var(ncids(i), 'CS_H2SO4', NF90_FLOAT, [nz_id,dt_id], cs_id),'def_var' )
                        call check( nf90_def_var(ncids(i), 'conc_pp', NF90_FLOAT, [nz_id,nspecp_id, nbins_id, dt_id], concpp_id),'def_var' )
                        call check( nf90_def_var(ncids(i), 'pH', NF90_FLOAT, [nz_id, nbins_id, dt_id], ph_id),'def_var' )

                        call check( nf90_def_var_deflate(ncids(i), conc_id, shuffle=1, deflate=1,deflate_level=9),'compression' )
                                           
                        call check( nf90_put_att(ncids(i), NF90_GLOBAL, 'Information', 'This netcdf file contains number and particle phase concentrations saved along the trajectory'), 'att1')
                        call check( nf90_put_att(ncids(i), nconc_id, 'units', '[#/m^3]'), 'conc')
                        call check( nf90_put_att(ncids(i), jnucl_n_id, 'units', '[#/m^3]'), 'jnucl_n')
                        call check( nf90_put_att(ncids(i), jnucl_d_id, 'units', '[#/m^3]'), 'jnucl_d')
                        call check( nf90_put_att(ncids(i), jnucl_i_id, 'units', '[#/m^3]'), 'jnucl_i')
                        call check( nf90_put_att(ncids(i), jnucl_m_id, 'units', '[#/m^3]'), 'jnucl_m')
                        call check( nf90_put_att(ncids(i), jnucl_and_id, 'units', '[#/m^3]'), 'jnucl_and')
                        call check( nf90_put_att(ncids(i), cs_id, 'units', '[s-1]'), 'CS_H2SO4')
                        call check( nf90_put_att(ncids(i), concpp_id, 'units', '[#/cm^3]'), 'CS_H2SO4')
                        call check( nf90_put_att(ncids(i), pH_id, 'units', '[ ]'), 'pH')

                        call check( nf90_enddef(ncids(i)) )  

                        general_ids(7:15)=(/nconc_id,jnucl_n_id,jnucl_d_id, jnucl_i_id, jnucl_m_id, jnucl_and_id, cs_id, concpp_id, ph_id/)
                        
                    elseif (writefname) THEN
                            ! write(*,*) 'L54:', ncids(i), general_ids(1), openfname, writefname,i
                            call check( nf90_put_var(ncids(i), general_ids(7),N_bins,start=(/1,1,increment/), count=(/nz, nr_bins,1/) ),'write-nconc' )
                            call check( nf90_put_var(ncids(i), general_ids(8),Jnucl_N,start=(/1,increment/), count=(/nz/) )         ,'write-jnucl_n' )
                            call check( nf90_put_var(ncids(i), general_ids(9),Jnucl_D,start=(/1,increment/), count=(/nz/) )         ,'write-jnucl_d' )
                            call check( nf90_put_var(ncids(i), general_ids(10),Jnucl_I,start=(/1,increment/), count=(/nz/) )        ,'write-jnucl_i' )
                            call check( nf90_put_var(ncids(i), general_ids(11),Jnucl_M,start=(/1,increment/), count=(/nz/) )        ,'write-jnucl_m' )
                            call check( nf90_put_var(ncids(i), general_ids(12),Jnucl_AND,start=(/1,increment/), count=(/nz/) )        ,'write-jnucl_AND' )
                            call check( nf90_put_var(ncids(i), general_ids(13),CS_H2SO4,start=(/1,increment/), count=(/nz/) )       ,'write-CS_H2SO4' )
                            call check( nf90_put_var(ncids(i), general_ids(14),conc_pp,start=(/1,1,1, increment/), count=(/nz, NSPEC_P, nr_bins,1/) ),'write-conc_pp' )
                            call check( nf90_put_var(ncids(i), general_ids(15),pH,start=(/1,1, increment/), count=(/nz, nr_bins,1/) ) ,'write-pH' )
                    end if    
                end if

        end do

        if (closefile) THEN
            if (I==1) call check( nf90_close(ncids(i)),'close') 
            if (I==2) call check( nf90_close(ncids(i)),'close') 
            if (I==3) call check( nf90_close(ncids(i)),'close') 
        end if    
        
    end subroutine open_write_nc_files



    subroutine check(status, operation)
        use netcdf
        implicit none
        integer, intent(in) :: status
        character(len=*),intent(in),optional :: operation
        
        if (status == NF90_NOERR) return
        print *, "Error encountered during ", operation 
        print *, nf90_strerror(status)
        STOP 1
    
    end subroutine check

END MODULE output     
