module acdc_system_3

implicit none

character(len=*), parameter :: tag = 'test'			! tag for the simulation system

logical, parameter :: small_set_mode = .true.			! small set of clusters (default mode)

integer, parameter :: nclust = 24						! number of clusters, molecules and ions
integer, parameter :: neq = 80							! number of equations
integer, parameter :: nclust_max = 24, neq_max = 80	! same as nclust and neq in the small set mode

integer, parameter :: n_variables = 1
logical, parameter :: variable_temp = .true.
logical, parameter :: variable_cs = .false.
logical, parameter :: variable_ipr = .false.
real(kind(1.d0)), parameter :: rh = 0.d0					! RH in %
real(kind(1.d0)), parameter :: pw = 0.d0					! partial pressure of water in Pa

integer, parameter :: n_monomer_types = 2
integer, parameter :: n_charges = 1			! number of charging states
integer, parameter :: n1Ii = 1, n1Io = 5			! cluster indices for monomers and ions

integer, parameter :: nout_all(1) = (/31/), nout_neu = 31			! indices for outgoing fluxes

integer, parameter :: nclust_coag_saved = 24			! number of saved scavenged clusters
integer, parameter :: nclust_out_saved = 24				! number of saved outgrown clusters
integer, parameter :: nrange_coag_saved(2) = (/33, 56/)	! index range of saved scavenged clusters
integer, parameter :: nrange_out_saved(2) = (/57, 80/)	! index range of saved outgrown clusters

integer, parameter :: n_mol_types = 2
integer, parameter :: nmolIi = 1, nmolIo = 2			! molecule indices for the used species

integer, parameter :: n_1Ii_clusters = 5				! number molecules and clusters containing 1 Ii molecule

integer, parameter :: n_neutrals = 24			! number of neutral molecules and clusters
integer, parameter :: n_negatives = 0			! negative
integer, parameter :: n_positives = 0			! positive
integer, parameter :: nclust_nogen = 24			! number of clusters and molecules excluding generic ions

integer, parameter :: n_neutral_monomers = 2, neutral_monomers(2) = (/1, 5/)			! number and indices of neutral monomers
integer, parameter :: n_negative_monomers = 0			! negative
integer, parameter :: n_positive_monomers = 0			! positive

integer, parameter :: n_neutral_clusters = 22			! number of neutral clusters
integer, parameter :: n_negative_clusters = 0			! negative
integer, parameter :: n_positive_clusters = 0			! positive

real(kind(1.d0)), parameter :: mass_max = 1343.2800
real(kind(1.d0)), parameter :: diameter_max = 0.9726
real(kind(1.d0)), parameter :: mob_diameter_max = 1.2726
integer, parameter :: ij_ind_max(2) = (/4, 4/)		! maximum molecular content
integer, parameter :: n_bound = 9		! number of clusters at the system boundary
integer, parameter :: n_comp_out = 24		! number of clusters growing out in different compositions

integer, parameter :: nmols_out_neutral(2, 2) = reshape((/5, 4, 4, 5/),(/2, 2/))			! criteria for outgrowing neutrals


contains

subroutine n_Ii_in_clusters_3(n_Ii)
	implicit none
	integer :: n_Ii(24)

	n_Ii = (/1, 2, 3, 4, 0, 1, 2, 3, 4, 0, &
		&1, 2, 3, 4, 0, 1, 2, 3, 4, 0, &
		&1, 2, 3, 4/)

end subroutine n_Ii_in_clusters_3

subroutine clusters_with_1_Ii_3(cluster_numbers)
	implicit none
	integer :: cluster_numbers(5)

	cluster_numbers = (/1, 6, 11, 16, 21/)

end subroutine clusters_with_1_Ii_3

subroutine cluster_arrays_3(neutral_clusters)
	implicit none
	integer :: neutral_clusters(22)

	neutral_clusters = (/2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, &
		&21, 22, 23, 24/)

end subroutine cluster_arrays_3

subroutine get_charging_state_3(charging_state)
	implicit none
	integer :: charging_state(24)

	charging_state = (/1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
		&1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
		&1, 1, 1, 1/)

end subroutine get_charging_state_3

subroutine get_mass_3(mass)
	implicit none
	real(kind(1.d0)) :: mass(24)

	mass = (/175.9100, 351.8200, 527.7300, 703.6400, 159.9100, 335.8200, 511.7300, 687.6400, 863.5500, 319.8200, &
		&495.7300, 671.6400, 847.5500, 1023.4600, 479.7300, 655.6400, 831.5500, 1007.4600, 1183.3700, 639.6400, &
		&815.5500, 991.4600, 1167.3700, 1343.2800/)

end subroutine get_mass_3

subroutine get_diameter_3(diameter)
	implicit none

	real(kind(1.d0)) :: diameter(24)

	 diameter = (/0.4939, 0.6223, 0.7124, 0.7840, 0.4785, 0.6127, 0.7051, 0.7781, 0.8394, 0.6028, &
		&0.6977, 0.7720, 0.8342, 0.8883, 0.6901, 0.7658, 0.8289, 0.8837, 0.9324, 0.7595, &
		&0.8236, 0.8790, 0.9282, 0.9726/)	! dry value

end subroutine get_diameter_3

subroutine get_mob_diameter_3(mob_diameter)
	implicit none

	real(kind(1.d0)) :: mob_diameter(24)

	 mob_diameter = (/0.7939, 0.9223, 1.0124, 1.0840, 0.7785, 0.9127, 1.0051, 1.0781, 1.1394, 0.9028, &
		&0.9977, 1.0720, 1.1342, 1.1883, 0.9901, 1.0658, 1.1289, 1.1837, 1.2324, 1.0595, &
		&1.1236, 1.1790, 1.2282, 1.2726/)	! dry value

end subroutine get_mob_diameter_3

subroutine cluster_names_3(clust)
	implicit none
	character(len=11), dimension(32) :: clust

	clust(1)(:) = '1Ii'
	clust(2)(:) = '2Ii'
	clust(3)(:) = '3Ii'
	clust(4)(:) = '4Ii'
	clust(5)(:) = '1Io'
	clust(6)(:) = '1Ii1Io'
	clust(7)(:) = '2Ii1Io'
	clust(8)(:) = '3Ii1Io'
	clust(9)(:) = '4Ii1Io'
	clust(10)(:) = '2Io'
	clust(11)(:) = '1Ii2Io'
	clust(12)(:) = '2Ii2Io'
	clust(13)(:) = '3Ii2Io'
	clust(14)(:) = '4Ii2Io'
	clust(15)(:) = '3Io'
	clust(16)(:) = '1Ii3Io'
	clust(17)(:) = '2Ii3Io'
	clust(18)(:) = '3Ii3Io'
	clust(19)(:) = '4Ii3Io'
	clust(20)(:) = '4Io'
	clust(21)(:) = '1Ii4Io'
	clust(22)(:) = '2Ii4Io'
	clust(23)(:) = '3Ii4Io'
	clust(24)(:) = '4Ii4Io'
	clust(25)(:) = 'source'
	clust(26)(:) = 'coag'
	clust(27)(:) = 'wall'
	clust(28)(:) = 'dil'
	clust(29)(:) = 'insink'
	clust(30)(:) = 'rec'
	clust(31)(:) = 'out_neu'
	clust(32)(:) = 'bound'

	! Additional flux elements not listed here:
	! Elements from 33 to 56: coagulation loss of each cluster
	! Elements from 57 to 80: each outgrown cluster composition

end subroutine cluster_names_3

subroutine monomer_names_3(clust_mon)
	implicit none
	character(len=11), dimension(2) :: clust_mon

	clust_mon(1)(:) = '1Ii'
	clust_mon(2)(:) = '1Io'

end subroutine monomer_names_3

subroutine molecule_names_3(labels)
	implicit none
	character(len=11), dimension(2) :: labels

	labels(1)(:) = 'Ii'
	labels(2)(:) = 'Io'

end subroutine molecule_names_3

subroutine monomer_indices_3(n_monomers)
	implicit none
	integer :: n_monomers(2)

	n_monomers = (/1, 5/)

end subroutine monomer_indices_3

subroutine get_bound_3(bound_clusters,nmols_bound)
	implicit none
	integer :: bound_clusters(9), nmols_bound(9,2)

	nmols_bound(1,:) = (/4, 0/)
	nmols_bound(2,:) = (/4, 1/)
	nmols_bound(3,:) = (/4, 2/)
	nmols_bound(4,:) = (/4, 3/)
	nmols_bound(5,:) = (/0, 4/)
	nmols_bound(6,:) = (/1, 4/)
	nmols_bound(7,:) = (/2, 4/)
	nmols_bound(8,:) = (/3, 4/)
	nmols_bound(9,:) = (/4, 4/)

	bound_clusters = (/4, 9, 14, 19, 20, 21, 22, 23, 24/)

end subroutine get_bound_3

subroutine get_diameter_out_3(diameter_out)
	implicit none
	real(kind(1.d0)) :: diameter_out(24)

	diameter_out = (/1.0134, 1.0510, 1.0862, 1.1192, 1.0098, 1.0477, 1.0831, 1.1163, 1.1477, 1.0444, &
		&1.0800, 1.1134, 1.1449, 1.1747, 1.0768, 1.1104, 1.1421, 1.1721, 1.2006, 1.1075, &
		&1.1393, 1.1694, 1.1981, 1.2254/)	! dry value

end subroutine get_diameter_out_3

subroutine get_nmols_out_3(nmols_out)
	implicit none
	integer :: nmols_out(24,2)

	nmols_out(1,:) = (/5, 4/)
	nmols_out(2,:) = (/6, 4/)
	nmols_out(3,:) = (/7, 4/)
	nmols_out(4,:) = (/8, 4/)
	nmols_out(5,:) = (/4, 5/)
	nmols_out(6,:) = (/5, 5/)
	nmols_out(7,:) = (/6, 5/)
	nmols_out(8,:) = (/7, 5/)
	nmols_out(9,:) = (/8, 5/)
	nmols_out(10,:) = (/4, 6/)
	nmols_out(11,:) = (/5, 6/)
	nmols_out(12,:) = (/6, 6/)
	nmols_out(13,:) = (/7, 6/)
	nmols_out(14,:) = (/8, 6/)
	nmols_out(15,:) = (/4, 7/)
	nmols_out(16,:) = (/5, 7/)
	nmols_out(17,:) = (/6, 7/)
	nmols_out(18,:) = (/7, 7/)
	nmols_out(19,:) = (/8, 7/)
	nmols_out(20,:) = (/4, 8/)
	nmols_out(21,:) = (/5, 8/)
	nmols_out(22,:) = (/6, 8/)
	nmols_out(23,:) = (/7, 8/)
	nmols_out(24,:) = (/8, 8/)

end subroutine get_nmols_out_3

subroutine molecule_names_nocharge_3(labels_nocharge)
	implicit none
	character(len=11), dimension(2) :: labels_nocharge

	labels_nocharge(1)(:) = 'Ii'
	labels_nocharge(2)(:) = 'Io'

end subroutine molecule_names_nocharge_3

subroutine get_nmols_nocharge_3(nmols_nocharge_clust)
	implicit none
	integer :: nmols_nocharge_clust(24,2)

	nmols_nocharge_clust = 0

	nmols_nocharge_clust(1,:) = (/1, 0/)
	nmols_nocharge_clust(2,:) = (/2, 0/)
	nmols_nocharge_clust(3,:) = (/3, 0/)
	nmols_nocharge_clust(4,:) = (/4, 0/)
	nmols_nocharge_clust(5,:) = (/0, 1/)
	nmols_nocharge_clust(6,:) = (/1, 1/)
	nmols_nocharge_clust(7,:) = (/2, 1/)
	nmols_nocharge_clust(8,:) = (/3, 1/)
	nmols_nocharge_clust(9,:) = (/4, 1/)
	nmols_nocharge_clust(10,:) = (/0, 2/)
	nmols_nocharge_clust(11,:) = (/1, 2/)
	nmols_nocharge_clust(12,:) = (/2, 2/)
	nmols_nocharge_clust(13,:) = (/3, 2/)
	nmols_nocharge_clust(14,:) = (/4, 2/)
	nmols_nocharge_clust(15,:) = (/0, 3/)
	nmols_nocharge_clust(16,:) = (/1, 3/)
	nmols_nocharge_clust(17,:) = (/2, 3/)
	nmols_nocharge_clust(18,:) = (/3, 3/)
	nmols_nocharge_clust(19,:) = (/4, 3/)
	nmols_nocharge_clust(20,:) = (/0, 4/)
	nmols_nocharge_clust(21,:) = (/1, 4/)
	nmols_nocharge_clust(22,:) = (/2, 4/)
	nmols_nocharge_clust(23,:) = (/3, 4/)
	nmols_nocharge_clust(24,:) = (/4, 4/)

end subroutine get_nmols_nocharge_3

subroutine get_nmols_nocharge_out_3(nmols_nocharge_out)
	implicit none
	integer :: nmols_nocharge_out(24,2)

	nmols_nocharge_out = 0

	nmols_nocharge_out(1,:) = (/5, 4/)
	nmols_nocharge_out(2,:) = (/6, 4/)
	nmols_nocharge_out(3,:) = (/7, 4/)
	nmols_nocharge_out(4,:) = (/8, 4/)
	nmols_nocharge_out(5,:) = (/4, 5/)
	nmols_nocharge_out(6,:) = (/5, 5/)
	nmols_nocharge_out(7,:) = (/6, 5/)
	nmols_nocharge_out(8,:) = (/7, 5/)
	nmols_nocharge_out(9,:) = (/8, 5/)
	nmols_nocharge_out(10,:) = (/4, 6/)
	nmols_nocharge_out(11,:) = (/5, 6/)
	nmols_nocharge_out(12,:) = (/6, 6/)
	nmols_nocharge_out(13,:) = (/7, 6/)
	nmols_nocharge_out(14,:) = (/8, 6/)
	nmols_nocharge_out(15,:) = (/4, 7/)
	nmols_nocharge_out(16,:) = (/5, 7/)
	nmols_nocharge_out(17,:) = (/6, 7/)
	nmols_nocharge_out(18,:) = (/7, 7/)
	nmols_nocharge_out(19,:) = (/8, 7/)
	nmols_nocharge_out(20,:) = (/4, 8/)
	nmols_nocharge_out(21,:) = (/5, 8/)
	nmols_nocharge_out(22,:) = (/6, 8/)
	nmols_nocharge_out(23,:) = (/7, 8/)
	nmols_nocharge_out(24,:) = (/8, 8/)

end subroutine get_nmols_nocharge_out_3


end module acdc_system_3

