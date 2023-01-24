module acdc_system_1

implicit none

character(len=*), parameter :: tag = 'test'			! tag for the simulation system

logical, parameter :: small_set_mode = .true.			! small set of clusters (default mode)

integer, parameter :: nclust = 54						! number of clusters, molecules and ions
integer, parameter :: neq = 198							! number of equations
integer, parameter :: nclust_max = 54, neq_max = 198	! same as nclust and neq in the small set mode

integer, parameter :: n_variables = 3
logical, parameter :: variable_temp = .true.
logical, parameter :: variable_cs = .false.
logical, parameter :: variable_ipr = .true.
real(kind(1.d0)), parameter :: rh = 0.d0					! RH in %
real(kind(1.d0)), parameter :: pw = 0.d0					! partial pressure of water in Pa

integer, parameter :: n_monomer_types = 4
integer, parameter :: n_charges = 3			! number of charging states
integer, parameter :: n1A = 1, n1N = 3, n1B = 17, n1N1P = 35, nneg = 53, npos = 54			! cluster indices for monomers and ions

integer, parameter :: nout_all(3) = (/61, 62, 63/), nout_neu = 61, nout_neg = 62, nout_pos = 63			! indices for outgoing fluxes

integer, parameter :: nclust_coag_saved = 54			! number of saved scavenged clusters
integer, parameter :: nclust_out_saved = 80				! number of saved outgrown clusters
integer, parameter :: nrange_coag_saved(2) = (/65, 118/)	! index range of saved scavenged clusters
integer, parameter :: nrange_out_saved(2) = (/119, 198/)	! index range of saved outgrown clusters

integer, parameter :: n_mol_types = 4
integer, parameter :: nmolA = 1, nmolB = 2, nmolN = 3, nmolP = 4			! molecule indices for the used species

integer, parameter :: n_1A_clusters = 2				! number molecules and clusters containing 1 A molecule

integer, parameter :: n_neutrals = 16			! number of neutral molecules and clusters
integer, parameter :: n_negatives = 19			! negative
integer, parameter :: n_positives = 19			! positive
integer, parameter :: nclust_nogen = 52			! number of clusters and molecules excluding generic ions

integer, parameter :: n_neutral_monomers = 2, neutral_monomers(2) = (/1, 3/)			! number and indices of neutral monomers
integer, parameter :: n_negative_monomers = 1, negative_monomers(1) = (/17/)			! negative
integer, parameter :: n_positive_monomers = 1, positive_monomers(1) = (/35/)			! positive

integer, parameter :: n_neutral_clusters = 14			! number of neutral clusters
integer, parameter :: n_negative_clusters = 17			! negative
integer, parameter :: n_positive_clusters = 17			! positive

real(kind(1.d0)), parameter :: mass_max = 576.60
real(kind(1.d0)), parameter :: diameter_max = 1.07
real(kind(1.d0)), parameter :: mob_diameter_max = 1.37
integer, parameter :: ij_ind_max(4) = (/5, 1, 5, 1/)		! maximum molecular content
integer, parameter :: n_bound = 12		! number of clusters at the system boundary
integer, parameter :: n_comp_out = 80		! number of clusters growing out in different compositions

integer, parameter :: nmols_out_neutral(1, 4) = reshape((/6, 0, 5, 0/),(/1, 4/))			! criteria for outgrowing neutrals
integer, parameter :: nmols_out_negative(1, 4) = reshape((/5, 1, 3, 0/),(/1, 4/))			! negatives
integer, parameter :: nmols_out_positive(1, 4) = reshape((/5, 0, 6, 1/),(/1, 4/))			! positives


contains

subroutine n_A_in_clusters_1(n_A)
	implicit none
	integer :: n_A(54)

	n_A = (/1, 2, 0, 1, 2, 3, 2, 3, 4, 3, &
		&4, 5, 4, 5, 4, 5, 1, 2, 3, 4, &
		&5, 1, 2, 3, 4, 5, 1, 2, 3, 4, &
		&5, 4, 5, 5, 0, 1, 2, 0, 1, 2, &
		&3, 0, 1, 2, 3, 4, 2, 3, 4, 3, &
		&4, 5, 0, 0/)

end subroutine n_A_in_clusters_1

subroutine clusters_with_1_A_1(cluster_numbers)
	implicit none
	integer :: cluster_numbers(2)

	cluster_numbers = (/1, 4/)

end subroutine clusters_with_1_A_1

subroutine arrays_1(neutrals, negatives, positives)
	implicit none
	integer :: neutrals(16), negatives(19), positives(19)

	neutrals = (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16/)
	negatives = (/17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, &
		&32, 33, 34, 53/)
	positives = (/35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, &
		&50, 51, 52, 54/)

end subroutine arrays_1

subroutine cluster_arrays_1(neutral_clusters, negative_clusters, positive_clusters)
	implicit none
	integer :: neutral_clusters(14), negative_clusters(17), positive_clusters(17)

	neutral_clusters = (/2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16/)
	negative_clusters = (/18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, &
		&32, 33, 34/)
	positive_clusters = (/36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, &
		&50, 51, 52/)

end subroutine cluster_arrays_1

subroutine get_charging_state_1(charging_state)
	implicit none
	integer :: charging_state(54)

	charging_state = (/1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
		&1, 1, 1, 1, 1, 1, 2, 2, 2, 2, &
		&2, 2, 2, 2, 2, 2, 2, 2, 2, 2, &
		&2, 2, 2, 2, 3, 3, 3, 3, 3, 3, &
		&3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
		&3, 3, 2, 3/)

end subroutine get_charging_state_1

subroutine get_mass_1(mass)
	implicit none
	real(kind(1.d0)) :: mass(54)

	mass = (/98.08, 196.16, 17.04, 115.12, 213.20, 311.28, 230.24, 328.32, 426.40, 345.36, &
		&443.44, 541.52, 460.48, 558.56, 477.52, 575.60, 97.08, 195.16, 293.24, 391.32, &
		&489.40, 114.12, 212.20, 310.28, 408.36, 506.44, 131.16, 229.24, 327.32, 425.40, &
		&523.48, 442.44, 540.52, 557.56, 18.04, 116.12, 214.20, 35.08, 133.16, 231.24, &
		&329.32, 52.12, 150.20, 248.28, 346.36, 444.44, 265.32, 363.40, 461.48, 380.44, &
		&478.52, 576.60, 32.00, 19.02/)

end subroutine get_mass_1

subroutine get_diameter_1(diameter)
	implicit none

	real(kind(1.d0)) :: diameter(54)

	 diameter = (/0.55, 0.70, 0.43, 0.63, 0.75, 0.84, 0.79, 0.87, 0.94, 0.91, &
		&0.97, 1.03, 1.00, 1.05, 1.02, 1.07, 0.55, 0.70, 0.80, 0.88, &
		&0.95, 0.63, 0.75, 0.84, 0.91, 0.97, 0.69, 0.79, 0.87, 0.94, &
		&1.00, 0.97, 1.03, 1.05, 0.43, 0.63, 0.75, 0.54, 0.69, 0.79, &
		&0.87, 0.62, 0.74, 0.83, 0.91, 0.97, 0.87, 0.94, 1.00, 0.96, &
		&1.02, 1.07, 0.45, 0.39/)	! dry value

end subroutine get_diameter_1

subroutine get_mob_diameter_1(mob_diameter)
	implicit none

	real(kind(1.d0)) :: mob_diameter(54)

	 mob_diameter = (/0.85, 1.00, 0.73, 0.93, 1.05, 1.14, 1.09, 1.17, 1.24, 1.21, &
		&1.27, 1.33, 1.30, 1.35, 1.32, 1.37, 0.85, 1.00, 1.10, 1.18, &
		&1.25, 0.93, 1.05, 1.14, 1.21, 1.27, 0.99, 1.09, 1.17, 1.24, &
		&1.30, 1.27, 1.33, 1.35, 0.73, 0.93, 1.05, 0.84, 0.99, 1.09, &
		&1.17, 0.92, 1.04, 1.13, 1.21, 1.27, 1.17, 1.24, 1.30, 1.26, &
		&1.32, 1.37, 0.75, 0.69/)	! dry value

end subroutine get_mob_diameter_1

subroutine cluster_names_1(clust)
	implicit none
	character(len=11), dimension(64) :: clust

	clust(1)(:) = '1A'
	clust(2)(:) = '2A'
	clust(3)(:) = '1N'
	clust(4)(:) = '1A1N'
	clust(5)(:) = '2A1N'
	clust(6)(:) = '3A1N'
	clust(7)(:) = '2A2N'
	clust(8)(:) = '3A2N'
	clust(9)(:) = '4A2N'
	clust(10)(:) = '3A3N'
	clust(11)(:) = '4A3N'
	clust(12)(:) = '5A3N'
	clust(13)(:) = '4A4N'
	clust(14)(:) = '5A4N'
	clust(15)(:) = '4A5N'
	clust(16)(:) = '5A5N'
	clust(17)(:) = '1B'
	clust(18)(:) = '1A1B'
	clust(19)(:) = '2A1B'
	clust(20)(:) = '3A1B'
	clust(21)(:) = '4A1B'
	clust(22)(:) = '1B1N'
	clust(23)(:) = '1A1B1N'
	clust(24)(:) = '2A1B1N'
	clust(25)(:) = '3A1B1N'
	clust(26)(:) = '4A1B1N'
	clust(27)(:) = '1B2N'
	clust(28)(:) = '1A1B2N'
	clust(29)(:) = '2A1B2N'
	clust(30)(:) = '3A1B2N'
	clust(31)(:) = '4A1B2N'
	clust(32)(:) = '3A1B3N'
	clust(33)(:) = '4A1B3N'
	clust(34)(:) = '4A1B4N'
	clust(35)(:) = '1N1P'
	clust(36)(:) = '1A1N1P'
	clust(37)(:) = '2A1N1P'
	clust(38)(:) = '2N1P'
	clust(39)(:) = '1A2N1P'
	clust(40)(:) = '2A2N1P'
	clust(41)(:) = '3A2N1P'
	clust(42)(:) = '3N1P'
	clust(43)(:) = '1A3N1P'
	clust(44)(:) = '2A3N1P'
	clust(45)(:) = '3A3N1P'
	clust(46)(:) = '4A3N1P'
	clust(47)(:) = '2A4N1P'
	clust(48)(:) = '3A4N1P'
	clust(49)(:) = '4A4N1P'
	clust(50)(:) = '3A5N1P'
	clust(51)(:) = '4A5N1P'
	clust(52)(:) = '5A5N1P'
	clust(53)(:) = 'neg'
	clust(54)(:) = 'pos'
	clust(55)(:) = 'source'
	clust(56)(:) = 'coag'
	clust(57)(:) = 'wall'
	clust(58)(:) = 'dil'
	clust(59)(:) = 'insink'
	clust(60)(:) = 'rec'
	clust(61)(:) = 'out_neu'
	clust(62)(:) = 'out_neg'
	clust(63)(:) = 'out_pos'
	clust(64)(:) = 'bound'

	! Additional flux elements not listed here:
	! Elements from 65 to 118: coagulation loss of each cluster
	! Elements from 119 to 198: each outgrown cluster composition

end subroutine cluster_names_1

subroutine monomer_names_1(clust_mon)
	implicit none
	character(len=11), dimension(4) :: clust_mon

	clust_mon(1)(:) = '1A'
	clust_mon(2)(:) = '1N'
	clust_mon(3)(:) = '1B'
	clust_mon(4)(:) = '1N1P'

end subroutine monomer_names_1

subroutine molecule_names_1(labels)
	implicit none
	character(len=11), dimension(4) :: labels

	labels(1)(:) = 'A'
	labels(2)(:) = 'B'
	labels(3)(:) = 'N'
	labels(4)(:) = 'P'

end subroutine molecule_names_1

subroutine monomer_indices_1(n_monomers)
	implicit none
	integer :: n_monomers(4)

	n_monomers = (/1, 3, 17, 35/)

end subroutine monomer_indices_1

subroutine get_bound_1(bound_clusters,nmols_bound)
	implicit none
	integer :: bound_clusters(12), nmols_bound(12,4)

	nmols_bound(1,:) = (/5, 0, 3, 0/)
	nmols_bound(2,:) = (/5, 0, 4, 0/)
	nmols_bound(3,:) = (/4, 0, 5, 0/)
	nmols_bound(4,:) = (/5, 0, 5, 0/)
	nmols_bound(5,:) = (/4, 1, 0, 0/)
	nmols_bound(6,:) = (/4, 1, 1, 0/)
	nmols_bound(7,:) = (/4, 1, 2, 0/)
	nmols_bound(8,:) = (/4, 1, 3, 0/)
	nmols_bound(9,:) = (/4, 1, 4, 0/)
	nmols_bound(10,:) = (/3, 0, 5, 1/)
	nmols_bound(11,:) = (/4, 0, 5, 1/)
	nmols_bound(12,:) = (/5, 0, 5, 1/)

	bound_clusters = (/12, 14, 15, 16, 21, 26, 31, 33, 34, 50, 51, 52/)

end subroutine get_bound_1

subroutine get_diameter_out_1(diameter_out)
	implicit none
	real(kind(1.d0)) :: diameter_out(80)

	diameter_out = (/1.12, 1.08, 1.10, 1.16, 1.12, 1.14, 1.10, 1.14, 1.12, 1.14, &
		&1.18, 1.16, 1.18, 1.20, 1.22, 1.17, 1.19, 1.20, 1.22, 1.16, &
		&1.20, 1.18, 1.12, 1.16, 1.20, 1.24, 1.22, 1.24, 1.24, 1.26, &
		&1.28, 1.21, 1.23, 1.24, 1.26, 1.28, 1.22, 1.26, 1.24, 1.18, &
		&1.22, 1.26, 1.29, 1.27, 1.29, 1.29, 1.31, 1.32, 1.25, 1.26, &
		&1.28, 1.29, 1.31, 1.26, 1.32, 1.27, 1.31, 1.14, 1.29, 1.24, &
		&1.27, 1.31, 1.34, 1.32, 1.34, 1.29, 1.32, 1.16, 1.20, 1.26, &
		&1.31, 1.14, 1.20, 1.25, 1.29, 1.32, 1.35, 1.34, 1.35, 1.28/)	! dry value

end subroutine get_diameter_out_1

subroutine get_nmols_out_1(nmols_out)
	implicit none
	integer :: nmols_out(80,4)

	nmols_out(1,:) = (/6, 0, 5, 0/)
	nmols_out(2,:) = (/5, 1, 3, 0/)
	nmols_out(3,:) = (/5, 1, 4, 0/)
	nmols_out(4,:) = (/7, 0, 5, 0/)
	nmols_out(5,:) = (/6, 1, 3, 0/)
	nmols_out(6,:) = (/6, 1, 4, 0/)
	nmols_out(7,:) = (/5, 0, 6, 1/)
	nmols_out(8,:) = (/6, 0, 6, 0/)
	nmols_out(9,:) = (/5, 1, 5, 0/)
	nmols_out(10,:) = (/6, 0, 6, 1/)
	nmols_out(11,:) = (/7, 0, 6, 0/)
	nmols_out(12,:) = (/6, 1, 5, 0/)
	nmols_out(13,:) = (/7, 0, 6, 1/)
	nmols_out(14,:) = (/8, 0, 5, 0/)
	nmols_out(15,:) = (/8, 0, 6, 0/)
	nmols_out(16,:) = (/7, 1, 3, 0/)
	nmols_out(17,:) = (/7, 1, 4, 0/)
	nmols_out(18,:) = (/7, 1, 5, 0/)
	nmols_out(19,:) = (/8, 0, 6, 1/)
	nmols_out(20,:) = (/6, 0, 7, 0/)
	nmols_out(21,:) = (/7, 0, 7, 0/)
	nmols_out(22,:) = (/6, 1, 6, 0/)
	nmols_out(23,:) = (/5, 0, 7, 1/)
	nmols_out(24,:) = (/6, 0, 7, 1/)
	nmols_out(25,:) = (/7, 0, 7, 1/)
	nmols_out(26,:) = (/8, 0, 7, 0/)
	nmols_out(27,:) = (/7, 1, 6, 0/)
	nmols_out(28,:) = (/8, 0, 7, 1/)
	nmols_out(29,:) = (/9, 0, 5, 0/)
	nmols_out(30,:) = (/9, 0, 6, 0/)
	nmols_out(31,:) = (/9, 0, 7, 0/)
	nmols_out(32,:) = (/8, 1, 3, 0/)
	nmols_out(33,:) = (/8, 1, 4, 0/)
	nmols_out(34,:) = (/8, 1, 5, 0/)
	nmols_out(35,:) = (/8, 1, 6, 0/)
	nmols_out(36,:) = (/9, 0, 7, 1/)
	nmols_out(37,:) = (/7, 0, 8, 0/)
	nmols_out(38,:) = (/8, 0, 8, 0/)
	nmols_out(39,:) = (/7, 1, 7, 0/)
	nmols_out(40,:) = (/6, 0, 8, 1/)
	nmols_out(41,:) = (/7, 0, 8, 1/)
	nmols_out(42,:) = (/8, 0, 8, 1/)
	nmols_out(43,:) = (/9, 0, 8, 0/)
	nmols_out(44,:) = (/8, 1, 7, 0/)
	nmols_out(45,:) = (/9, 0, 8, 1/)
	nmols_out(46,:) = (/10, 0, 6, 0/)
	nmols_out(47,:) = (/10, 0, 7, 0/)
	nmols_out(48,:) = (/10, 0, 8, 0/)
	nmols_out(49,:) = (/9, 1, 3, 0/)
	nmols_out(50,:) = (/9, 1, 4, 0/)
	nmols_out(51,:) = (/9, 1, 5, 0/)
	nmols_out(52,:) = (/9, 1, 6, 0/)
	nmols_out(53,:) = (/9, 1, 7, 0/)
	nmols_out(54,:) = (/9, 0, 6, 1/)
	nmols_out(55,:) = (/10, 0, 8, 1/)
	nmols_out(56,:) = (/8, 0, 9, 0/)
	nmols_out(57,:) = (/9, 0, 9, 0/)
	nmols_out(58,:) = (/5, 1, 6, 0/)
	nmols_out(59,:) = (/8, 1, 8, 0/)
	nmols_out(60,:) = (/7, 0, 9, 1/)
	nmols_out(61,:) = (/8, 0, 9, 1/)
	nmols_out(62,:) = (/9, 0, 9, 1/)
	nmols_out(63,:) = (/10, 0, 9, 0/)
	nmols_out(64,:) = (/9, 1, 8, 0/)
	nmols_out(65,:) = (/10, 0, 9, 1/)
	nmols_out(66,:) = (/8, 0, 10, 0/)
	nmols_out(67,:) = (/9, 0, 10, 0/)
	nmols_out(68,:) = (/5, 1, 7, 0/)
	nmols_out(69,:) = (/6, 1, 7, 0/)
	nmols_out(70,:) = (/7, 1, 8, 0/)
	nmols_out(71,:) = (/8, 1, 9, 0/)
	nmols_out(72,:) = (/5, 0, 8, 1/)
	nmols_out(73,:) = (/6, 0, 9, 1/)
	nmols_out(74,:) = (/7, 0, 10, 1/)
	nmols_out(75,:) = (/8, 0, 10, 1/)
	nmols_out(76,:) = (/9, 0, 10, 1/)
	nmols_out(77,:) = (/10, 0, 10, 0/)
	nmols_out(78,:) = (/9, 1, 9, 0/)
	nmols_out(79,:) = (/10, 0, 10, 1/)
	nmols_out(80,:) = (/10, 0, 5, 0/)

end subroutine get_nmols_out_1

subroutine molecule_names_nocharge_1(labels_nocharge)
	implicit none
	character(len=11), dimension(2) :: labels_nocharge

	labels_nocharge(1)(:) = 'A'
	labels_nocharge(2)(:) = 'N'

end subroutine molecule_names_nocharge_1

subroutine get_nmols_nocharge_1(nmols_nocharge_clust)
	implicit none
	integer :: nmols_nocharge_clust(54,2)

	nmols_nocharge_clust = 0

	nmols_nocharge_clust(1,:) = (/1, 0/)
	nmols_nocharge_clust(2,:) = (/2, 0/)
	nmols_nocharge_clust(3,:) = (/0, 1/)
	nmols_nocharge_clust(4,:) = (/1, 1/)
	nmols_nocharge_clust(5,:) = (/2, 1/)
	nmols_nocharge_clust(6,:) = (/3, 1/)
	nmols_nocharge_clust(7,:) = (/2, 2/)
	nmols_nocharge_clust(8,:) = (/3, 2/)
	nmols_nocharge_clust(9,:) = (/4, 2/)
	nmols_nocharge_clust(10,:) = (/3, 3/)
	nmols_nocharge_clust(11,:) = (/4, 3/)
	nmols_nocharge_clust(12,:) = (/5, 3/)
	nmols_nocharge_clust(13,:) = (/4, 4/)
	nmols_nocharge_clust(14,:) = (/5, 4/)
	nmols_nocharge_clust(15,:) = (/4, 5/)
	nmols_nocharge_clust(16,:) = (/5, 5/)
	nmols_nocharge_clust(17,:) = (/1, 0/)
	nmols_nocharge_clust(18,:) = (/2, 0/)
	nmols_nocharge_clust(19,:) = (/3, 0/)
	nmols_nocharge_clust(20,:) = (/4, 0/)
	nmols_nocharge_clust(21,:) = (/5, 0/)
	nmols_nocharge_clust(22,:) = (/1, 1/)
	nmols_nocharge_clust(23,:) = (/2, 1/)
	nmols_nocharge_clust(24,:) = (/3, 1/)
	nmols_nocharge_clust(25,:) = (/4, 1/)
	nmols_nocharge_clust(26,:) = (/5, 1/)
	nmols_nocharge_clust(27,:) = (/1, 2/)
	nmols_nocharge_clust(28,:) = (/2, 2/)
	nmols_nocharge_clust(29,:) = (/3, 2/)
	nmols_nocharge_clust(30,:) = (/4, 2/)
	nmols_nocharge_clust(31,:) = (/5, 2/)
	nmols_nocharge_clust(32,:) = (/4, 3/)
	nmols_nocharge_clust(33,:) = (/5, 3/)
	nmols_nocharge_clust(34,:) = (/5, 4/)
	nmols_nocharge_clust(35,:) = (/0, 1/)
	nmols_nocharge_clust(36,:) = (/1, 1/)
	nmols_nocharge_clust(37,:) = (/2, 1/)
	nmols_nocharge_clust(38,:) = (/0, 2/)
	nmols_nocharge_clust(39,:) = (/1, 2/)
	nmols_nocharge_clust(40,:) = (/2, 2/)
	nmols_nocharge_clust(41,:) = (/3, 2/)
	nmols_nocharge_clust(42,:) = (/0, 3/)
	nmols_nocharge_clust(43,:) = (/1, 3/)
	nmols_nocharge_clust(44,:) = (/2, 3/)
	nmols_nocharge_clust(45,:) = (/3, 3/)
	nmols_nocharge_clust(46,:) = (/4, 3/)
	nmols_nocharge_clust(47,:) = (/2, 4/)
	nmols_nocharge_clust(48,:) = (/3, 4/)
	nmols_nocharge_clust(49,:) = (/4, 4/)
	nmols_nocharge_clust(50,:) = (/3, 5/)
	nmols_nocharge_clust(51,:) = (/4, 5/)
	nmols_nocharge_clust(52,:) = (/5, 5/)

end subroutine get_nmols_nocharge_1

subroutine get_nmols_nocharge_out_1(nmols_nocharge_out)
	implicit none
	integer :: nmols_nocharge_out(80,2)

	nmols_nocharge_out = 0

	nmols_nocharge_out(1,:) = (/6, 5/)
	nmols_nocharge_out(2,:) = (/6, 3/)
	nmols_nocharge_out(3,:) = (/6, 4/)
	nmols_nocharge_out(4,:) = (/7, 5/)
	nmols_nocharge_out(5,:) = (/7, 3/)
	nmols_nocharge_out(6,:) = (/7, 4/)
	nmols_nocharge_out(7,:) = (/5, 6/)
	nmols_nocharge_out(8,:) = (/6, 6/)
	nmols_nocharge_out(9,:) = (/6, 5/)
	nmols_nocharge_out(10,:) = (/6, 6/)
	nmols_nocharge_out(11,:) = (/7, 6/)
	nmols_nocharge_out(12,:) = (/7, 5/)
	nmols_nocharge_out(13,:) = (/7, 6/)
	nmols_nocharge_out(14,:) = (/8, 5/)
	nmols_nocharge_out(15,:) = (/8, 6/)
	nmols_nocharge_out(16,:) = (/8, 3/)
	nmols_nocharge_out(17,:) = (/8, 4/)
	nmols_nocharge_out(18,:) = (/8, 5/)
	nmols_nocharge_out(19,:) = (/8, 6/)
	nmols_nocharge_out(20,:) = (/6, 7/)
	nmols_nocharge_out(21,:) = (/7, 7/)
	nmols_nocharge_out(22,:) = (/7, 6/)
	nmols_nocharge_out(23,:) = (/5, 7/)
	nmols_nocharge_out(24,:) = (/6, 7/)
	nmols_nocharge_out(25,:) = (/7, 7/)
	nmols_nocharge_out(26,:) = (/8, 7/)
	nmols_nocharge_out(27,:) = (/8, 6/)
	nmols_nocharge_out(28,:) = (/8, 7/)
	nmols_nocharge_out(29,:) = (/9, 5/)
	nmols_nocharge_out(30,:) = (/9, 6/)
	nmols_nocharge_out(31,:) = (/9, 7/)
	nmols_nocharge_out(32,:) = (/9, 3/)
	nmols_nocharge_out(33,:) = (/9, 4/)
	nmols_nocharge_out(34,:) = (/9, 5/)
	nmols_nocharge_out(35,:) = (/9, 6/)
	nmols_nocharge_out(36,:) = (/9, 7/)
	nmols_nocharge_out(37,:) = (/7, 8/)
	nmols_nocharge_out(38,:) = (/8, 8/)
	nmols_nocharge_out(39,:) = (/8, 7/)
	nmols_nocharge_out(40,:) = (/6, 8/)
	nmols_nocharge_out(41,:) = (/7, 8/)
	nmols_nocharge_out(42,:) = (/8, 8/)
	nmols_nocharge_out(43,:) = (/9, 8/)
	nmols_nocharge_out(44,:) = (/9, 7/)
	nmols_nocharge_out(45,:) = (/9, 8/)
	nmols_nocharge_out(46,:) = (/10, 6/)
	nmols_nocharge_out(47,:) = (/10, 7/)
	nmols_nocharge_out(48,:) = (/10, 8/)
	nmols_nocharge_out(49,:) = (/10, 3/)
	nmols_nocharge_out(50,:) = (/10, 4/)
	nmols_nocharge_out(51,:) = (/10, 5/)
	nmols_nocharge_out(52,:) = (/10, 6/)
	nmols_nocharge_out(53,:) = (/10, 7/)
	nmols_nocharge_out(54,:) = (/9, 6/)
	nmols_nocharge_out(55,:) = (/10, 8/)
	nmols_nocharge_out(56,:) = (/8, 9/)
	nmols_nocharge_out(57,:) = (/9, 9/)
	nmols_nocharge_out(58,:) = (/6, 6/)
	nmols_nocharge_out(59,:) = (/9, 8/)
	nmols_nocharge_out(60,:) = (/7, 9/)
	nmols_nocharge_out(61,:) = (/8, 9/)
	nmols_nocharge_out(62,:) = (/9, 9/)
	nmols_nocharge_out(63,:) = (/10, 9/)
	nmols_nocharge_out(64,:) = (/10, 8/)
	nmols_nocharge_out(65,:) = (/10, 9/)
	nmols_nocharge_out(66,:) = (/8, 10/)
	nmols_nocharge_out(67,:) = (/9, 10/)
	nmols_nocharge_out(68,:) = (/6, 7/)
	nmols_nocharge_out(69,:) = (/7, 7/)
	nmols_nocharge_out(70,:) = (/8, 8/)
	nmols_nocharge_out(71,:) = (/9, 9/)
	nmols_nocharge_out(72,:) = (/5, 8/)
	nmols_nocharge_out(73,:) = (/6, 9/)
	nmols_nocharge_out(74,:) = (/7, 10/)
	nmols_nocharge_out(75,:) = (/8, 10/)
	nmols_nocharge_out(76,:) = (/9, 10/)
	nmols_nocharge_out(77,:) = (/10, 10/)
	nmols_nocharge_out(78,:) = (/10, 9/)
	nmols_nocharge_out(79,:) = (/10, 10/)
	nmols_nocharge_out(80,:) = (/10, 5/)

end subroutine get_nmols_nocharge_out_1


end module acdc_system_1

