module acdc_system_2

implicit none

character(len=*), parameter :: tag = 'test'			! tag for the simulation system

logical, parameter :: small_set_mode = .true.			! small set of clusters (default mode)

integer, parameter :: nclust = 55						! number of clusters, molecules and ions
integer, parameter :: neq = 196							! number of equations
integer, parameter :: nclust_max = 55, neq_max = 196	! same as nclust and neq in the small set mode

integer, parameter :: n_variables = 3
logical, parameter :: variable_temp = .true.
logical, parameter :: variable_cs = .false.
logical, parameter :: variable_ipr = .true.
real(kind(1.d0)), parameter :: rh = 0.d0					! RH in %
real(kind(1.d0)), parameter :: pw = 0.d0					! partial pressure of water in Pa

integer, parameter :: n_monomer_types = 4
integer, parameter :: n_charges = 3			! number of charging states
integer, parameter :: n1A = 1, n1D = 5, n1B = 25, n1D1P = 40, nneg = 54, npos = 55			! cluster indices for monomers and ions

integer, parameter :: nout_all(3) = (/62, 63, 64/), nout_neu = 62, nout_neg = 63, nout_pos = 64			! indices for outgoing fluxes

integer, parameter :: nclust_coag_saved = 55			! number of saved scavenged clusters
integer, parameter :: nclust_out_saved = 76				! number of saved outgrown clusters
integer, parameter :: nrange_coag_saved(2) = (/66, 120/)	! index range of saved scavenged clusters
integer, parameter :: nrange_out_saved(2) = (/121, 196/)	! index range of saved outgrown clusters

integer, parameter :: n_mol_types = 4
integer, parameter :: nmolA = 1, nmolB = 2, nmolD = 3, nmolP = 4			! molecule indices for the used species

integer, parameter :: n_1A_clusters = 5				! number molecules and clusters containing 1 A molecule

integer, parameter :: n_neutrals = 24			! number of neutral molecules and clusters
integer, parameter :: n_negatives = 16			! negative
integer, parameter :: n_positives = 15			! positive
integer, parameter :: nclust_nogen = 53			! number of clusters and molecules excluding generic ions

integer, parameter :: n_neutral_monomers = 2, neutral_monomers(2) = (/1, 5/)			! number and indices of neutral monomers
integer, parameter :: n_negative_monomers = 1, negative_monomers(1) = (/25/)			! negative
integer, parameter :: n_positive_monomers = 1, positive_monomers(1) = (/40/)			! positive

integer, parameter :: n_neutral_clusters = 22			! number of neutral clusters
integer, parameter :: n_negative_clusters = 14			! negative
integer, parameter :: n_positive_clusters = 13			! positive

real(kind(1.d0)), parameter :: mass_max = 573.64
real(kind(1.d0)), parameter :: diameter_max = 1.15
real(kind(1.d0)), parameter :: mob_diameter_max = 1.45
integer, parameter :: ij_ind_max(4) = (/4, 1, 4, 1/)		! maximum molecular content
integer, parameter :: n_bound = 18		! number of clusters at the system boundary
integer, parameter :: n_comp_out = 76		! number of clusters growing out in different compositions

integer, parameter :: nmols_out_neutral(2, 4) = reshape((/5, 4, 0, 0, 4, 5, 0, 0/),(/2, 4/))			! criteria for outgrowing neutrals
integer, parameter :: nmols_out_negative(1, 4) = reshape((/4, 1, 2, 0/),(/1, 4/))			! negatives
integer, parameter :: nmols_out_positive(1, 4) = reshape((/3, 0, 5, 1/),(/1, 4/))			! positives


contains

subroutine n_A_in_clusters_2(n_A)
	implicit none
	integer :: n_A(55)

	n_A = (/1, 2, 3, 4, 0, 1, 2, 3, 4, 0, &
		&1, 2, 3, 4, 0, 1, 2, 3, 4, 0, &
		&1, 2, 3, 4, 1, 2, 3, 4, 1, 2, &
		&3, 4, 1, 2, 3, 4, 3, 4, 4, 0, &
		&1, 2, 0, 1, 2, 0, 1, 2, 3, 4, &
		&2, 3, 4, 0, 0/)

end subroutine n_A_in_clusters_2

subroutine clusters_with_1_A_2(cluster_numbers)
	implicit none
	integer :: cluster_numbers(5)

	cluster_numbers = (/1, 6, 11, 16, 21/)

end subroutine clusters_with_1_A_2

subroutine arrays_2(neutrals, negatives, positives)
	implicit none
	integer :: neutrals(24), negatives(16), positives(15)

	neutrals = (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, &
		&21, 22, 23, 24/)
	negatives = (/25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 54/)
	positives = (/40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 55/)

end subroutine arrays_2

subroutine cluster_arrays_2(neutral_clusters, negative_clusters, positive_clusters)
	implicit none
	integer :: neutral_clusters(22), negative_clusters(14), positive_clusters(13)

	neutral_clusters = (/2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, &
		&21, 22, 23, 24/)
	negative_clusters = (/26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39/)
	positive_clusters = (/41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53/)

end subroutine cluster_arrays_2

subroutine get_charging_state_2(charging_state)
	implicit none
	integer :: charging_state(55)

	charging_state = (/1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
		&1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
		&1, 1, 1, 1, 2, 2, 2, 2, 2, 2, &
		&2, 2, 2, 2, 2, 2, 2, 2, 2, 3, &
		&3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
		&3, 3, 3, 2, 3/)

end subroutine get_charging_state_2

subroutine get_mass_2(mass)
	implicit none
	real(kind(1.d0)) :: mass(55)

	mass = (/98.08, 196.16, 294.24, 392.32, 45.08, 143.16, 241.24, 339.32, 437.40, 90.16, &
		&188.24, 286.32, 384.40, 482.48, 135.24, 233.32, 331.40, 429.48, 527.56, 180.32, &
		&278.40, 376.48, 474.56, 572.64, 97.08, 195.16, 293.24, 391.32, 142.16, 240.24, &
		&338.32, 436.40, 187.24, 285.32, 383.40, 481.48, 428.48, 526.56, 571.64, 46.08, &
		&144.16, 242.24, 91.16, 189.24, 287.32, 136.24, 234.32, 332.40, 430.48, 528.56, &
		&377.48, 475.56, 573.64, 32.00, 19.02/)

end subroutine get_mass_2

subroutine get_diameter_2(diameter)
	implicit none

	real(kind(1.d0)) :: diameter(55)

	 diameter = (/0.55, 0.70, 0.80, 0.88, 0.59, 0.72, 0.82, 0.90, 0.96, 0.75, &
		&0.84, 0.91, 0.98, 1.03, 0.86, 0.93, 0.99, 1.04, 1.09, 0.94, &
		&1.00, 1.06, 1.11, 1.15, 0.55, 0.70, 0.80, 0.88, 0.72, 0.82, &
		&0.90, 0.96, 0.84, 0.91, 0.98, 1.03, 1.04, 1.09, 1.15, 0.59, &
		&0.72, 0.82, 0.75, 0.84, 0.91, 0.86, 0.93, 0.99, 1.04, 1.09, &
		&1.06, 1.11, 1.15, 0.45, 0.39/)	! dry value

end subroutine get_diameter_2

subroutine get_mob_diameter_2(mob_diameter)
	implicit none

	real(kind(1.d0)) :: mob_diameter(55)

	 mob_diameter = (/0.85, 1.00, 1.10, 1.18, 0.89, 1.02, 1.12, 1.20, 1.26, 1.05, &
		&1.14, 1.21, 1.28, 1.33, 1.16, 1.23, 1.29, 1.34, 1.39, 1.24, &
		&1.30, 1.36, 1.41, 1.45, 0.85, 1.00, 1.10, 1.18, 1.02, 1.12, &
		&1.20, 1.26, 1.14, 1.21, 1.28, 1.33, 1.34, 1.39, 1.45, 0.89, &
		&1.02, 1.12, 1.05, 1.14, 1.21, 1.16, 1.23, 1.29, 1.34, 1.39, &
		&1.36, 1.41, 1.45, 0.75, 0.69/)	! dry value

end subroutine get_mob_diameter_2

subroutine cluster_names_2(clust)
	implicit none
	character(len=11), dimension(65) :: clust

	clust(1)(:) = '1A'
	clust(2)(:) = '2A'
	clust(3)(:) = '3A'
	clust(4)(:) = '4A'
	clust(5)(:) = '1D'
	clust(6)(:) = '1A1D'
	clust(7)(:) = '2A1D'
	clust(8)(:) = '3A1D'
	clust(9)(:) = '4A1D'
	clust(10)(:) = '2D'
	clust(11)(:) = '1A2D'
	clust(12)(:) = '2A2D'
	clust(13)(:) = '3A2D'
	clust(14)(:) = '4A2D'
	clust(15)(:) = '3D'
	clust(16)(:) = '1A3D'
	clust(17)(:) = '2A3D'
	clust(18)(:) = '3A3D'
	clust(19)(:) = '4A3D'
	clust(20)(:) = '4D'
	clust(21)(:) = '1A4D'
	clust(22)(:) = '2A4D'
	clust(23)(:) = '3A4D'
	clust(24)(:) = '4A4D'
	clust(25)(:) = '1B'
	clust(26)(:) = '1A1B'
	clust(27)(:) = '2A1B'
	clust(28)(:) = '3A1B'
	clust(29)(:) = '1B1D'
	clust(30)(:) = '1A1B1D'
	clust(31)(:) = '2A1B1D'
	clust(32)(:) = '3A1B1D'
	clust(33)(:) = '1B2D'
	clust(34)(:) = '1A1B2D'
	clust(35)(:) = '2A1B2D'
	clust(36)(:) = '3A1B2D'
	clust(37)(:) = '2A1B3D'
	clust(38)(:) = '3A1B3D'
	clust(39)(:) = '3A1B4D'
	clust(40)(:) = '1D1P'
	clust(41)(:) = '1A1D1P'
	clust(42)(:) = '2A1D1P'
	clust(43)(:) = '2D1P'
	clust(44)(:) = '1A2D1P'
	clust(45)(:) = '2A2D1P'
	clust(46)(:) = '3D1P'
	clust(47)(:) = '1A3D1P'
	clust(48)(:) = '2A3D1P'
	clust(49)(:) = '3A3D1P'
	clust(50)(:) = '4A3D1P'
	clust(51)(:) = '2A4D1P'
	clust(52)(:) = '3A4D1P'
	clust(53)(:) = '4A4D1P'
	clust(54)(:) = 'neg'
	clust(55)(:) = 'pos'
	clust(56)(:) = 'source'
	clust(57)(:) = 'coag'
	clust(58)(:) = 'wall'
	clust(59)(:) = 'dil'
	clust(60)(:) = 'insink'
	clust(61)(:) = 'rec'
	clust(62)(:) = 'out_neu'
	clust(63)(:) = 'out_neg'
	clust(64)(:) = 'out_pos'
	clust(65)(:) = 'bound'

	! Additional flux elements not listed here:
	! Elements from 66 to 120: coagulation loss of each cluster
	! Elements from 121 to 196: each outgrown cluster composition

end subroutine cluster_names_2

subroutine monomer_names_2(clust_mon)
	implicit none
	character(len=11), dimension(4) :: clust_mon

	clust_mon(1)(:) = '1A'
	clust_mon(2)(:) = '1D'
	clust_mon(3)(:) = '1B'
	clust_mon(4)(:) = '1D1P'

end subroutine monomer_names_2

subroutine molecule_names_2(labels)
	implicit none
	character(len=11), dimension(4) :: labels

	labels(1)(:) = 'A'
	labels(2)(:) = 'B'
	labels(3)(:) = 'D'
	labels(4)(:) = 'P'

end subroutine molecule_names_2

subroutine monomer_indices_2(n_monomers)
	implicit none
	integer :: n_monomers(4)

	n_monomers = (/1, 5, 25, 40/)

end subroutine monomer_indices_2

subroutine get_bound_2(bound_clusters,nmols_bound)
	implicit none
	integer :: bound_clusters(18), nmols_bound(18,4)

	nmols_bound(1,:) = (/4, 0, 0, 0/)
	nmols_bound(2,:) = (/4, 0, 1, 0/)
	nmols_bound(3,:) = (/4, 0, 2, 0/)
	nmols_bound(4,:) = (/4, 0, 3, 0/)
	nmols_bound(5,:) = (/0, 0, 4, 0/)
	nmols_bound(6,:) = (/1, 0, 4, 0/)
	nmols_bound(7,:) = (/2, 0, 4, 0/)
	nmols_bound(8,:) = (/3, 0, 4, 0/)
	nmols_bound(9,:) = (/4, 0, 4, 0/)
	nmols_bound(10,:) = (/3, 1, 0, 0/)
	nmols_bound(11,:) = (/3, 1, 1, 0/)
	nmols_bound(12,:) = (/3, 1, 2, 0/)
	nmols_bound(13,:) = (/3, 1, 3, 0/)
	nmols_bound(14,:) = (/3, 1, 4, 0/)
	nmols_bound(15,:) = (/4, 0, 3, 1/)
	nmols_bound(16,:) = (/2, 0, 4, 1/)
	nmols_bound(17,:) = (/3, 0, 4, 1/)
	nmols_bound(18,:) = (/4, 0, 4, 1/)

	bound_clusters = (/4, 9, 14, 19, 20, 21, 22, 23, 24, 28, 32, 36, 38, 39, 50, 51, 52, 53/)

end subroutine get_bound_2

subroutine get_diameter_out_2(diameter_out)
	implicit none
	real(kind(1.d0)) :: diameter_out(76)

	diameter_out = (/1.19, 1.08, 1.14, 1.19, 1.23, 1.13, 1.18, 1.23, 1.27, 1.17, &
		&1.22, 1.27, 1.30, 1.21, 1.26, 1.30, 1.20, 1.16, 1.20, 1.24, &
		&1.24, 1.24, 1.27, 1.27, 1.27, 1.31, 1.31, 1.31, 1.34, 1.34, &
		&1.34, 1.25, 1.21, 1.25, 1.28, 1.28, 1.28, 1.32, 1.32, 1.32, &
		&1.35, 1.35, 1.35, 1.38, 1.38, 1.38, 1.29, 1.26, 1.29, 1.32, &
		&1.32, 1.32, 1.36, 1.36, 1.36, 1.39, 1.39, 1.39, 1.41, 1.41, &
		&1.41, 1.33, 1.30, 1.33, 1.36, 1.36, 1.36, 1.39, 1.39, 1.39, &
		&1.42, 1.42, 1.42, 1.45, 1.45, 1.45/)	! dry value

end subroutine get_diameter_out_2

subroutine get_nmols_out_2(nmols_out)
	implicit none
	integer :: nmols_out(76,4)

	nmols_out(1,:) = (/5, 0, 4, 0/)
	nmols_out(2,:) = (/4, 1, 2, 0/)
	nmols_out(3,:) = (/4, 1, 3, 0/)
	nmols_out(4,:) = (/4, 1, 4, 0/)
	nmols_out(5,:) = (/6, 0, 4, 0/)
	nmols_out(6,:) = (/5, 1, 2, 0/)
	nmols_out(7,:) = (/5, 1, 3, 0/)
	nmols_out(8,:) = (/5, 1, 4, 0/)
	nmols_out(9,:) = (/7, 0, 4, 0/)
	nmols_out(10,:) = (/6, 1, 2, 0/)
	nmols_out(11,:) = (/6, 1, 3, 0/)
	nmols_out(12,:) = (/6, 1, 4, 0/)
	nmols_out(13,:) = (/8, 0, 4, 0/)
	nmols_out(14,:) = (/7, 1, 2, 0/)
	nmols_out(15,:) = (/7, 1, 3, 0/)
	nmols_out(16,:) = (/7, 1, 4, 0/)
	nmols_out(17,:) = (/4, 0, 5, 0/)
	nmols_out(18,:) = (/3, 0, 5, 1/)
	nmols_out(19,:) = (/4, 0, 5, 1/)
	nmols_out(20,:) = (/5, 0, 5, 0/)
	nmols_out(21,:) = (/4, 1, 5, 0/)
	nmols_out(22,:) = (/5, 0, 5, 1/)
	nmols_out(23,:) = (/6, 0, 5, 0/)
	nmols_out(24,:) = (/5, 1, 5, 0/)
	nmols_out(25,:) = (/6, 0, 5, 1/)
	nmols_out(26,:) = (/7, 0, 5, 0/)
	nmols_out(27,:) = (/6, 1, 5, 0/)
	nmols_out(28,:) = (/7, 0, 5, 1/)
	nmols_out(29,:) = (/8, 0, 5, 0/)
	nmols_out(30,:) = (/7, 1, 5, 0/)
	nmols_out(31,:) = (/8, 0, 5, 1/)
	nmols_out(32,:) = (/4, 0, 6, 0/)
	nmols_out(33,:) = (/3, 0, 6, 1/)
	nmols_out(34,:) = (/4, 0, 6, 1/)
	nmols_out(35,:) = (/5, 0, 6, 0/)
	nmols_out(36,:) = (/4, 1, 6, 0/)
	nmols_out(37,:) = (/5, 0, 6, 1/)
	nmols_out(38,:) = (/6, 0, 6, 0/)
	nmols_out(39,:) = (/5, 1, 6, 0/)
	nmols_out(40,:) = (/6, 0, 6, 1/)
	nmols_out(41,:) = (/7, 0, 6, 0/)
	nmols_out(42,:) = (/6, 1, 6, 0/)
	nmols_out(43,:) = (/7, 0, 6, 1/)
	nmols_out(44,:) = (/8, 0, 6, 0/)
	nmols_out(45,:) = (/7, 1, 6, 0/)
	nmols_out(46,:) = (/8, 0, 6, 1/)
	nmols_out(47,:) = (/4, 0, 7, 0/)
	nmols_out(48,:) = (/3, 0, 7, 1/)
	nmols_out(49,:) = (/4, 0, 7, 1/)
	nmols_out(50,:) = (/5, 0, 7, 0/)
	nmols_out(51,:) = (/4, 1, 7, 0/)
	nmols_out(52,:) = (/5, 0, 7, 1/)
	nmols_out(53,:) = (/6, 0, 7, 0/)
	nmols_out(54,:) = (/5, 1, 7, 0/)
	nmols_out(55,:) = (/6, 0, 7, 1/)
	nmols_out(56,:) = (/7, 0, 7, 0/)
	nmols_out(57,:) = (/6, 1, 7, 0/)
	nmols_out(58,:) = (/7, 0, 7, 1/)
	nmols_out(59,:) = (/8, 0, 7, 0/)
	nmols_out(60,:) = (/7, 1, 7, 0/)
	nmols_out(61,:) = (/8, 0, 7, 1/)
	nmols_out(62,:) = (/4, 0, 8, 0/)
	nmols_out(63,:) = (/3, 0, 8, 1/)
	nmols_out(64,:) = (/4, 0, 8, 1/)
	nmols_out(65,:) = (/5, 0, 8, 0/)
	nmols_out(66,:) = (/4, 1, 8, 0/)
	nmols_out(67,:) = (/5, 0, 8, 1/)
	nmols_out(68,:) = (/6, 0, 8, 0/)
	nmols_out(69,:) = (/5, 1, 8, 0/)
	nmols_out(70,:) = (/6, 0, 8, 1/)
	nmols_out(71,:) = (/7, 0, 8, 0/)
	nmols_out(72,:) = (/6, 1, 8, 0/)
	nmols_out(73,:) = (/7, 0, 8, 1/)
	nmols_out(74,:) = (/8, 0, 8, 0/)
	nmols_out(75,:) = (/7, 1, 8, 0/)
	nmols_out(76,:) = (/8, 0, 8, 1/)

end subroutine get_nmols_out_2

subroutine molecule_names_nocharge_2(labels_nocharge)
	implicit none
	character(len=11), dimension(2) :: labels_nocharge

	labels_nocharge(1)(:) = 'A'
	labels_nocharge(2)(:) = 'D'

end subroutine molecule_names_nocharge_2

subroutine get_nmols_nocharge_2(nmols_nocharge_clust)
	implicit none
	integer :: nmols_nocharge_clust(55,2)

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
	nmols_nocharge_clust(25,:) = (/1, 0/)
	nmols_nocharge_clust(26,:) = (/2, 0/)
	nmols_nocharge_clust(27,:) = (/3, 0/)
	nmols_nocharge_clust(28,:) = (/4, 0/)
	nmols_nocharge_clust(29,:) = (/1, 1/)
	nmols_nocharge_clust(30,:) = (/2, 1/)
	nmols_nocharge_clust(31,:) = (/3, 1/)
	nmols_nocharge_clust(32,:) = (/4, 1/)
	nmols_nocharge_clust(33,:) = (/1, 2/)
	nmols_nocharge_clust(34,:) = (/2, 2/)
	nmols_nocharge_clust(35,:) = (/3, 2/)
	nmols_nocharge_clust(36,:) = (/4, 2/)
	nmols_nocharge_clust(37,:) = (/3, 3/)
	nmols_nocharge_clust(38,:) = (/4, 3/)
	nmols_nocharge_clust(39,:) = (/4, 4/)
	nmols_nocharge_clust(40,:) = (/0, 1/)
	nmols_nocharge_clust(41,:) = (/1, 1/)
	nmols_nocharge_clust(42,:) = (/2, 1/)
	nmols_nocharge_clust(43,:) = (/0, 2/)
	nmols_nocharge_clust(44,:) = (/1, 2/)
	nmols_nocharge_clust(45,:) = (/2, 2/)
	nmols_nocharge_clust(46,:) = (/0, 3/)
	nmols_nocharge_clust(47,:) = (/1, 3/)
	nmols_nocharge_clust(48,:) = (/2, 3/)
	nmols_nocharge_clust(49,:) = (/3, 3/)
	nmols_nocharge_clust(50,:) = (/4, 3/)
	nmols_nocharge_clust(51,:) = (/2, 4/)
	nmols_nocharge_clust(52,:) = (/3, 4/)
	nmols_nocharge_clust(53,:) = (/4, 4/)

end subroutine get_nmols_nocharge_2

subroutine get_nmols_nocharge_out_2(nmols_nocharge_out)
	implicit none
	integer :: nmols_nocharge_out(76,2)

	nmols_nocharge_out = 0

	nmols_nocharge_out(1,:) = (/5, 4/)
	nmols_nocharge_out(2,:) = (/5, 2/)
	nmols_nocharge_out(3,:) = (/5, 3/)
	nmols_nocharge_out(4,:) = (/5, 4/)
	nmols_nocharge_out(5,:) = (/6, 4/)
	nmols_nocharge_out(6,:) = (/6, 2/)
	nmols_nocharge_out(7,:) = (/6, 3/)
	nmols_nocharge_out(8,:) = (/6, 4/)
	nmols_nocharge_out(9,:) = (/7, 4/)
	nmols_nocharge_out(10,:) = (/7, 2/)
	nmols_nocharge_out(11,:) = (/7, 3/)
	nmols_nocharge_out(12,:) = (/7, 4/)
	nmols_nocharge_out(13,:) = (/8, 4/)
	nmols_nocharge_out(14,:) = (/8, 2/)
	nmols_nocharge_out(15,:) = (/8, 3/)
	nmols_nocharge_out(16,:) = (/8, 4/)
	nmols_nocharge_out(17,:) = (/4, 5/)
	nmols_nocharge_out(18,:) = (/3, 5/)
	nmols_nocharge_out(19,:) = (/4, 5/)
	nmols_nocharge_out(20,:) = (/5, 5/)
	nmols_nocharge_out(21,:) = (/5, 5/)
	nmols_nocharge_out(22,:) = (/5, 5/)
	nmols_nocharge_out(23,:) = (/6, 5/)
	nmols_nocharge_out(24,:) = (/6, 5/)
	nmols_nocharge_out(25,:) = (/6, 5/)
	nmols_nocharge_out(26,:) = (/7, 5/)
	nmols_nocharge_out(27,:) = (/7, 5/)
	nmols_nocharge_out(28,:) = (/7, 5/)
	nmols_nocharge_out(29,:) = (/8, 5/)
	nmols_nocharge_out(30,:) = (/8, 5/)
	nmols_nocharge_out(31,:) = (/8, 5/)
	nmols_nocharge_out(32,:) = (/4, 6/)
	nmols_nocharge_out(33,:) = (/3, 6/)
	nmols_nocharge_out(34,:) = (/4, 6/)
	nmols_nocharge_out(35,:) = (/5, 6/)
	nmols_nocharge_out(36,:) = (/5, 6/)
	nmols_nocharge_out(37,:) = (/5, 6/)
	nmols_nocharge_out(38,:) = (/6, 6/)
	nmols_nocharge_out(39,:) = (/6, 6/)
	nmols_nocharge_out(40,:) = (/6, 6/)
	nmols_nocharge_out(41,:) = (/7, 6/)
	nmols_nocharge_out(42,:) = (/7, 6/)
	nmols_nocharge_out(43,:) = (/7, 6/)
	nmols_nocharge_out(44,:) = (/8, 6/)
	nmols_nocharge_out(45,:) = (/8, 6/)
	nmols_nocharge_out(46,:) = (/8, 6/)
	nmols_nocharge_out(47,:) = (/4, 7/)
	nmols_nocharge_out(48,:) = (/3, 7/)
	nmols_nocharge_out(49,:) = (/4, 7/)
	nmols_nocharge_out(50,:) = (/5, 7/)
	nmols_nocharge_out(51,:) = (/5, 7/)
	nmols_nocharge_out(52,:) = (/5, 7/)
	nmols_nocharge_out(53,:) = (/6, 7/)
	nmols_nocharge_out(54,:) = (/6, 7/)
	nmols_nocharge_out(55,:) = (/6, 7/)
	nmols_nocharge_out(56,:) = (/7, 7/)
	nmols_nocharge_out(57,:) = (/7, 7/)
	nmols_nocharge_out(58,:) = (/7, 7/)
	nmols_nocharge_out(59,:) = (/8, 7/)
	nmols_nocharge_out(60,:) = (/8, 7/)
	nmols_nocharge_out(61,:) = (/8, 7/)
	nmols_nocharge_out(62,:) = (/4, 8/)
	nmols_nocharge_out(63,:) = (/3, 8/)
	nmols_nocharge_out(64,:) = (/4, 8/)
	nmols_nocharge_out(65,:) = (/5, 8/)
	nmols_nocharge_out(66,:) = (/5, 8/)
	nmols_nocharge_out(67,:) = (/5, 8/)
	nmols_nocharge_out(68,:) = (/6, 8/)
	nmols_nocharge_out(69,:) = (/6, 8/)
	nmols_nocharge_out(70,:) = (/6, 8/)
	nmols_nocharge_out(71,:) = (/7, 8/)
	nmols_nocharge_out(72,:) = (/7, 8/)
	nmols_nocharge_out(73,:) = (/7, 8/)
	nmols_nocharge_out(74,:) = (/8, 8/)
	nmols_nocharge_out(75,:) = (/8, 8/)
	nmols_nocharge_out(76,:) = (/8, 8/)

end subroutine get_nmols_nocharge_out_2


end module acdc_system_2

