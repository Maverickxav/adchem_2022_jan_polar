module acdc_system_2

implicit none

character(len=*), parameter :: tag = 'test'			! tag for the simulation system

logical, parameter :: small_set_mode = .true.			! small set of clusters (default mode)

integer, parameter :: nclust = 53						! number of clusters, molecules and ions
integer, parameter :: neq = 188							! number of equations
integer, parameter :: nclust_max = 53, neq_max = 188	! same as nclust and neq in the small set mode

integer, parameter :: n_variables = 3
logical, parameter :: variable_temp = .true.
logical, parameter :: variable_cs = .false.
logical, parameter :: variable_ipr = .true.
real(kind(1.d0)), parameter :: rh = 0.d0					! RH in %
real(kind(1.d0)), parameter :: pw = 0.d0					! partial pressure of water in Pa

integer, parameter :: n_monomer_types = 4
integer, parameter :: n_charges = 3			! number of charging states
integer, parameter :: n1A = 1, n1D = 5, n1B = 25, n1D1P = 40, nneg = 52, npos = 53			! cluster indices for monomers and ions

integer, parameter :: nout_all(3) = (/60, 61, 62/), nout_neu = 60, nout_neg = 61, nout_pos = 62			! indices for outgoing fluxes

integer, parameter :: nclust_coag_saved = 53			! number of saved scavenged clusters
integer, parameter :: nclust_out_saved = 72				! number of saved outgrown clusters
integer, parameter :: nrange_coag_saved(2) = (/64, 116/)	! index range of saved scavenged clusters
integer, parameter :: nrange_out_saved(2) = (/117, 188/)	! index range of saved outgrown clusters

integer, parameter :: n_mol_types = 4
integer, parameter :: nmolA = 1, nmolB = 2, nmolD = 3, nmolP = 4			! molecule indices for the used species

integer, parameter :: n_1A_clusters = 5				! number molecules and clusters containing 1 A molecule

integer, parameter :: n_neutrals = 24			! number of neutral molecules and clusters
integer, parameter :: n_negatives = 16			! negative
integer, parameter :: n_positives = 13			! positive
integer, parameter :: nclust_nogen = 51			! number of clusters and molecules excluding generic ions

integer, parameter :: n_neutral_monomers = 2, neutral_monomers(2) = (/1, 5/)			! number and indices of neutral monomers
integer, parameter :: n_negative_monomers = 1, negative_monomers(1) = (/25/)			! negative
integer, parameter :: n_positive_monomers = 1, positive_monomers(1) = (/40/)			! positive

integer, parameter :: n_neutral_clusters = 22			! number of neutral clusters
integer, parameter :: n_negative_clusters = 14			! negative
integer, parameter :: n_positive_clusters = 11			! positive

real(kind(1.d0)), parameter :: mass_max = 572.6400
real(kind(1.d0)), parameter :: diameter_max = 1.1500
real(kind(1.d0)), parameter :: mob_diameter_max = 1.4500
integer, parameter :: ij_ind_max(4) = (/4, 1, 4, 1/)		! maximum molecular content
integer, parameter :: n_bound = 17		! number of clusters at the system boundary
integer, parameter :: n_comp_out = 72		! number of clusters growing out in different compositions

integer, parameter :: nmols_out_neutral(2, 4) = reshape((/5, 4, 0, 0, 4, 5, 0, 0/),(/2, 4/))			! criteria for outgrowing neutrals
integer, parameter :: nmols_out_negative(1, 4) = reshape((/4, 1, 2, 0/),(/1, 4/))			! negatives
integer, parameter :: nmols_out_positive(1, 4) = reshape((/3, 0, 5, 1/),(/1, 4/))			! positives


contains

subroutine n_A_in_clusters_2(n_A)
	implicit none
	integer :: n_A(53)

	n_A = (/1, 2, 3, 4, 0, 1, 2, 3, 4, 0, &
		&1, 2, 3, 4, 0, 1, 2, 3, 4, 0, &
		&1, 2, 3, 4, 1, 2, 3, 4, 1, 2, &
		&3, 4, 1, 2, 3, 4, 3, 4, 4, 0, &
		&1, 2, 0, 1, 2, 0, 1, 2, 3, 2, &
		&3, 0, 0/)

end subroutine n_A_in_clusters_2

subroutine clusters_with_1_A_2(cluster_numbers)
	implicit none
	integer :: cluster_numbers(5)

	cluster_numbers = (/1, 6, 11, 16, 21/)

end subroutine clusters_with_1_A_2

subroutine arrays_2(neutrals, negatives, positives)
	implicit none
	integer :: neutrals(24), negatives(16), positives(13)

	neutrals = (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, &
		&21, 22, 23, 24/)
	negatives = (/25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 52/)
	positives = (/40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 53/)

end subroutine arrays_2

subroutine cluster_arrays_2(neutral_clusters, negative_clusters, positive_clusters)
	implicit none
	integer :: neutral_clusters(22), negative_clusters(14), positive_clusters(11)

	neutral_clusters = (/2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, &
		&21, 22, 23, 24/)
	negative_clusters = (/26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39/)
	positive_clusters = (/41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51/)

end subroutine cluster_arrays_2

subroutine get_charging_state_2(charging_state)
	implicit none
	integer :: charging_state(53)

	charging_state = (/1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
		&1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
		&1, 1, 1, 1, 2, 2, 2, 2, 2, 2, &
		&2, 2, 2, 2, 2, 2, 2, 2, 2, 3, &
		&3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
		&3, 2, 3/)

end subroutine get_charging_state_2

subroutine get_mass_2(mass)
	implicit none
	real(kind(1.d0)) :: mass(53)

	mass = (/98.0800, 196.1600, 294.2400, 392.3200, 45.0800, 143.1600, 241.2400, 339.3200, 437.4000, 90.1600, &
		&188.2400, 286.3200, 384.4000, 482.4800, 135.2400, 233.3200, 331.4000, 429.4800, 527.5600, 180.3200, &
		&278.4000, 376.4800, 474.5600, 572.6400, 97.0800, 195.1600, 293.2400, 391.3200, 142.1600, 240.2400, &
		&338.3200, 436.4000, 187.2400, 285.3200, 383.4000, 481.4800, 428.4800, 526.5600, 571.6400, 46.0800, &
		&144.1600, 242.2400, 91.1600, 189.2400, 287.3200, 136.2400, 234.3200, 332.4000, 430.4800, 377.4800, &
		&475.5600, 32.0000, 19.0200/)

end subroutine get_mass_2

subroutine get_diameter_2(diameter)
	implicit none

	real(kind(1.d0)) :: diameter(53)

	 diameter = (/0.5539, 0.6979, 0.7989, 0.8793, 0.5946, 0.7245, 0.8194, 0.8963, 0.9619, 0.7492, &
		&0.8389, 0.9128, 0.9762, 1.0324, 0.8576, 0.9286, 0.9901, 1.0448, 1.0944, 0.9439, &
		&1.0036, 1.0570, 1.1055, 1.1500, 0.5520, 0.6967, 0.7980, 0.8786, 0.7234, 0.8186, &
		&0.8956, 0.9613, 0.8381, 0.9121, 0.9756, 1.0319, 1.0443, 1.0939, 1.1496, 0.5946, &
		&0.7245, 0.8194, 0.7492, 0.8389, 0.9128, 0.8576, 0.9286, 0.9901, 1.0448, 1.0570, &
		&1.1055, 0.4464, 0.3926/)	! dry value

end subroutine get_diameter_2

subroutine get_mob_diameter_2(mob_diameter)
	implicit none

	real(kind(1.d0)) :: mob_diameter(53)

	 mob_diameter = (/0.8539, 0.9979, 1.0989, 1.1793, 0.8946, 1.0245, 1.1194, 1.1963, 1.2619, 1.0492, &
		&1.1389, 1.2128, 1.2762, 1.3324, 1.1576, 1.2286, 1.2901, 1.3448, 1.3944, 1.2439, &
		&1.3036, 1.3570, 1.4055, 1.4500, 0.8520, 0.9967, 1.0980, 1.1786, 1.0234, 1.1186, &
		&1.1956, 1.2613, 1.1381, 1.2121, 1.2756, 1.3319, 1.3443, 1.3939, 1.4496, 0.8946, &
		&1.0245, 1.1194, 1.0492, 1.1389, 1.2128, 1.1576, 1.2286, 1.2901, 1.3448, 1.3570, &
		&1.4055, 0.7464, 0.6926/)	! dry value

end subroutine get_mob_diameter_2

subroutine cluster_names_2(clust)
	implicit none
	character(len=11), dimension(63) :: clust

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
	clust(50)(:) = '2A4D1P'
	clust(51)(:) = '3A4D1P'
	clust(52)(:) = 'neg'
	clust(53)(:) = 'pos'
	clust(54)(:) = 'source'
	clust(55)(:) = 'coag'
	clust(56)(:) = 'wall'
	clust(57)(:) = 'dil'
	clust(58)(:) = 'insink'
	clust(59)(:) = 'rec'
	clust(60)(:) = 'out_neu'
	clust(61)(:) = 'out_neg'
	clust(62)(:) = 'out_pos'
	clust(63)(:) = 'bound'

	! Additional flux elements not listed here:
	! Elements from 64 to 116: coagulation loss of each cluster
	! Elements from 117 to 188: each outgrown cluster composition

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
	integer :: bound_clusters(17), nmols_bound(17,4)

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
	nmols_bound(15,:) = (/3, 0, 3, 1/)
	nmols_bound(16,:) = (/2, 0, 4, 1/)
	nmols_bound(17,:) = (/3, 0, 4, 1/)

	bound_clusters = (/4, 9, 14, 19, 20, 21, 22, 23, 24, 28, 32, 36, 38, 39, 49, 50, 51/)

end subroutine get_bound_2

subroutine get_diameter_out_2(diameter_out)
	implicit none
	real(kind(1.d0)) :: diameter_out(72)

	diameter_out = (/1.1913, 1.0825, 1.1393, 1.1909, 1.2300, 1.1289, 1.1814, 1.2296, 1.2664, 1.1717, &
		&1.2207, 1.2660, 1.3007, 1.2116, 1.2576, 1.3004, 1.2007, 1.1601, 1.2388, 1.2384, &
		&1.2007, 1.2747, 1.2743, 1.2388, 1.3086, 1.3083, 1.2747, 1.3409, 1.3406, 1.3086, &
		&1.2475, 1.2100, 1.2829, 1.2825, 1.2475, 1.3164, 1.3161, 1.2829, 1.3483, 1.3480, &
		&1.3164, 1.3788, 1.3785, 1.3483, 1.2910, 1.2561, 1.3241, 1.3238, 1.2910, 1.3557, &
		&1.3554, 1.3241, 1.3858, 1.3855, 1.3557, 1.4147, 1.4144, 1.3858, 1.3317, 1.2990, &
		&1.3629, 1.3626, 1.3317, 1.3928, 1.3925, 1.3629, 1.4214, 1.4211, 1.3928, 1.4489, &
		&1.4486, 1.4214/)	! dry value

end subroutine get_diameter_out_2

subroutine get_nmols_out_2(nmols_out)
	implicit none
	integer :: nmols_out(72,4)

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
	nmols_out(19,:) = (/5, 0, 5, 0/)
	nmols_out(20,:) = (/4, 1, 5, 0/)
	nmols_out(21,:) = (/4, 0, 5, 1/)
	nmols_out(22,:) = (/6, 0, 5, 0/)
	nmols_out(23,:) = (/5, 1, 5, 0/)
	nmols_out(24,:) = (/5, 0, 5, 1/)
	nmols_out(25,:) = (/7, 0, 5, 0/)
	nmols_out(26,:) = (/6, 1, 5, 0/)
	nmols_out(27,:) = (/6, 0, 5, 1/)
	nmols_out(28,:) = (/8, 0, 5, 0/)
	nmols_out(29,:) = (/7, 1, 5, 0/)
	nmols_out(30,:) = (/7, 0, 5, 1/)
	nmols_out(31,:) = (/4, 0, 6, 0/)
	nmols_out(32,:) = (/3, 0, 6, 1/)
	nmols_out(33,:) = (/5, 0, 6, 0/)
	nmols_out(34,:) = (/4, 1, 6, 0/)
	nmols_out(35,:) = (/4, 0, 6, 1/)
	nmols_out(36,:) = (/6, 0, 6, 0/)
	nmols_out(37,:) = (/5, 1, 6, 0/)
	nmols_out(38,:) = (/5, 0, 6, 1/)
	nmols_out(39,:) = (/7, 0, 6, 0/)
	nmols_out(40,:) = (/6, 1, 6, 0/)
	nmols_out(41,:) = (/6, 0, 6, 1/)
	nmols_out(42,:) = (/8, 0, 6, 0/)
	nmols_out(43,:) = (/7, 1, 6, 0/)
	nmols_out(44,:) = (/7, 0, 6, 1/)
	nmols_out(45,:) = (/4, 0, 7, 0/)
	nmols_out(46,:) = (/3, 0, 7, 1/)
	nmols_out(47,:) = (/5, 0, 7, 0/)
	nmols_out(48,:) = (/4, 1, 7, 0/)
	nmols_out(49,:) = (/4, 0, 7, 1/)
	nmols_out(50,:) = (/6, 0, 7, 0/)
	nmols_out(51,:) = (/5, 1, 7, 0/)
	nmols_out(52,:) = (/5, 0, 7, 1/)
	nmols_out(53,:) = (/7, 0, 7, 0/)
	nmols_out(54,:) = (/6, 1, 7, 0/)
	nmols_out(55,:) = (/6, 0, 7, 1/)
	nmols_out(56,:) = (/8, 0, 7, 0/)
	nmols_out(57,:) = (/7, 1, 7, 0/)
	nmols_out(58,:) = (/7, 0, 7, 1/)
	nmols_out(59,:) = (/4, 0, 8, 0/)
	nmols_out(60,:) = (/3, 0, 8, 1/)
	nmols_out(61,:) = (/5, 0, 8, 0/)
	nmols_out(62,:) = (/4, 1, 8, 0/)
	nmols_out(63,:) = (/4, 0, 8, 1/)
	nmols_out(64,:) = (/6, 0, 8, 0/)
	nmols_out(65,:) = (/5, 1, 8, 0/)
	nmols_out(66,:) = (/5, 0, 8, 1/)
	nmols_out(67,:) = (/7, 0, 8, 0/)
	nmols_out(68,:) = (/6, 1, 8, 0/)
	nmols_out(69,:) = (/6, 0, 8, 1/)
	nmols_out(70,:) = (/8, 0, 8, 0/)
	nmols_out(71,:) = (/7, 1, 8, 0/)
	nmols_out(72,:) = (/7, 0, 8, 1/)

end subroutine get_nmols_out_2

subroutine molecule_names_nocharge_2(labels_nocharge)
	implicit none
	character(len=11), dimension(2) :: labels_nocharge

	labels_nocharge(1)(:) = 'A'
	labels_nocharge(2)(:) = 'D'

end subroutine molecule_names_nocharge_2

subroutine get_nmols_nocharge_2(nmols_nocharge_clust)
	implicit none
	integer :: nmols_nocharge_clust(53,2)

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
	nmols_nocharge_clust(50,:) = (/2, 4/)
	nmols_nocharge_clust(51,:) = (/3, 4/)

end subroutine get_nmols_nocharge_2

subroutine get_nmols_nocharge_out_2(nmols_nocharge_out)
	implicit none
	integer :: nmols_nocharge_out(72,2)

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
	nmols_nocharge_out(19,:) = (/5, 5/)
	nmols_nocharge_out(20,:) = (/5, 5/)
	nmols_nocharge_out(21,:) = (/4, 5/)
	nmols_nocharge_out(22,:) = (/6, 5/)
	nmols_nocharge_out(23,:) = (/6, 5/)
	nmols_nocharge_out(24,:) = (/5, 5/)
	nmols_nocharge_out(25,:) = (/7, 5/)
	nmols_nocharge_out(26,:) = (/7, 5/)
	nmols_nocharge_out(27,:) = (/6, 5/)
	nmols_nocharge_out(28,:) = (/8, 5/)
	nmols_nocharge_out(29,:) = (/8, 5/)
	nmols_nocharge_out(30,:) = (/7, 5/)
	nmols_nocharge_out(31,:) = (/4, 6/)
	nmols_nocharge_out(32,:) = (/3, 6/)
	nmols_nocharge_out(33,:) = (/5, 6/)
	nmols_nocharge_out(34,:) = (/5, 6/)
	nmols_nocharge_out(35,:) = (/4, 6/)
	nmols_nocharge_out(36,:) = (/6, 6/)
	nmols_nocharge_out(37,:) = (/6, 6/)
	nmols_nocharge_out(38,:) = (/5, 6/)
	nmols_nocharge_out(39,:) = (/7, 6/)
	nmols_nocharge_out(40,:) = (/7, 6/)
	nmols_nocharge_out(41,:) = (/6, 6/)
	nmols_nocharge_out(42,:) = (/8, 6/)
	nmols_nocharge_out(43,:) = (/8, 6/)
	nmols_nocharge_out(44,:) = (/7, 6/)
	nmols_nocharge_out(45,:) = (/4, 7/)
	nmols_nocharge_out(46,:) = (/3, 7/)
	nmols_nocharge_out(47,:) = (/5, 7/)
	nmols_nocharge_out(48,:) = (/5, 7/)
	nmols_nocharge_out(49,:) = (/4, 7/)
	nmols_nocharge_out(50,:) = (/6, 7/)
	nmols_nocharge_out(51,:) = (/6, 7/)
	nmols_nocharge_out(52,:) = (/5, 7/)
	nmols_nocharge_out(53,:) = (/7, 7/)
	nmols_nocharge_out(54,:) = (/7, 7/)
	nmols_nocharge_out(55,:) = (/6, 7/)
	nmols_nocharge_out(56,:) = (/8, 7/)
	nmols_nocharge_out(57,:) = (/8, 7/)
	nmols_nocharge_out(58,:) = (/7, 7/)
	nmols_nocharge_out(59,:) = (/4, 8/)
	nmols_nocharge_out(60,:) = (/3, 8/)
	nmols_nocharge_out(61,:) = (/5, 8/)
	nmols_nocharge_out(62,:) = (/5, 8/)
	nmols_nocharge_out(63,:) = (/4, 8/)
	nmols_nocharge_out(64,:) = (/6, 8/)
	nmols_nocharge_out(65,:) = (/6, 8/)
	nmols_nocharge_out(66,:) = (/5, 8/)
	nmols_nocharge_out(67,:) = (/7, 8/)
	nmols_nocharge_out(68,:) = (/7, 8/)
	nmols_nocharge_out(69,:) = (/6, 8/)
	nmols_nocharge_out(70,:) = (/8, 8/)
	nmols_nocharge_out(71,:) = (/8, 8/)
	nmols_nocharge_out(72,:) = (/7, 8/)

end subroutine get_nmols_nocharge_out_2


end module acdc_system_2

