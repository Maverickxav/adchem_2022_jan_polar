module acdc_system_1

implicit none

character(len=*), parameter :: tag = 'test'			! tag for the simulation system

logical, parameter :: small_set_mode = .true.			! small set of clusters (default mode)

integer, parameter :: nclust = 65						! number of clusters, molecules and ions
integer, parameter :: neq = 245							! number of equations
integer, parameter :: nclust_max = 65, neq_max = 245	! same as nclust and neq in the small set mode

integer, parameter :: n_variables = 3
logical, parameter :: variable_temp = .true.
logical, parameter :: variable_cs = .false.
logical, parameter :: variable_ipr = .true.
real(kind(1.d0)), parameter :: rh = 0.d0					! RH in %
real(kind(1.d0)), parameter :: pw = 0.d0					! partial pressure of water in Pa

integer, parameter :: n_monomer_types = 4
integer, parameter :: n_charges = 3			! number of charging states
integer, parameter :: n1A = 1, n1N = 4, n1B = 25, n1P1N = 47, nneg = 64, npos = 65			! cluster indices for monomers and ions

integer, parameter :: nout_all(3) = (/72, 73, 74/), nout_neu = 72, nout_neg = 73, nout_pos = 74			! indices for outgoing fluxes

integer, parameter :: nclust_coag_saved = 65			! number of saved scavenged clusters
integer, parameter :: nclust_out_saved = 105				! number of saved outgrown clusters
integer, parameter :: nrange_coag_saved(2) = (/76, 140/)	! index range of saved scavenged clusters
integer, parameter :: nrange_out_saved(2) = (/141, 245/)	! index range of saved outgrown clusters

integer, parameter :: n_mol_types = 4
integer, parameter :: nmolA = 1, nmolB = 2, nmolP = 3, nmolN = 4			! molecule indices for the used species

integer, parameter :: n_1A_clusters = 3				! number molecules and clusters containing 1 A molecule

integer, parameter :: n_neutrals = 24			! number of neutral molecules and clusters
integer, parameter :: n_negatives = 23			! negative
integer, parameter :: n_positives = 18			! positive
integer, parameter :: nclust_nogen = 63			! number of clusters and molecules excluding generic ions

integer, parameter :: n_neutral_monomers = 2, neutral_monomers(2) = (/1, 4/)			! number and indices of neutral monomers
integer, parameter :: n_negative_monomers = 1, negative_monomers(1) = (/25/)			! negative
integer, parameter :: n_positive_monomers = 1, positive_monomers(1) = (/47/)			! positive

integer, parameter :: n_neutral_clusters = 22			! number of neutral clusters
integer, parameter :: n_negative_clusters = 21			! negative
integer, parameter :: n_positive_clusters = 16			! positive

real(kind(1.d0)), parameter :: mass_max = 691.7200
real(kind(1.d0)), parameter :: diameter_max = 1.1411
real(kind(1.d0)), parameter :: mob_diameter_max = 1.4411
integer, parameter :: ij_ind_max(4) = (/6, 1, 1, 6/)		! maximum molecular content
integer, parameter :: n_bound = 11		! number of clusters at the system boundary
integer, parameter :: n_comp_out = 105		! number of clusters growing out in different compositions

integer, parameter :: nmols_out_neutral(1, 4) = reshape((/7, 0, 0, 6/),(/1, 4/))			! criteria for outgrowing neutrals
integer, parameter :: nmols_out_negative(1, 4) = reshape((/6, 1, 0, 3/),(/1, 4/))			! negatives
integer, parameter :: nmols_out_positive(1, 4) = reshape((/6, 0, 1, 7/),(/1, 4/))			! positives


contains

subroutine n_A_in_clusters_1(n_A)
	implicit none
	integer :: n_A(65)

	n_A = (/1, 2, 3, 0, 1, 2, 3, 1, 2, 3, &
		&4, 2, 3, 4, 5, 3, 4, 5, 6, 4, &
		&5, 6, 5, 6, 1, 2, 3, 4, 1, 2, &
		&3, 4, 2, 3, 4, 5, 3, 4, 5, 6, &
		&4, 5, 6, 5, 6, 6, 0, 1, 0, 1, &
		&2, 1, 2, 3, 2, 3, 4, 3, 4, 5, &
		&4, 5, 6, 0, 0/)

end subroutine n_A_in_clusters_1

subroutine clusters_with_1_A_1(cluster_numbers)
	implicit none
	integer :: cluster_numbers(3)

	cluster_numbers = (/1, 5, 8/)

end subroutine clusters_with_1_A_1

subroutine arrays_1(neutrals, negatives, positives)
	implicit none
	integer :: neutrals(24), negatives(23), positives(18)

	neutrals = (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, &
		&21, 22, 23, 24/)
	negatives = (/25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, &
		&40, 41, 42, 43, 44, 45, 46, 64/)
	positives = (/47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, &
		&62, 63, 65/)

end subroutine arrays_1

subroutine cluster_arrays_1(neutral_clusters, negative_clusters, positive_clusters)
	implicit none
	integer :: neutral_clusters(22), negative_clusters(21), positive_clusters(16)

	neutral_clusters = (/2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, &
		&21, 22, 23, 24/)
	negative_clusters = (/26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, &
		&40, 41, 42, 43, 44, 45, 46/)
	positive_clusters = (/48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, &
		&62, 63/)

end subroutine cluster_arrays_1

subroutine get_charging_state_1(charging_state)
	implicit none
	integer :: charging_state(65)

	charging_state = (/1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
		&1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
		&1, 1, 1, 1, 2, 2, 2, 2, 2, 2, &
		&2, 2, 2, 2, 2, 2, 2, 2, 2, 2, &
		&2, 2, 2, 2, 2, 2, 3, 3, 3, 3, &
		&3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
		&3, 3, 3, 2, 3/)

end subroutine get_charging_state_1

subroutine get_mass_1(mass)
	implicit none
	real(kind(1.d0)) :: mass(65)

	mass = (/98.0800, 196.1600, 294.2400, 17.0400, 115.1200, 213.2000, 311.2800, 132.1600, 230.2400, 328.3200, &
		&426.4000, 247.2800, 345.3600, 443.4400, 541.5200, 362.4000, 460.4800, 558.5600, 656.6400, 477.5200, &
		&575.6000, 673.6800, 592.6400, 690.7200, 97.0800, 195.1600, 293.2400, 391.3200, 114.1200, 212.2000, &
		&310.2800, 408.3600, 229.2400, 327.3200, 425.4000, 523.4800, 344.3600, 442.4400, 540.5200, 638.6000, &
		&459.4800, 557.5600, 655.6400, 574.6000, 672.6800, 689.7200, 18.0400, 116.1200, 35.0800, 133.1600, &
		&231.2400, 150.2000, 248.2800, 346.3600, 265.3200, 363.4000, 461.4800, 380.4400, 478.5200, 576.6000, &
		&495.5600, 593.6400, 691.7200, 32.0000, 19.0200/)

end subroutine get_mass_1

subroutine get_diameter_1(diameter)
	implicit none

	real(kind(1.d0)) :: diameter(65)

	 diameter = (/0.5539, 0.6979, 0.7989, 0.4266, 0.6280, 0.7475, 0.8376, 0.6877, 0.7912, 0.8729, &
		&0.9417, 0.8305, 0.9057, 0.9701, 1.0269, 0.9362, 0.9968, 1.0509, 1.0998, 1.0222, &
		&1.0738, 1.1208, 1.0958, 1.1411, 0.5520, 0.6967, 0.7980, 0.8786, 0.6265, 0.7464, &
		&0.8367, 0.9109, 0.7902, 0.8722, 0.9411, 1.0011, 0.9050, 0.9694, 1.0263, 1.0775, &
		&0.9962, 1.0503, 1.0994, 1.0733, 1.1204, 1.1406, 0.4266, 0.6280, 0.5375, 0.6877, &
		&0.7912, 0.7386, 0.8305, 0.9057, 0.8665, 0.9362, 0.9968, 0.9648, 1.0222, 1.0738, &
		&1.0464, 1.0958, 1.1411, 0.4464, 0.3926/)	! dry value

end subroutine get_diameter_1

subroutine get_mob_diameter_1(mob_diameter)
	implicit none

	real(kind(1.d0)) :: mob_diameter(65)

	 mob_diameter = (/0.8539, 0.9979, 1.0989, 0.7266, 0.9280, 1.0475, 1.1376, 0.9877, 1.0912, 1.1729, &
		&1.2417, 1.1305, 1.2057, 1.2701, 1.3269, 1.2362, 1.2968, 1.3509, 1.3998, 1.3222, &
		&1.3738, 1.4208, 1.3958, 1.4411, 0.8520, 0.9967, 1.0980, 1.1786, 0.9265, 1.0464, &
		&1.1367, 1.2109, 1.0902, 1.1722, 1.2411, 1.3011, 1.2050, 1.2694, 1.3263, 1.3775, &
		&1.2962, 1.3503, 1.3994, 1.3733, 1.4204, 1.4406, 0.7266, 0.9280, 0.8375, 0.9877, &
		&1.0912, 1.0386, 1.1305, 1.2057, 1.1665, 1.2362, 1.2968, 1.2648, 1.3222, 1.3738, &
		&1.3464, 1.3958, 1.4411, 0.7464, 0.6926/)	! dry value

end subroutine get_mob_diameter_1

subroutine cluster_names_1(clust)
	implicit none
	character(len=11), dimension(75) :: clust

	clust(1)(:) = '1A'
	clust(2)(:) = '2A'
	clust(3)(:) = '3A'
	clust(4)(:) = '1N'
	clust(5)(:) = '1A1N'
	clust(6)(:) = '2A1N'
	clust(7)(:) = '3A1N'
	clust(8)(:) = '1A2N'
	clust(9)(:) = '2A2N'
	clust(10)(:) = '3A2N'
	clust(11)(:) = '4A2N'
	clust(12)(:) = '2A3N'
	clust(13)(:) = '3A3N'
	clust(14)(:) = '4A3N'
	clust(15)(:) = '5A3N'
	clust(16)(:) = '3A4N'
	clust(17)(:) = '4A4N'
	clust(18)(:) = '5A4N'
	clust(19)(:) = '6A4N'
	clust(20)(:) = '4A5N'
	clust(21)(:) = '5A5N'
	clust(22)(:) = '6A5N'
	clust(23)(:) = '5A6N'
	clust(24)(:) = '6A6N'
	clust(25)(:) = '1B'
	clust(26)(:) = '1A1B'
	clust(27)(:) = '2A1B'
	clust(28)(:) = '3A1B'
	clust(29)(:) = '1B1N'
	clust(30)(:) = '1A1B1N'
	clust(31)(:) = '2A1B1N'
	clust(32)(:) = '3A1B1N'
	clust(33)(:) = '1A1B2N'
	clust(34)(:) = '2A1B2N'
	clust(35)(:) = '3A1B2N'
	clust(36)(:) = '4A1B2N'
	clust(37)(:) = '2A1B3N'
	clust(38)(:) = '3A1B3N'
	clust(39)(:) = '4A1B3N'
	clust(40)(:) = '5A1B3N'
	clust(41)(:) = '3A1B4N'
	clust(42)(:) = '4A1B4N'
	clust(43)(:) = '5A1B4N'
	clust(44)(:) = '4A1B5N'
	clust(45)(:) = '5A1B5N'
	clust(46)(:) = '5A1B6N'
	clust(47)(:) = '1P1N'
	clust(48)(:) = '1A1P1N'
	clust(49)(:) = '1P2N'
	clust(50)(:) = '1A1P2N'
	clust(51)(:) = '2A1P2N'
	clust(52)(:) = '1A1P3N'
	clust(53)(:) = '2A1P3N'
	clust(54)(:) = '3A1P3N'
	clust(55)(:) = '2A1P4N'
	clust(56)(:) = '3A1P4N'
	clust(57)(:) = '4A1P4N'
	clust(58)(:) = '3A1P5N'
	clust(59)(:) = '4A1P5N'
	clust(60)(:) = '5A1P5N'
	clust(61)(:) = '4A1P6N'
	clust(62)(:) = '5A1P6N'
	clust(63)(:) = '6A1P6N'
	clust(64)(:) = 'neg'
	clust(65)(:) = 'pos'
	clust(66)(:) = 'source'
	clust(67)(:) = 'coag'
	clust(68)(:) = 'wall'
	clust(69)(:) = 'dil'
	clust(70)(:) = 'insink'
	clust(71)(:) = 'rec'
	clust(72)(:) = 'out_neu'
	clust(73)(:) = 'out_neg'
	clust(74)(:) = 'out_pos'
	clust(75)(:) = 'bound'

	! Additional flux elements not listed here:
	! Elements from 76 to 140: coagulation loss of each cluster
	! Elements from 141 to 245: each outgrown cluster composition

end subroutine cluster_names_1

subroutine monomer_names_1(clust_mon)
	implicit none
	character(len=11), dimension(4) :: clust_mon

	clust_mon(1)(:) = '1A'
	clust_mon(2)(:) = '1N'
	clust_mon(3)(:) = '1B'
	clust_mon(4)(:) = '1P1N'

end subroutine monomer_names_1

subroutine molecule_names_1(labels)
	implicit none
	character(len=11), dimension(4) :: labels

	labels(1)(:) = 'A'
	labels(2)(:) = 'B'
	labels(3)(:) = 'P'
	labels(4)(:) = 'N'

end subroutine molecule_names_1

subroutine monomer_indices_1(n_monomers)
	implicit none
	integer :: n_monomers(4)

	n_monomers = (/1, 4, 25, 47/)

end subroutine monomer_indices_1

subroutine get_bound_1(bound_clusters,nmols_bound)
	implicit none
	integer :: bound_clusters(11), nmols_bound(11,4)

	nmols_bound(1,:) = (/6, 0, 0, 4/)
	nmols_bound(2,:) = (/6, 0, 0, 5/)
	nmols_bound(3,:) = (/5, 0, 0, 6/)
	nmols_bound(4,:) = (/6, 0, 0, 6/)
	nmols_bound(5,:) = (/5, 1, 0, 3/)
	nmols_bound(6,:) = (/5, 1, 0, 4/)
	nmols_bound(7,:) = (/5, 1, 0, 5/)
	nmols_bound(8,:) = (/5, 1, 0, 6/)
	nmols_bound(9,:) = (/4, 0, 1, 6/)
	nmols_bound(10,:) = (/5, 0, 1, 6/)
	nmols_bound(11,:) = (/6, 0, 1, 6/)

	bound_clusters = (/19, 22, 23, 24, 40, 43, 45, 46, 61, 62, 63/)

end subroutine get_bound_1

subroutine get_diameter_out_1(diameter_out)
	implicit none
	real(kind(1.d0)) :: diameter_out(105)

	diameter_out = (/1.1830, 1.1243, 1.1444, 1.1638, 1.1826, 1.2222, 1.1674, 1.1861, 1.2042, 1.2218, &
		&1.2590, 1.2076, 1.2251, 1.2421, 1.2586, 1.1606, 1.2012, 1.2008, 1.2012, 1.2393, &
		&1.2389, 1.2393, 1.2751, 1.2748, 1.2751, 1.2189, 1.2185, 1.1795, 1.2189, 1.2559, &
		&1.2555, 1.2559, 1.2908, 1.2905, 1.2908, 1.2938, 1.3091, 1.3240, 1.2778, 1.2934, &
		&1.3087, 1.3237, 1.3240, 1.2361, 1.2721, 1.2717, 1.1978, 1.2361, 1.2721, 1.3062, &
		&1.3059, 1.3062, 1.3386, 1.3383, 1.3386, 1.3413, 1.3556, 1.3695, 1.3264, 1.3410, &
		&1.3552, 1.3692, 1.3695, 1.2879, 1.3212, 1.3209, 1.2528, 1.2879, 1.3212, 1.3529, &
		&1.3526, 1.3529, 1.3832, 1.3829, 1.3832, 1.3857, 1.3991, 1.4122, 1.2617, 1.3718, &
		&1.3854, 1.3988, 1.4119, 1.4122, 1.3359, 1.3669, 1.3666, 1.3033, 1.3359, 1.3669, &
		&1.3966, 1.3963, 1.3966, 1.4250, 1.4247, 1.4250, 1.3806, 1.4097, 1.4094, 1.3502, &
		&1.3806, 1.4097, 1.4377, 1.4374, 1.4377/)	! dry value

end subroutine get_diameter_out_1

subroutine get_nmols_out_1(nmols_out)
	implicit none
	integer :: nmols_out(105,4)

	nmols_out(1,:) = (/7, 0, 0, 6/)
	nmols_out(2,:) = (/6, 1, 0, 3/)
	nmols_out(3,:) = (/6, 1, 0, 4/)
	nmols_out(4,:) = (/6, 1, 0, 5/)
	nmols_out(5,:) = (/6, 1, 0, 6/)
	nmols_out(6,:) = (/8, 0, 0, 6/)
	nmols_out(7,:) = (/7, 1, 0, 3/)
	nmols_out(8,:) = (/7, 1, 0, 4/)
	nmols_out(9,:) = (/7, 1, 0, 5/)
	nmols_out(10,:) = (/7, 1, 0, 6/)
	nmols_out(11,:) = (/9, 0, 0, 6/)
	nmols_out(12,:) = (/8, 1, 0, 3/)
	nmols_out(13,:) = (/8, 1, 0, 4/)
	nmols_out(14,:) = (/8, 1, 0, 5/)
	nmols_out(15,:) = (/8, 1, 0, 6/)
	nmols_out(16,:) = (/6, 0, 1, 7/)
	nmols_out(17,:) = (/7, 0, 0, 7/)
	nmols_out(18,:) = (/6, 1, 0, 7/)
	nmols_out(19,:) = (/7, 0, 1, 7/)
	nmols_out(20,:) = (/8, 0, 0, 7/)
	nmols_out(21,:) = (/7, 1, 0, 7/)
	nmols_out(22,:) = (/8, 0, 1, 7/)
	nmols_out(23,:) = (/9, 0, 0, 7/)
	nmols_out(24,:) = (/8, 1, 0, 7/)
	nmols_out(25,:) = (/9, 0, 1, 7/)
	nmols_out(26,:) = (/7, 0, 0, 8/)
	nmols_out(27,:) = (/6, 1, 0, 8/)
	nmols_out(28,:) = (/6, 0, 1, 8/)
	nmols_out(29,:) = (/7, 0, 1, 8/)
	nmols_out(30,:) = (/8, 0, 0, 8/)
	nmols_out(31,:) = (/7, 1, 0, 8/)
	nmols_out(32,:) = (/8, 0, 1, 8/)
	nmols_out(33,:) = (/9, 0, 0, 8/)
	nmols_out(34,:) = (/8, 1, 0, 8/)
	nmols_out(35,:) = (/9, 0, 1, 8/)
	nmols_out(36,:) = (/10, 0, 0, 6/)
	nmols_out(37,:) = (/10, 0, 0, 7/)
	nmols_out(38,:) = (/10, 0, 0, 8/)
	nmols_out(39,:) = (/9, 1, 0, 5/)
	nmols_out(40,:) = (/9, 1, 0, 6/)
	nmols_out(41,:) = (/9, 1, 0, 7/)
	nmols_out(42,:) = (/9, 1, 0, 8/)
	nmols_out(43,:) = (/10, 0, 1, 8/)
	nmols_out(44,:) = (/7, 0, 0, 9/)
	nmols_out(45,:) = (/8, 0, 0, 9/)
	nmols_out(46,:) = (/7, 1, 0, 9/)
	nmols_out(47,:) = (/6, 0, 1, 9/)
	nmols_out(48,:) = (/7, 0, 1, 9/)
	nmols_out(49,:) = (/8, 0, 1, 9/)
	nmols_out(50,:) = (/9, 0, 0, 9/)
	nmols_out(51,:) = (/8, 1, 0, 9/)
	nmols_out(52,:) = (/9, 0, 1, 9/)
	nmols_out(53,:) = (/10, 0, 0, 9/)
	nmols_out(54,:) = (/9, 1, 0, 9/)
	nmols_out(55,:) = (/10, 0, 1, 9/)
	nmols_out(56,:) = (/11, 0, 0, 7/)
	nmols_out(57,:) = (/11, 0, 0, 8/)
	nmols_out(58,:) = (/11, 0, 0, 9/)
	nmols_out(59,:) = (/10, 1, 0, 6/)
	nmols_out(60,:) = (/10, 1, 0, 7/)
	nmols_out(61,:) = (/10, 1, 0, 8/)
	nmols_out(62,:) = (/10, 1, 0, 9/)
	nmols_out(63,:) = (/11, 0, 1, 9/)
	nmols_out(64,:) = (/8, 0, 0, 10/)
	nmols_out(65,:) = (/9, 0, 0, 10/)
	nmols_out(66,:) = (/8, 1, 0, 10/)
	nmols_out(67,:) = (/7, 0, 1, 10/)
	nmols_out(68,:) = (/8, 0, 1, 10/)
	nmols_out(69,:) = (/9, 0, 1, 10/)
	nmols_out(70,:) = (/10, 0, 0, 10/)
	nmols_out(71,:) = (/9, 1, 0, 10/)
	nmols_out(72,:) = (/10, 0, 1, 10/)
	nmols_out(73,:) = (/11, 0, 0, 10/)
	nmols_out(74,:) = (/10, 1, 0, 10/)
	nmols_out(75,:) = (/11, 0, 1, 10/)
	nmols_out(76,:) = (/12, 0, 0, 8/)
	nmols_out(77,:) = (/12, 0, 0, 9/)
	nmols_out(78,:) = (/12, 0, 0, 10/)
	nmols_out(79,:) = (/9, 1, 0, 4/)
	nmols_out(80,:) = (/11, 1, 0, 7/)
	nmols_out(81,:) = (/11, 1, 0, 8/)
	nmols_out(82,:) = (/11, 1, 0, 9/)
	nmols_out(83,:) = (/11, 1, 0, 10/)
	nmols_out(84,:) = (/12, 0, 1, 10/)
	nmols_out(85,:) = (/9, 0, 0, 11/)
	nmols_out(86,:) = (/10, 0, 0, 11/)
	nmols_out(87,:) = (/9, 1, 0, 11/)
	nmols_out(88,:) = (/8, 0, 1, 11/)
	nmols_out(89,:) = (/9, 0, 1, 11/)
	nmols_out(90,:) = (/10, 0, 1, 11/)
	nmols_out(91,:) = (/11, 0, 0, 11/)
	nmols_out(92,:) = (/10, 1, 0, 11/)
	nmols_out(93,:) = (/11, 0, 1, 11/)
	nmols_out(94,:) = (/12, 0, 0, 11/)
	nmols_out(95,:) = (/11, 1, 0, 11/)
	nmols_out(96,:) = (/12, 0, 1, 11/)
	nmols_out(97,:) = (/10, 0, 0, 12/)
	nmols_out(98,:) = (/11, 0, 0, 12/)
	nmols_out(99,:) = (/10, 1, 0, 12/)
	nmols_out(100,:) = (/9, 0, 1, 12/)
	nmols_out(101,:) = (/10, 0, 1, 12/)
	nmols_out(102,:) = (/11, 0, 1, 12/)
	nmols_out(103,:) = (/12, 0, 0, 12/)
	nmols_out(104,:) = (/11, 1, 0, 12/)
	nmols_out(105,:) = (/12, 0, 1, 12/)

end subroutine get_nmols_out_1

subroutine molecule_names_nocharge_1(labels_nocharge)
	implicit none
	character(len=11), dimension(2) :: labels_nocharge

	labels_nocharge(1)(:) = 'A'
	labels_nocharge(2)(:) = 'N'

end subroutine molecule_names_nocharge_1

subroutine get_nmols_nocharge_1(nmols_nocharge_clust)
	implicit none
	integer :: nmols_nocharge_clust(65,2)

	nmols_nocharge_clust = 0

	nmols_nocharge_clust(1,:) = (/1, 0/)
	nmols_nocharge_clust(2,:) = (/2, 0/)
	nmols_nocharge_clust(3,:) = (/3, 0/)
	nmols_nocharge_clust(4,:) = (/0, 1/)
	nmols_nocharge_clust(5,:) = (/1, 1/)
	nmols_nocharge_clust(6,:) = (/2, 1/)
	nmols_nocharge_clust(7,:) = (/3, 1/)
	nmols_nocharge_clust(8,:) = (/1, 2/)
	nmols_nocharge_clust(9,:) = (/2, 2/)
	nmols_nocharge_clust(10,:) = (/3, 2/)
	nmols_nocharge_clust(11,:) = (/4, 2/)
	nmols_nocharge_clust(12,:) = (/2, 3/)
	nmols_nocharge_clust(13,:) = (/3, 3/)
	nmols_nocharge_clust(14,:) = (/4, 3/)
	nmols_nocharge_clust(15,:) = (/5, 3/)
	nmols_nocharge_clust(16,:) = (/3, 4/)
	nmols_nocharge_clust(17,:) = (/4, 4/)
	nmols_nocharge_clust(18,:) = (/5, 4/)
	nmols_nocharge_clust(19,:) = (/6, 4/)
	nmols_nocharge_clust(20,:) = (/4, 5/)
	nmols_nocharge_clust(21,:) = (/5, 5/)
	nmols_nocharge_clust(22,:) = (/6, 5/)
	nmols_nocharge_clust(23,:) = (/5, 6/)
	nmols_nocharge_clust(24,:) = (/6, 6/)
	nmols_nocharge_clust(25,:) = (/1, 0/)
	nmols_nocharge_clust(26,:) = (/2, 0/)
	nmols_nocharge_clust(27,:) = (/3, 0/)
	nmols_nocharge_clust(28,:) = (/4, 0/)
	nmols_nocharge_clust(29,:) = (/1, 1/)
	nmols_nocharge_clust(30,:) = (/2, 1/)
	nmols_nocharge_clust(31,:) = (/3, 1/)
	nmols_nocharge_clust(32,:) = (/4, 1/)
	nmols_nocharge_clust(33,:) = (/2, 2/)
	nmols_nocharge_clust(34,:) = (/3, 2/)
	nmols_nocharge_clust(35,:) = (/4, 2/)
	nmols_nocharge_clust(36,:) = (/5, 2/)
	nmols_nocharge_clust(37,:) = (/3, 3/)
	nmols_nocharge_clust(38,:) = (/4, 3/)
	nmols_nocharge_clust(39,:) = (/5, 3/)
	nmols_nocharge_clust(40,:) = (/6, 3/)
	nmols_nocharge_clust(41,:) = (/4, 4/)
	nmols_nocharge_clust(42,:) = (/5, 4/)
	nmols_nocharge_clust(43,:) = (/6, 4/)
	nmols_nocharge_clust(44,:) = (/5, 5/)
	nmols_nocharge_clust(45,:) = (/6, 5/)
	nmols_nocharge_clust(46,:) = (/6, 6/)
	nmols_nocharge_clust(47,:) = (/0, 1/)
	nmols_nocharge_clust(48,:) = (/1, 1/)
	nmols_nocharge_clust(49,:) = (/0, 2/)
	nmols_nocharge_clust(50,:) = (/1, 2/)
	nmols_nocharge_clust(51,:) = (/2, 2/)
	nmols_nocharge_clust(52,:) = (/1, 3/)
	nmols_nocharge_clust(53,:) = (/2, 3/)
	nmols_nocharge_clust(54,:) = (/3, 3/)
	nmols_nocharge_clust(55,:) = (/2, 4/)
	nmols_nocharge_clust(56,:) = (/3, 4/)
	nmols_nocharge_clust(57,:) = (/4, 4/)
	nmols_nocharge_clust(58,:) = (/3, 5/)
	nmols_nocharge_clust(59,:) = (/4, 5/)
	nmols_nocharge_clust(60,:) = (/5, 5/)
	nmols_nocharge_clust(61,:) = (/4, 6/)
	nmols_nocharge_clust(62,:) = (/5, 6/)
	nmols_nocharge_clust(63,:) = (/6, 6/)

end subroutine get_nmols_nocharge_1

subroutine get_nmols_nocharge_out_1(nmols_nocharge_out)
	implicit none
	integer :: nmols_nocharge_out(105,2)

	nmols_nocharge_out = 0

	nmols_nocharge_out(1,:) = (/7, 6/)
	nmols_nocharge_out(2,:) = (/7, 3/)
	nmols_nocharge_out(3,:) = (/7, 4/)
	nmols_nocharge_out(4,:) = (/7, 5/)
	nmols_nocharge_out(5,:) = (/7, 6/)
	nmols_nocharge_out(6,:) = (/8, 6/)
	nmols_nocharge_out(7,:) = (/8, 3/)
	nmols_nocharge_out(8,:) = (/8, 4/)
	nmols_nocharge_out(9,:) = (/8, 5/)
	nmols_nocharge_out(10,:) = (/8, 6/)
	nmols_nocharge_out(11,:) = (/9, 6/)
	nmols_nocharge_out(12,:) = (/9, 3/)
	nmols_nocharge_out(13,:) = (/9, 4/)
	nmols_nocharge_out(14,:) = (/9, 5/)
	nmols_nocharge_out(15,:) = (/9, 6/)
	nmols_nocharge_out(16,:) = (/6, 7/)
	nmols_nocharge_out(17,:) = (/7, 7/)
	nmols_nocharge_out(18,:) = (/7, 7/)
	nmols_nocharge_out(19,:) = (/7, 7/)
	nmols_nocharge_out(20,:) = (/8, 7/)
	nmols_nocharge_out(21,:) = (/8, 7/)
	nmols_nocharge_out(22,:) = (/8, 7/)
	nmols_nocharge_out(23,:) = (/9, 7/)
	nmols_nocharge_out(24,:) = (/9, 7/)
	nmols_nocharge_out(25,:) = (/9, 7/)
	nmols_nocharge_out(26,:) = (/7, 8/)
	nmols_nocharge_out(27,:) = (/7, 8/)
	nmols_nocharge_out(28,:) = (/6, 8/)
	nmols_nocharge_out(29,:) = (/7, 8/)
	nmols_nocharge_out(30,:) = (/8, 8/)
	nmols_nocharge_out(31,:) = (/8, 8/)
	nmols_nocharge_out(32,:) = (/8, 8/)
	nmols_nocharge_out(33,:) = (/9, 8/)
	nmols_nocharge_out(34,:) = (/9, 8/)
	nmols_nocharge_out(35,:) = (/9, 8/)
	nmols_nocharge_out(36,:) = (/10, 6/)
	nmols_nocharge_out(37,:) = (/10, 7/)
	nmols_nocharge_out(38,:) = (/10, 8/)
	nmols_nocharge_out(39,:) = (/10, 5/)
	nmols_nocharge_out(40,:) = (/10, 6/)
	nmols_nocharge_out(41,:) = (/10, 7/)
	nmols_nocharge_out(42,:) = (/10, 8/)
	nmols_nocharge_out(43,:) = (/10, 8/)
	nmols_nocharge_out(44,:) = (/7, 9/)
	nmols_nocharge_out(45,:) = (/8, 9/)
	nmols_nocharge_out(46,:) = (/8, 9/)
	nmols_nocharge_out(47,:) = (/6, 9/)
	nmols_nocharge_out(48,:) = (/7, 9/)
	nmols_nocharge_out(49,:) = (/8, 9/)
	nmols_nocharge_out(50,:) = (/9, 9/)
	nmols_nocharge_out(51,:) = (/9, 9/)
	nmols_nocharge_out(52,:) = (/9, 9/)
	nmols_nocharge_out(53,:) = (/10, 9/)
	nmols_nocharge_out(54,:) = (/10, 9/)
	nmols_nocharge_out(55,:) = (/10, 9/)
	nmols_nocharge_out(56,:) = (/11, 7/)
	nmols_nocharge_out(57,:) = (/11, 8/)
	nmols_nocharge_out(58,:) = (/11, 9/)
	nmols_nocharge_out(59,:) = (/11, 6/)
	nmols_nocharge_out(60,:) = (/11, 7/)
	nmols_nocharge_out(61,:) = (/11, 8/)
	nmols_nocharge_out(62,:) = (/11, 9/)
	nmols_nocharge_out(63,:) = (/11, 9/)
	nmols_nocharge_out(64,:) = (/8, 10/)
	nmols_nocharge_out(65,:) = (/9, 10/)
	nmols_nocharge_out(66,:) = (/9, 10/)
	nmols_nocharge_out(67,:) = (/7, 10/)
	nmols_nocharge_out(68,:) = (/8, 10/)
	nmols_nocharge_out(69,:) = (/9, 10/)
	nmols_nocharge_out(70,:) = (/10, 10/)
	nmols_nocharge_out(71,:) = (/10, 10/)
	nmols_nocharge_out(72,:) = (/10, 10/)
	nmols_nocharge_out(73,:) = (/11, 10/)
	nmols_nocharge_out(74,:) = (/11, 10/)
	nmols_nocharge_out(75,:) = (/11, 10/)
	nmols_nocharge_out(76,:) = (/12, 8/)
	nmols_nocharge_out(77,:) = (/12, 9/)
	nmols_nocharge_out(78,:) = (/12, 10/)
	nmols_nocharge_out(79,:) = (/10, 4/)
	nmols_nocharge_out(80,:) = (/12, 7/)
	nmols_nocharge_out(81,:) = (/12, 8/)
	nmols_nocharge_out(82,:) = (/12, 9/)
	nmols_nocharge_out(83,:) = (/12, 10/)
	nmols_nocharge_out(84,:) = (/12, 10/)
	nmols_nocharge_out(85,:) = (/9, 11/)
	nmols_nocharge_out(86,:) = (/10, 11/)
	nmols_nocharge_out(87,:) = (/10, 11/)
	nmols_nocharge_out(88,:) = (/8, 11/)
	nmols_nocharge_out(89,:) = (/9, 11/)
	nmols_nocharge_out(90,:) = (/10, 11/)
	nmols_nocharge_out(91,:) = (/11, 11/)
	nmols_nocharge_out(92,:) = (/11, 11/)
	nmols_nocharge_out(93,:) = (/11, 11/)
	nmols_nocharge_out(94,:) = (/12, 11/)
	nmols_nocharge_out(95,:) = (/12, 11/)
	nmols_nocharge_out(96,:) = (/12, 11/)
	nmols_nocharge_out(97,:) = (/10, 12/)
	nmols_nocharge_out(98,:) = (/11, 12/)
	nmols_nocharge_out(99,:) = (/11, 12/)
	nmols_nocharge_out(100,:) = (/9, 12/)
	nmols_nocharge_out(101,:) = (/10, 12/)
	nmols_nocharge_out(102,:) = (/11, 12/)
	nmols_nocharge_out(103,:) = (/12, 12/)
	nmols_nocharge_out(104,:) = (/12, 12/)
	nmols_nocharge_out(105,:) = (/12, 12/)

end subroutine get_nmols_nocharge_out_1


end module acdc_system_1

