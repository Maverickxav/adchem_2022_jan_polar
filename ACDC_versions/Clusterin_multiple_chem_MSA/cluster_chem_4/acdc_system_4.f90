module acdc_system_4

implicit none

character(len=*), parameter :: tag = 'test'			! tag for the simulation system

logical, parameter :: small_set_mode = .true.			! small set of clusters (default mode)

integer, parameter :: nclust = 58						! number of clusters, molecules and ions
integer, parameter :: neq = 222							! number of equations
integer, parameter :: nclust_max = 58, neq_max = 222	! same as nclust and neq in the small set mode

integer, parameter :: n_variables = 1
logical, parameter :: variable_temp = .true.
logical, parameter :: variable_cs = .false.
logical, parameter :: variable_ipr = .false.
real(kind(1.d0)), parameter :: rh = 0.d0					! RH in %
real(kind(1.d0)), parameter :: pw = 0.d0					! partial pressure of water in Pa

integer, parameter :: n_monomer_types = 3
integer, parameter :: n_charges = 1			! number of charging states
integer, parameter :: n1A = 1, n1D = 5, n1Ms = 25			! cluster indices for monomers and ions

integer, parameter :: nout_all(1) = (/65/), nout_neu = 65			! indices for outgoing fluxes

integer, parameter :: nclust_coag_saved = 58			! number of saved scavenged clusters
integer, parameter :: nclust_out_saved = 98				! number of saved outgrown clusters
integer, parameter :: nrange_coag_saved(2) = (/67, 124/)	! index range of saved scavenged clusters
integer, parameter :: nrange_out_saved(2) = (/125, 222/)	! index range of saved outgrown clusters

integer, parameter :: n_mol_types = 3
integer, parameter :: nmolA = 1, nmolMs = 2, nmolD = 3			! molecule indices for the used species

integer, parameter :: n_1A_clusters = 18				! number molecules and clusters containing 1 A molecule

integer, parameter :: n_neutrals = 58			! number of neutral molecules and clusters
integer, parameter :: n_negatives = 0			! negative
integer, parameter :: n_positives = 0			! positive
integer, parameter :: nclust_nogen = 58			! number of clusters and molecules excluding generic ions

integer, parameter :: n_neutral_monomers = 3, neutral_monomers(3) = (/1, 5, 25/)			! number and indices of neutral monomers
integer, parameter :: n_negative_monomers = 0			! negative
integer, parameter :: n_positive_monomers = 0			! positive

integer, parameter :: n_neutral_clusters = 55			! number of neutral clusters
integer, parameter :: n_negative_clusters = 0			! negative
integer, parameter :: n_positive_clusters = 0			! positive

real(kind(1.d0)), parameter :: mass_max = 582.5700
real(kind(1.d0)), parameter :: diameter_max = 1.1765
real(kind(1.d0)), parameter :: mob_diameter_max = 1.4765
integer, parameter :: ij_ind_max(3) = (/4, 4, 4/)		! maximum molecular content
integer, parameter :: n_bound = 13		! number of clusters at the system boundary
integer, parameter :: n_comp_out = 98		! number of clusters growing out in different compositions

integer, parameter :: nmols_out_neutral(3, 3) = reshape((/5, 4, 3, 0, 1, 2, 4, 4, 4/),(/3, 3/))			! criteria for outgrowing neutrals


contains

subroutine n_A_in_clusters_4(n_A)
	implicit none
	integer :: n_A(58)

	n_A = (/1, 2, 3, 4, 0, 1, 2, 3, 4, 0, &
		&1, 2, 3, 4, 0, 1, 2, 3, 4, 0, &
		&1, 2, 3, 4, 0, 1, 2, 3, 0, 1, &
		&2, 3, 0, 1, 2, 3, 0, 1, 2, 3, &
		&1, 2, 1, 1, 2, 3, 1, 2, 1, 1, &
		&2, 3, 1, 2, 1, 3, 2, 1/)

end subroutine n_A_in_clusters_4

subroutine clusters_with_1_A_4(cluster_numbers)
	implicit none
	integer :: cluster_numbers(18)

	cluster_numbers = (/1, 6, 11, 16, 21, 26, 30, 34, 38, 41, 43, 44, 47, 49, 50, &
		&53, 55, 58/)

end subroutine clusters_with_1_A_4

subroutine cluster_arrays_4(neutral_clusters)
	implicit none
	integer :: neutral_clusters(55)

	neutral_clusters = (/2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, &
		&21, 22, 23, 24, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, &
		&41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58/)

end subroutine cluster_arrays_4

subroutine get_charging_state_4(charging_state)
	implicit none
	integer :: charging_state(58)

	charging_state = (/1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
		&1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
		&1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
		&1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
		&1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
		&1, 1, 1, 1, 1, 1, 1, 1/)

end subroutine get_charging_state_4

subroutine get_mass_4(mass)
	implicit none
	real(kind(1.d0)) :: mass(58)

	mass = (/98.0800, 196.1600, 294.2400, 392.3200, 45.0800, 143.1600, 241.2400, 339.3200, 437.4000, 90.1600, &
		&188.2400, 286.3200, 384.4000, 482.4800, 135.2400, 233.3200, 331.4000, 429.4800, 527.5600, 180.3200, &
		&278.4000, 376.4800, 474.5600, 572.6400, 96.1100, 194.1900, 292.2700, 390.3500, 192.2200, 290.3000, &
		&388.3800, 486.4600, 288.3300, 386.4100, 484.4900, 582.5700, 384.4400, 239.2700, 337.3500, 435.4300, &
		&335.3800, 433.4600, 431.4900, 284.3500, 382.4300, 480.5100, 380.4600, 478.5400, 476.5700, 329.4300, &
		&427.5100, 525.5900, 425.5400, 523.6200, 521.6500, 570.6700, 568.7000, 566.7300/)

end subroutine get_mass_4

subroutine get_diameter_4(diameter)
	implicit none

	real(kind(1.d0)) :: diameter(58)

	 diameter = (/0.5539, 0.6979, 0.7989, 0.8793, 0.5946, 0.7245, 0.8194, 0.8963, 0.9619, 0.7492, &
		&0.8389, 0.9128, 0.9762, 1.0324, 0.8576, 0.9286, 0.9901, 1.0448, 1.0944, 0.9439, &
		&1.0036, 1.0570, 1.1055, 1.1500, 0.5904, 0.7216, 0.8172, 0.8945, 0.7439, 0.8347, &
		&0.9092, 0.9731, 0.8515, 0.9234, 0.9856, 1.0408, 0.9372, 0.8368, 0.9110, 0.9747, &
		&0.9251, 0.9871, 0.9992, 0.9269, 0.9886, 1.0435, 1.0007, 1.0543, 1.0650, 1.0022, &
		&1.0557, 1.1042, 1.0663, 1.1140, 1.1235, 1.1590, 1.1678, 1.1765/)	! dry value

end subroutine get_diameter_4

subroutine get_mob_diameter_4(mob_diameter)
	implicit none

	real(kind(1.d0)) :: mob_diameter(58)

	 mob_diameter = (/0.8539, 0.9979, 1.0989, 1.1793, 0.8946, 1.0245, 1.1194, 1.1963, 1.2619, 1.0492, &
		&1.1389, 1.2128, 1.2762, 1.3324, 1.1576, 1.2286, 1.2901, 1.3448, 1.3944, 1.2439, &
		&1.3036, 1.3570, 1.4055, 1.4500, 0.8904, 1.0216, 1.1172, 1.1945, 1.0439, 1.1347, &
		&1.2092, 1.2731, 1.1515, 1.2234, 1.2856, 1.3408, 1.2372, 1.1368, 1.2110, 1.2747, &
		&1.2251, 1.2871, 1.2992, 1.2269, 1.2886, 1.3435, 1.3007, 1.3543, 1.3650, 1.3022, &
		&1.3557, 1.4042, 1.3663, 1.4140, 1.4235, 1.4590, 1.4678, 1.4765/)	! dry value

end subroutine get_mob_diameter_4

subroutine cluster_names_4(clust)
	implicit none
	character(len=11), dimension(66) :: clust

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
	clust(25)(:) = '1Ms'
	clust(26)(:) = '1A1Ms'
	clust(27)(:) = '2A1Ms'
	clust(28)(:) = '3A1Ms'
	clust(29)(:) = '2Ms'
	clust(30)(:) = '1A2Ms'
	clust(31)(:) = '2A2Ms'
	clust(32)(:) = '3A2Ms'
	clust(33)(:) = '3Ms'
	clust(34)(:) = '1A3Ms'
	clust(35)(:) = '2A3Ms'
	clust(36)(:) = '3A3Ms'
	clust(37)(:) = '4Ms'
	clust(38)(:) = '1A1Ms1D'
	clust(39)(:) = '2A1Ms1D'
	clust(40)(:) = '3A1Ms1D'
	clust(41)(:) = '1A2Ms1D'
	clust(42)(:) = '2A2Ms1D'
	clust(43)(:) = '1A3Ms1D'
	clust(44)(:) = '1A1Ms2D'
	clust(45)(:) = '2A1Ms2D'
	clust(46)(:) = '3A1Ms2D'
	clust(47)(:) = '1A2Ms2D'
	clust(48)(:) = '2A2Ms2D'
	clust(49)(:) = '1A3Ms2D'
	clust(50)(:) = '1A1Ms3D'
	clust(51)(:) = '2A1Ms3D'
	clust(52)(:) = '3A1Ms3D'
	clust(53)(:) = '1A2Ms3D'
	clust(54)(:) = '2A2Ms3D'
	clust(55)(:) = '1A3Ms3D'
	clust(56)(:) = '3A1Ms4D'
	clust(57)(:) = '2A2Ms4D'
	clust(58)(:) = '1A3Ms4D'
	clust(59)(:) = 'source'
	clust(60)(:) = 'coag'
	clust(61)(:) = 'wall'
	clust(62)(:) = 'dil'
	clust(63)(:) = 'insink'
	clust(64)(:) = 'rec'
	clust(65)(:) = 'out_neu'
	clust(66)(:) = 'bound'

	! Additional flux elements not listed here:
	! Elements from 67 to 124: coagulation loss of each cluster
	! Elements from 125 to 222: each outgrown cluster composition

end subroutine cluster_names_4

subroutine monomer_names_4(clust_mon)
	implicit none
	character(len=11), dimension(3) :: clust_mon

	clust_mon(1)(:) = '1A'
	clust_mon(2)(:) = '1D'
	clust_mon(3)(:) = '1Ms'

end subroutine monomer_names_4

subroutine molecule_names_4(labels)
	implicit none
	character(len=11), dimension(3) :: labels

	labels(1)(:) = 'A'
	labels(2)(:) = 'Ms'
	labels(3)(:) = 'D'

end subroutine molecule_names_4

subroutine monomer_indices_4(n_monomers)
	implicit none
	integer :: n_monomers(3)

	n_monomers = (/1, 5, 25/)

end subroutine monomer_indices_4

subroutine get_bound_4(bound_clusters,nmols_bound)
	implicit none
	integer :: bound_clusters(13), nmols_bound(13,3)

	nmols_bound(1,:) = (/4, 0, 0/)
	nmols_bound(2,:) = (/4, 0, 1/)
	nmols_bound(3,:) = (/4, 0, 2/)
	nmols_bound(4,:) = (/4, 0, 3/)
	nmols_bound(5,:) = (/0, 0, 4/)
	nmols_bound(6,:) = (/1, 0, 4/)
	nmols_bound(7,:) = (/2, 0, 4/)
	nmols_bound(8,:) = (/3, 0, 4/)
	nmols_bound(9,:) = (/4, 0, 4/)
	nmols_bound(10,:) = (/0, 4, 0/)
	nmols_bound(11,:) = (/3, 1, 4/)
	nmols_bound(12,:) = (/2, 2, 4/)
	nmols_bound(13,:) = (/1, 3, 4/)

	bound_clusters = (/4, 9, 14, 19, 20, 21, 22, 23, 24, 37, 56, 57, 58/)

end subroutine get_bound_4

subroutine get_diameter_out_4(diameter_out)
	implicit none
	real(kind(1.d0)) :: diameter_out(98)

	diameter_out = (/1.1913, 1.1997, 1.2079, 1.2300, 1.2378, 1.2456, 1.2532, 1.2664, 1.2738, 1.2811, &
		&1.2883, 1.3007, 1.3078, 1.3147, 1.3216, 1.2388, 1.2465, 1.2542, 1.2747, 1.2820, &
		&1.2892, 1.2964, 1.3086, 1.3156, 1.3224, 1.3292, 1.3409, 1.3475, 1.3541, 1.3606, &
		&1.2829, 1.2901, 1.2972, 1.3164, 1.3233, 1.3301, 1.3368, 1.3483, 1.3549, 1.3613, &
		&1.3678, 1.3788, 1.3851, 1.3913, 1.3974, 1.3241, 1.3309, 1.3376, 1.3557, 1.3621, &
		&1.3686, 1.3749, 1.3858, 1.3920, 1.3982, 1.4042, 1.4147, 1.4207, 1.4266, 1.4324, &
		&1.3629, 1.3693, 1.3757, 1.3928, 1.3989, 1.4050, 1.4110, 1.4214, 1.3533, 1.2955, &
		&1.4273, 1.4331, 1.4389, 1.4489, 1.3467, 1.3835, 1.3284, 1.4546, 1.4602, 1.4658, &
		&1.3351, 1.3598, 1.3662, 1.3726, 1.3897, 1.3959, 1.4020, 1.3360, 1.3670, 1.3733, &
		&1.3741, 1.4035, 1.4095, 1.4103, 1.4382, 1.4440, 1.4713, 1.4768/)	! dry value

end subroutine get_diameter_out_4

subroutine get_nmols_out_4(nmols_out)
	implicit none
	integer :: nmols_out(98,3)

	nmols_out(1,:) = (/5, 0, 4/)
	nmols_out(2,:) = (/4, 1, 4/)
	nmols_out(3,:) = (/3, 2, 4/)
	nmols_out(4,:) = (/6, 0, 4/)
	nmols_out(5,:) = (/5, 1, 4/)
	nmols_out(6,:) = (/4, 2, 4/)
	nmols_out(7,:) = (/3, 3, 4/)
	nmols_out(8,:) = (/7, 0, 4/)
	nmols_out(9,:) = (/6, 1, 4/)
	nmols_out(10,:) = (/5, 2, 4/)
	nmols_out(11,:) = (/4, 3, 4/)
	nmols_out(12,:) = (/8, 0, 4/)
	nmols_out(13,:) = (/7, 1, 4/)
	nmols_out(14,:) = (/6, 2, 4/)
	nmols_out(15,:) = (/5, 3, 4/)
	nmols_out(16,:) = (/5, 0, 5/)
	nmols_out(17,:) = (/4, 1, 5/)
	nmols_out(18,:) = (/3, 2, 5/)
	nmols_out(19,:) = (/6, 0, 5/)
	nmols_out(20,:) = (/5, 1, 5/)
	nmols_out(21,:) = (/4, 2, 5/)
	nmols_out(22,:) = (/3, 3, 5/)
	nmols_out(23,:) = (/7, 0, 5/)
	nmols_out(24,:) = (/6, 1, 5/)
	nmols_out(25,:) = (/5, 2, 5/)
	nmols_out(26,:) = (/4, 3, 5/)
	nmols_out(27,:) = (/8, 0, 5/)
	nmols_out(28,:) = (/7, 1, 5/)
	nmols_out(29,:) = (/6, 2, 5/)
	nmols_out(30,:) = (/5, 3, 5/)
	nmols_out(31,:) = (/5, 0, 6/)
	nmols_out(32,:) = (/4, 1, 6/)
	nmols_out(33,:) = (/3, 2, 6/)
	nmols_out(34,:) = (/6, 0, 6/)
	nmols_out(35,:) = (/5, 1, 6/)
	nmols_out(36,:) = (/4, 2, 6/)
	nmols_out(37,:) = (/3, 3, 6/)
	nmols_out(38,:) = (/7, 0, 6/)
	nmols_out(39,:) = (/6, 1, 6/)
	nmols_out(40,:) = (/5, 2, 6/)
	nmols_out(41,:) = (/4, 3, 6/)
	nmols_out(42,:) = (/8, 0, 6/)
	nmols_out(43,:) = (/7, 1, 6/)
	nmols_out(44,:) = (/6, 2, 6/)
	nmols_out(45,:) = (/5, 3, 6/)
	nmols_out(46,:) = (/5, 0, 7/)
	nmols_out(47,:) = (/4, 1, 7/)
	nmols_out(48,:) = (/3, 2, 7/)
	nmols_out(49,:) = (/6, 0, 7/)
	nmols_out(50,:) = (/5, 1, 7/)
	nmols_out(51,:) = (/4, 2, 7/)
	nmols_out(52,:) = (/3, 3, 7/)
	nmols_out(53,:) = (/7, 0, 7/)
	nmols_out(54,:) = (/6, 1, 7/)
	nmols_out(55,:) = (/5, 2, 7/)
	nmols_out(56,:) = (/4, 3, 7/)
	nmols_out(57,:) = (/8, 0, 7/)
	nmols_out(58,:) = (/7, 1, 7/)
	nmols_out(59,:) = (/6, 2, 7/)
	nmols_out(60,:) = (/5, 3, 7/)
	nmols_out(61,:) = (/5, 0, 8/)
	nmols_out(62,:) = (/4, 1, 8/)
	nmols_out(63,:) = (/3, 2, 8/)
	nmols_out(64,:) = (/6, 0, 8/)
	nmols_out(65,:) = (/5, 1, 8/)
	nmols_out(66,:) = (/4, 2, 8/)
	nmols_out(67,:) = (/3, 3, 8/)
	nmols_out(68,:) = (/7, 0, 8/)
	nmols_out(69,:) = (/6, 3, 4/)
	nmols_out(70,:) = (/3, 4, 4/)
	nmols_out(71,:) = (/6, 1, 8/)
	nmols_out(72,:) = (/5, 2, 8/)
	nmols_out(73,:) = (/4, 3, 8/)
	nmols_out(74,:) = (/8, 0, 8/)
	nmols_out(75,:) = (/7, 2, 4/)
	nmols_out(76,:) = (/7, 3, 4/)
	nmols_out(77,:) = (/4, 4, 4/)
	nmols_out(78,:) = (/7, 1, 8/)
	nmols_out(79,:) = (/6, 2, 8/)
	nmols_out(80,:) = (/5, 3, 8/)
	nmols_out(81,:) = (/3, 5, 4/)
	nmols_out(82,:) = (/5, 4, 4/)
	nmols_out(83,:) = (/4, 5, 4/)
	nmols_out(84,:) = (/3, 6, 4/)
	nmols_out(85,:) = (/6, 4, 4/)
	nmols_out(86,:) = (/5, 5, 4/)
	nmols_out(87,:) = (/4, 6, 4/)
	nmols_out(88,:) = (/3, 4, 5/)
	nmols_out(89,:) = (/4, 4, 5/)
	nmols_out(90,:) = (/3, 5, 5/)
	nmols_out(91,:) = (/3, 4, 6/)
	nmols_out(92,:) = (/4, 4, 6/)
	nmols_out(93,:) = (/3, 5, 6/)
	nmols_out(94,:) = (/3, 4, 7/)
	nmols_out(95,:) = (/4, 4, 7/)
	nmols_out(96,:) = (/3, 5, 7/)
	nmols_out(97,:) = (/4, 4, 8/)
	nmols_out(98,:) = (/3, 5, 8/)

end subroutine get_nmols_out_4

subroutine molecule_names_nocharge_4(labels_nocharge)
	implicit none
	character(len=11), dimension(3) :: labels_nocharge

	labels_nocharge(1)(:) = 'A'
	labels_nocharge(2)(:) = 'Ms'
	labels_nocharge(3)(:) = 'D'

end subroutine molecule_names_nocharge_4

subroutine get_nmols_nocharge_4(nmols_nocharge_clust)
	implicit none
	integer :: nmols_nocharge_clust(58,3)

	nmols_nocharge_clust = 0

	nmols_nocharge_clust(1,:) = (/1, 0, 0/)
	nmols_nocharge_clust(2,:) = (/2, 0, 0/)
	nmols_nocharge_clust(3,:) = (/3, 0, 0/)
	nmols_nocharge_clust(4,:) = (/4, 0, 0/)
	nmols_nocharge_clust(5,:) = (/0, 0, 1/)
	nmols_nocharge_clust(6,:) = (/1, 0, 1/)
	nmols_nocharge_clust(7,:) = (/2, 0, 1/)
	nmols_nocharge_clust(8,:) = (/3, 0, 1/)
	nmols_nocharge_clust(9,:) = (/4, 0, 1/)
	nmols_nocharge_clust(10,:) = (/0, 0, 2/)
	nmols_nocharge_clust(11,:) = (/1, 0, 2/)
	nmols_nocharge_clust(12,:) = (/2, 0, 2/)
	nmols_nocharge_clust(13,:) = (/3, 0, 2/)
	nmols_nocharge_clust(14,:) = (/4, 0, 2/)
	nmols_nocharge_clust(15,:) = (/0, 0, 3/)
	nmols_nocharge_clust(16,:) = (/1, 0, 3/)
	nmols_nocharge_clust(17,:) = (/2, 0, 3/)
	nmols_nocharge_clust(18,:) = (/3, 0, 3/)
	nmols_nocharge_clust(19,:) = (/4, 0, 3/)
	nmols_nocharge_clust(20,:) = (/0, 0, 4/)
	nmols_nocharge_clust(21,:) = (/1, 0, 4/)
	nmols_nocharge_clust(22,:) = (/2, 0, 4/)
	nmols_nocharge_clust(23,:) = (/3, 0, 4/)
	nmols_nocharge_clust(24,:) = (/4, 0, 4/)
	nmols_nocharge_clust(25,:) = (/0, 1, 0/)
	nmols_nocharge_clust(26,:) = (/1, 1, 0/)
	nmols_nocharge_clust(27,:) = (/2, 1, 0/)
	nmols_nocharge_clust(28,:) = (/3, 1, 0/)
	nmols_nocharge_clust(29,:) = (/0, 2, 0/)
	nmols_nocharge_clust(30,:) = (/1, 2, 0/)
	nmols_nocharge_clust(31,:) = (/2, 2, 0/)
	nmols_nocharge_clust(32,:) = (/3, 2, 0/)
	nmols_nocharge_clust(33,:) = (/0, 3, 0/)
	nmols_nocharge_clust(34,:) = (/1, 3, 0/)
	nmols_nocharge_clust(35,:) = (/2, 3, 0/)
	nmols_nocharge_clust(36,:) = (/3, 3, 0/)
	nmols_nocharge_clust(37,:) = (/0, 4, 0/)
	nmols_nocharge_clust(38,:) = (/1, 1, 1/)
	nmols_nocharge_clust(39,:) = (/2, 1, 1/)
	nmols_nocharge_clust(40,:) = (/3, 1, 1/)
	nmols_nocharge_clust(41,:) = (/1, 2, 1/)
	nmols_nocharge_clust(42,:) = (/2, 2, 1/)
	nmols_nocharge_clust(43,:) = (/1, 3, 1/)
	nmols_nocharge_clust(44,:) = (/1, 1, 2/)
	nmols_nocharge_clust(45,:) = (/2, 1, 2/)
	nmols_nocharge_clust(46,:) = (/3, 1, 2/)
	nmols_nocharge_clust(47,:) = (/1, 2, 2/)
	nmols_nocharge_clust(48,:) = (/2, 2, 2/)
	nmols_nocharge_clust(49,:) = (/1, 3, 2/)
	nmols_nocharge_clust(50,:) = (/1, 1, 3/)
	nmols_nocharge_clust(51,:) = (/2, 1, 3/)
	nmols_nocharge_clust(52,:) = (/3, 1, 3/)
	nmols_nocharge_clust(53,:) = (/1, 2, 3/)
	nmols_nocharge_clust(54,:) = (/2, 2, 3/)
	nmols_nocharge_clust(55,:) = (/1, 3, 3/)
	nmols_nocharge_clust(56,:) = (/3, 1, 4/)
	nmols_nocharge_clust(57,:) = (/2, 2, 4/)
	nmols_nocharge_clust(58,:) = (/1, 3, 4/)

end subroutine get_nmols_nocharge_4

subroutine get_nmols_nocharge_out_4(nmols_nocharge_out)
	implicit none
	integer :: nmols_nocharge_out(98,3)

	nmols_nocharge_out = 0

	nmols_nocharge_out(1,:) = (/5, 0, 4/)
	nmols_nocharge_out(2,:) = (/4, 1, 4/)
	nmols_nocharge_out(3,:) = (/3, 2, 4/)
	nmols_nocharge_out(4,:) = (/6, 0, 4/)
	nmols_nocharge_out(5,:) = (/5, 1, 4/)
	nmols_nocharge_out(6,:) = (/4, 2, 4/)
	nmols_nocharge_out(7,:) = (/3, 3, 4/)
	nmols_nocharge_out(8,:) = (/7, 0, 4/)
	nmols_nocharge_out(9,:) = (/6, 1, 4/)
	nmols_nocharge_out(10,:) = (/5, 2, 4/)
	nmols_nocharge_out(11,:) = (/4, 3, 4/)
	nmols_nocharge_out(12,:) = (/8, 0, 4/)
	nmols_nocharge_out(13,:) = (/7, 1, 4/)
	nmols_nocharge_out(14,:) = (/6, 2, 4/)
	nmols_nocharge_out(15,:) = (/5, 3, 4/)
	nmols_nocharge_out(16,:) = (/5, 0, 5/)
	nmols_nocharge_out(17,:) = (/4, 1, 5/)
	nmols_nocharge_out(18,:) = (/3, 2, 5/)
	nmols_nocharge_out(19,:) = (/6, 0, 5/)
	nmols_nocharge_out(20,:) = (/5, 1, 5/)
	nmols_nocharge_out(21,:) = (/4, 2, 5/)
	nmols_nocharge_out(22,:) = (/3, 3, 5/)
	nmols_nocharge_out(23,:) = (/7, 0, 5/)
	nmols_nocharge_out(24,:) = (/6, 1, 5/)
	nmols_nocharge_out(25,:) = (/5, 2, 5/)
	nmols_nocharge_out(26,:) = (/4, 3, 5/)
	nmols_nocharge_out(27,:) = (/8, 0, 5/)
	nmols_nocharge_out(28,:) = (/7, 1, 5/)
	nmols_nocharge_out(29,:) = (/6, 2, 5/)
	nmols_nocharge_out(30,:) = (/5, 3, 5/)
	nmols_nocharge_out(31,:) = (/5, 0, 6/)
	nmols_nocharge_out(32,:) = (/4, 1, 6/)
	nmols_nocharge_out(33,:) = (/3, 2, 6/)
	nmols_nocharge_out(34,:) = (/6, 0, 6/)
	nmols_nocharge_out(35,:) = (/5, 1, 6/)
	nmols_nocharge_out(36,:) = (/4, 2, 6/)
	nmols_nocharge_out(37,:) = (/3, 3, 6/)
	nmols_nocharge_out(38,:) = (/7, 0, 6/)
	nmols_nocharge_out(39,:) = (/6, 1, 6/)
	nmols_nocharge_out(40,:) = (/5, 2, 6/)
	nmols_nocharge_out(41,:) = (/4, 3, 6/)
	nmols_nocharge_out(42,:) = (/8, 0, 6/)
	nmols_nocharge_out(43,:) = (/7, 1, 6/)
	nmols_nocharge_out(44,:) = (/6, 2, 6/)
	nmols_nocharge_out(45,:) = (/5, 3, 6/)
	nmols_nocharge_out(46,:) = (/5, 0, 7/)
	nmols_nocharge_out(47,:) = (/4, 1, 7/)
	nmols_nocharge_out(48,:) = (/3, 2, 7/)
	nmols_nocharge_out(49,:) = (/6, 0, 7/)
	nmols_nocharge_out(50,:) = (/5, 1, 7/)
	nmols_nocharge_out(51,:) = (/4, 2, 7/)
	nmols_nocharge_out(52,:) = (/3, 3, 7/)
	nmols_nocharge_out(53,:) = (/7, 0, 7/)
	nmols_nocharge_out(54,:) = (/6, 1, 7/)
	nmols_nocharge_out(55,:) = (/5, 2, 7/)
	nmols_nocharge_out(56,:) = (/4, 3, 7/)
	nmols_nocharge_out(57,:) = (/8, 0, 7/)
	nmols_nocharge_out(58,:) = (/7, 1, 7/)
	nmols_nocharge_out(59,:) = (/6, 2, 7/)
	nmols_nocharge_out(60,:) = (/5, 3, 7/)
	nmols_nocharge_out(61,:) = (/5, 0, 8/)
	nmols_nocharge_out(62,:) = (/4, 1, 8/)
	nmols_nocharge_out(63,:) = (/3, 2, 8/)
	nmols_nocharge_out(64,:) = (/6, 0, 8/)
	nmols_nocharge_out(65,:) = (/5, 1, 8/)
	nmols_nocharge_out(66,:) = (/4, 2, 8/)
	nmols_nocharge_out(67,:) = (/3, 3, 8/)
	nmols_nocharge_out(68,:) = (/7, 0, 8/)
	nmols_nocharge_out(69,:) = (/6, 3, 4/)
	nmols_nocharge_out(70,:) = (/3, 4, 4/)
	nmols_nocharge_out(71,:) = (/6, 1, 8/)
	nmols_nocharge_out(72,:) = (/5, 2, 8/)
	nmols_nocharge_out(73,:) = (/4, 3, 8/)
	nmols_nocharge_out(74,:) = (/8, 0, 8/)
	nmols_nocharge_out(75,:) = (/7, 2, 4/)
	nmols_nocharge_out(76,:) = (/7, 3, 4/)
	nmols_nocharge_out(77,:) = (/4, 4, 4/)
	nmols_nocharge_out(78,:) = (/7, 1, 8/)
	nmols_nocharge_out(79,:) = (/6, 2, 8/)
	nmols_nocharge_out(80,:) = (/5, 3, 8/)
	nmols_nocharge_out(81,:) = (/3, 5, 4/)
	nmols_nocharge_out(82,:) = (/5, 4, 4/)
	nmols_nocharge_out(83,:) = (/4, 5, 4/)
	nmols_nocharge_out(84,:) = (/3, 6, 4/)
	nmols_nocharge_out(85,:) = (/6, 4, 4/)
	nmols_nocharge_out(86,:) = (/5, 5, 4/)
	nmols_nocharge_out(87,:) = (/4, 6, 4/)
	nmols_nocharge_out(88,:) = (/3, 4, 5/)
	nmols_nocharge_out(89,:) = (/4, 4, 5/)
	nmols_nocharge_out(90,:) = (/3, 5, 5/)
	nmols_nocharge_out(91,:) = (/3, 4, 6/)
	nmols_nocharge_out(92,:) = (/4, 4, 6/)
	nmols_nocharge_out(93,:) = (/3, 5, 6/)
	nmols_nocharge_out(94,:) = (/3, 4, 7/)
	nmols_nocharge_out(95,:) = (/4, 4, 7/)
	nmols_nocharge_out(96,:) = (/3, 5, 7/)
	nmols_nocharge_out(97,:) = (/4, 4, 8/)
	nmols_nocharge_out(98,:) = (/3, 5, 8/)

end subroutine get_nmols_nocharge_out_4


end module acdc_system_4

