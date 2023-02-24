module acdc_system_AMsD

implicit none

integer, parameter :: nclust = 58						! number of clusters, molecules and ions
integer, parameter :: neq = 65							! number of equations

real(kind(1.d0)), parameter :: temp = 298.15d0				! temperature in K
real(kind(1.d0)), parameter :: rh = 0.d0					! RH in %

real(kind(1.d0)), parameter :: cs_exponent_default = -1.60000000000000d+00			! parameters for the exponential loss
real(kind(1.d0)), parameter :: cs_coefficient_default = 2.60000000000000d-03

integer, parameter :: n_monomer_types = 3
integer, parameter :: n1A = 1, n1D = 5, n1Ms = 25			! cluster indices for monomers and ions
integer, parameter :: nout_neu = 64

integer, parameter :: n_mol_types = 3
integer, parameter :: nmolA = 1, nmolMs = 2, nmolD = 3			! molecule indices for the used species

integer, parameter :: n_1A_clusters = 18				! number molecules and clusters containing 1 A molecule

integer, parameter :: n_neutrals = 58			! number of neutral molecules and clusters
integer, parameter :: n_negatives = 0			! number of negative molecules and clusters
integer, parameter :: n_positives = 0			! number of positive molecules and clusters

integer, parameter :: n_neutral_monomers = 3			! number of neutral monomers
integer, parameter :: n_negative_monomers = 0			! number of negative monomers
integer, parameter :: n_positive_monomers = 0			! number of positive monomers

integer, parameter :: n_neutral_clusters = 55			! number of neutral clusters
integer, parameter :: n_negative_clusters = 0			! number of negative clusters
integer, parameter :: n_positive_clusters = 0			! number of positive clusters

real(kind(1.d0)), parameter :: mass_max = 582.57
real(kind(1.d0)), parameter :: diameter_max = 1.18
real(kind(1.d0)), parameter :: mob_diameter_max = 1.48

integer, parameter :: nmols_out_neutral(3, 3) = reshape((/5, 4, 3, 0, 1, 2, 4, 4, 4/),(/3, 3/))			! criteria for outgrowing neutrals


contains

subroutine n_A_in_clusters(n_A)
	implicit none
	integer :: n_A(58)

	n_A = (/1, 2, 3, 4, 0, 1, 2, 3, 4, 0, &
		&1, 2, 3, 4, 0, 1, 2, 3, 4, 0, &
		&1, 2, 3, 4, 0, 1, 2, 3, 0, 1, &
		&2, 3, 0, 1, 2, 3, 0, 1, 2, 3, &
		&1, 2, 1, 1, 2, 3, 1, 2, 1, 1, &
		&2, 3, 1, 2, 1, 3, 2, 1/)

end subroutine n_A_in_clusters

subroutine clusters_with_1_A(cluster_numbers)
	implicit none
	integer :: cluster_numbers(18)

	cluster_numbers = (/1, 6, 11, 16, 21, 26, 30, 34, 38, 41, 43, 44, 47, 49, 50, &
		&53, 55, 58/)

end subroutine clusters_with_1_A

subroutine monomer_arrays(neutral_monomers)
	implicit none
	integer :: neutral_monomers(3)

	neutral_monomers = (/1, 5, 25/)

end subroutine monomer_arrays

subroutine cluster_arrays(neutral_clusters)
	implicit none
	integer :: neutral_clusters(55)

	neutral_clusters = (/2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, &
		&21, 22, 23, 24, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, &
		&41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58/)

end subroutine cluster_arrays

subroutine get_mass(mass)
	implicit none
	real(kind(1.d0)) :: mass(58)

	mass = (/98.08, 196.16, 294.24, 392.32, 45.08, 143.16, 241.24, 339.32, 437.40, 90.16, &
		&188.24, 286.32, 384.40, 482.48, 135.24, 233.32, 331.40, 429.48, 527.56, 180.32, &
		&278.40, 376.48, 474.56, 572.64, 96.11, 194.19, 292.27, 390.35, 192.22, 290.30, &
		&388.38, 486.46, 288.33, 386.41, 484.49, 582.57, 384.44, 239.27, 337.35, 435.43, &
		&335.38, 433.46, 431.49, 284.35, 382.43, 480.51, 380.46, 478.54, 476.57, 329.43, &
		&427.51, 525.59, 425.54, 523.62, 521.65, 570.67, 568.70, 566.73/)

end subroutine get_mass

subroutine get_diameter(diameter)
	implicit none

	real(kind(1.d0)) :: diameter(58)

	 diameter = (/0.55, 0.70, 0.80, 0.88, 0.59, 0.72, 0.82, 0.90, 0.96, 0.75, &
		&0.84, 0.91, 0.98, 1.03, 0.86, 0.93, 0.99, 1.04, 1.09, 0.94, &
		&1.00, 1.06, 1.11, 1.15, 0.59, 0.72, 0.82, 0.89, 0.74, 0.83, &
		&0.91, 0.97, 0.85, 0.92, 0.99, 1.04, 0.94, 0.84, 0.91, 0.97, &
		&0.93, 0.99, 1.00, 0.93, 0.99, 1.04, 1.00, 1.05, 1.06, 1.00, &
		&1.06, 1.10, 1.07, 1.11, 1.12, 1.16, 1.17, 1.18/)	! dry value

end subroutine get_diameter

subroutine get_mob_diameter(mob_diameter)
	implicit none

	real(kind(1.d0)) :: mob_diameter(58)

	 mob_diameter = (/0.85, 1.00, 1.10, 1.18, 0.89, 1.02, 1.12, 1.20, 1.26, 1.05, &
		&1.14, 1.21, 1.28, 1.33, 1.16, 1.23, 1.29, 1.34, 1.39, 1.24, &
		&1.30, 1.36, 1.41, 1.45, 0.89, 1.02, 1.12, 1.19, 1.04, 1.13, &
		&1.21, 1.27, 1.15, 1.22, 1.29, 1.34, 1.24, 1.14, 1.21, 1.27, &
		&1.23, 1.29, 1.30, 1.23, 1.29, 1.34, 1.30, 1.35, 1.36, 1.30, &
		&1.36, 1.40, 1.37, 1.41, 1.42, 1.46, 1.47, 1.48/)	! dry value

end subroutine get_mob_diameter

subroutine cluster_names(clust)
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
	clust(63)(:) = 'rec'
	clust(64)(:) = 'out_neu'
	clust(65)(:) = 'bound'

end subroutine cluster_names

subroutine molecule_names(labels)
	implicit none
	character(len=11), dimension(3) :: labels

	labels(1)(:) = 'A'
	labels(2)(:) = 'Ms'
	labels(3)(:) = 'D'

end subroutine molecule_names

subroutine monomer_indices(n_monomers)
	implicit none
	integer :: n_monomers(3)

	n_monomers = (/1, 5, 25/)

end subroutine monomer_indices


end module acdc_system_AMsD

