!!!!!!!!!!!!!!!!!!!!!!!! KEY !!!!!!!!!!!!!!!!!!!!!!!!
! Cluster 1: 1Ii
! Cluster 2: 2Ii
! Cluster 3: 3Ii
! Cluster 4: 4Ii
! Cluster 5: 1Io
! Cluster 6: 1Ii1Io
! Cluster 7: 2Ii1Io
! Cluster 8: 3Ii1Io
! Cluster 9: 4Ii1Io
! Cluster 10: 2Io
! Cluster 11: 1Ii2Io
! Cluster 12: 2Ii2Io
! Cluster 13: 3Ii2Io
! Cluster 14: 4Ii2Io
! Cluster 15: 3Io
! Cluster 16: 1Ii3Io
! Cluster 17: 2Ii3Io
! Cluster 18: 3Ii3Io
! Cluster 19: 4Ii3Io
! Cluster 20: 4Io
! Cluster 21: 1Ii4Io
! Cluster 22: 2Ii4Io
! Cluster 23: 3Ii4Io
! Cluster 24: 4Ii4Io
! 25 is for source fluxes
! 26 is for coagulation losses
! 27 is for wall losses
! 28 is for dilution losses
! 29 is for user-defined losses
! 30 is for recombination of positive and negative charger ions with each other
! 31 is for collisions that lead succesfully out of the system as neutrals
! 32 is for collisions that lead out of the system, but the resulting cluster is brought back
! 33-56 are for coagulation losses of each cluster
! 57-80 are for collisions that lead out of the system to different cluster compositions
!!!!!!!!!!!!!!!!!!!!!! END KEY !!!!!!!!!!!!!!!!!!!!!!

! differential equations: f = dc/dt
subroutine feval_3(neqn,t,c,f,coef,ipar)
	implicit none
	integer, parameter :: nclust = 24
	integer :: neqn, ipar(4), i, j, k, n
	real(kind(1.d0)) :: t, c(neqn), f(neqn), coef(1)
	logical, save :: isconst(80) = .false.
	real(kind(1.d0)), save :: coef_quad(nclust,nclust,32)=0.d0,coef_lin(32,32,nclust)=0.d0,source(80)=0.d0
	integer, save :: ind_quad_loss(32,0:50)=0,ind_quad_form(32,0:208)=0
	integer, save :: ind_quad_loss_extra(32,0:48)=0,ind_quad_form_extra(32,0:144)=0
	integer, save :: ind_lin_loss(32,0:26)=0,ind_lin_form(32,0:48)=0,fitted(80,0:80)=0
	integer, save :: ind_out_clust(nclust,nclust)=0
	real(kind(1.d0)), save :: n_quad_form_extra(32,0:48)=0.d0
	real(kind(1.d0)), parameter :: pi=4.d0*atan(1.d0)

	! these parameters are read at the very beginning and eg. if source terms or collision rates have changed
	if (ipar(1) .eq. 0) then
		! after this the same values are used until some other routine tells otherwise
		ipar(1) = 1
		call initialize_parameters_3(coef_quad,coef_lin,ind_quad_loss,ind_quad_form,&
		&	ind_lin_loss,ind_lin_form,ind_quad_loss_extra,ind_quad_form_extra,&
		&	ind_out_clust,&
		&	source,isconst,fitted,coef,ipar)
		n_quad_form_extra(:,1:48) = real(ind_quad_form_extra(:,3:144:3),kind=kind(1.d0))
	end if

	! override the isconst settings, if needed
	if (ipar(4) .eq. 1) then
		isconst = .false.
	end if

	f = 0.d0
	do i = 1, neqn-48 ! loop over all the clusters + other fluxes
		if (.not. isconst(i)) then
			! first calculate coefficients for all loss terms
			do n = 1, 2*ind_quad_loss(i,0)-1, 2 ! loop over all collisions removing this cluster
				f(i) = f(i)-coef_quad(i,ind_quad_loss(i,n),ind_quad_loss(i,n+1))*c(ind_quad_loss(i,n))
			end do
			do n = 1, 2*ind_lin_loss(i,0)-1, 2 ! loop over all evaporations + wall losses etc. removing this cluster
				f(i) = f(i)-coef_lin(ind_lin_loss(i,n),ind_lin_loss(i,n+1),i)
			end do
			f(i) = f(i)*c(i) ! multiplying loss coefficients with current concentration
			! then add all terms that form this cluster
			do n = 1, 2*ind_quad_form(i,0)-1, 2 ! loop over all collisions forming this cluster
				f(i) = f(i)+coef_quad(ind_quad_form(i,n),ind_quad_form(i,n+1),i)*&
				&      c(ind_quad_form(i,n))*c(ind_quad_form(i,n+1))
			end do
			do n = 1, 3*ind_quad_form_extra(i,0)-2, 3 ! loop over all collisions forming this cluster as an extra product
				f(i) = f(i)+coef_quad(ind_quad_form_extra(i,n),ind_quad_form_extra(i,n+1),i)*&
				&      n_quad_form_extra(i,(n+2)/3)*c(ind_quad_form_extra(i,n))*c(ind_quad_form_extra(i,n+1))
			end do
			do n = 1, 2*ind_lin_form(i,0)-1, 2 ! loop over all evaporations forming this cluster
				f(i) = f(i)+coef_lin(i,ind_lin_form(i,n),ind_lin_form(i,n+1))*c(ind_lin_form(i,n+1))
			end do
			! finally, add possible external sources
			f(i) = f(i) + source(i)
		end if
	end do
	! coagulation for each cluster
	i = 26
	do n = 1, 2*ind_lin_form(i,0)-1, 2 ! loop over all coagulations
		j = 32+ind_lin_form(i,n+1)
		f(j) = f(j)+coef_lin(i,ind_lin_form(i,n),ind_lin_form(i,n+1))*c(ind_lin_form(i,n+1))
	end do
	! flux out for each outgrown cluster
	do i = 31, 31 ! loop over all charging states
		do n = 1, 2*ind_quad_form(i,0)-1, 2 ! loop over all collisions going out through this channel
			j = 56+ind_out_clust(ind_quad_form(i,n),ind_quad_form(i,n+1))
			f(j) = f(j)+coef_quad(ind_quad_form(i,n),ind_quad_form(i,n+1),i)*&
			&      c(ind_quad_form(i,n))*c(ind_quad_form(i,n+1))
		end do
	end do
	! loop over the clusters whose concentrations are fitted (e.g. [1A]+[1A1D]=const.)
	do n = 1, fitted(1,0)
		f(fitted(n+1,0)) = 0.d0 ! don't use value from birth-death equations
		do i = 1, fitted(n+1,1) ! loop over the clusters used for the fitting
			f(fitted(n+1,0)) = f(fitted(n+1,0))-f(fitted(n+1,i+1)) ! (e.g. d[1A]/dt = -d[1A1D]/dt)
		end do
	end do

end subroutine feval_3

!-----------------------------------------------------------

! jacobian of the differential equations: dfdc(i,j) = df(i)/dc(j)
! not using this, since the solver works faster using a numerical jacobian...
subroutine jeval_3(neqn,t,c,ml,mu,dfdc,ldim,coef,ipar)
	implicit none
	integer :: ldim,neqn,ierr,ipar(4),i,j,k,n
	real(kind(1.d0)) :: t,c(neqn),dfdc(ldim,neqn),coef(1),ml,mu

end subroutine jeval_3

!-----------------------------------------------------------

subroutine formation_3(neqn,c,j_tot,j_by_charge,j_by_cluster,j_all,ipar,coef)
	implicit none
	integer, parameter :: nclust = 24
	integer :: neqn, i, n, ipar(4)
	real(kind(1.d0)) :: c(neqn), j_add, j_tot,j_by_charge(4),j_by_cluster(neqn),j_all(neqn,4),coef(1)
	real(kind(1.d0)), save :: coef_quad(nclust,nclust,32)=0.d0,coef_lin(32,32,nclust)=0.d0,source(80)=0.d0
	integer, save :: charge(nclust)=0,ind_quad_loss(32,0:50)=0,ind_quad_form(32,0:208)=0
	integer, save :: ind_quad_loss_extra(32,0:48)=0,ind_quad_form_extra(32,0:144)=0
	integer, save :: ind_lin_loss(32,0:26)=0,ind_lin_form(32,0:48)=0,fitted(80,0:80)=0
	integer, save :: ind_out_clust(nclust,nclust)=0
	logical, save :: isconst(80)=.false.

	if (ipar(3) .eq. 0) then
		ipar(3) = 1
		call initialize_parameters_3(coef_quad,coef_lin,ind_quad_loss,ind_quad_form,&
			&ind_lin_loss,ind_lin_form,ind_quad_loss_extra,ind_quad_form_extra,&
		&	ind_out_clust,&
		&	source,isconst,fitted,coef,ipar)
		! cluster charges
		charge = (/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
	end if
	j_tot = 0.d0			! total formation rate
	j_by_charge = 0.d0		! contributions of neutrals, negatives, positives and recombinations
	j_by_cluster = 0.d0		! contributions of different clusters
	j_all = 0.d0			! contributions of different clusters and different charging states
	do n = 1, 2*ind_quad_form(31,0)-1, 2 ! loop over all collisions leading out of the system
		j_add = coef_quad(ind_quad_form(31,n),ind_quad_form(31,n+1),31)*c(ind_quad_form(31,n))*c(ind_quad_form(31,n+1))
		j_tot = j_tot + j_add
		j_by_charge(1) = j_by_charge(1) + j_add
		j_by_cluster(ind_quad_form(31,n)) = j_by_cluster(ind_quad_form(31,n)) + j_add
		j_by_cluster(ind_quad_form(31,n+1)) = j_by_cluster(ind_quad_form(31,n+1)) + j_add
		j_all(ind_quad_form(31,n),1) = j_all(ind_quad_form(31,n),1) + j_add
		j_all(ind_quad_form(31,n+1),1) = j_all(ind_quad_form(31,n+1),1) + j_add
	end do

end subroutine formation_3

!-----------------------------------------------------------

subroutine initialize_parameters_3(coef_quad,coef_lin,ind_quad_loss,ind_quad_form,&
			&ind_lin_loss,ind_lin_form,ind_quad_loss_extra,ind_quad_form_extra,&
		&	ind_out_clust,&
		&	source,isconst,fitted,coef,ipar)

	use acdc_simulation_setup_3, only : sources_and_constants_3

	implicit none
	integer, parameter :: nclust=24, neqn=80
	logical :: isconst(80)
	real(kind(1.d0)) :: coef_quad(nclust,nclust,32),coef_lin(32,32,nclust)
	real(kind(1.d0)) :: source(80)
	integer :: ind_quad_loss(32,0:50),ind_quad_form(32,0:208),ind_quad_loss_extra(32,0:48),ind_quad_form_extra(32,0:144)
	integer :: ind_lin_loss(32,0:26),ind_lin_form(32,0:48)
	integer :: fitted(80,0:80)
	integer :: ind_out_clust(nclust,nclust)
	real(kind(1.d0)) :: coef(1)
	integer :: ipar(4)

	call sources_and_constants_3(80,source,isconst,fitted)

	ind_quad_loss = 0
	ind_quad_form = 0
	ind_lin_loss = 0
	ind_lin_form = 0
	ind_quad_loss_extra = 0
	ind_quad_form_extra = 0
	ind_out_clust = 0

	! ind_quad_loss(i,0:2n) = (/n,j1,k1,j2,k2,...,jn,kn/) gives the cluster numbers for all the collisions
	!  i + j1 -> k1 etc. through which cluster i is lost

	! ind_quad_form(k,0:2n) = (/n,i1,j1,i2,j2,...,in,jn/) gives the cluster numbers for all the collisions
	!  i1 + j1 -> k etc. through which cluster k is formed

	! ind_lin_loss(k,0:2n) = (/n,i1,j1,i2,j2,...,lossn,lossn/) gives the cluster numbers for all the evaporations
	!  k -> i1 + j1 etc. and other losses k -> wall etc. through which cluster k is lost

	! ind_lin_form(i,0:2n) = (/n,j1,k1,j2,k2,...,jn,kn/) gives the cluster numbers for all the evaporations
	!  k1 -> i + j1 etc. through which cluster i is formed

	! ind_quad_loss_extra(i,0:3n) = (/n,j1,k1,c1,...,jn,kn,cn/) gives the cluster numbers and coefficients 
	!  i + j1 -> c1*k1 etc. for additional collision products k (e.g. monomers from the boundary)
	!  in collisions where cluster i is lost

	! ind_quad_form_extra(k,0:2n) = (/n,i1,j1,c1,...,in,jn,cn/) gives the cluster numbers and coefficients 
	!  i1 + j1 -> c1*k etc. for additional ways of forming k (e.g. as a monomer from the boundary)


	! Cluster 1: 1Ii
	ind_quad_loss(1,0:42) = (/ 21, 1,2, 1,2, 2,3, 3,4, 5,6, 6,7, 7,8, 8,9, 10,11&
			&, 11,12, 12,13, 13,14, 15,16, 16,17, 17,18, 18,19, 20,21, 21,22, 22,23&
			&, 23,24, 24,31 /)
	ind_quad_form_extra(1,0:144) = (/ 48, 2,3,1, 2,4,2, 2,8,1, 2,9,2, 2,13,1, 2,14,2, 2,18,1, 2,19,2, 3,3,2&
			&, 3,4,3, 3,7,1, 3,8,2, 3,9,3, 3,12,1, 3,13,2, 3,14,3, 3,17,1, 3,18,2, 3,19,3&
			&, 4,4,4, 4,6,1, 4,7,2, 4,8,3, 4,9,4, 4,11,1, 4,12,2, 4,13,3, 4,14,4, 4,16,1&
			&, 4,17,2, 4,18,3, 4,19,4, 6,9,1, 6,14,1, 7,8,1, 7,9,2, 7,13,1, 7,14,2, 8,8,2&
			&, 8,9,3, 8,12,1, 8,13,2, 8,14,3, 9,9,4, 9,11,1, 9,12,2, 9,13,3, 9,14,4 /)
	ind_lin_loss(1,0:2) = (/ 1, 26,26 /)
	ind_lin_form(1,0:40) = (/ 20, 1,2, 1,2, 2,3, 3,4, 5,6, 6,7, 7,8, 8,9, 10,11&
			&, 11,12, 12,13, 13,14, 15,16, 16,17, 17,18, 18,19, 20,21, 21,22, 22,23&
			&, 23,24 /)

	! Cluster 2: 2Ii
	ind_quad_loss(2,0:50) = (/ 25, 1,3, 2,4, 2,4, 3,4, 4,4, 5,7, 6,8, 7,9, 8,9&
			&, 9,9, 10,12, 11,13, 12,14, 13,14, 14,14, 15,17, 16,18, 17,19, 18,19&
			&, 19,19, 20,22, 21,23, 22,24, 23,31, 24,31 /)
	ind_quad_loss_extra(2,0:24) = (/ 8, 3,1,1, 4,1,2, 8,1,1, 9,1,2, 13,1,1, 14,1,2, 18,1,1, 19,1,2 /)
	ind_quad_form(2,0:2) = (/ 1, 1,1 /)
	ind_lin_loss(2,0:4) = (/ 2, 1,1, 26,26 /)
	ind_lin_form(2,0:30) = (/ 15, 1,3, 2,4, 2,4, 5,7, 6,8, 7,9, 10,12, 11,13, 12,14&
			&, 15,17, 16,18, 17,19, 20,22, 21,23, 22,24 /)

	! Cluster 3: 3Ii
	ind_quad_loss(3,0:50) = (/ 25, 1,4, 2,4, 3,4, 3,4, 4,4, 5,8, 6,9, 7,9, 8,9&
			&, 9,9, 10,13, 11,14, 12,14, 13,14, 14,14, 15,18, 16,19, 17,19, 18,19&
			&, 19,19, 20,23, 21,24, 22,31, 23,31, 24,31 /)
	ind_quad_loss_extra(3,0:39) = (/ 13, 2,1,1, 3,1,2, 3,1,2, 4,1,3, 7,1,1, 8,1,2, 9,1,3, 12,1,1, 13,1,2&
			&, 14,1,3, 17,1,1, 18,1,2, 19,1,3 /)
	ind_quad_form(3,0:2) = (/ 1, 1,2 /)
	ind_lin_loss(3,0:4) = (/ 2, 1,2, 26,26 /)
	ind_lin_form(3,0:18) = (/ 9, 1,4, 5,8, 6,9, 10,13, 11,14, 15,18, 16,19, 20,23, 21,24 /)

	! Cluster 4: 4Ii
	ind_quad_loss(4,0:48) = (/ 24, 2,4, 3,4, 4,4, 4,4, 5,9, 6,9, 7,9, 8,9, 9,9&
			&, 10,14, 11,14, 12,14, 13,14, 14,14, 15,19, 16,19, 17,19, 18,19, 19,19&
			&, 20,24, 21,31, 22,31, 23,31, 24,31 /)
	ind_quad_loss_extra(4,0:48) = (/ 16, 2,1,2, 3,1,3, 4,1,4, 4,1,4, 6,1,1, 7,1,2, 8,1,3, 9,1,4, 11,1,1&
			&, 12,1,2, 13,1,3, 14,1,4, 16,1,1, 17,1,2, 18,1,3, 19,1,4 /)
	ind_quad_form(4,0:14) = (/ 7, 1,3, 2,2, 2,3, 2,4, 3,3, 3,4, 4,4 /)
	ind_lin_loss(4,0:6) = (/ 3, 1,3, 2,2, 26,26 /)
	ind_lin_form(4,0:8) = (/ 4, 5,9, 10,14, 15,19, 20,24 /)

	! Cluster 5: 1Io
	ind_quad_loss(5,0:42) = (/ 21, 1,6, 2,7, 3,8, 4,9, 5,10, 5,10, 6,11, 7,12, 8,13&
			&, 9,14, 10,15, 11,16, 12,17, 13,18, 14,19, 15,20, 16,21, 17,22, 18,23&
			&, 19,24, 24,31 /)
	ind_quad_form_extra(5,0:144) = (/ 48, 6,20,1, 6,21,1, 6,22,1, 7,20,1, 7,21,1, 8,20,1, 10,15,1, 10,16,1, 10,17,1&
			&, 10,18,1, 10,20,2, 10,21,2, 10,22,2, 10,23,2, 11,15,1, 11,16,1, 11,17,1, 11,20,2, 11,21,2&
			&, 11,22,2, 12,15,1, 12,16,1, 12,20,2, 12,21,2, 13,15,1, 13,20,2, 15,15,2, 15,16,2, 15,17,2&
			&, 15,18,2, 15,20,3, 15,21,3, 15,22,3, 15,23,3, 16,16,2, 16,17,2, 16,20,3, 16,21,3, 16,22,3&
			&, 17,20,3, 17,21,3, 18,20,3, 20,20,4, 20,21,4, 20,22,4, 20,23,4, 21,21,4, 21,22,4 /)
	ind_lin_loss(5,0:2) = (/ 1, 26,26 /)
	ind_lin_form(5,0:40) = (/ 20, 1,6, 2,7, 3,8, 4,9, 5,10, 5,10, 6,11, 7,12, 8,13&
			&, 9,14, 10,15, 11,16, 12,17, 13,18, 14,19, 15,20, 16,21, 17,22, 18,23&
			&, 19,24 /)

	! Cluster 6: 1Ii1Io
	ind_quad_loss(6,0:50) = (/ 25, 1,7, 2,8, 3,9, 4,9, 5,11, 6,12, 6,12, 7,13, 8,14&
			&, 9,14, 10,16, 11,17, 12,18, 13,19, 14,19, 15,21, 16,22, 17,23, 18,24&
			&, 19,31, 20,21, 21,22, 22,23, 23,31, 24,31 /)
	ind_quad_loss_extra(6,0:18) = (/ 6, 4,1,1, 9,1,1, 14,1,1, 20,5,1, 21,5,1, 22,5,1 /)
	ind_quad_form(6,0:2) = (/ 1, 1,5 /)
	ind_lin_loss(6,0:4) = (/ 2, 1,5, 26,26 /)
	ind_lin_form(6,0:32) = (/ 16, 1,7, 2,8, 3,9, 5,11, 6,12, 6,12, 7,13, 8,14, 10,16&
			&, 11,17, 12,18, 13,19, 15,21, 16,22, 17,23, 18,24 /)

	! Cluster 7: 2Ii1Io
	ind_quad_loss(7,0:50) = (/ 25, 1,8, 2,9, 3,9, 4,9, 5,12, 6,13, 7,14, 7,14, 8,14&
			&, 9,14, 10,17, 11,18, 12,19, 13,19, 14,19, 15,22, 16,23, 17,24, 18,31&
			&, 19,31, 20,22, 21,23, 22,31, 23,31, 24,31 /)
	ind_quad_loss_extra(7,0:24) = (/ 8, 3,1,1, 4,1,2, 8,1,1, 9,1,2, 13,1,1, 14,1,2, 20,5,1, 21,5,1 /)
	ind_quad_form(7,0:4) = (/ 2, 1,6, 2,5 /)
	ind_lin_loss(7,0:6) = (/ 3, 1,6, 2,5, 26,26 /)
	ind_lin_form(7,0:24) = (/ 12, 1,8, 2,9, 5,12, 6,13, 7,14, 7,14, 10,17, 11,18, 12,19&
			&, 15,22, 16,23, 17,24 /)

	! Cluster 8: 3Ii1Io
	ind_quad_loss(8,0:50) = (/ 25, 1,9, 2,9, 3,9, 4,9, 5,13, 6,14, 7,14, 8,14, 8,14&
			&, 9,14, 10,18, 11,19, 12,19, 13,19, 14,19, 15,23, 16,24, 17,31, 18,31&
			&, 19,31, 20,23, 21,31, 22,31, 23,31, 24,31 /)
	ind_quad_loss_extra(8,0:33) = (/ 11, 2,1,1, 3,1,2, 4,1,3, 7,1,1, 8,1,2, 8,1,2, 9,1,3, 12,1,1, 13,1,2&
			&, 14,1,3, 20,5,1 /)
	ind_quad_form(8,0:6) = (/ 3, 1,7, 2,6, 3,5 /)
	ind_lin_loss(8,0:8) = (/ 4, 1,7, 2,6, 3,5, 26,26 /)
	ind_lin_form(8,0:14) = (/ 7, 1,9, 5,13, 6,14, 10,18, 11,19, 15,23, 16,24 /)

	! Cluster 9: 4Ii1Io
	ind_quad_loss(9,0:48) = (/ 24, 2,9, 3,9, 4,9, 5,14, 6,14, 7,14, 8,14, 9,14, 9,14&
			&, 10,19, 11,19, 12,19, 13,19, 14,19, 15,24, 16,31, 17,31, 18,31, 19,31&
			&, 20,31, 21,31, 22,31, 23,31, 24,31 /)
	ind_quad_loss_extra(9,0:36) = (/ 12, 2,1,2, 3,1,3, 4,1,4, 6,1,1, 7,1,2, 8,1,3, 9,1,4, 9,1,4, 11,1,1&
			&, 12,1,2, 13,1,3, 14,1,4 /)
	ind_quad_form(9,0:26) = (/ 13, 1,8, 2,7, 2,8, 2,9, 3,6, 3,7, 3,8, 3,9, 4,5&
			&, 4,6, 4,7, 4,8, 4,9 /)
	ind_lin_loss(9,0:10) = (/ 5, 1,8, 2,7, 3,6, 4,5, 26,26 /)
	ind_lin_form(9,0:6) = (/ 3, 5,14, 10,19, 15,24 /)

	! Cluster 10: 2Io
	ind_quad_loss(10,0:50) = (/ 25, 1,11, 2,12, 3,13, 4,14, 5,15, 6,16, 7,17, 8,18, 9,19&
			&, 10,20, 10,20, 11,21, 12,22, 13,23, 14,24, 15,20, 16,21, 17,22, 18,23&
			&, 19,31, 20,20, 21,21, 22,22, 23,23, 24,31 /)
	ind_quad_loss_extra(10,0:24) = (/ 8, 15,5,1, 16,5,1, 17,5,1, 18,5,1, 20,5,2, 21,5,2, 22,5,2, 23,5,2 /)
	ind_quad_form(10,0:2) = (/ 1, 5,5 /)
	ind_lin_loss(10,0:4) = (/ 2, 5,5, 26,26 /)
	ind_lin_form(10,0:30) = (/ 15, 1,11, 2,12, 3,13, 4,14, 5,15, 6,16, 7,17, 8,18, 9,19&
			&, 10,20, 10,20, 11,21, 12,22, 13,23, 14,24 /)

	! Cluster 11: 1Ii2Io
	ind_quad_loss(11,0:50) = (/ 25, 1,12, 2,13, 3,14, 4,14, 5,16, 6,17, 7,18, 8,19, 9,19&
			&, 10,21, 11,22, 11,22, 12,23, 13,24, 14,31, 15,21, 16,22, 17,23, 18,31&
			&, 19,31, 20,21, 21,22, 22,23, 23,31, 24,31 /)
	ind_quad_loss_extra(11,0:24) = (/ 8, 4,1,1, 9,1,1, 15,5,1, 16,5,1, 17,5,1, 20,5,2, 21,5,2, 22,5,2 /)
	ind_quad_form(11,0:4) = (/ 2, 1,10, 5,6 /)
	ind_lin_loss(11,0:6) = (/ 3, 1,10, 5,6, 26,26 /)
	ind_lin_form(11,0:24) = (/ 12, 1,12, 2,13, 3,14, 5,16, 6,17, 7,18, 8,19, 10,21, 11,22&
			&, 11,22, 12,23, 13,24 /)

	! Cluster 12: 2Ii2Io
	ind_quad_loss(12,0:50) = (/ 25, 1,13, 2,14, 3,14, 4,14, 5,17, 6,18, 7,19, 8,19, 9,19&
			&, 10,22, 11,23, 12,24, 12,24, 13,31, 14,31, 15,22, 16,23, 17,31, 18,31&
			&, 19,31, 20,22, 21,23, 22,31, 23,31, 24,31 /)
	ind_quad_loss_extra(12,0:24) = (/ 8, 3,1,1, 4,1,2, 8,1,1, 9,1,2, 15,5,1, 16,5,1, 20,5,2, 21,5,2 /)
	ind_quad_form(12,0:8) = (/ 4, 1,11, 2,10, 5,7, 6,6 /)
	ind_lin_loss(12,0:10) = (/ 5, 1,11, 2,10, 5,7, 6,6, 26,26 /)
	ind_lin_form(12,0:18) = (/ 9, 1,13, 2,14, 5,17, 6,18, 7,19, 10,22, 11,23, 12,24, 12,24 /)

	! Cluster 13: 3Ii2Io
	ind_quad_loss(13,0:50) = (/ 25, 1,14, 2,14, 3,14, 4,14, 5,18, 6,19, 7,19, 8,19, 9,19&
			&, 10,23, 11,24, 12,31, 13,31, 13,31, 14,31, 15,23, 16,31, 17,31, 18,31&
			&, 19,31, 20,23, 21,31, 22,31, 23,31, 24,31 /)
	ind_quad_loss_extra(13,0:24) = (/ 8, 2,1,1, 3,1,2, 4,1,3, 7,1,1, 8,1,2, 9,1,3, 15,5,1, 20,5,2 /)
	ind_quad_form(13,0:10) = (/ 5, 1,12, 2,11, 3,10, 5,8, 6,7 /)
	ind_lin_loss(13,0:12) = (/ 6, 1,12, 2,11, 3,10, 5,8, 6,7, 26,26 /)
	ind_lin_form(13,0:10) = (/ 5, 1,14, 5,18, 6,19, 10,23, 11,24 /)

	! Cluster 14: 4Ii2Io
	ind_quad_loss(14,0:48) = (/ 24, 2,14, 3,14, 4,14, 5,19, 6,19, 7,19, 8,19, 9,19, 10,24&
			&, 11,31, 12,31, 13,31, 14,31, 14,31, 15,31, 16,31, 17,31, 18,31, 19,31&
			&, 20,31, 21,31, 22,31, 23,31, 24,31 /)
	ind_quad_loss_extra(14,0:21) = (/ 7, 2,1,2, 3,1,3, 4,1,4, 6,1,1, 7,1,2, 8,1,3, 9,1,4 /)
	ind_quad_form(14,0:44) = (/ 22, 1,13, 2,12, 2,13, 2,14, 3,11, 3,12, 3,13, 3,14, 4,10&
			&, 4,11, 4,12, 4,13, 4,14, 5,9, 6,8, 6,9, 7,7, 7,8, 7,9&
			&, 8,8, 8,9, 9,9 /)
	ind_lin_loss(14,0:16) = (/ 8, 1,13, 2,12, 3,11, 4,10, 5,9, 6,8, 7,7, 26,26 /)
	ind_lin_form(14,0:4) = (/ 2, 5,19, 10,24 /)

	! Cluster 15: 3Io
	ind_quad_loss(15,0:50) = (/ 25, 1,16, 2,17, 3,18, 4,19, 5,20, 6,21, 7,22, 8,23, 9,24&
			&, 10,20, 11,21, 12,22, 13,23, 14,31, 15,20, 15,20, 16,21, 17,22, 18,23&
			&, 19,31, 20,20, 21,21, 22,22, 23,23, 24,31 /)
	ind_quad_loss_extra(15,0:39) = (/ 13, 10,5,1, 11,5,1, 12,5,1, 13,5,1, 15,5,2, 15,5,2, 16,5,2, 17,5,2, 18,5,2&
			&, 20,5,3, 21,5,3, 22,5,3, 23,5,3 /)
	ind_quad_form(15,0:2) = (/ 1, 5,10 /)
	ind_lin_loss(15,0:4) = (/ 2, 5,10, 26,26 /)
	ind_lin_form(15,0:18) = (/ 9, 1,16, 2,17, 3,18, 4,19, 5,20, 6,21, 7,22, 8,23, 9,24 /)

	! Cluster 16: 1Ii3Io
	ind_quad_loss(16,0:50) = (/ 25, 1,17, 2,18, 3,19, 4,19, 5,21, 6,22, 7,23, 8,24, 9,31&
			&, 10,21, 11,22, 12,23, 13,31, 14,31, 15,21, 16,22, 16,22, 17,23, 18,31&
			&, 19,31, 20,21, 21,22, 22,23, 23,31, 24,31 /)
	ind_quad_loss_extra(16,0:33) = (/ 11, 4,1,1, 10,5,1, 11,5,1, 12,5,1, 15,5,2, 16,5,2, 16,5,2, 17,5,2, 20,5,3&
			&, 21,5,3, 22,5,3 /)
	ind_quad_form(16,0:6) = (/ 3, 1,15, 5,11, 6,10 /)
	ind_lin_loss(16,0:8) = (/ 4, 1,15, 5,11, 6,10, 26,26 /)
	ind_lin_form(16,0:14) = (/ 7, 1,17, 2,18, 3,19, 5,21, 6,22, 7,23, 8,24 /)

	! Cluster 17: 2Ii3Io
	ind_quad_loss(17,0:50) = (/ 25, 1,18, 2,19, 3,19, 4,19, 5,22, 6,23, 7,24, 8,31, 9,31&
			&, 10,22, 11,23, 12,31, 13,31, 14,31, 15,22, 16,23, 17,31, 17,31, 18,31&
			&, 19,31, 20,22, 21,23, 22,31, 23,31, 24,31 /)
	ind_quad_loss_extra(17,0:24) = (/ 8, 3,1,1, 4,1,2, 10,5,1, 11,5,1, 15,5,2, 16,5,2, 20,5,3, 21,5,3 /)
	ind_quad_form(17,0:10) = (/ 5, 1,16, 2,15, 5,12, 6,11, 7,10 /)
	ind_lin_loss(17,0:12) = (/ 6, 1,16, 2,15, 5,12, 6,11, 7,10, 26,26 /)
	ind_lin_form(17,0:10) = (/ 5, 1,18, 2,19, 5,22, 6,23, 7,24 /)

	! Cluster 18: 3Ii3Io
	ind_quad_loss(18,0:50) = (/ 25, 1,19, 2,19, 3,19, 4,19, 5,23, 6,24, 7,31, 8,31, 9,31&
			&, 10,23, 11,31, 12,31, 13,31, 14,31, 15,23, 16,31, 17,31, 18,31, 18,31&
			&, 19,31, 20,23, 21,31, 22,31, 23,31, 24,31 /)
	ind_quad_loss_extra(18,0:18) = (/ 6, 2,1,1, 3,1,2, 4,1,3, 10,5,1, 15,5,2, 20,5,3 /)
	ind_quad_form(18,0:14) = (/ 7, 1,17, 2,16, 3,15, 5,13, 6,12, 7,11, 8,10 /)
	ind_lin_loss(18,0:16) = (/ 8, 1,17, 2,16, 3,15, 5,13, 6,12, 7,11, 8,10, 26,26 /)
	ind_lin_form(18,0:6) = (/ 3, 1,19, 5,23, 6,24 /)

	! Cluster 19: 4Ii3Io
	ind_quad_loss(19,0:48) = (/ 24, 2,19, 3,19, 4,19, 5,24, 6,31, 7,31, 8,31, 9,31, 10,31&
			&, 11,31, 12,31, 13,31, 14,31, 15,31, 16,31, 17,31, 18,31, 19,31, 19,31&
			&, 20,31, 21,31, 22,31, 23,31, 24,31 /)
	ind_quad_loss_extra(19,0:9) = (/ 3, 2,1,2, 3,1,3, 4,1,4 /)
	ind_quad_form(19,0:56) = (/ 28, 1,18, 2,17, 2,18, 2,19, 3,16, 3,17, 3,18, 3,19, 4,15&
			&, 4,16, 4,17, 4,18, 4,19, 5,14, 6,13, 6,14, 7,12, 7,13, 7,14&
			&, 8,11, 8,12, 8,13, 8,14, 9,10, 9,11, 9,12, 9,13, 9,14 /)
	ind_lin_loss(19,0:20) = (/ 10, 1,18, 2,17, 3,16, 4,15, 5,14, 6,13, 7,12, 8,11, 9,10&
			&, 26,26 /)
	ind_lin_form(19,0:2) = (/ 1, 5,24 /)

	! Cluster 20: 4Io
	ind_quad_loss(20,0:48) = (/ 24, 1,21, 2,22, 3,23, 4,24, 6,21, 7,22, 8,23, 9,31, 10,20&
			&, 11,21, 12,22, 13,23, 14,31, 15,20, 16,21, 17,22, 18,23, 19,31, 20,20&
			&, 20,20, 21,21, 22,22, 23,23, 24,31 /)
	ind_quad_loss_extra(20,0:48) = (/ 16, 6,5,1, 7,5,1, 8,5,1, 10,5,2, 11,5,2, 12,5,2, 13,5,2, 15,5,3, 16,5,3&
			&, 17,5,3, 18,5,3, 20,5,4, 20,5,4, 21,5,4, 22,5,4, 23,5,4 /)
	ind_quad_form(20,0:14) = (/ 7, 5,15, 10,10, 10,15, 10,20, 15,15, 15,20, 20,20 /)
	ind_lin_loss(20,0:6) = (/ 3, 5,15, 10,10, 26,26 /)
	ind_lin_form(20,0:8) = (/ 4, 1,21, 2,22, 3,23, 4,24 /)

	! Cluster 21: 1Ii4Io
	ind_quad_loss(21,0:48) = (/ 24, 1,22, 2,23, 3,24, 4,31, 6,22, 7,23, 8,31, 9,31, 10,21&
			&, 11,22, 12,23, 13,31, 14,31, 15,21, 16,22, 17,23, 18,31, 19,31, 20,21&
			&, 21,22, 21,22, 22,23, 23,31, 24,31 /)
	ind_quad_loss_extra(21,0:36) = (/ 12, 6,5,1, 7,5,1, 10,5,2, 11,5,2, 12,5,2, 15,5,3, 16,5,3, 17,5,3, 20,5,4&
			&, 21,5,4, 21,5,4, 22,5,4 /)
	ind_quad_form(21,0:26) = (/ 13, 1,20, 5,16, 6,15, 6,20, 10,11, 10,16, 10,21, 11,15, 11,20&
			&, 15,16, 15,21, 16,20, 20,21 /)
	ind_lin_loss(21,0:10) = (/ 5, 1,20, 5,16, 6,15, 10,11, 26,26 /)
	ind_lin_form(21,0:6) = (/ 3, 1,22, 2,23, 3,24 /)

	! Cluster 22: 2Ii4Io
	ind_quad_loss(22,0:48) = (/ 24, 1,23, 2,24, 3,31, 4,31, 6,23, 7,31, 8,31, 9,31, 10,22&
			&, 11,23, 12,31, 13,31, 14,31, 15,22, 16,23, 17,31, 18,31, 19,31, 20,22&
			&, 21,23, 22,31, 22,31, 23,31, 24,31 /)
	ind_quad_loss_extra(22,0:21) = (/ 7, 6,5,1, 10,5,2, 11,5,2, 15,5,3, 16,5,3, 20,5,4, 21,5,4 /)
	ind_quad_form(22,0:44) = (/ 22, 1,21, 2,20, 5,17, 6,16, 6,21, 7,15, 7,20, 10,12, 10,17&
			&, 10,22, 11,11, 11,16, 11,21, 12,15, 12,20, 15,17, 15,22, 16,16, 16,21&
			&, 17,20, 20,22, 21,21 /)
	ind_lin_loss(22,0:16) = (/ 8, 1,21, 2,20, 5,17, 6,16, 7,15, 10,12, 11,11, 26,26 /)
	ind_lin_form(22,0:4) = (/ 2, 1,23, 2,24 /)

	! Cluster 23: 3Ii4Io
	ind_quad_loss(23,0:48) = (/ 24, 1,24, 2,31, 3,31, 4,31, 6,31, 7,31, 8,31, 9,31, 10,23&
			&, 11,31, 12,31, 13,31, 14,31, 15,23, 16,31, 17,31, 18,31, 19,31, 20,23&
			&, 21,31, 22,31, 23,31, 23,31, 24,31 /)
	ind_quad_loss_extra(23,0:9) = (/ 3, 10,5,2, 15,5,3, 20,5,4 /)
	ind_quad_form(23,0:56) = (/ 28, 1,22, 2,21, 3,20, 5,18, 6,17, 6,22, 7,16, 7,21, 8,15&
			&, 8,20, 10,13, 10,18, 10,23, 11,12, 11,17, 11,22, 12,16, 12,21, 13,15&
			&, 13,20, 15,18, 15,23, 16,17, 16,22, 17,21, 18,20, 20,23, 21,22 /)
	ind_lin_loss(23,0:20) = (/ 10, 1,22, 2,21, 3,20, 5,18, 6,17, 7,16, 8,15, 10,13, 11,12&
			&, 26,26 /)
	ind_lin_form(23,0:2) = (/ 1, 1,24 /)

	! Cluster 24: 4Ii4Io
	ind_quad_loss(24,0:50) = (/ 25, 1,31, 2,31, 3,31, 4,31, 5,31, 6,31, 7,31, 8,31, 9,31&
			&, 10,31, 11,31, 12,31, 13,31, 14,31, 15,31, 16,31, 17,31, 18,31, 19,31&
			&, 20,31, 21,31, 22,31, 23,31, 24,31, 24,31 /)
	ind_quad_form(24,0:24) = (/ 12, 1,23, 2,22, 3,21, 4,20, 5,19, 6,18, 7,17, 8,16, 9,15&
			&, 10,14, 11,13, 12,12 /)
	ind_lin_loss(24,0:26) = (/ 13, 1,23, 2,22, 3,21, 4,20, 5,19, 6,18, 7,17, 8,16, 9,15&
			&, 10,14, 11,13, 12,12, 26,26 /)

	! Cluster 26: coag
	ind_lin_form(26,0:48) = (/ 24, 26,1, 26,2, 26,3, 26,4, 26,5, 26,6, 26,7, 26,8, 26,9&
			&, 26,10, 26,11, 26,12, 26,13, 26,14, 26,15, 26,16, 26,17, 26,18, 26,19&
			&, 26,20, 26,21, 26,22, 26,23, 26,24 /)

	! Cluster 31: out_neu
	ind_quad_form(31,0:208) = (/ 104, 1,24, 2,23, 2,24, 3,22, 3,23, 3,24, 4,21, 4,22, 4,23&
			&, 4,24, 5,24, 6,19, 6,23, 6,24, 7,18, 7,19, 7,22, 7,23, 7,24&
			&, 8,17, 8,18, 8,19, 8,21, 8,22, 8,23, 8,24, 9,16, 9,17, 9,18&
			&, 9,19, 9,20, 9,21, 9,22, 9,23, 9,24, 10,19, 10,24, 11,14, 11,18&
			&, 11,19, 11,23, 11,24, 12,13, 12,14, 12,17, 12,18, 12,19, 12,22, 12,23&
			&, 12,24, 13,13, 13,14, 13,16, 13,17, 13,18, 13,19, 13,21, 13,22, 13,23&
			&, 13,24, 14,14, 14,15, 14,16, 14,17, 14,18, 14,19, 14,20, 14,21, 14,22&
			&, 14,23, 14,24, 15,19, 15,24, 16,18, 16,19, 16,23, 16,24, 17,17, 17,18&
			&, 17,19, 17,22, 17,23, 17,24, 18,18, 18,19, 18,21, 18,22, 18,23, 18,24&
			&, 19,19, 19,20, 19,21, 19,22, 19,23, 19,24, 20,24, 21,23, 21,24, 22,22&
			&, 22,23, 22,24, 23,23, 23,24, 24,24 /)

	! Index of the outgrowing cluster for all collisions out
	call get_ind_out_clust_3(ind_out_clust)

	call get_rate_coefs_3(coef_quad,coef_lin,coef)

end subroutine initialize_parameters_3

!-----------------------------------------------------------

subroutine get_ind_out_clust_3(ind_out_clust)
	implicit none
	integer :: ind_out_clust(24,24)

	ind_out_clust = 0

	ind_out_clust(1,24) = 1
	ind_out_clust(24,1) = 1
	ind_out_clust(2,23) = 1
	ind_out_clust(23,2) = 1
	ind_out_clust(2,24) = 2
	ind_out_clust(24,2) = 2
	ind_out_clust(3,22) = 1
	ind_out_clust(22,3) = 1
	ind_out_clust(3,23) = 2
	ind_out_clust(23,3) = 2
	ind_out_clust(3,24) = 3
	ind_out_clust(24,3) = 3
	ind_out_clust(4,21) = 1
	ind_out_clust(21,4) = 1
	ind_out_clust(4,22) = 2
	ind_out_clust(22,4) = 2
	ind_out_clust(4,23) = 3
	ind_out_clust(23,4) = 3
	ind_out_clust(4,24) = 4
	ind_out_clust(24,4) = 4
	ind_out_clust(5,24) = 5
	ind_out_clust(24,5) = 5
	ind_out_clust(6,19) = 1
	ind_out_clust(19,6) = 1
	ind_out_clust(6,23) = 5
	ind_out_clust(23,6) = 5
	ind_out_clust(6,24) = 6
	ind_out_clust(24,6) = 6
	ind_out_clust(7,18) = 1
	ind_out_clust(18,7) = 1
	ind_out_clust(7,19) = 2
	ind_out_clust(19,7) = 2
	ind_out_clust(7,22) = 5
	ind_out_clust(22,7) = 5
	ind_out_clust(7,23) = 6
	ind_out_clust(23,7) = 6
	ind_out_clust(7,24) = 7
	ind_out_clust(24,7) = 7
	ind_out_clust(8,17) = 1
	ind_out_clust(17,8) = 1
	ind_out_clust(8,18) = 2
	ind_out_clust(18,8) = 2
	ind_out_clust(8,19) = 3
	ind_out_clust(19,8) = 3
	ind_out_clust(8,21) = 5
	ind_out_clust(21,8) = 5
	ind_out_clust(8,22) = 6
	ind_out_clust(22,8) = 6
	ind_out_clust(8,23) = 7
	ind_out_clust(23,8) = 7
	ind_out_clust(8,24) = 8
	ind_out_clust(24,8) = 8
	ind_out_clust(9,16) = 1
	ind_out_clust(16,9) = 1
	ind_out_clust(9,17) = 2
	ind_out_clust(17,9) = 2
	ind_out_clust(9,18) = 3
	ind_out_clust(18,9) = 3
	ind_out_clust(9,19) = 4
	ind_out_clust(19,9) = 4
	ind_out_clust(9,20) = 5
	ind_out_clust(20,9) = 5
	ind_out_clust(9,21) = 6
	ind_out_clust(21,9) = 6
	ind_out_clust(9,22) = 7
	ind_out_clust(22,9) = 7
	ind_out_clust(9,23) = 8
	ind_out_clust(23,9) = 8
	ind_out_clust(9,24) = 9
	ind_out_clust(24,9) = 9
	ind_out_clust(10,19) = 5
	ind_out_clust(19,10) = 5
	ind_out_clust(10,24) = 10
	ind_out_clust(24,10) = 10
	ind_out_clust(11,14) = 1
	ind_out_clust(14,11) = 1
	ind_out_clust(11,18) = 5
	ind_out_clust(18,11) = 5
	ind_out_clust(11,19) = 6
	ind_out_clust(19,11) = 6
	ind_out_clust(11,23) = 10
	ind_out_clust(23,11) = 10
	ind_out_clust(11,24) = 11
	ind_out_clust(24,11) = 11
	ind_out_clust(12,13) = 1
	ind_out_clust(13,12) = 1
	ind_out_clust(12,14) = 2
	ind_out_clust(14,12) = 2
	ind_out_clust(12,17) = 5
	ind_out_clust(17,12) = 5
	ind_out_clust(12,18) = 6
	ind_out_clust(18,12) = 6
	ind_out_clust(12,19) = 7
	ind_out_clust(19,12) = 7
	ind_out_clust(12,22) = 10
	ind_out_clust(22,12) = 10
	ind_out_clust(12,23) = 11
	ind_out_clust(23,12) = 11
	ind_out_clust(12,24) = 12
	ind_out_clust(24,12) = 12
	ind_out_clust(13,13) = 2
	ind_out_clust(13,13) = 2
	ind_out_clust(13,14) = 3
	ind_out_clust(14,13) = 3
	ind_out_clust(13,16) = 5
	ind_out_clust(16,13) = 5
	ind_out_clust(13,17) = 6
	ind_out_clust(17,13) = 6
	ind_out_clust(13,18) = 7
	ind_out_clust(18,13) = 7
	ind_out_clust(13,19) = 8
	ind_out_clust(19,13) = 8
	ind_out_clust(13,21) = 10
	ind_out_clust(21,13) = 10
	ind_out_clust(13,22) = 11
	ind_out_clust(22,13) = 11
	ind_out_clust(13,23) = 12
	ind_out_clust(23,13) = 12
	ind_out_clust(13,24) = 13
	ind_out_clust(24,13) = 13
	ind_out_clust(14,14) = 4
	ind_out_clust(14,14) = 4
	ind_out_clust(14,15) = 5
	ind_out_clust(15,14) = 5
	ind_out_clust(14,16) = 6
	ind_out_clust(16,14) = 6
	ind_out_clust(14,17) = 7
	ind_out_clust(17,14) = 7
	ind_out_clust(14,18) = 8
	ind_out_clust(18,14) = 8
	ind_out_clust(14,19) = 9
	ind_out_clust(19,14) = 9
	ind_out_clust(14,20) = 10
	ind_out_clust(20,14) = 10
	ind_out_clust(14,21) = 11
	ind_out_clust(21,14) = 11
	ind_out_clust(14,22) = 12
	ind_out_clust(22,14) = 12
	ind_out_clust(14,23) = 13
	ind_out_clust(23,14) = 13
	ind_out_clust(14,24) = 14
	ind_out_clust(24,14) = 14
	ind_out_clust(15,19) = 10
	ind_out_clust(19,15) = 10
	ind_out_clust(15,24) = 15
	ind_out_clust(24,15) = 15
	ind_out_clust(16,18) = 10
	ind_out_clust(18,16) = 10
	ind_out_clust(16,19) = 11
	ind_out_clust(19,16) = 11
	ind_out_clust(16,23) = 15
	ind_out_clust(23,16) = 15
	ind_out_clust(16,24) = 16
	ind_out_clust(24,16) = 16
	ind_out_clust(17,17) = 10
	ind_out_clust(17,17) = 10
	ind_out_clust(17,18) = 11
	ind_out_clust(18,17) = 11
	ind_out_clust(17,19) = 12
	ind_out_clust(19,17) = 12
	ind_out_clust(17,22) = 15
	ind_out_clust(22,17) = 15
	ind_out_clust(17,23) = 16
	ind_out_clust(23,17) = 16
	ind_out_clust(17,24) = 17
	ind_out_clust(24,17) = 17
	ind_out_clust(18,18) = 12
	ind_out_clust(18,18) = 12
	ind_out_clust(18,19) = 13
	ind_out_clust(19,18) = 13
	ind_out_clust(18,21) = 15
	ind_out_clust(21,18) = 15
	ind_out_clust(18,22) = 16
	ind_out_clust(22,18) = 16
	ind_out_clust(18,23) = 17
	ind_out_clust(23,18) = 17
	ind_out_clust(18,24) = 18
	ind_out_clust(24,18) = 18
	ind_out_clust(19,19) = 14
	ind_out_clust(19,19) = 14
	ind_out_clust(19,20) = 15
	ind_out_clust(20,19) = 15
	ind_out_clust(19,21) = 16
	ind_out_clust(21,19) = 16
	ind_out_clust(19,22) = 17
	ind_out_clust(22,19) = 17
	ind_out_clust(19,23) = 18
	ind_out_clust(23,19) = 18
	ind_out_clust(19,24) = 19
	ind_out_clust(24,19) = 19
	ind_out_clust(20,24) = 20
	ind_out_clust(24,20) = 20
	ind_out_clust(21,23) = 20
	ind_out_clust(23,21) = 20
	ind_out_clust(21,24) = 21
	ind_out_clust(24,21) = 21
	ind_out_clust(22,22) = 20
	ind_out_clust(22,22) = 20
	ind_out_clust(22,23) = 21
	ind_out_clust(23,22) = 21
	ind_out_clust(22,24) = 22
	ind_out_clust(24,22) = 22
	ind_out_clust(23,23) = 22
	ind_out_clust(23,23) = 22
	ind_out_clust(23,24) = 23
	ind_out_clust(24,23) = 23
	ind_out_clust(24,24) = 24
	ind_out_clust(24,24) = 24

end subroutine get_ind_out_clust_3

!-----------------------------------------------------------

subroutine get_rate_coefs_3(coef_quad,coef_lin,coef)
	implicit none
	real(kind(1.d0)) :: coef_quad(24,24,32),coef_lin(32,32,24)
	real(kind(1.d0)) :: K(24,24),E(24,24),cs(24)
	real(kind(1.d0)) :: coef(1)

	coef_quad = 0.d0
	coef_lin = 0.d0
	call get_coll_3(K,coef(1))
	call get_evap_3(E,K,coef(1))
	call get_losses_3(cs,coef(1))


	coef_quad(1,1,2) = 0.5d0*K(1,1)	! 1Ii + 1Ii -> 2Ii
	coef_lin(1,1,2) = E(1,1)	! 2Ii -> 1Ii + 1Ii

	coef_quad(2,1,3) = K(1,2)	! 2Ii + 1Ii -> 3Ii
	coef_quad(1,2,3) = K(1,2)
	coef_lin(2,1,3) = E(1,2)	! 3Ii -> 2Ii + 1Ii
	coef_lin(1,2,3) = E(1,2)
	coef_quad(2,2,4) = 0.5d0*K(2,2)	! 2Ii + 2Ii -> 4Ii
	coef_lin(2,2,4) = E(2,2)	! 4Ii -> 2Ii + 2Ii

	coef_quad(3,1,4) = K(1,3)	! 3Ii + 1Ii -> 4Ii
	coef_quad(1,3,4) = K(1,3)
	coef_lin(3,1,4) = E(1,3)	! 4Ii -> 3Ii + 1Ii
	coef_lin(1,3,4) = E(1,3)
	coef_quad(3,2,4) = K(2,3)	! 3Ii + 2Ii -> boundary -> 4Ii
	coef_quad(2,3,4) = K(2,3)
	coef_quad(3,2,1) = K(2,3)	! 3Ii + 2Ii -> boundary -> 1Ii
	coef_quad(2,3,1) = K(2,3)
	coef_quad(3,3,4) = 0.5d0*K(3,3)	! 3Ii + 3Ii -> boundary -> 4Ii
	coef_quad(3,3,1) = 0.5d0*K(3,3)	! 3Ii + 3Ii -> boundary -> 1Ii

	coef_quad(4,2,4) = K(2,4)	! 4Ii + 2Ii -> boundary -> 4Ii
	coef_quad(2,4,4) = K(2,4)
	coef_quad(4,2,1) = K(2,4)	! 4Ii + 2Ii -> boundary -> 1Ii
	coef_quad(2,4,1) = K(2,4)
	coef_quad(4,3,4) = K(3,4)	! 4Ii + 3Ii -> boundary -> 4Ii
	coef_quad(3,4,4) = K(3,4)
	coef_quad(4,3,1) = K(3,4)	! 4Ii + 3Ii -> boundary -> 1Ii
	coef_quad(3,4,1) = K(3,4)
	coef_quad(4,4,4) = 0.5d0*K(4,4)	! 4Ii + 4Ii -> boundary -> 4Ii
	coef_quad(4,4,1) = 0.5d0*K(4,4)	! 4Ii + 4Ii -> boundary -> 1Ii

	coef_quad(5,1,6) = K(1,5)	! 1Io + 1Ii -> 1Ii1Io
	coef_quad(1,5,6) = K(1,5)
	coef_lin(5,1,6) = E(1,5)	! 1Ii1Io -> 1Io + 1Ii
	coef_lin(1,5,6) = E(1,5)
	coef_quad(5,2,7) = K(2,5)	! 1Io + 2Ii -> 2Ii1Io
	coef_quad(2,5,7) = K(2,5)
	coef_lin(5,2,7) = E(2,5)	! 2Ii1Io -> 1Io + 2Ii
	coef_lin(2,5,7) = E(2,5)
	coef_quad(5,3,8) = K(3,5)	! 1Io + 3Ii -> 3Ii1Io
	coef_quad(3,5,8) = K(3,5)
	coef_lin(5,3,8) = E(3,5)	! 3Ii1Io -> 1Io + 3Ii
	coef_lin(3,5,8) = E(3,5)
	coef_quad(5,4,9) = K(4,5)	! 1Io + 4Ii -> 4Ii1Io
	coef_quad(4,5,9) = K(4,5)
	coef_lin(5,4,9) = E(4,5)	! 4Ii1Io -> 1Io + 4Ii
	coef_lin(4,5,9) = E(4,5)
	coef_quad(5,5,10) = 0.5d0*K(5,5)	! 1Io + 1Io -> 2Io
	coef_lin(5,5,10) = E(5,5)	! 2Io -> 1Io + 1Io

	coef_quad(6,1,7) = K(1,6)	! 1Ii1Io + 1Ii -> 2Ii1Io
	coef_quad(1,6,7) = K(1,6)
	coef_lin(6,1,7) = E(1,6)	! 2Ii1Io -> 1Ii1Io + 1Ii
	coef_lin(1,6,7) = E(1,6)
	coef_quad(6,2,8) = K(2,6)	! 1Ii1Io + 2Ii -> 3Ii1Io
	coef_quad(2,6,8) = K(2,6)
	coef_lin(6,2,8) = E(2,6)	! 3Ii1Io -> 1Ii1Io + 2Ii
	coef_lin(2,6,8) = E(2,6)
	coef_quad(6,3,9) = K(3,6)	! 1Ii1Io + 3Ii -> 4Ii1Io
	coef_quad(3,6,9) = K(3,6)
	coef_lin(6,3,9) = E(3,6)	! 4Ii1Io -> 1Ii1Io + 3Ii
	coef_lin(3,6,9) = E(3,6)
	coef_quad(6,4,9) = K(4,6)	! 1Ii1Io + 4Ii -> boundary -> 4Ii1Io
	coef_quad(4,6,9) = K(4,6)
	coef_quad(6,4,1) = K(4,6)	! 1Ii1Io + 4Ii -> boundary -> 1Ii
	coef_quad(4,6,1) = K(4,6)
	coef_quad(6,5,11) = K(5,6)	! 1Ii1Io + 1Io -> 1Ii2Io
	coef_quad(5,6,11) = K(5,6)
	coef_lin(6,5,11) = E(5,6)	! 1Ii2Io -> 1Ii1Io + 1Io
	coef_lin(5,6,11) = E(5,6)
	coef_quad(6,6,12) = 0.5d0*K(6,6)	! 1Ii1Io + 1Ii1Io -> 2Ii2Io
	coef_lin(6,6,12) = E(6,6)	! 2Ii2Io -> 1Ii1Io + 1Ii1Io

	coef_quad(7,1,8) = K(1,7)	! 2Ii1Io + 1Ii -> 3Ii1Io
	coef_quad(1,7,8) = K(1,7)
	coef_lin(7,1,8) = E(1,7)	! 3Ii1Io -> 2Ii1Io + 1Ii
	coef_lin(1,7,8) = E(1,7)
	coef_quad(7,2,9) = K(2,7)	! 2Ii1Io + 2Ii -> 4Ii1Io
	coef_quad(2,7,9) = K(2,7)
	coef_lin(7,2,9) = E(2,7)	! 4Ii1Io -> 2Ii1Io + 2Ii
	coef_lin(2,7,9) = E(2,7)
	coef_quad(7,3,9) = K(3,7)	! 2Ii1Io + 3Ii -> boundary -> 4Ii1Io
	coef_quad(3,7,9) = K(3,7)
	coef_quad(7,3,1) = K(3,7)	! 2Ii1Io + 3Ii -> boundary -> 1Ii
	coef_quad(3,7,1) = K(3,7)
	coef_quad(7,4,9) = K(4,7)	! 2Ii1Io + 4Ii -> boundary -> 4Ii1Io
	coef_quad(4,7,9) = K(4,7)
	coef_quad(7,4,1) = K(4,7)	! 2Ii1Io + 4Ii -> boundary -> 1Ii
	coef_quad(4,7,1) = K(4,7)
	coef_quad(7,5,12) = K(5,7)	! 2Ii1Io + 1Io -> 2Ii2Io
	coef_quad(5,7,12) = K(5,7)
	coef_lin(7,5,12) = E(5,7)	! 2Ii2Io -> 2Ii1Io + 1Io
	coef_lin(5,7,12) = E(5,7)
	coef_quad(7,6,13) = K(6,7)	! 2Ii1Io + 1Ii1Io -> 3Ii2Io
	coef_quad(6,7,13) = K(6,7)
	coef_lin(7,6,13) = E(6,7)	! 3Ii2Io -> 2Ii1Io + 1Ii1Io
	coef_lin(6,7,13) = E(6,7)
	coef_quad(7,7,14) = 0.5d0*K(7,7)	! 2Ii1Io + 2Ii1Io -> 4Ii2Io
	coef_lin(7,7,14) = E(7,7)	! 4Ii2Io -> 2Ii1Io + 2Ii1Io

	coef_quad(8,1,9) = K(1,8)	! 3Ii1Io + 1Ii -> 4Ii1Io
	coef_quad(1,8,9) = K(1,8)
	coef_lin(8,1,9) = E(1,8)	! 4Ii1Io -> 3Ii1Io + 1Ii
	coef_lin(1,8,9) = E(1,8)
	coef_quad(8,2,9) = K(2,8)	! 3Ii1Io + 2Ii -> boundary -> 4Ii1Io
	coef_quad(2,8,9) = K(2,8)
	coef_quad(8,2,1) = K(2,8)	! 3Ii1Io + 2Ii -> boundary -> 1Ii
	coef_quad(2,8,1) = K(2,8)
	coef_quad(8,3,9) = K(3,8)	! 3Ii1Io + 3Ii -> boundary -> 4Ii1Io
	coef_quad(3,8,9) = K(3,8)
	coef_quad(8,3,1) = K(3,8)	! 3Ii1Io + 3Ii -> boundary -> 1Ii
	coef_quad(3,8,1) = K(3,8)
	coef_quad(8,4,9) = K(4,8)	! 3Ii1Io + 4Ii -> boundary -> 4Ii1Io
	coef_quad(4,8,9) = K(4,8)
	coef_quad(8,4,1) = K(4,8)	! 3Ii1Io + 4Ii -> boundary -> 1Ii
	coef_quad(4,8,1) = K(4,8)
	coef_quad(8,5,13) = K(5,8)	! 3Ii1Io + 1Io -> 3Ii2Io
	coef_quad(5,8,13) = K(5,8)
	coef_lin(8,5,13) = E(5,8)	! 3Ii2Io -> 3Ii1Io + 1Io
	coef_lin(5,8,13) = E(5,8)
	coef_quad(8,6,14) = K(6,8)	! 3Ii1Io + 1Ii1Io -> 4Ii2Io
	coef_quad(6,8,14) = K(6,8)
	coef_lin(8,6,14) = E(6,8)	! 4Ii2Io -> 3Ii1Io + 1Ii1Io
	coef_lin(6,8,14) = E(6,8)
	coef_quad(8,7,14) = K(7,8)	! 3Ii1Io + 2Ii1Io -> boundary -> 4Ii2Io
	coef_quad(7,8,14) = K(7,8)
	coef_quad(8,7,1) = K(7,8)	! 3Ii1Io + 2Ii1Io -> boundary -> 1Ii
	coef_quad(7,8,1) = K(7,8)
	coef_quad(8,8,14) = 0.5d0*K(8,8)	! 3Ii1Io + 3Ii1Io -> boundary -> 4Ii2Io
	coef_quad(8,8,1) = 0.5d0*K(8,8)	! 3Ii1Io + 3Ii1Io -> boundary -> 1Ii

	coef_quad(9,2,9) = K(2,9)	! 4Ii1Io + 2Ii -> boundary -> 4Ii1Io
	coef_quad(2,9,9) = K(2,9)
	coef_quad(9,2,1) = K(2,9)	! 4Ii1Io + 2Ii -> boundary -> 1Ii
	coef_quad(2,9,1) = K(2,9)
	coef_quad(9,3,9) = K(3,9)	! 4Ii1Io + 3Ii -> boundary -> 4Ii1Io
	coef_quad(3,9,9) = K(3,9)
	coef_quad(9,3,1) = K(3,9)	! 4Ii1Io + 3Ii -> boundary -> 1Ii
	coef_quad(3,9,1) = K(3,9)
	coef_quad(9,4,9) = K(4,9)	! 4Ii1Io + 4Ii -> boundary -> 4Ii1Io
	coef_quad(4,9,9) = K(4,9)
	coef_quad(9,4,1) = K(4,9)	! 4Ii1Io + 4Ii -> boundary -> 1Ii
	coef_quad(4,9,1) = K(4,9)
	coef_quad(9,5,14) = K(5,9)	! 4Ii1Io + 1Io -> 4Ii2Io
	coef_quad(5,9,14) = K(5,9)
	coef_lin(9,5,14) = E(5,9)	! 4Ii2Io -> 4Ii1Io + 1Io
	coef_lin(5,9,14) = E(5,9)
	coef_quad(9,6,14) = K(6,9)	! 4Ii1Io + 1Ii1Io -> boundary -> 4Ii2Io
	coef_quad(6,9,14) = K(6,9)
	coef_quad(9,6,1) = K(6,9)	! 4Ii1Io + 1Ii1Io -> boundary -> 1Ii
	coef_quad(6,9,1) = K(6,9)
	coef_quad(9,7,14) = K(7,9)	! 4Ii1Io + 2Ii1Io -> boundary -> 4Ii2Io
	coef_quad(7,9,14) = K(7,9)
	coef_quad(9,7,1) = K(7,9)	! 4Ii1Io + 2Ii1Io -> boundary -> 1Ii
	coef_quad(7,9,1) = K(7,9)
	coef_quad(9,8,14) = K(8,9)	! 4Ii1Io + 3Ii1Io -> boundary -> 4Ii2Io
	coef_quad(8,9,14) = K(8,9)
	coef_quad(9,8,1) = K(8,9)	! 4Ii1Io + 3Ii1Io -> boundary -> 1Ii
	coef_quad(8,9,1) = K(8,9)
	coef_quad(9,9,14) = 0.5d0*K(9,9)	! 4Ii1Io + 4Ii1Io -> boundary -> 4Ii2Io
	coef_quad(9,9,1) = 0.5d0*K(9,9)	! 4Ii1Io + 4Ii1Io -> boundary -> 1Ii

	coef_quad(10,1,11) = K(1,10)	! 2Io + 1Ii -> 1Ii2Io
	coef_quad(1,10,11) = K(1,10)
	coef_lin(10,1,11) = E(1,10)	! 1Ii2Io -> 2Io + 1Ii
	coef_lin(1,10,11) = E(1,10)
	coef_quad(10,2,12) = K(2,10)	! 2Io + 2Ii -> 2Ii2Io
	coef_quad(2,10,12) = K(2,10)
	coef_lin(10,2,12) = E(2,10)	! 2Ii2Io -> 2Io + 2Ii
	coef_lin(2,10,12) = E(2,10)
	coef_quad(10,3,13) = K(3,10)	! 2Io + 3Ii -> 3Ii2Io
	coef_quad(3,10,13) = K(3,10)
	coef_lin(10,3,13) = E(3,10)	! 3Ii2Io -> 2Io + 3Ii
	coef_lin(3,10,13) = E(3,10)
	coef_quad(10,4,14) = K(4,10)	! 2Io + 4Ii -> 4Ii2Io
	coef_quad(4,10,14) = K(4,10)
	coef_lin(10,4,14) = E(4,10)	! 4Ii2Io -> 2Io + 4Ii
	coef_lin(4,10,14) = E(4,10)
	coef_quad(10,5,15) = K(5,10)	! 2Io + 1Io -> 3Io
	coef_quad(5,10,15) = K(5,10)
	coef_lin(10,5,15) = E(5,10)	! 3Io -> 2Io + 1Io
	coef_lin(5,10,15) = E(5,10)
	coef_quad(10,6,16) = K(6,10)	! 2Io + 1Ii1Io -> 1Ii3Io
	coef_quad(6,10,16) = K(6,10)
	coef_lin(10,6,16) = E(6,10)	! 1Ii3Io -> 2Io + 1Ii1Io
	coef_lin(6,10,16) = E(6,10)
	coef_quad(10,7,17) = K(7,10)	! 2Io + 2Ii1Io -> 2Ii3Io
	coef_quad(7,10,17) = K(7,10)
	coef_lin(10,7,17) = E(7,10)	! 2Ii3Io -> 2Io + 2Ii1Io
	coef_lin(7,10,17) = E(7,10)
	coef_quad(10,8,18) = K(8,10)	! 2Io + 3Ii1Io -> 3Ii3Io
	coef_quad(8,10,18) = K(8,10)
	coef_lin(10,8,18) = E(8,10)	! 3Ii3Io -> 2Io + 3Ii1Io
	coef_lin(8,10,18) = E(8,10)
	coef_quad(10,9,19) = K(9,10)	! 2Io + 4Ii1Io -> 4Ii3Io
	coef_quad(9,10,19) = K(9,10)
	coef_lin(10,9,19) = E(9,10)	! 4Ii3Io -> 2Io + 4Ii1Io
	coef_lin(9,10,19) = E(9,10)
	coef_quad(10,10,20) = 0.5d0*K(10,10)	! 2Io + 2Io -> 4Io
	coef_lin(10,10,20) = E(10,10)	! 4Io -> 2Io + 2Io

	coef_quad(11,1,12) = K(1,11)	! 1Ii2Io + 1Ii -> 2Ii2Io
	coef_quad(1,11,12) = K(1,11)
	coef_lin(11,1,12) = E(1,11)	! 2Ii2Io -> 1Ii2Io + 1Ii
	coef_lin(1,11,12) = E(1,11)
	coef_quad(11,2,13) = K(2,11)	! 1Ii2Io + 2Ii -> 3Ii2Io
	coef_quad(2,11,13) = K(2,11)
	coef_lin(11,2,13) = E(2,11)	! 3Ii2Io -> 1Ii2Io + 2Ii
	coef_lin(2,11,13) = E(2,11)
	coef_quad(11,3,14) = K(3,11)	! 1Ii2Io + 3Ii -> 4Ii2Io
	coef_quad(3,11,14) = K(3,11)
	coef_lin(11,3,14) = E(3,11)	! 4Ii2Io -> 1Ii2Io + 3Ii
	coef_lin(3,11,14) = E(3,11)
	coef_quad(11,4,14) = K(4,11)	! 1Ii2Io + 4Ii -> boundary -> 4Ii2Io
	coef_quad(4,11,14) = K(4,11)
	coef_quad(11,4,1) = K(4,11)	! 1Ii2Io + 4Ii -> boundary -> 1Ii
	coef_quad(4,11,1) = K(4,11)
	coef_quad(11,5,16) = K(5,11)	! 1Ii2Io + 1Io -> 1Ii3Io
	coef_quad(5,11,16) = K(5,11)
	coef_lin(11,5,16) = E(5,11)	! 1Ii3Io -> 1Ii2Io + 1Io
	coef_lin(5,11,16) = E(5,11)
	coef_quad(11,6,17) = K(6,11)	! 1Ii2Io + 1Ii1Io -> 2Ii3Io
	coef_quad(6,11,17) = K(6,11)
	coef_lin(11,6,17) = E(6,11)	! 2Ii3Io -> 1Ii2Io + 1Ii1Io
	coef_lin(6,11,17) = E(6,11)
	coef_quad(11,7,18) = K(7,11)	! 1Ii2Io + 2Ii1Io -> 3Ii3Io
	coef_quad(7,11,18) = K(7,11)
	coef_lin(11,7,18) = E(7,11)	! 3Ii3Io -> 1Ii2Io + 2Ii1Io
	coef_lin(7,11,18) = E(7,11)
	coef_quad(11,8,19) = K(8,11)	! 1Ii2Io + 3Ii1Io -> 4Ii3Io
	coef_quad(8,11,19) = K(8,11)
	coef_lin(11,8,19) = E(8,11)	! 4Ii3Io -> 1Ii2Io + 3Ii1Io
	coef_lin(8,11,19) = E(8,11)
	coef_quad(11,9,19) = K(9,11)	! 1Ii2Io + 4Ii1Io -> boundary -> 4Ii3Io
	coef_quad(9,11,19) = K(9,11)
	coef_quad(11,9,1) = K(9,11)	! 1Ii2Io + 4Ii1Io -> boundary -> 1Ii
	coef_quad(9,11,1) = K(9,11)
	coef_quad(11,10,21) = K(10,11)	! 1Ii2Io + 2Io -> 1Ii4Io
	coef_quad(10,11,21) = K(10,11)
	coef_lin(11,10,21) = E(10,11)	! 1Ii4Io -> 1Ii2Io + 2Io
	coef_lin(10,11,21) = E(10,11)
	coef_quad(11,11,22) = 0.5d0*K(11,11)	! 1Ii2Io + 1Ii2Io -> 2Ii4Io
	coef_lin(11,11,22) = E(11,11)	! 2Ii4Io -> 1Ii2Io + 1Ii2Io

	coef_quad(12,1,13) = K(1,12)	! 2Ii2Io + 1Ii -> 3Ii2Io
	coef_quad(1,12,13) = K(1,12)
	coef_lin(12,1,13) = E(1,12)	! 3Ii2Io -> 2Ii2Io + 1Ii
	coef_lin(1,12,13) = E(1,12)
	coef_quad(12,2,14) = K(2,12)	! 2Ii2Io + 2Ii -> 4Ii2Io
	coef_quad(2,12,14) = K(2,12)
	coef_lin(12,2,14) = E(2,12)	! 4Ii2Io -> 2Ii2Io + 2Ii
	coef_lin(2,12,14) = E(2,12)
	coef_quad(12,3,14) = K(3,12)	! 2Ii2Io + 3Ii -> boundary -> 4Ii2Io
	coef_quad(3,12,14) = K(3,12)
	coef_quad(12,3,1) = K(3,12)	! 2Ii2Io + 3Ii -> boundary -> 1Ii
	coef_quad(3,12,1) = K(3,12)
	coef_quad(12,4,14) = K(4,12)	! 2Ii2Io + 4Ii -> boundary -> 4Ii2Io
	coef_quad(4,12,14) = K(4,12)
	coef_quad(12,4,1) = K(4,12)	! 2Ii2Io + 4Ii -> boundary -> 1Ii
	coef_quad(4,12,1) = K(4,12)
	coef_quad(12,5,17) = K(5,12)	! 2Ii2Io + 1Io -> 2Ii3Io
	coef_quad(5,12,17) = K(5,12)
	coef_lin(12,5,17) = E(5,12)	! 2Ii3Io -> 2Ii2Io + 1Io
	coef_lin(5,12,17) = E(5,12)
	coef_quad(12,6,18) = K(6,12)	! 2Ii2Io + 1Ii1Io -> 3Ii3Io
	coef_quad(6,12,18) = K(6,12)
	coef_lin(12,6,18) = E(6,12)	! 3Ii3Io -> 2Ii2Io + 1Ii1Io
	coef_lin(6,12,18) = E(6,12)
	coef_quad(12,7,19) = K(7,12)	! 2Ii2Io + 2Ii1Io -> 4Ii3Io
	coef_quad(7,12,19) = K(7,12)
	coef_lin(12,7,19) = E(7,12)	! 4Ii3Io -> 2Ii2Io + 2Ii1Io
	coef_lin(7,12,19) = E(7,12)
	coef_quad(12,8,19) = K(8,12)	! 2Ii2Io + 3Ii1Io -> boundary -> 4Ii3Io
	coef_quad(8,12,19) = K(8,12)
	coef_quad(12,8,1) = K(8,12)	! 2Ii2Io + 3Ii1Io -> boundary -> 1Ii
	coef_quad(8,12,1) = K(8,12)
	coef_quad(12,9,19) = K(9,12)	! 2Ii2Io + 4Ii1Io -> boundary -> 4Ii3Io
	coef_quad(9,12,19) = K(9,12)
	coef_quad(12,9,1) = K(9,12)	! 2Ii2Io + 4Ii1Io -> boundary -> 1Ii
	coef_quad(9,12,1) = K(9,12)
	coef_quad(12,10,22) = K(10,12)	! 2Ii2Io + 2Io -> 2Ii4Io
	coef_quad(10,12,22) = K(10,12)
	coef_lin(12,10,22) = E(10,12)	! 2Ii4Io -> 2Ii2Io + 2Io
	coef_lin(10,12,22) = E(10,12)
	coef_quad(12,11,23) = K(11,12)	! 2Ii2Io + 1Ii2Io -> 3Ii4Io
	coef_quad(11,12,23) = K(11,12)
	coef_lin(12,11,23) = E(11,12)	! 3Ii4Io -> 2Ii2Io + 1Ii2Io
	coef_lin(11,12,23) = E(11,12)
	coef_quad(12,12,24) = 0.5d0*K(12,12)	! 2Ii2Io + 2Ii2Io -> 4Ii4Io
	coef_lin(12,12,24) = E(12,12)	! 4Ii4Io -> 2Ii2Io + 2Ii2Io

	coef_quad(13,1,14) = K(1,13)	! 3Ii2Io + 1Ii -> 4Ii2Io
	coef_quad(1,13,14) = K(1,13)
	coef_lin(13,1,14) = E(1,13)	! 4Ii2Io -> 3Ii2Io + 1Ii
	coef_lin(1,13,14) = E(1,13)
	coef_quad(13,2,14) = K(2,13)	! 3Ii2Io + 2Ii -> boundary -> 4Ii2Io
	coef_quad(2,13,14) = K(2,13)
	coef_quad(13,2,1) = K(2,13)	! 3Ii2Io + 2Ii -> boundary -> 1Ii
	coef_quad(2,13,1) = K(2,13)
	coef_quad(13,3,14) = K(3,13)	! 3Ii2Io + 3Ii -> boundary -> 4Ii2Io
	coef_quad(3,13,14) = K(3,13)
	coef_quad(13,3,1) = K(3,13)	! 3Ii2Io + 3Ii -> boundary -> 1Ii
	coef_quad(3,13,1) = K(3,13)
	coef_quad(13,4,14) = K(4,13)	! 3Ii2Io + 4Ii -> boundary -> 4Ii2Io
	coef_quad(4,13,14) = K(4,13)
	coef_quad(13,4,1) = K(4,13)	! 3Ii2Io + 4Ii -> boundary -> 1Ii
	coef_quad(4,13,1) = K(4,13)
	coef_quad(13,5,18) = K(5,13)	! 3Ii2Io + 1Io -> 3Ii3Io
	coef_quad(5,13,18) = K(5,13)
	coef_lin(13,5,18) = E(5,13)	! 3Ii3Io -> 3Ii2Io + 1Io
	coef_lin(5,13,18) = E(5,13)
	coef_quad(13,6,19) = K(6,13)	! 3Ii2Io + 1Ii1Io -> 4Ii3Io
	coef_quad(6,13,19) = K(6,13)
	coef_lin(13,6,19) = E(6,13)	! 4Ii3Io -> 3Ii2Io + 1Ii1Io
	coef_lin(6,13,19) = E(6,13)
	coef_quad(13,7,19) = K(7,13)	! 3Ii2Io + 2Ii1Io -> boundary -> 4Ii3Io
	coef_quad(7,13,19) = K(7,13)
	coef_quad(13,7,1) = K(7,13)	! 3Ii2Io + 2Ii1Io -> boundary -> 1Ii
	coef_quad(7,13,1) = K(7,13)
	coef_quad(13,8,19) = K(8,13)	! 3Ii2Io + 3Ii1Io -> boundary -> 4Ii3Io
	coef_quad(8,13,19) = K(8,13)
	coef_quad(13,8,1) = K(8,13)	! 3Ii2Io + 3Ii1Io -> boundary -> 1Ii
	coef_quad(8,13,1) = K(8,13)
	coef_quad(13,9,19) = K(9,13)	! 3Ii2Io + 4Ii1Io -> boundary -> 4Ii3Io
	coef_quad(9,13,19) = K(9,13)
	coef_quad(13,9,1) = K(9,13)	! 3Ii2Io + 4Ii1Io -> boundary -> 1Ii
	coef_quad(9,13,1) = K(9,13)
	coef_quad(13,10,23) = K(10,13)	! 3Ii2Io + 2Io -> 3Ii4Io
	coef_quad(10,13,23) = K(10,13)
	coef_lin(13,10,23) = E(10,13)	! 3Ii4Io -> 3Ii2Io + 2Io
	coef_lin(10,13,23) = E(10,13)
	coef_quad(13,11,24) = K(11,13)	! 3Ii2Io + 1Ii2Io -> 4Ii4Io
	coef_quad(11,13,24) = K(11,13)
	coef_lin(13,11,24) = E(11,13)	! 4Ii4Io -> 3Ii2Io + 1Ii2Io
	coef_lin(11,13,24) = E(11,13)
	coef_quad(13,12,31) = K(12,13)	! 3Ii2Io + 2Ii2Io -> out_neu
	coef_quad(12,13,31) = K(12,13)
	coef_quad(13,13,31) = 0.5d0*K(13,13)	! 3Ii2Io + 3Ii2Io -> out_neu

	coef_quad(14,2,14) = K(2,14)	! 4Ii2Io + 2Ii -> boundary -> 4Ii2Io
	coef_quad(2,14,14) = K(2,14)
	coef_quad(14,2,1) = K(2,14)	! 4Ii2Io + 2Ii -> boundary -> 1Ii
	coef_quad(2,14,1) = K(2,14)
	coef_quad(14,3,14) = K(3,14)	! 4Ii2Io + 3Ii -> boundary -> 4Ii2Io
	coef_quad(3,14,14) = K(3,14)
	coef_quad(14,3,1) = K(3,14)	! 4Ii2Io + 3Ii -> boundary -> 1Ii
	coef_quad(3,14,1) = K(3,14)
	coef_quad(14,4,14) = K(4,14)	! 4Ii2Io + 4Ii -> boundary -> 4Ii2Io
	coef_quad(4,14,14) = K(4,14)
	coef_quad(14,4,1) = K(4,14)	! 4Ii2Io + 4Ii -> boundary -> 1Ii
	coef_quad(4,14,1) = K(4,14)
	coef_quad(14,5,19) = K(5,14)	! 4Ii2Io + 1Io -> 4Ii3Io
	coef_quad(5,14,19) = K(5,14)
	coef_lin(14,5,19) = E(5,14)	! 4Ii3Io -> 4Ii2Io + 1Io
	coef_lin(5,14,19) = E(5,14)
	coef_quad(14,6,19) = K(6,14)	! 4Ii2Io + 1Ii1Io -> boundary -> 4Ii3Io
	coef_quad(6,14,19) = K(6,14)
	coef_quad(14,6,1) = K(6,14)	! 4Ii2Io + 1Ii1Io -> boundary -> 1Ii
	coef_quad(6,14,1) = K(6,14)
	coef_quad(14,7,19) = K(7,14)	! 4Ii2Io + 2Ii1Io -> boundary -> 4Ii3Io
	coef_quad(7,14,19) = K(7,14)
	coef_quad(14,7,1) = K(7,14)	! 4Ii2Io + 2Ii1Io -> boundary -> 1Ii
	coef_quad(7,14,1) = K(7,14)
	coef_quad(14,8,19) = K(8,14)	! 4Ii2Io + 3Ii1Io -> boundary -> 4Ii3Io
	coef_quad(8,14,19) = K(8,14)
	coef_quad(14,8,1) = K(8,14)	! 4Ii2Io + 3Ii1Io -> boundary -> 1Ii
	coef_quad(8,14,1) = K(8,14)
	coef_quad(14,9,19) = K(9,14)	! 4Ii2Io + 4Ii1Io -> boundary -> 4Ii3Io
	coef_quad(9,14,19) = K(9,14)
	coef_quad(14,9,1) = K(9,14)	! 4Ii2Io + 4Ii1Io -> boundary -> 1Ii
	coef_quad(9,14,1) = K(9,14)
	coef_quad(14,10,24) = K(10,14)	! 4Ii2Io + 2Io -> 4Ii4Io
	coef_quad(10,14,24) = K(10,14)
	coef_lin(14,10,24) = E(10,14)	! 4Ii4Io -> 4Ii2Io + 2Io
	coef_lin(10,14,24) = E(10,14)
	coef_quad(14,11,31) = K(11,14)	! 4Ii2Io + 1Ii2Io -> out_neu
	coef_quad(11,14,31) = K(11,14)
	coef_quad(14,12,31) = K(12,14)	! 4Ii2Io + 2Ii2Io -> out_neu
	coef_quad(12,14,31) = K(12,14)
	coef_quad(14,13,31) = K(13,14)	! 4Ii2Io + 3Ii2Io -> out_neu
	coef_quad(13,14,31) = K(13,14)
	coef_quad(14,14,31) = 0.5d0*K(14,14)	! 4Ii2Io + 4Ii2Io -> out_neu

	coef_quad(15,1,16) = K(1,15)	! 3Io + 1Ii -> 1Ii3Io
	coef_quad(1,15,16) = K(1,15)
	coef_lin(15,1,16) = E(1,15)	! 1Ii3Io -> 3Io + 1Ii
	coef_lin(1,15,16) = E(1,15)
	coef_quad(15,2,17) = K(2,15)	! 3Io + 2Ii -> 2Ii3Io
	coef_quad(2,15,17) = K(2,15)
	coef_lin(15,2,17) = E(2,15)	! 2Ii3Io -> 3Io + 2Ii
	coef_lin(2,15,17) = E(2,15)
	coef_quad(15,3,18) = K(3,15)	! 3Io + 3Ii -> 3Ii3Io
	coef_quad(3,15,18) = K(3,15)
	coef_lin(15,3,18) = E(3,15)	! 3Ii3Io -> 3Io + 3Ii
	coef_lin(3,15,18) = E(3,15)
	coef_quad(15,4,19) = K(4,15)	! 3Io + 4Ii -> 4Ii3Io
	coef_quad(4,15,19) = K(4,15)
	coef_lin(15,4,19) = E(4,15)	! 4Ii3Io -> 3Io + 4Ii
	coef_lin(4,15,19) = E(4,15)
	coef_quad(15,5,20) = K(5,15)	! 3Io + 1Io -> 4Io
	coef_quad(5,15,20) = K(5,15)
	coef_lin(15,5,20) = E(5,15)	! 4Io -> 3Io + 1Io
	coef_lin(5,15,20) = E(5,15)
	coef_quad(15,6,21) = K(6,15)	! 3Io + 1Ii1Io -> 1Ii4Io
	coef_quad(6,15,21) = K(6,15)
	coef_lin(15,6,21) = E(6,15)	! 1Ii4Io -> 3Io + 1Ii1Io
	coef_lin(6,15,21) = E(6,15)
	coef_quad(15,7,22) = K(7,15)	! 3Io + 2Ii1Io -> 2Ii4Io
	coef_quad(7,15,22) = K(7,15)
	coef_lin(15,7,22) = E(7,15)	! 2Ii4Io -> 3Io + 2Ii1Io
	coef_lin(7,15,22) = E(7,15)
	coef_quad(15,8,23) = K(8,15)	! 3Io + 3Ii1Io -> 3Ii4Io
	coef_quad(8,15,23) = K(8,15)
	coef_lin(15,8,23) = E(8,15)	! 3Ii4Io -> 3Io + 3Ii1Io
	coef_lin(8,15,23) = E(8,15)
	coef_quad(15,9,24) = K(9,15)	! 3Io + 4Ii1Io -> 4Ii4Io
	coef_quad(9,15,24) = K(9,15)
	coef_lin(15,9,24) = E(9,15)	! 4Ii4Io -> 3Io + 4Ii1Io
	coef_lin(9,15,24) = E(9,15)
	coef_quad(15,10,20) = K(10,15)	! 3Io + 2Io -> boundary -> 4Io
	coef_quad(10,15,20) = K(10,15)
	coef_quad(15,10,5) = K(10,15)	! 3Io + 2Io -> boundary -> 1Io
	coef_quad(10,15,5) = K(10,15)
	coef_quad(15,11,21) = K(11,15)	! 3Io + 1Ii2Io -> boundary -> 1Ii4Io
	coef_quad(11,15,21) = K(11,15)
	coef_quad(15,11,5) = K(11,15)	! 3Io + 1Ii2Io -> boundary -> 1Io
	coef_quad(11,15,5) = K(11,15)
	coef_quad(15,12,22) = K(12,15)	! 3Io + 2Ii2Io -> boundary -> 2Ii4Io
	coef_quad(12,15,22) = K(12,15)
	coef_quad(15,12,5) = K(12,15)	! 3Io + 2Ii2Io -> boundary -> 1Io
	coef_quad(12,15,5) = K(12,15)
	coef_quad(15,13,23) = K(13,15)	! 3Io + 3Ii2Io -> boundary -> 3Ii4Io
	coef_quad(13,15,23) = K(13,15)
	coef_quad(15,13,5) = K(13,15)	! 3Io + 3Ii2Io -> boundary -> 1Io
	coef_quad(13,15,5) = K(13,15)
	coef_quad(15,14,31) = K(14,15)	! 3Io + 4Ii2Io -> out_neu
	coef_quad(14,15,31) = K(14,15)
	coef_quad(15,15,20) = 0.5d0*K(15,15)	! 3Io + 3Io -> boundary -> 4Io
	coef_quad(15,15,5) = 0.5d0*K(15,15)	! 3Io + 3Io -> boundary -> 1Io

	coef_quad(16,1,17) = K(1,16)	! 1Ii3Io + 1Ii -> 2Ii3Io
	coef_quad(1,16,17) = K(1,16)
	coef_lin(16,1,17) = E(1,16)	! 2Ii3Io -> 1Ii3Io + 1Ii
	coef_lin(1,16,17) = E(1,16)
	coef_quad(16,2,18) = K(2,16)	! 1Ii3Io + 2Ii -> 3Ii3Io
	coef_quad(2,16,18) = K(2,16)
	coef_lin(16,2,18) = E(2,16)	! 3Ii3Io -> 1Ii3Io + 2Ii
	coef_lin(2,16,18) = E(2,16)
	coef_quad(16,3,19) = K(3,16)	! 1Ii3Io + 3Ii -> 4Ii3Io
	coef_quad(3,16,19) = K(3,16)
	coef_lin(16,3,19) = E(3,16)	! 4Ii3Io -> 1Ii3Io + 3Ii
	coef_lin(3,16,19) = E(3,16)
	coef_quad(16,4,19) = K(4,16)	! 1Ii3Io + 4Ii -> boundary -> 4Ii3Io
	coef_quad(4,16,19) = K(4,16)
	coef_quad(16,4,1) = K(4,16)	! 1Ii3Io + 4Ii -> boundary -> 1Ii
	coef_quad(4,16,1) = K(4,16)
	coef_quad(16,5,21) = K(5,16)	! 1Ii3Io + 1Io -> 1Ii4Io
	coef_quad(5,16,21) = K(5,16)
	coef_lin(16,5,21) = E(5,16)	! 1Ii4Io -> 1Ii3Io + 1Io
	coef_lin(5,16,21) = E(5,16)
	coef_quad(16,6,22) = K(6,16)	! 1Ii3Io + 1Ii1Io -> 2Ii4Io
	coef_quad(6,16,22) = K(6,16)
	coef_lin(16,6,22) = E(6,16)	! 2Ii4Io -> 1Ii3Io + 1Ii1Io
	coef_lin(6,16,22) = E(6,16)
	coef_quad(16,7,23) = K(7,16)	! 1Ii3Io + 2Ii1Io -> 3Ii4Io
	coef_quad(7,16,23) = K(7,16)
	coef_lin(16,7,23) = E(7,16)	! 3Ii4Io -> 1Ii3Io + 2Ii1Io
	coef_lin(7,16,23) = E(7,16)
	coef_quad(16,8,24) = K(8,16)	! 1Ii3Io + 3Ii1Io -> 4Ii4Io
	coef_quad(8,16,24) = K(8,16)
	coef_lin(16,8,24) = E(8,16)	! 4Ii4Io -> 1Ii3Io + 3Ii1Io
	coef_lin(8,16,24) = E(8,16)
	coef_quad(16,9,31) = K(9,16)	! 1Ii3Io + 4Ii1Io -> out_neu
	coef_quad(9,16,31) = K(9,16)
	coef_quad(16,10,21) = K(10,16)	! 1Ii3Io + 2Io -> boundary -> 1Ii4Io
	coef_quad(10,16,21) = K(10,16)
	coef_quad(16,10,5) = K(10,16)	! 1Ii3Io + 2Io -> boundary -> 1Io
	coef_quad(10,16,5) = K(10,16)
	coef_quad(16,11,22) = K(11,16)	! 1Ii3Io + 1Ii2Io -> boundary -> 2Ii4Io
	coef_quad(11,16,22) = K(11,16)
	coef_quad(16,11,5) = K(11,16)	! 1Ii3Io + 1Ii2Io -> boundary -> 1Io
	coef_quad(11,16,5) = K(11,16)
	coef_quad(16,12,23) = K(12,16)	! 1Ii3Io + 2Ii2Io -> boundary -> 3Ii4Io
	coef_quad(12,16,23) = K(12,16)
	coef_quad(16,12,5) = K(12,16)	! 1Ii3Io + 2Ii2Io -> boundary -> 1Io
	coef_quad(12,16,5) = K(12,16)
	coef_quad(16,13,31) = K(13,16)	! 1Ii3Io + 3Ii2Io -> out_neu
	coef_quad(13,16,31) = K(13,16)
	coef_quad(16,14,31) = K(14,16)	! 1Ii3Io + 4Ii2Io -> out_neu
	coef_quad(14,16,31) = K(14,16)
	coef_quad(16,15,21) = K(15,16)	! 1Ii3Io + 3Io -> boundary -> 1Ii4Io
	coef_quad(15,16,21) = K(15,16)
	coef_quad(16,15,5) = K(15,16)	! 1Ii3Io + 3Io -> boundary -> 1Io
	coef_quad(15,16,5) = K(15,16)
	coef_quad(16,16,22) = 0.5d0*K(16,16)	! 1Ii3Io + 1Ii3Io -> boundary -> 2Ii4Io
	coef_quad(16,16,5) = 0.5d0*K(16,16)	! 1Ii3Io + 1Ii3Io -> boundary -> 1Io

	coef_quad(17,1,18) = K(1,17)	! 2Ii3Io + 1Ii -> 3Ii3Io
	coef_quad(1,17,18) = K(1,17)
	coef_lin(17,1,18) = E(1,17)	! 3Ii3Io -> 2Ii3Io + 1Ii
	coef_lin(1,17,18) = E(1,17)
	coef_quad(17,2,19) = K(2,17)	! 2Ii3Io + 2Ii -> 4Ii3Io
	coef_quad(2,17,19) = K(2,17)
	coef_lin(17,2,19) = E(2,17)	! 4Ii3Io -> 2Ii3Io + 2Ii
	coef_lin(2,17,19) = E(2,17)
	coef_quad(17,3,19) = K(3,17)	! 2Ii3Io + 3Ii -> boundary -> 4Ii3Io
	coef_quad(3,17,19) = K(3,17)
	coef_quad(17,3,1) = K(3,17)	! 2Ii3Io + 3Ii -> boundary -> 1Ii
	coef_quad(3,17,1) = K(3,17)
	coef_quad(17,4,19) = K(4,17)	! 2Ii3Io + 4Ii -> boundary -> 4Ii3Io
	coef_quad(4,17,19) = K(4,17)
	coef_quad(17,4,1) = K(4,17)	! 2Ii3Io + 4Ii -> boundary -> 1Ii
	coef_quad(4,17,1) = K(4,17)
	coef_quad(17,5,22) = K(5,17)	! 2Ii3Io + 1Io -> 2Ii4Io
	coef_quad(5,17,22) = K(5,17)
	coef_lin(17,5,22) = E(5,17)	! 2Ii4Io -> 2Ii3Io + 1Io
	coef_lin(5,17,22) = E(5,17)
	coef_quad(17,6,23) = K(6,17)	! 2Ii3Io + 1Ii1Io -> 3Ii4Io
	coef_quad(6,17,23) = K(6,17)
	coef_lin(17,6,23) = E(6,17)	! 3Ii4Io -> 2Ii3Io + 1Ii1Io
	coef_lin(6,17,23) = E(6,17)
	coef_quad(17,7,24) = K(7,17)	! 2Ii3Io + 2Ii1Io -> 4Ii4Io
	coef_quad(7,17,24) = K(7,17)
	coef_lin(17,7,24) = E(7,17)	! 4Ii4Io -> 2Ii3Io + 2Ii1Io
	coef_lin(7,17,24) = E(7,17)
	coef_quad(17,8,31) = K(8,17)	! 2Ii3Io + 3Ii1Io -> out_neu
	coef_quad(8,17,31) = K(8,17)
	coef_quad(17,9,31) = K(9,17)	! 2Ii3Io + 4Ii1Io -> out_neu
	coef_quad(9,17,31) = K(9,17)
	coef_quad(17,10,22) = K(10,17)	! 2Ii3Io + 2Io -> boundary -> 2Ii4Io
	coef_quad(10,17,22) = K(10,17)
	coef_quad(17,10,5) = K(10,17)	! 2Ii3Io + 2Io -> boundary -> 1Io
	coef_quad(10,17,5) = K(10,17)
	coef_quad(17,11,23) = K(11,17)	! 2Ii3Io + 1Ii2Io -> boundary -> 3Ii4Io
	coef_quad(11,17,23) = K(11,17)
	coef_quad(17,11,5) = K(11,17)	! 2Ii3Io + 1Ii2Io -> boundary -> 1Io
	coef_quad(11,17,5) = K(11,17)
	coef_quad(17,12,31) = K(12,17)	! 2Ii3Io + 2Ii2Io -> out_neu
	coef_quad(12,17,31) = K(12,17)
	coef_quad(17,13,31) = K(13,17)	! 2Ii3Io + 3Ii2Io -> out_neu
	coef_quad(13,17,31) = K(13,17)
	coef_quad(17,14,31) = K(14,17)	! 2Ii3Io + 4Ii2Io -> out_neu
	coef_quad(14,17,31) = K(14,17)
	coef_quad(17,15,22) = K(15,17)	! 2Ii3Io + 3Io -> boundary -> 2Ii4Io
	coef_quad(15,17,22) = K(15,17)
	coef_quad(17,15,5) = K(15,17)	! 2Ii3Io + 3Io -> boundary -> 1Io
	coef_quad(15,17,5) = K(15,17)
	coef_quad(17,16,23) = K(16,17)	! 2Ii3Io + 1Ii3Io -> boundary -> 3Ii4Io
	coef_quad(16,17,23) = K(16,17)
	coef_quad(17,16,5) = K(16,17)	! 2Ii3Io + 1Ii3Io -> boundary -> 1Io
	coef_quad(16,17,5) = K(16,17)
	coef_quad(17,17,31) = 0.5d0*K(17,17)	! 2Ii3Io + 2Ii3Io -> out_neu

	coef_quad(18,1,19) = K(1,18)	! 3Ii3Io + 1Ii -> 4Ii3Io
	coef_quad(1,18,19) = K(1,18)
	coef_lin(18,1,19) = E(1,18)	! 4Ii3Io -> 3Ii3Io + 1Ii
	coef_lin(1,18,19) = E(1,18)
	coef_quad(18,2,19) = K(2,18)	! 3Ii3Io + 2Ii -> boundary -> 4Ii3Io
	coef_quad(2,18,19) = K(2,18)
	coef_quad(18,2,1) = K(2,18)	! 3Ii3Io + 2Ii -> boundary -> 1Ii
	coef_quad(2,18,1) = K(2,18)
	coef_quad(18,3,19) = K(3,18)	! 3Ii3Io + 3Ii -> boundary -> 4Ii3Io
	coef_quad(3,18,19) = K(3,18)
	coef_quad(18,3,1) = K(3,18)	! 3Ii3Io + 3Ii -> boundary -> 1Ii
	coef_quad(3,18,1) = K(3,18)
	coef_quad(18,4,19) = K(4,18)	! 3Ii3Io + 4Ii -> boundary -> 4Ii3Io
	coef_quad(4,18,19) = K(4,18)
	coef_quad(18,4,1) = K(4,18)	! 3Ii3Io + 4Ii -> boundary -> 1Ii
	coef_quad(4,18,1) = K(4,18)
	coef_quad(18,5,23) = K(5,18)	! 3Ii3Io + 1Io -> 3Ii4Io
	coef_quad(5,18,23) = K(5,18)
	coef_lin(18,5,23) = E(5,18)	! 3Ii4Io -> 3Ii3Io + 1Io
	coef_lin(5,18,23) = E(5,18)
	coef_quad(18,6,24) = K(6,18)	! 3Ii3Io + 1Ii1Io -> 4Ii4Io
	coef_quad(6,18,24) = K(6,18)
	coef_lin(18,6,24) = E(6,18)	! 4Ii4Io -> 3Ii3Io + 1Ii1Io
	coef_lin(6,18,24) = E(6,18)
	coef_quad(18,7,31) = K(7,18)	! 3Ii3Io + 2Ii1Io -> out_neu
	coef_quad(7,18,31) = K(7,18)
	coef_quad(18,8,31) = K(8,18)	! 3Ii3Io + 3Ii1Io -> out_neu
	coef_quad(8,18,31) = K(8,18)
	coef_quad(18,9,31) = K(9,18)	! 3Ii3Io + 4Ii1Io -> out_neu
	coef_quad(9,18,31) = K(9,18)
	coef_quad(18,10,23) = K(10,18)	! 3Ii3Io + 2Io -> boundary -> 3Ii4Io
	coef_quad(10,18,23) = K(10,18)
	coef_quad(18,10,5) = K(10,18)	! 3Ii3Io + 2Io -> boundary -> 1Io
	coef_quad(10,18,5) = K(10,18)
	coef_quad(18,11,31) = K(11,18)	! 3Ii3Io + 1Ii2Io -> out_neu
	coef_quad(11,18,31) = K(11,18)
	coef_quad(18,12,31) = K(12,18)	! 3Ii3Io + 2Ii2Io -> out_neu
	coef_quad(12,18,31) = K(12,18)
	coef_quad(18,13,31) = K(13,18)	! 3Ii3Io + 3Ii2Io -> out_neu
	coef_quad(13,18,31) = K(13,18)
	coef_quad(18,14,31) = K(14,18)	! 3Ii3Io + 4Ii2Io -> out_neu
	coef_quad(14,18,31) = K(14,18)
	coef_quad(18,15,23) = K(15,18)	! 3Ii3Io + 3Io -> boundary -> 3Ii4Io
	coef_quad(15,18,23) = K(15,18)
	coef_quad(18,15,5) = K(15,18)	! 3Ii3Io + 3Io -> boundary -> 1Io
	coef_quad(15,18,5) = K(15,18)
	coef_quad(18,16,31) = K(16,18)	! 3Ii3Io + 1Ii3Io -> out_neu
	coef_quad(16,18,31) = K(16,18)
	coef_quad(18,17,31) = K(17,18)	! 3Ii3Io + 2Ii3Io -> out_neu
	coef_quad(17,18,31) = K(17,18)
	coef_quad(18,18,31) = 0.5d0*K(18,18)	! 3Ii3Io + 3Ii3Io -> out_neu

	coef_quad(19,2,19) = K(2,19)	! 4Ii3Io + 2Ii -> boundary -> 4Ii3Io
	coef_quad(2,19,19) = K(2,19)
	coef_quad(19,2,1) = K(2,19)	! 4Ii3Io + 2Ii -> boundary -> 1Ii
	coef_quad(2,19,1) = K(2,19)
	coef_quad(19,3,19) = K(3,19)	! 4Ii3Io + 3Ii -> boundary -> 4Ii3Io
	coef_quad(3,19,19) = K(3,19)
	coef_quad(19,3,1) = K(3,19)	! 4Ii3Io + 3Ii -> boundary -> 1Ii
	coef_quad(3,19,1) = K(3,19)
	coef_quad(19,4,19) = K(4,19)	! 4Ii3Io + 4Ii -> boundary -> 4Ii3Io
	coef_quad(4,19,19) = K(4,19)
	coef_quad(19,4,1) = K(4,19)	! 4Ii3Io + 4Ii -> boundary -> 1Ii
	coef_quad(4,19,1) = K(4,19)
	coef_quad(19,5,24) = K(5,19)	! 4Ii3Io + 1Io -> 4Ii4Io
	coef_quad(5,19,24) = K(5,19)
	coef_lin(19,5,24) = E(5,19)	! 4Ii4Io -> 4Ii3Io + 1Io
	coef_lin(5,19,24) = E(5,19)
	coef_quad(19,6,31) = K(6,19)	! 4Ii3Io + 1Ii1Io -> out_neu
	coef_quad(6,19,31) = K(6,19)
	coef_quad(19,7,31) = K(7,19)	! 4Ii3Io + 2Ii1Io -> out_neu
	coef_quad(7,19,31) = K(7,19)
	coef_quad(19,8,31) = K(8,19)	! 4Ii3Io + 3Ii1Io -> out_neu
	coef_quad(8,19,31) = K(8,19)
	coef_quad(19,9,31) = K(9,19)	! 4Ii3Io + 4Ii1Io -> out_neu
	coef_quad(9,19,31) = K(9,19)
	coef_quad(19,10,31) = K(10,19)	! 4Ii3Io + 2Io -> out_neu
	coef_quad(10,19,31) = K(10,19)
	coef_quad(19,11,31) = K(11,19)	! 4Ii3Io + 1Ii2Io -> out_neu
	coef_quad(11,19,31) = K(11,19)
	coef_quad(19,12,31) = K(12,19)	! 4Ii3Io + 2Ii2Io -> out_neu
	coef_quad(12,19,31) = K(12,19)
	coef_quad(19,13,31) = K(13,19)	! 4Ii3Io + 3Ii2Io -> out_neu
	coef_quad(13,19,31) = K(13,19)
	coef_quad(19,14,31) = K(14,19)	! 4Ii3Io + 4Ii2Io -> out_neu
	coef_quad(14,19,31) = K(14,19)
	coef_quad(19,15,31) = K(15,19)	! 4Ii3Io + 3Io -> out_neu
	coef_quad(15,19,31) = K(15,19)
	coef_quad(19,16,31) = K(16,19)	! 4Ii3Io + 1Ii3Io -> out_neu
	coef_quad(16,19,31) = K(16,19)
	coef_quad(19,17,31) = K(17,19)	! 4Ii3Io + 2Ii3Io -> out_neu
	coef_quad(17,19,31) = K(17,19)
	coef_quad(19,18,31) = K(18,19)	! 4Ii3Io + 3Ii3Io -> out_neu
	coef_quad(18,19,31) = K(18,19)
	coef_quad(19,19,31) = 0.5d0*K(19,19)	! 4Ii3Io + 4Ii3Io -> out_neu

	coef_quad(20,1,21) = K(1,20)	! 4Io + 1Ii -> 1Ii4Io
	coef_quad(1,20,21) = K(1,20)
	coef_lin(20,1,21) = E(1,20)	! 1Ii4Io -> 4Io + 1Ii
	coef_lin(1,20,21) = E(1,20)
	coef_quad(20,2,22) = K(2,20)	! 4Io + 2Ii -> 2Ii4Io
	coef_quad(2,20,22) = K(2,20)
	coef_lin(20,2,22) = E(2,20)	! 2Ii4Io -> 4Io + 2Ii
	coef_lin(2,20,22) = E(2,20)
	coef_quad(20,3,23) = K(3,20)	! 4Io + 3Ii -> 3Ii4Io
	coef_quad(3,20,23) = K(3,20)
	coef_lin(20,3,23) = E(3,20)	! 3Ii4Io -> 4Io + 3Ii
	coef_lin(3,20,23) = E(3,20)
	coef_quad(20,4,24) = K(4,20)	! 4Io + 4Ii -> 4Ii4Io
	coef_quad(4,20,24) = K(4,20)
	coef_lin(20,4,24) = E(4,20)	! 4Ii4Io -> 4Io + 4Ii
	coef_lin(4,20,24) = E(4,20)
	coef_quad(20,6,21) = K(6,20)	! 4Io + 1Ii1Io -> boundary -> 1Ii4Io
	coef_quad(6,20,21) = K(6,20)
	coef_quad(20,6,5) = K(6,20)	! 4Io + 1Ii1Io -> boundary -> 1Io
	coef_quad(6,20,5) = K(6,20)
	coef_quad(20,7,22) = K(7,20)	! 4Io + 2Ii1Io -> boundary -> 2Ii4Io
	coef_quad(7,20,22) = K(7,20)
	coef_quad(20,7,5) = K(7,20)	! 4Io + 2Ii1Io -> boundary -> 1Io
	coef_quad(7,20,5) = K(7,20)
	coef_quad(20,8,23) = K(8,20)	! 4Io + 3Ii1Io -> boundary -> 3Ii4Io
	coef_quad(8,20,23) = K(8,20)
	coef_quad(20,8,5) = K(8,20)	! 4Io + 3Ii1Io -> boundary -> 1Io
	coef_quad(8,20,5) = K(8,20)
	coef_quad(20,9,31) = K(9,20)	! 4Io + 4Ii1Io -> out_neu
	coef_quad(9,20,31) = K(9,20)
	coef_quad(20,10,20) = K(10,20)	! 4Io + 2Io -> boundary -> 4Io
	coef_quad(10,20,20) = K(10,20)
	coef_quad(20,10,5) = K(10,20)	! 4Io + 2Io -> boundary -> 1Io
	coef_quad(10,20,5) = K(10,20)
	coef_quad(20,11,21) = K(11,20)	! 4Io + 1Ii2Io -> boundary -> 1Ii4Io
	coef_quad(11,20,21) = K(11,20)
	coef_quad(20,11,5) = K(11,20)	! 4Io + 1Ii2Io -> boundary -> 1Io
	coef_quad(11,20,5) = K(11,20)
	coef_quad(20,12,22) = K(12,20)	! 4Io + 2Ii2Io -> boundary -> 2Ii4Io
	coef_quad(12,20,22) = K(12,20)
	coef_quad(20,12,5) = K(12,20)	! 4Io + 2Ii2Io -> boundary -> 1Io
	coef_quad(12,20,5) = K(12,20)
	coef_quad(20,13,23) = K(13,20)	! 4Io + 3Ii2Io -> boundary -> 3Ii4Io
	coef_quad(13,20,23) = K(13,20)
	coef_quad(20,13,5) = K(13,20)	! 4Io + 3Ii2Io -> boundary -> 1Io
	coef_quad(13,20,5) = K(13,20)
	coef_quad(20,14,31) = K(14,20)	! 4Io + 4Ii2Io -> out_neu
	coef_quad(14,20,31) = K(14,20)
	coef_quad(20,15,20) = K(15,20)	! 4Io + 3Io -> boundary -> 4Io
	coef_quad(15,20,20) = K(15,20)
	coef_quad(20,15,5) = K(15,20)	! 4Io + 3Io -> boundary -> 1Io
	coef_quad(15,20,5) = K(15,20)
	coef_quad(20,16,21) = K(16,20)	! 4Io + 1Ii3Io -> boundary -> 1Ii4Io
	coef_quad(16,20,21) = K(16,20)
	coef_quad(20,16,5) = K(16,20)	! 4Io + 1Ii3Io -> boundary -> 1Io
	coef_quad(16,20,5) = K(16,20)
	coef_quad(20,17,22) = K(17,20)	! 4Io + 2Ii3Io -> boundary -> 2Ii4Io
	coef_quad(17,20,22) = K(17,20)
	coef_quad(20,17,5) = K(17,20)	! 4Io + 2Ii3Io -> boundary -> 1Io
	coef_quad(17,20,5) = K(17,20)
	coef_quad(20,18,23) = K(18,20)	! 4Io + 3Ii3Io -> boundary -> 3Ii4Io
	coef_quad(18,20,23) = K(18,20)
	coef_quad(20,18,5) = K(18,20)	! 4Io + 3Ii3Io -> boundary -> 1Io
	coef_quad(18,20,5) = K(18,20)
	coef_quad(20,19,31) = K(19,20)	! 4Io + 4Ii3Io -> out_neu
	coef_quad(19,20,31) = K(19,20)
	coef_quad(20,20,20) = 0.5d0*K(20,20)	! 4Io + 4Io -> boundary -> 4Io
	coef_quad(20,20,5) = 0.5d0*K(20,20)	! 4Io + 4Io -> boundary -> 1Io

	coef_quad(21,1,22) = K(1,21)	! 1Ii4Io + 1Ii -> 2Ii4Io
	coef_quad(1,21,22) = K(1,21)
	coef_lin(21,1,22) = E(1,21)	! 2Ii4Io -> 1Ii4Io + 1Ii
	coef_lin(1,21,22) = E(1,21)
	coef_quad(21,2,23) = K(2,21)	! 1Ii4Io + 2Ii -> 3Ii4Io
	coef_quad(2,21,23) = K(2,21)
	coef_lin(21,2,23) = E(2,21)	! 3Ii4Io -> 1Ii4Io + 2Ii
	coef_lin(2,21,23) = E(2,21)
	coef_quad(21,3,24) = K(3,21)	! 1Ii4Io + 3Ii -> 4Ii4Io
	coef_quad(3,21,24) = K(3,21)
	coef_lin(21,3,24) = E(3,21)	! 4Ii4Io -> 1Ii4Io + 3Ii
	coef_lin(3,21,24) = E(3,21)
	coef_quad(21,4,31) = K(4,21)	! 1Ii4Io + 4Ii -> out_neu
	coef_quad(4,21,31) = K(4,21)
	coef_quad(21,6,22) = K(6,21)	! 1Ii4Io + 1Ii1Io -> boundary -> 2Ii4Io
	coef_quad(6,21,22) = K(6,21)
	coef_quad(21,6,5) = K(6,21)	! 1Ii4Io + 1Ii1Io -> boundary -> 1Io
	coef_quad(6,21,5) = K(6,21)
	coef_quad(21,7,23) = K(7,21)	! 1Ii4Io + 2Ii1Io -> boundary -> 3Ii4Io
	coef_quad(7,21,23) = K(7,21)
	coef_quad(21,7,5) = K(7,21)	! 1Ii4Io + 2Ii1Io -> boundary -> 1Io
	coef_quad(7,21,5) = K(7,21)
	coef_quad(21,8,31) = K(8,21)	! 1Ii4Io + 3Ii1Io -> out_neu
	coef_quad(8,21,31) = K(8,21)
	coef_quad(21,9,31) = K(9,21)	! 1Ii4Io + 4Ii1Io -> out_neu
	coef_quad(9,21,31) = K(9,21)
	coef_quad(21,10,21) = K(10,21)	! 1Ii4Io + 2Io -> boundary -> 1Ii4Io
	coef_quad(10,21,21) = K(10,21)
	coef_quad(21,10,5) = K(10,21)	! 1Ii4Io + 2Io -> boundary -> 1Io
	coef_quad(10,21,5) = K(10,21)
	coef_quad(21,11,22) = K(11,21)	! 1Ii4Io + 1Ii2Io -> boundary -> 2Ii4Io
	coef_quad(11,21,22) = K(11,21)
	coef_quad(21,11,5) = K(11,21)	! 1Ii4Io + 1Ii2Io -> boundary -> 1Io
	coef_quad(11,21,5) = K(11,21)
	coef_quad(21,12,23) = K(12,21)	! 1Ii4Io + 2Ii2Io -> boundary -> 3Ii4Io
	coef_quad(12,21,23) = K(12,21)
	coef_quad(21,12,5) = K(12,21)	! 1Ii4Io + 2Ii2Io -> boundary -> 1Io
	coef_quad(12,21,5) = K(12,21)
	coef_quad(21,13,31) = K(13,21)	! 1Ii4Io + 3Ii2Io -> out_neu
	coef_quad(13,21,31) = K(13,21)
	coef_quad(21,14,31) = K(14,21)	! 1Ii4Io + 4Ii2Io -> out_neu
	coef_quad(14,21,31) = K(14,21)
	coef_quad(21,15,21) = K(15,21)	! 1Ii4Io + 3Io -> boundary -> 1Ii4Io
	coef_quad(15,21,21) = K(15,21)
	coef_quad(21,15,5) = K(15,21)	! 1Ii4Io + 3Io -> boundary -> 1Io
	coef_quad(15,21,5) = K(15,21)
	coef_quad(21,16,22) = K(16,21)	! 1Ii4Io + 1Ii3Io -> boundary -> 2Ii4Io
	coef_quad(16,21,22) = K(16,21)
	coef_quad(21,16,5) = K(16,21)	! 1Ii4Io + 1Ii3Io -> boundary -> 1Io
	coef_quad(16,21,5) = K(16,21)
	coef_quad(21,17,23) = K(17,21)	! 1Ii4Io + 2Ii3Io -> boundary -> 3Ii4Io
	coef_quad(17,21,23) = K(17,21)
	coef_quad(21,17,5) = K(17,21)	! 1Ii4Io + 2Ii3Io -> boundary -> 1Io
	coef_quad(17,21,5) = K(17,21)
	coef_quad(21,18,31) = K(18,21)	! 1Ii4Io + 3Ii3Io -> out_neu
	coef_quad(18,21,31) = K(18,21)
	coef_quad(21,19,31) = K(19,21)	! 1Ii4Io + 4Ii3Io -> out_neu
	coef_quad(19,21,31) = K(19,21)
	coef_quad(21,20,21) = K(20,21)	! 1Ii4Io + 4Io -> boundary -> 1Ii4Io
	coef_quad(20,21,21) = K(20,21)
	coef_quad(21,20,5) = K(20,21)	! 1Ii4Io + 4Io -> boundary -> 1Io
	coef_quad(20,21,5) = K(20,21)
	coef_quad(21,21,22) = 0.5d0*K(21,21)	! 1Ii4Io + 1Ii4Io -> boundary -> 2Ii4Io
	coef_quad(21,21,5) = 0.5d0*K(21,21)	! 1Ii4Io + 1Ii4Io -> boundary -> 1Io

	coef_quad(22,1,23) = K(1,22)	! 2Ii4Io + 1Ii -> 3Ii4Io
	coef_quad(1,22,23) = K(1,22)
	coef_lin(22,1,23) = E(1,22)	! 3Ii4Io -> 2Ii4Io + 1Ii
	coef_lin(1,22,23) = E(1,22)
	coef_quad(22,2,24) = K(2,22)	! 2Ii4Io + 2Ii -> 4Ii4Io
	coef_quad(2,22,24) = K(2,22)
	coef_lin(22,2,24) = E(2,22)	! 4Ii4Io -> 2Ii4Io + 2Ii
	coef_lin(2,22,24) = E(2,22)
	coef_quad(22,3,31) = K(3,22)	! 2Ii4Io + 3Ii -> out_neu
	coef_quad(3,22,31) = K(3,22)
	coef_quad(22,4,31) = K(4,22)	! 2Ii4Io + 4Ii -> out_neu
	coef_quad(4,22,31) = K(4,22)
	coef_quad(22,6,23) = K(6,22)	! 2Ii4Io + 1Ii1Io -> boundary -> 3Ii4Io
	coef_quad(6,22,23) = K(6,22)
	coef_quad(22,6,5) = K(6,22)	! 2Ii4Io + 1Ii1Io -> boundary -> 1Io
	coef_quad(6,22,5) = K(6,22)
	coef_quad(22,7,31) = K(7,22)	! 2Ii4Io + 2Ii1Io -> out_neu
	coef_quad(7,22,31) = K(7,22)
	coef_quad(22,8,31) = K(8,22)	! 2Ii4Io + 3Ii1Io -> out_neu
	coef_quad(8,22,31) = K(8,22)
	coef_quad(22,9,31) = K(9,22)	! 2Ii4Io + 4Ii1Io -> out_neu
	coef_quad(9,22,31) = K(9,22)
	coef_quad(22,10,22) = K(10,22)	! 2Ii4Io + 2Io -> boundary -> 2Ii4Io
	coef_quad(10,22,22) = K(10,22)
	coef_quad(22,10,5) = K(10,22)	! 2Ii4Io + 2Io -> boundary -> 1Io
	coef_quad(10,22,5) = K(10,22)
	coef_quad(22,11,23) = K(11,22)	! 2Ii4Io + 1Ii2Io -> boundary -> 3Ii4Io
	coef_quad(11,22,23) = K(11,22)
	coef_quad(22,11,5) = K(11,22)	! 2Ii4Io + 1Ii2Io -> boundary -> 1Io
	coef_quad(11,22,5) = K(11,22)
	coef_quad(22,12,31) = K(12,22)	! 2Ii4Io + 2Ii2Io -> out_neu
	coef_quad(12,22,31) = K(12,22)
	coef_quad(22,13,31) = K(13,22)	! 2Ii4Io + 3Ii2Io -> out_neu
	coef_quad(13,22,31) = K(13,22)
	coef_quad(22,14,31) = K(14,22)	! 2Ii4Io + 4Ii2Io -> out_neu
	coef_quad(14,22,31) = K(14,22)
	coef_quad(22,15,22) = K(15,22)	! 2Ii4Io + 3Io -> boundary -> 2Ii4Io
	coef_quad(15,22,22) = K(15,22)
	coef_quad(22,15,5) = K(15,22)	! 2Ii4Io + 3Io -> boundary -> 1Io
	coef_quad(15,22,5) = K(15,22)
	coef_quad(22,16,23) = K(16,22)	! 2Ii4Io + 1Ii3Io -> boundary -> 3Ii4Io
	coef_quad(16,22,23) = K(16,22)
	coef_quad(22,16,5) = K(16,22)	! 2Ii4Io + 1Ii3Io -> boundary -> 1Io
	coef_quad(16,22,5) = K(16,22)
	coef_quad(22,17,31) = K(17,22)	! 2Ii4Io + 2Ii3Io -> out_neu
	coef_quad(17,22,31) = K(17,22)
	coef_quad(22,18,31) = K(18,22)	! 2Ii4Io + 3Ii3Io -> out_neu
	coef_quad(18,22,31) = K(18,22)
	coef_quad(22,19,31) = K(19,22)	! 2Ii4Io + 4Ii3Io -> out_neu
	coef_quad(19,22,31) = K(19,22)
	coef_quad(22,20,22) = K(20,22)	! 2Ii4Io + 4Io -> boundary -> 2Ii4Io
	coef_quad(20,22,22) = K(20,22)
	coef_quad(22,20,5) = K(20,22)	! 2Ii4Io + 4Io -> boundary -> 1Io
	coef_quad(20,22,5) = K(20,22)
	coef_quad(22,21,23) = K(21,22)	! 2Ii4Io + 1Ii4Io -> boundary -> 3Ii4Io
	coef_quad(21,22,23) = K(21,22)
	coef_quad(22,21,5) = K(21,22)	! 2Ii4Io + 1Ii4Io -> boundary -> 1Io
	coef_quad(21,22,5) = K(21,22)
	coef_quad(22,22,31) = 0.5d0*K(22,22)	! 2Ii4Io + 2Ii4Io -> out_neu

	coef_quad(23,1,24) = K(1,23)	! 3Ii4Io + 1Ii -> 4Ii4Io
	coef_quad(1,23,24) = K(1,23)
	coef_lin(23,1,24) = E(1,23)	! 4Ii4Io -> 3Ii4Io + 1Ii
	coef_lin(1,23,24) = E(1,23)
	coef_quad(23,2,31) = K(2,23)	! 3Ii4Io + 2Ii -> out_neu
	coef_quad(2,23,31) = K(2,23)
	coef_quad(23,3,31) = K(3,23)	! 3Ii4Io + 3Ii -> out_neu
	coef_quad(3,23,31) = K(3,23)
	coef_quad(23,4,31) = K(4,23)	! 3Ii4Io + 4Ii -> out_neu
	coef_quad(4,23,31) = K(4,23)
	coef_quad(23,6,31) = K(6,23)	! 3Ii4Io + 1Ii1Io -> out_neu
	coef_quad(6,23,31) = K(6,23)
	coef_quad(23,7,31) = K(7,23)	! 3Ii4Io + 2Ii1Io -> out_neu
	coef_quad(7,23,31) = K(7,23)
	coef_quad(23,8,31) = K(8,23)	! 3Ii4Io + 3Ii1Io -> out_neu
	coef_quad(8,23,31) = K(8,23)
	coef_quad(23,9,31) = K(9,23)	! 3Ii4Io + 4Ii1Io -> out_neu
	coef_quad(9,23,31) = K(9,23)
	coef_quad(23,10,23) = K(10,23)	! 3Ii4Io + 2Io -> boundary -> 3Ii4Io
	coef_quad(10,23,23) = K(10,23)
	coef_quad(23,10,5) = K(10,23)	! 3Ii4Io + 2Io -> boundary -> 1Io
	coef_quad(10,23,5) = K(10,23)
	coef_quad(23,11,31) = K(11,23)	! 3Ii4Io + 1Ii2Io -> out_neu
	coef_quad(11,23,31) = K(11,23)
	coef_quad(23,12,31) = K(12,23)	! 3Ii4Io + 2Ii2Io -> out_neu
	coef_quad(12,23,31) = K(12,23)
	coef_quad(23,13,31) = K(13,23)	! 3Ii4Io + 3Ii2Io -> out_neu
	coef_quad(13,23,31) = K(13,23)
	coef_quad(23,14,31) = K(14,23)	! 3Ii4Io + 4Ii2Io -> out_neu
	coef_quad(14,23,31) = K(14,23)
	coef_quad(23,15,23) = K(15,23)	! 3Ii4Io + 3Io -> boundary -> 3Ii4Io
	coef_quad(15,23,23) = K(15,23)
	coef_quad(23,15,5) = K(15,23)	! 3Ii4Io + 3Io -> boundary -> 1Io
	coef_quad(15,23,5) = K(15,23)
	coef_quad(23,16,31) = K(16,23)	! 3Ii4Io + 1Ii3Io -> out_neu
	coef_quad(16,23,31) = K(16,23)
	coef_quad(23,17,31) = K(17,23)	! 3Ii4Io + 2Ii3Io -> out_neu
	coef_quad(17,23,31) = K(17,23)
	coef_quad(23,18,31) = K(18,23)	! 3Ii4Io + 3Ii3Io -> out_neu
	coef_quad(18,23,31) = K(18,23)
	coef_quad(23,19,31) = K(19,23)	! 3Ii4Io + 4Ii3Io -> out_neu
	coef_quad(19,23,31) = K(19,23)
	coef_quad(23,20,23) = K(20,23)	! 3Ii4Io + 4Io -> boundary -> 3Ii4Io
	coef_quad(20,23,23) = K(20,23)
	coef_quad(23,20,5) = K(20,23)	! 3Ii4Io + 4Io -> boundary -> 1Io
	coef_quad(20,23,5) = K(20,23)
	coef_quad(23,21,31) = K(21,23)	! 3Ii4Io + 1Ii4Io -> out_neu
	coef_quad(21,23,31) = K(21,23)
	coef_quad(23,22,31) = K(22,23)	! 3Ii4Io + 2Ii4Io -> out_neu
	coef_quad(22,23,31) = K(22,23)
	coef_quad(23,23,31) = 0.5d0*K(23,23)	! 3Ii4Io + 3Ii4Io -> out_neu

	coef_quad(24,1,31) = K(1,24)	! 4Ii4Io + 1Ii -> out_neu
	coef_quad(1,24,31) = K(1,24)
	coef_quad(24,2,31) = K(2,24)	! 4Ii4Io + 2Ii -> out_neu
	coef_quad(2,24,31) = K(2,24)
	coef_quad(24,3,31) = K(3,24)	! 4Ii4Io + 3Ii -> out_neu
	coef_quad(3,24,31) = K(3,24)
	coef_quad(24,4,31) = K(4,24)	! 4Ii4Io + 4Ii -> out_neu
	coef_quad(4,24,31) = K(4,24)
	coef_quad(24,5,31) = K(5,24)	! 4Ii4Io + 1Io -> out_neu
	coef_quad(5,24,31) = K(5,24)
	coef_quad(24,6,31) = K(6,24)	! 4Ii4Io + 1Ii1Io -> out_neu
	coef_quad(6,24,31) = K(6,24)
	coef_quad(24,7,31) = K(7,24)	! 4Ii4Io + 2Ii1Io -> out_neu
	coef_quad(7,24,31) = K(7,24)
	coef_quad(24,8,31) = K(8,24)	! 4Ii4Io + 3Ii1Io -> out_neu
	coef_quad(8,24,31) = K(8,24)
	coef_quad(24,9,31) = K(9,24)	! 4Ii4Io + 4Ii1Io -> out_neu
	coef_quad(9,24,31) = K(9,24)
	coef_quad(24,10,31) = K(10,24)	! 4Ii4Io + 2Io -> out_neu
	coef_quad(10,24,31) = K(10,24)
	coef_quad(24,11,31) = K(11,24)	! 4Ii4Io + 1Ii2Io -> out_neu
	coef_quad(11,24,31) = K(11,24)
	coef_quad(24,12,31) = K(12,24)	! 4Ii4Io + 2Ii2Io -> out_neu
	coef_quad(12,24,31) = K(12,24)
	coef_quad(24,13,31) = K(13,24)	! 4Ii4Io + 3Ii2Io -> out_neu
	coef_quad(13,24,31) = K(13,24)
	coef_quad(24,14,31) = K(14,24)	! 4Ii4Io + 4Ii2Io -> out_neu
	coef_quad(14,24,31) = K(14,24)
	coef_quad(24,15,31) = K(15,24)	! 4Ii4Io + 3Io -> out_neu
	coef_quad(15,24,31) = K(15,24)
	coef_quad(24,16,31) = K(16,24)	! 4Ii4Io + 1Ii3Io -> out_neu
	coef_quad(16,24,31) = K(16,24)
	coef_quad(24,17,31) = K(17,24)	! 4Ii4Io + 2Ii3Io -> out_neu
	coef_quad(17,24,31) = K(17,24)
	coef_quad(24,18,31) = K(18,24)	! 4Ii4Io + 3Ii3Io -> out_neu
	coef_quad(18,24,31) = K(18,24)
	coef_quad(24,19,31) = K(19,24)	! 4Ii4Io + 4Ii3Io -> out_neu
	coef_quad(19,24,31) = K(19,24)
	coef_quad(24,20,31) = K(20,24)	! 4Ii4Io + 4Io -> out_neu
	coef_quad(20,24,31) = K(20,24)
	coef_quad(24,21,31) = K(21,24)	! 4Ii4Io + 1Ii4Io -> out_neu
	coef_quad(21,24,31) = K(21,24)
	coef_quad(24,22,31) = K(22,24)	! 4Ii4Io + 2Ii4Io -> out_neu
	coef_quad(22,24,31) = K(22,24)
	coef_quad(24,23,31) = K(23,24)	! 4Ii4Io + 3Ii4Io -> out_neu
	coef_quad(23,24,31) = K(23,24)
	coef_quad(24,24,31) = 0.5d0*K(24,24)	! 4Ii4Io + 4Ii4Io -> out_neu

	coef_lin(26,26,24) = cs(24)	! 4Ii4Io -> coag
	coef_lin(26,26,23) = cs(23)	! 3Ii4Io -> coag
	coef_lin(26,26,22) = cs(22)	! 2Ii4Io -> coag
	coef_lin(26,26,21) = cs(21)	! 1Ii4Io -> coag
	coef_lin(26,26,20) = cs(20)	! 4Io -> coag
	coef_lin(26,26,19) = cs(19)	! 4Ii3Io -> coag
	coef_lin(26,26,18) = cs(18)	! 3Ii3Io -> coag
	coef_lin(26,26,17) = cs(17)	! 2Ii3Io -> coag
	coef_lin(26,26,16) = cs(16)	! 1Ii3Io -> coag
	coef_lin(26,26,15) = cs(15)	! 3Io -> coag
	coef_lin(26,26,14) = cs(14)	! 4Ii2Io -> coag
	coef_lin(26,26,13) = cs(13)	! 3Ii2Io -> coag
	coef_lin(26,26,12) = cs(12)	! 2Ii2Io -> coag
	coef_lin(26,26,11) = cs(11)	! 1Ii2Io -> coag
	coef_lin(26,26,10) = cs(10)	! 2Io -> coag
	coef_lin(26,26,9) = cs(9)	! 4Ii1Io -> coag
	coef_lin(26,26,8) = cs(8)	! 3Ii1Io -> coag
	coef_lin(26,26,7) = cs(7)	! 2Ii1Io -> coag
	coef_lin(26,26,6) = cs(6)	! 1Ii1Io -> coag
	coef_lin(26,26,5) = cs(5)	! 1Io -> coag
	coef_lin(26,26,4) = cs(4)	! 4Ii -> coag
	coef_lin(26,26,3) = cs(3)	! 3Ii -> coag
	coef_lin(26,26,2) = cs(2)	! 2Ii -> coag
	coef_lin(26,26,1) = cs(1)	! 1Ii -> coag

end subroutine get_rate_coefs_3

!-----------------------------------------------------------

subroutine get_losses_3(cs,temperature)

	use acdc_simulation_setup_3, only : get_scav_parameters_3

	implicit none
	real(kind(1.d0)) :: cs(24),temperature

	! coagulation sink

	! cs obtained from an external routine
	call get_scav_parameters_3(temperature=temperature,cs_per_clust=cs)

	cs(1) = 0.00000000000000d+00	! coagulation loss of 1Ii
	cs(5) = 0.00000000000000d+00	! coagulation loss of 1Io

end subroutine get_losses_3

!-----------------------------------------------------------

subroutine get_coll_3(K,temperature)
	implicit none
	integer, parameter :: nclust = 24
	real(kind(1.d0)) :: K(nclust,nclust), temperature

	! collision coefficients

	K = 0.d0
	K(1,1) = 1.18908702156147d-17*sqrt(temperature)	! 1Ii + 1Ii
	K(2,1) = 1.31483366159986d-17*sqrt(temperature)	! 2Ii + 1Ii
	K(1,2) = K(2,1)
	K(2,2) = 1.33470505383998d-17*sqrt(temperature)	! 2Ii + 2Ii
	K(3,1) = 1.44773175964888d-17*sqrt(temperature)	! 3Ii + 1Ii
	K(1,3) = K(3,1)
	K(3,2) = 1.40111428930154d-17*sqrt(temperature)	! 3Ii + 2Ii
	K(2,3) = K(3,2)
	K(3,3) = 1.42801854711333d-17*sqrt(temperature)	! 3Ii + 3Ii
	K(4,2) = 1.47585088498549d-17*sqrt(temperature)	! 4Ii + 2Ii
	K(2,4) = K(4,2)
	K(4,3) = 1.47360858963629d-17*sqrt(temperature)	! 4Ii + 3Ii
	K(3,4) = K(4,3)
	K(4,4) = 1.49815576862209d-17*sqrt(temperature)	! 4Ii + 4Ii
	K(5,1) = 1.18064366217448d-17*sqrt(temperature)	! 1Io + 1Ii
	K(1,5) = K(5,1)
	K(5,2) = 1.32063778466581d-17*sqrt(temperature)	! 1Io + 2Ii
	K(2,5) = K(5,2)
	K(5,3) = 1.46285619551422d-17*sqrt(temperature)	! 1Io + 3Ii
	K(3,5) = K(5,3)
	K(5,4) = 1.59578717110122d-17*sqrt(temperature)	! 1Io + 4Ii
	K(4,5) = K(5,4)
	K(5,5) = 1.17033759113648d-17*sqrt(temperature)	! 1Io + 1Io
	K(6,1) = 1.30258238670595d-17*sqrt(temperature)	! 1Ii1Io + 1Ii
	K(1,6) = K(6,1)
	K(6,2) = 1.32979819656409d-17*sqrt(temperature)	! 1Ii1Io + 2Ii
	K(2,6) = K(6,2)
	K(6,3) = 1.40067276930756d-17*sqrt(temperature)	! 1Ii1Io + 3Ii
	K(3,6) = K(6,3)
	K(6,4) = 1.47875153272340d-17*sqrt(temperature)	! 1Ii1Io + 4Ii
	K(4,6) = K(6,4)
	K(6,5) = 1.30737492072270d-17*sqrt(temperature)	! 1Ii1Io + 1Io
	K(5,6) = K(6,5)
	K(6,6) = 1.32439127328180d-17*sqrt(temperature)	! 1Ii1Io + 1Ii1Io
	K(7,1) = 1.43590556305513d-17*sqrt(temperature)	! 2Ii1Io + 1Ii
	K(1,7) = K(7,1)
	K(7,2) = 1.39452458943047d-17*sqrt(temperature)	! 2Ii1Io + 2Ii
	K(2,7) = K(7,2)
	K(7,3) = 1.42448135308741d-17*sqrt(temperature)	! 2Ii1Io + 3Ii
	K(3,7) = K(7,3)
	K(7,4) = 1.47229738866166d-17*sqrt(temperature)	! 2Ii1Io + 4Ii
	K(4,7) = K(7,4)
	K(7,5) = 1.45028527947434d-17*sqrt(temperature)	! 2Ii1Io + 1Io
	K(5,7) = K(7,5)
	K(7,6) = 1.39373540072696d-17*sqrt(temperature)	! 2Ii1Io + 1Ii1Io
	K(6,7) = K(7,6)
	K(7,7) = 1.42070975967104d-17*sqrt(temperature)	! 2Ii1Io + 2Ii1Io
	K(8,1) = 1.56224681431055d-17*sqrt(temperature)	! 3Ii1Io + 1Ii
	K(1,8) = K(8,1)
	K(8,2) = 1.46897264469969d-17*sqrt(temperature)	! 3Ii1Io + 2Ii
	K(2,8) = K(8,2)
	K(8,3) = 1.46910843137569d-17*sqrt(temperature)	! 3Ii1Io + 3Ii
	K(3,8) = K(8,3)
	K(8,4) = 1.49535829663980d-17*sqrt(temperature)	! 3Ii1Io + 4Ii
	K(4,8) = K(8,4)
	K(8,5) = 1.58408633387067d-17*sqrt(temperature)	! 3Ii1Io + 1Io
	K(5,8) = K(8,5)
	K(8,6) = 1.47159906668696d-17*sqrt(temperature)	! 3Ii1Io + 1Ii1Io
	K(6,8) = K(8,6)
	K(8,7) = 1.46761643291989d-17*sqrt(temperature)	! 3Ii1Io + 2Ii1Io
	K(7,8) = K(8,7)
	K(8,8) = 1.49242347905704d-17*sqrt(temperature)	! 3Ii1Io + 3Ii1Io
	K(9,2) = 1.54440247279563d-17*sqrt(temperature)	! 4Ii1Io + 2Ii
	K(2,9) = K(9,2)
	K(9,3) = 1.52053612336831d-17*sqrt(temperature)	! 4Ii1Io + 3Ii
	K(3,9) = K(9,3)
	K(9,4) = 1.52971500958417d-17*sqrt(temperature)	! 4Ii1Io + 4Ii
	K(4,9) = K(9,4)
	K(9,5) = 1.70880663616102d-17*sqrt(temperature)	! 4Ii1Io + 1Io
	K(5,9) = K(9,5)
	K(9,6) = 1.54980932785864d-17*sqrt(temperature)	! 4Ii1Io + 1Ii1Io
	K(6,9) = K(9,6)
	K(9,7) = 1.52086222491292d-17*sqrt(temperature)	! 4Ii1Io + 2Ii1Io
	K(7,9) = K(9,7)
	K(9,8) = 1.52815384106407d-17*sqrt(temperature)	! 4Ii1Io + 3Ii1Io
	K(8,9) = K(9,8)
	K(9,9) = 1.55017168155779d-17*sqrt(temperature)	! 4Ii1Io + 4Ii1Io
	K(10,1) = 1.29035923248779d-17*sqrt(temperature)	! 2Io + 1Ii
	K(1,10) = K(10,1)
	K(10,2) = 1.32522770336785d-17*sqrt(temperature)	! 2Io + 2Ii
	K(2,10) = K(10,2)
	K(10,3) = 1.40077717031022d-17*sqrt(temperature)	! 2Io + 3Ii
	K(3,10) = K(10,3)
	K(10,4) = 1.48236579285073d-17*sqrt(temperature)	! 2Io + 4Ii
	K(4,10) = K(10,4)
	K(10,5) = 1.29410146806684d-17*sqrt(temperature)	! 2Io + 1Io
	K(5,10) = K(10,5)
	K(10,6) = 1.31929827077980d-17*sqrt(temperature)	! 2Io + 1Ii1Io
	K(6,10) = K(10,6)
	K(10,7) = 1.39347522486206d-17*sqrt(temperature)	! 2Io + 2Ii1Io
	K(7,10) = K(10,7)
	K(10,8) = 1.47492498949607d-17*sqrt(temperature)	! 2Io + 3Ii1Io
	K(8,10) = K(10,8)
	K(10,9) = 1.55606167286879d-17*sqrt(temperature)	! 2Io + 4Ii1Io
	K(9,10) = K(10,9)
	K(10,10) = 1.31365952976051d-17*sqrt(temperature)	! 2Io + 2Io
	K(11,1) = 1.42401575533924d-17*sqrt(temperature)	! 1Ii2Io + 1Ii
	K(1,11) = K(11,1)
	K(11,2) = 1.38800079229097d-17*sqrt(temperature)	! 1Ii2Io + 2Ii
	K(2,11) = K(11,2)
	K(11,3) = 1.42109582216599d-17*sqrt(temperature)	! 1Ii2Io + 3Ii
	K(3,11) = K(11,3)
	K(11,4) = 1.47120484705213d-17*sqrt(temperature)	! 1Ii2Io + 4Ii
	K(4,11) = K(11,4)
	K(11,5) = 1.43763422683783d-17*sqrt(temperature)	! 1Ii2Io + 1Io
	K(5,11) = K(11,5)
	K(11,6) = 1.38685462884339d-17*sqrt(temperature)	! 1Ii2Io + 1Ii1Io
	K(6,11) = K(11,6)
	K(11,7) = 1.41708299125391d-17*sqrt(temperature)	! 1Ii2Io + 2Ii1Io
	K(7,11) = K(11,7)
	K(11,8) = 1.46633751350798d-17*sqrt(temperature)	! 1Ii2Io + 3Ii1Io
	K(8,11) = K(11,8)
	K(11,9) = 1.52145864415270d-17*sqrt(temperature)	! 1Ii2Io + 4Ii1Io
	K(9,11) = K(11,9)
	K(11,10) = 1.38622022202634d-17*sqrt(temperature)	! 1Ii2Io + 2Io
	K(10,11) = K(11,10)
	K(11,11) = 1.41320799106120d-17*sqrt(temperature)	! 1Ii2Io + 1Ii2Io
	K(12,1) = 1.55109244812529d-17*sqrt(temperature)	! 2Ii2Io + 1Ii
	K(1,12) = K(12,1)
	K(12,2) = 1.46209929387367d-17*sqrt(temperature)	! 2Ii2Io + 2Ii
	K(2,12) = K(12,2)
	K(12,3) = 1.46465961541941d-17*sqrt(temperature)	! 2Ii2Io + 3Ii
	K(3,12) = K(12,3)
	K(12,4) = 1.49264800755343d-17*sqrt(temperature)	! 2Ii2Io + 4Ii
	K(4,12) = K(12,4)
	K(12,5) = 1.57230993187833d-17*sqrt(temperature)	! 2Ii2Io + 1Io
	K(5,12) = K(12,5)
	K(12,6) = 1.46444641617005d-17*sqrt(temperature)	! 2Ii2Io + 1Ii1Io
	K(6,12) = K(12,6)
	K(12,7) = 1.46298313746497d-17*sqrt(temperature)	! 2Ii2Io + 2Ii1Io
	K(7,12) = K(12,7)
	K(12,8) = 1.48957287626647d-17*sqrt(temperature)	! 2Ii2Io + 3Ii1Io
	K(8,12) = K(12,8)
	K(12,9) = 1.52670717343411d-17*sqrt(temperature)	! 2Ii2Io + 4Ii1Io
	K(9,12) = K(12,9)
	K(12,10) = 1.46747873142271d-17*sqrt(temperature)	! 2Ii2Io + 2Io
	K(10,12) = K(12,10)
	K(12,11) = 1.46151406991697d-17*sqrt(temperature)	! 2Ii2Io + 1Ii2Io
	K(11,12) = K(12,11)
	K(12,12) = 1.48657894137095d-17*sqrt(temperature)	! 2Ii2Io + 2Ii2Io
	K(13,1) = 1.67019470381656d-17*sqrt(temperature)	! 3Ii2Io + 1Ii
	K(1,13) = K(13,1)
	K(13,2) = 1.53758264059820d-17*sqrt(temperature)	! 3Ii2Io + 2Ii
	K(2,13) = K(13,2)
	K(13,3) = 1.51572774808914d-17*sqrt(temperature)	! 3Ii2Io + 3Ii
	K(3,13) = K(13,3)
	K(13,4) = 1.52632753509357d-17*sqrt(temperature)	! 3Ii2Io + 4Ii
	K(4,13) = K(13,4)
	K(13,5) = 1.69780356430948d-17*sqrt(temperature)	! 3Ii2Io + 1Io
	K(5,13) = K(13,5)
	K(13,6) = 1.54275540586532d-17*sqrt(temperature)	! 3Ii2Io + 1Ii1Io
	K(6,13) = K(13,6)
	K(13,7) = 1.51590205763308d-17*sqrt(temperature)	! 3Ii2Io + 2Ii1Io
	K(7,13) = K(13,7)
	K(13,8) = 1.52465231585432d-17*sqrt(temperature)	! 3Ii2Io + 3Ii1Io
	K(8,13) = K(13,8)
	K(13,9) = 1.54780649717345d-17*sqrt(temperature)	! 3Ii2Io + 4Ii1Io
	K(9,13) = K(13,9)
	K(13,10) = 1.54876112243261d-17*sqrt(temperature)	! 3Ii2Io + 2Io
	K(10,13) = K(13,10)
	K(13,11) = 1.51634181036549d-17*sqrt(temperature)	! 3Ii2Io + 1Ii2Io
	K(11,13) = K(13,11)
	K(13,12) = 1.52308903217961d-17*sqrt(temperature)	! 3Ii2Io + 2Ii2Io
	K(12,13) = K(13,12)
	K(13,13) = 1.54534732805773d-17*sqrt(temperature)	! 3Ii2Io + 3Ii2Io
	K(14,2) = 1.61174949948569d-17*sqrt(temperature)	! 4Ii2Io + 2Ii
	K(2,14) = K(14,2)
	K(14,3) = 1.56922118956600d-17*sqrt(temperature)	! 4Ii2Io + 3Ii
	K(3,14) = K(14,3)
	K(14,4) = 1.56530092706992d-17*sqrt(temperature)	! 4Ii2Io + 4Ii
	K(4,14) = K(14,4)
	K(14,5) = 1.81542074375935d-17*sqrt(temperature)	! 4Ii2Io + 1Io
	K(5,14) = K(14,5)
	K(14,6) = 1.61933760687658d-17*sqrt(temperature)	! 4Ii2Io + 1Ii1Io
	K(6,14) = K(14,6)
	K(14,7) = 1.57095063818453d-17*sqrt(temperature)	! 4Ii2Io + 2Ii1Io
	K(7,14) = K(14,7)
	K(14,8) = 1.56478874641713d-17*sqrt(temperature)	! 4Ii2Io + 3Ii1Io
	K(8,14) = K(14,8)
	K(14,9) = 1.57637334076994d-17*sqrt(temperature)	! 4Ii2Io + 4Ii1Io
	K(9,14) = K(14,9)
	K(14,10) = 1.62789018543170d-17*sqrt(temperature)	! 4Ii2Io + 2Io
	K(10,14) = K(14,10)
	K(14,11) = 1.57299630027022d-17*sqrt(temperature)	! 4Ii2Io + 1Ii2Io
	K(11,14) = K(14,11)
	K(14,12) = 1.56441509270127d-17*sqrt(temperature)	! 4Ii2Io + 2Ii2Io
	K(12,14) = K(14,12)
	K(14,13) = 1.57486954483769d-17*sqrt(temperature)	! 4Ii2Io + 3Ii2Io
	K(13,14) = K(14,13)
	K(14,14) = 1.59469278689348d-17*sqrt(temperature)	! 4Ii2Io + 4Ii2Io
	K(15,1) = 1.41206442197815d-17*sqrt(temperature)	! 3Io + 1Ii
	K(1,15) = K(15,1)
	K(15,2) = 1.38155365153126d-17*sqrt(temperature)	! 3Io + 2Ii
	K(2,15) = K(15,2)
	K(15,3) = 1.41787860479967d-17*sqrt(temperature)	! 3Io + 3Ii
	K(3,15) = K(15,3)
	K(15,4) = 1.47035233355283d-17*sqrt(temperature)	! 3Io + 4Ii
	K(4,15) = K(15,4)
	K(15,5) = 1.42490403938166d-17*sqrt(temperature)	! 3Io + 1Io
	K(5,15) = K(15,5)
	K(15,6) = 1.38004057547406d-17*sqrt(temperature)	! 3Io + 1Ii1Io
	K(6,15) = K(15,6)
	K(15,7) = 1.41361722264686d-17*sqrt(temperature)	! 3Io + 2Ii1Io
	K(7,15) = K(15,7)
	K(15,8) = 1.46529264466420d-17*sqrt(temperature)	! 3Io + 3Ii1Io
	K(8,15) = K(15,8)
	K(15,9) = 1.52235046487874d-17*sqrt(temperature)	! 3Io + 4Ii1Io
	K(9,15) = K(15,9)
	K(15,10) = 1.37902162963208d-17*sqrt(temperature)	! 3Io + 2Io
	K(10,15) = K(15,10)
	K(15,11) = 1.40948651095745d-17*sqrt(temperature)	! 3Io + 1Ii2Io
	K(11,15) = K(15,11)
	K(15,12) = 1.46027298159924d-17*sqrt(temperature)	! 3Io + 2Ii2Io
	K(12,15) = K(15,12)
	K(15,13) = 1.51707173659462d-17*sqrt(temperature)	! 3Io + 3Ii2Io
	K(13,15) = K(15,13)
	K(15,14) = 1.57538661722730d-17*sqrt(temperature)	! 3Io + 4Ii2Io
	K(14,15) = K(15,14)
	K(15,15) = 1.40550166322746d-17*sqrt(temperature)	! 3Io + 3Io
	K(16,1) = 1.53987119819895d-17*sqrt(temperature)	! 1Ii3Io + 1Ii
	K(1,16) = K(16,1)
	K(16,2) = 1.45523374868251d-17*sqrt(temperature)	! 1Ii3Io + 2Ii
	K(2,16) = K(16,2)
	K(16,3) = 1.46026736773777d-17*sqrt(temperature)	! 1Ii3Io + 3Ii
	K(3,16) = K(16,3)
	K(16,4) = 1.49003194209100d-17*sqrt(temperature)	! 1Ii3Io + 4Ii
	K(4,16) = K(16,4)
	K(16,5) = 1.56045699743831d-17*sqrt(temperature)	! 1Ii3Io + 1Io
	K(5,16) = K(16,5)
	K(16,6) = 1.45729624761553d-17*sqrt(temperature)	! 1Ii3Io + 1Ii1Io
	K(6,16) = K(16,6)
	K(16,7) = 1.45840254357293d-17*sqrt(temperature)	! 1Ii3Io + 2Ii1Io
	K(7,16) = K(16,7)
	K(16,8) = 1.48681337744890d-17*sqrt(temperature)	! 1Ii3Io + 3Ii1Io
	K(8,16) = K(16,8)
	K(16,9) = 1.52538345417531d-17*sqrt(temperature)	! 1Ii3Io + 4Ii1Io
	K(9,16) = K(16,9)
	K(16,10) = 1.46002942622132d-17*sqrt(temperature)	! 1Ii3Io + 2Io
	K(10,16) = K(16,10)
	K(16,11) = 1.45673936873684d-17*sqrt(temperature)	! 1Ii3Io + 1Ii2Io
	K(11,16) = K(16,11)
	K(16,12) = 1.48367293817953d-17*sqrt(temperature)	! 1Ii3Io + 2Ii2Io
	K(12,16) = K(16,12)
	K(16,13) = 1.52164599849689d-17*sqrt(temperature)	! 1Ii3Io + 3Ii2Io
	K(13,16) = K(16,13)
	K(16,14) = 1.56418966957763d-17*sqrt(temperature)	! 1Ii3Io + 4Ii2Io
	K(14,16) = K(16,14)
	K(16,15) = 1.45529800361419d-17*sqrt(temperature)	! 1Ii3Io + 3Io
	K(15,16) = K(16,15)
	K(16,16) = 1.48061720329670d-17*sqrt(temperature)	! 1Ii3Io + 1Ii3Io
	K(17,1) = 1.65967190554506d-17*sqrt(temperature)	! 2Ii3Io + 1Ii
	K(1,17) = K(17,1)
	K(17,2) = 1.53075085098254d-17*sqrt(temperature)	! 2Ii3Io + 2Ii
	K(2,17) = K(17,2)
	K(17,3) = 1.51093767164980d-17*sqrt(temperature)	! 2Ii3Io + 3Ii
	K(3,17) = K(17,3)
	K(17,4) = 1.52298153965683d-17*sqrt(temperature)	! 2Ii3Io + 4Ii
	K(4,17) = K(17,4)
	K(17,5) = 1.68673543634796d-17*sqrt(temperature)	! 2Ii3Io + 1Io
	K(5,17) = K(17,5)
	K(17,6) = 1.53568620758583d-17*sqrt(temperature)	! 2Ii3Io + 1Ii1Io
	K(6,17) = K(17,6)
	K(17,7) = 1.51095779930606d-17*sqrt(temperature)	! 2Ii3Io + 2Ii1Io
	K(7,17) = K(17,7)
	K(17,8) = 1.52119035616620d-17*sqrt(temperature)	! 2Ii3Io + 3Ii1Io
	K(8,17) = K(17,8)
	K(17,9) = 1.54550034551659d-17*sqrt(temperature)	! 2Ii3Io + 4Ii1Io
	K(9,17) = K(17,9)
	K(17,10) = 1.54144184524324d-17*sqrt(temperature)	! 2Ii3Io + 2Io
	K(10,17) = K(17,10)
	K(17,11) = 1.51123843679686d-17*sqrt(temperature)	! 2Ii3Io + 1Ii2Io
	K(11,17) = K(17,11)
	K(17,12) = 1.51950850928954d-17*sqrt(temperature)	! 2Ii3Io + 2Ii2Io
	K(12,17) = K(17,12)
	K(17,13) = 1.54294554793822d-17*sqrt(temperature)	! 2Ii3Io + 3Ii2Io
	K(13,17) = K(17,13)
	K(17,14) = 1.57344012636459d-17*sqrt(temperature)	! 2Ii3Io + 4Ii2Io
	K(14,17) = K(17,14)
	K(17,15) = 1.51180395730010d-17*sqrt(temperature)	! 2Ii3Io + 3Io
	K(15,17) = K(17,15)
	K(17,16) = 1.51794417886290d-17*sqrt(temperature)	! 2Ii3Io + 1Ii3Io
	K(16,17) = K(17,16)
	K(17,17) = 1.54044647306236d-17*sqrt(temperature)	! 2Ii3Io + 2Ii3Io
	K(18,1) = 1.77226022842801d-17*sqrt(temperature)	! 3Ii3Io + 1Ii
	K(1,18) = K(18,1)
	K(18,2) = 1.60508536388940d-17*sqrt(temperature)	! 3Ii3Io + 2Ii
	K(2,18) = K(18,2)
	K(18,3) = 1.56431932535653d-17*sqrt(temperature)	! 3Ii3Io + 3Ii
	K(3,18) = K(18,3)
	K(18,4) = 1.56162960552554d-17*sqrt(temperature)	! 3Ii3Io + 4Ii
	K(4,18) = K(18,4)
	K(18,5) = 1.80501361818480d-17*sqrt(temperature)	! 3Ii3Io + 1Io
	K(5,18) = K(18,5)
	K(18,6) = 1.61246676968525d-17*sqrt(temperature)	! 3Ii3Io + 1Ii1Io
	K(6,18) = K(18,6)
	K(18,7) = 1.56591663502189d-17*sqrt(temperature)	! 3Ii3Io + 2Ii1Io
	K(7,18) = K(18,7)
	K(18,8) = 1.56101906934527d-17*sqrt(temperature)	! 3Ii3Io + 3Ii1Io
	K(8,18) = K(18,8)
	K(18,9) = 1.57358050841694d-17*sqrt(temperature)	! 3Ii3Io + 4Ii1Io
	K(9,18) = K(18,9)
	K(18,10) = 1.62080120747642d-17*sqrt(temperature)	! 3Ii3Io + 2Io
	K(10,18) = K(18,10)
	K(18,11) = 1.56782577265358d-17*sqrt(temperature)	! 3Ii3Io + 1Ii2Io
	K(11,18) = K(18,11)
	K(18,12) = 1.56054477424591d-17*sqrt(temperature)	! 3Ii3Io + 2Ii2Io
	K(12,18) = K(18,12)
	K(18,13) = 1.57199617895623d-17*sqrt(temperature)	! 3Ii3Io + 3Ii2Io
	K(13,18) = K(18,13)
	K(18,14) = 1.59263592828332d-17*sqrt(temperature)	! 3Ii3Io + 4Ii2Io
	K(14,18) = K(18,14)
	K(18,15) = 1.57007485656619d-17*sqrt(temperature)	! 3Ii3Io + 3Io
	K(15,18) = K(18,15)
	K(18,16) = 1.56021630357829d-17*sqrt(temperature)	! 3Ii3Io + 1Ii3Io
	K(16,18) = K(18,16)
	K(18,17) = 1.57048477101259d-17*sqrt(temperature)	! 3Ii3Io + 2Ii3Io
	K(17,18) = K(18,17)
	K(18,18) = 1.59051042319672d-17*sqrt(temperature)	! 3Ii3Io + 3Ii3Io
	K(19,2) = 1.67734252181848d-17*sqrt(temperature)	! 4Ii3Io + 2Ii
	K(2,19) = K(19,2)
	K(19,3) = 1.61823841516934d-17*sqrt(temperature)	! 4Ii3Io + 3Ii
	K(3,19) = K(19,3)
	K(19,4) = 1.60280842071345d-17*sqrt(temperature)	! 4Ii3Io + 4Ii
	K(4,19) = K(19,4)
	K(19,5) = 1.91657022696002d-17*sqrt(temperature)	! 4Ii3Io + 1Io
	K(5,19) = K(19,5)
	K(19,6) = 1.68688213640403d-17*sqrt(temperature)	! 4Ii3Io + 1Ii1Io
	K(6,19) = K(19,6)
	K(19,7) = 1.62120691415508d-17*sqrt(temperature)	! 4Ii3Io + 2Ii1Io
	K(7,19) = K(19,7)
	K(19,8) = 1.60321431629445d-17*sqrt(temperature)	! 4Ii3Io + 3Ii1Io
	K(8,19) = K(19,8)
	K(19,9) = 1.60569474138282d-17*sqrt(temperature)	! 4Ii3Io + 4Ii1Io
	K(9,19) = K(19,9)
	K(19,10) = 1.69749586545773d-17*sqrt(temperature)	! 4Ii3Io + 2Io
	K(10,19) = K(19,10)
	K(19,11) = 1.62453339517400d-17*sqrt(temperature)	! 4Ii3Io + 1Ii2Io
	K(11,19) = K(19,11)
	K(19,12) = 1.60378040761808d-17*sqrt(temperature)	! 4Ii3Io + 2Ii2Io
	K(12,19) = K(19,12)
	K(19,13) = 1.60494035835451d-17*sqrt(temperature)	! 4Ii3Io + 3Ii2Io
	K(13,19) = K(19,13)
	K(19,14) = 1.61717489314726d-17*sqrt(temperature)	! 4Ii3Io + 4Ii2Io
	K(14,19) = K(19,14)
	K(19,15) = 1.62824939963151d-17*sqrt(temperature)	! 4Ii3Io + 3Io
	K(15,19) = K(19,15)
	K(19,16) = 1.60451754710442d-17*sqrt(temperature)	! 4Ii3Io + 1Ii3Io
	K(16,19) = K(19,16)
	K(19,17) = 1.60427411199568d-17*sqrt(temperature)	! 4Ii3Io + 2Ii3Io
	K(17,19) = K(19,17)
	K(19,18) = 1.61575524761397d-17*sqrt(temperature)	! 4Ii3Io + 3Ii3Io
	K(18,19) = K(19,18)
	K(19,19) = 1.63374890709404d-17*sqrt(temperature)	! 4Ii3Io + 4Ii3Io
	K(20,1) = 1.52858261149791d-17*sqrt(temperature)	! 4Io + 1Ii
	K(1,20) = K(20,1)
	K(20,2) = 1.44837926715315d-17*sqrt(temperature)	! 4Io + 2Ii
	K(2,20) = K(20,2)
	K(20,3) = 1.45593743374062d-17*sqrt(temperature)	! 4Io + 3Ii
	K(3,20) = K(20,3)
	K(20,4) = 1.48751780239860d-17*sqrt(temperature)	! 4Io + 4Ii
	K(4,20) = K(20,4)
	K(20,6) = 1.45015155027024d-17*sqrt(temperature)	! 4Io + 1Ii1Io
	K(6,20) = K(20,6)
	K(20,7) = 1.45388019764488d-17*sqrt(temperature)	! 4Io + 2Ii1Io
	K(7,20) = K(20,7)
	K(20,8) = 1.48415252124912d-17*sqrt(temperature)	! 4Io + 3Ii1Io
	K(8,20) = K(20,8)
	K(20,9) = 1.52419190403678d-17*sqrt(temperature)	! 4Io + 4Ii1Io
	K(9,20) = K(20,9)
	K(20,10) = 1.45257978456648d-17*sqrt(temperature)	! 4Io + 2Io
	K(10,20) = K(20,10)
	K(20,11) = 1.45201875293978d-17*sqrt(temperature)	! 4Io + 1Ii2Io
	K(11,20) = K(20,11)
	K(20,12) = 1.48086223935051d-17*sqrt(temperature)	! 4Io + 2Ii2Io
	K(12,20) = K(20,12)
	K(20,13) = 1.52033229185456d-17*sqrt(temperature)	! 4Io + 3Ii2Io
	K(13,20) = K(20,13)
	K(20,14) = 1.56412305516703d-17*sqrt(temperature)	! 4Io + 4Ii2Io
	K(14,20) = K(20,14)
	K(20,15) = 1.45037284555360d-17*sqrt(temperature)	! 4Io + 3Io
	K(15,20) = K(20,15)
	K(20,16) = 1.47765337686934d-17*sqrt(temperature)	! 4Io + 1Ii3Io
	K(16,20) = K(20,16)
	K(20,17) = 1.51650629658743d-17*sqrt(temperature)	! 4Io + 2Ii3Io
	K(17,20) = K(20,17)
	K(20,18) = 1.56004410541555d-17*sqrt(temperature)	! 4Io + 3Ii3Io
	K(18,20) = K(20,18)
	K(20,19) = 1.60543755518554d-17*sqrt(temperature)	! 4Io + 4Ii3Io
	K(19,20) = K(20,19)
	K(20,20) = 1.47453296655611d-17*sqrt(temperature)	! 4Io + 4Io
	K(21,1) = 1.64908938332254d-17*sqrt(temperature)	! 1Ii4Io + 1Ii
	K(1,21) = K(21,1)
	K(21,2) = 1.52390809511267d-17*sqrt(temperature)	! 1Ii4Io + 2Ii
	K(2,21) = K(21,2)
	K(21,3) = 1.50616808590993d-17*sqrt(temperature)	! 1Ii4Io + 3Ii
	K(3,21) = K(21,3)
	K(21,4) = 1.51968014962082d-17*sqrt(temperature)	! 1Ii4Io + 4Ii
	K(4,21) = K(21,4)
	K(21,6) = 1.52860259383652d-17*sqrt(temperature)	! 1Ii4Io + 1Ii1Io
	K(6,21) = K(21,6)
	K(21,7) = 1.50603154618510d-17*sqrt(temperature)	! 1Ii4Io + 2Ii1Io
	K(7,21) = K(21,7)
	K(21,8) = 1.51777101067943d-17*sqrt(temperature)	! 1Ii4Io + 3Ii1Io
	K(8,21) = K(21,8)
	K(21,9) = 1.54325707040825d-17*sqrt(temperature)	! 1Ii4Io + 4Ii1Io
	K(9,21) = K(21,9)
	K(21,10) = 1.53410456692660d-17*sqrt(temperature)	! 1Ii4Io + 2Io
	K(10,21) = K(21,10)
	K(21,11) = 1.50615052187473d-17*sqrt(temperature)	! 1Ii4Io + 1Ii2Io
	K(11,21) = K(21,11)
	K(21,12) = 1.51596857451440d-17*sqrt(temperature)	! 1Ii4Io + 2Ii2Io
	K(12,21) = K(21,12)
	K(21,13) = 1.54060493304800d-17*sqrt(temperature)	! 1Ii4Io + 3Ii2Io
	K(13,21) = K(21,13)
	K(21,14) = 1.57208956450324d-17*sqrt(temperature)	! 1Ii4Io + 4Ii2Io
	K(14,21) = K(21,14)
	K(21,15) = 1.50654902524234d-17*sqrt(temperature)	! 1Ii4Io + 3Io
	K(15,21) = K(21,15)
	K(21,16) = 1.51428088477377d-17*sqrt(temperature)	! 1Ii4Io + 1Ii3Io
	K(16,21) = K(21,16)
	K(21,17) = 1.53800682869494d-17*sqrt(temperature)	! 1Ii4Io + 2Ii3Io
	K(17,21) = K(21,17)
	K(21,18) = 1.56905070310608d-17*sqrt(temperature)	! 1Ii4Io + 3Ii3Io
	K(18,21) = K(21,18)
	K(21,19) = 1.60370105826641d-17*sqrt(temperature)	! 1Ii4Io + 4Ii3Io
	K(19,21) = K(21,19)
	K(21,20) = 1.51271672610203d-17*sqrt(temperature)	! 1Ii4Io + 4Io
	K(20,21) = K(21,20)
	K(21,21) = 1.53546639585670d-17*sqrt(temperature)	! 1Ii4Io + 1Ii4Io
	K(22,1) = 1.76229001732328d-17*sqrt(temperature)	! 2Ii4Io + 1Ii
	K(1,22) = K(22,1)
	K(22,2) = 1.59840364156290d-17*sqrt(temperature)	! 2Ii4Io + 2Ii
	K(2,22) = K(22,2)
	K(22,3) = 1.55942116534802d-17*sqrt(temperature)	! 2Ii4Io + 3Ii
	K(3,22) = K(22,3)
	K(22,4) = 1.55797821236995d-17*sqrt(temperature)	! 2Ii4Io + 4Ii
	K(4,22) = K(22,4)
	K(22,6) = 1.60557600144292d-17*sqrt(temperature)	! 2Ii4Io + 1Ii1Io
	K(6,22) = K(22,6)
	K(22,7) = 1.56088465909416d-17*sqrt(temperature)	! 2Ii4Io + 2Ii1Io
	K(7,22) = K(22,7)
	K(22,8) = 1.55726798429797d-17*sqrt(temperature)	! 2Ii4Io + 3Ii1Io
	K(8,22) = K(22,8)
	K(22,9) = 1.57081983597200d-17*sqrt(temperature)	! 2Ii4Io + 4Ii1Io
	K(9,22) = K(22,9)
	K(22,10) = 1.61368985897606d-17*sqrt(temperature)	! 2Ii4Io + 2Io
	K(10,22) = K(22,10)
	K(22,11) = 1.56265555320680d-17*sqrt(temperature)	! 2Ii4Io + 1Ii2Io
	K(11,22) = K(22,11)
	K(22,12) = 1.55669168739888d-17*sqrt(temperature)	! 2Ii4Io + 2Ii2Io
	K(12,22) = K(22,12)
	K(22,13) = 1.56915382983923d-17*sqrt(temperature)	! 2Ii4Io + 3Ii2Io
	K(13,22) = K(22,13)
	K(22,14) = 1.59062187698723d-17*sqrt(temperature)	! 2Ii4Io + 4Ii2Io
	K(14,22) = K(22,14)
	K(22,15) = 1.56476164028531d-17*sqrt(temperature)	! 2Ii4Io + 3Io
	K(15,22) = K(22,15)
	K(22,16) = 1.55625878362344d-17*sqrt(temperature)	! 2Ii4Io + 1Ii3Io
	K(16,22) = K(22,16)
	K(22,17) = 1.56755927310290d-17*sqrt(temperature)	! 2Ii4Io + 2Ii3Io
	K(17,22) = K(22,17)
	K(22,18) = 1.58842671712148d-17*sqrt(temperature)	! 2Ii4Io + 3Ii3Io
	K(18,22) = K(22,18)
	K(22,19) = 1.61438791691895d-17*sqrt(temperature)	! 2Ii4Io + 4Ii3Io
	K(19,22) = K(22,19)
	K(22,20) = 1.55597958982355d-17*sqrt(temperature)	! 2Ii4Io + 4Io
	K(20,22) = K(22,20)
	K(22,21) = 1.56604052305331d-17*sqrt(temperature)	! 2Ii4Io + 1Ii4Io
	K(21,22) = K(22,21)
	K(22,22) = 1.58627233633372d-17*sqrt(temperature)	! 2Ii4Io + 2Ii4Io
	K(23,1) = 1.86922201392806d-17*sqrt(temperature)	! 3Ii4Io + 1Ii
	K(1,23) = K(23,1)
	K(23,2) = 1.67086421166608d-17*sqrt(temperature)	! 3Ii4Io + 2Ii
	K(2,23) = K(23,2)
	K(23,3) = 1.61334171994947d-17*sqrt(temperature)	! 3Ii4Io + 3Ii
	K(3,23) = K(23,3)
	K(23,4) = 1.59900494823024d-17*sqrt(temperature)	! 3Ii4Io + 4Ii
	K(4,23) = K(23,4)
	K(23,6) = 1.68021705549531d-17*sqrt(temperature)	! 3Ii4Io + 1Ii1Io
	K(6,23) = K(23,6)
	K(23,7) = 1.61619228943298d-17*sqrt(temperature)	! 3Ii4Io + 2Ii1Io
	K(7,23) = K(23,7)
	K(23,8) = 1.59932378868562d-17*sqrt(temperature)	! 3Ii4Io + 3Ii1Io
	K(8,23) = K(23,8)
	K(23,9) = 1.60266642060422d-17*sqrt(temperature)	! 3Ii4Io + 4Ii1Io
	K(9,23) = K(23,9)
	K(23,10) = 1.69063339381306d-17*sqrt(temperature)	! 3Ii4Io + 2Io
	K(10,23) = K(23,10)
	K(23,11) = 1.61939681534634d-17*sqrt(temperature)	! 3Ii4Io + 1Ii2Io
	K(11,23) = K(23,11)
	K(23,12) = 1.59980074437675d-17*sqrt(temperature)	! 3Ii4Io + 2Ii2Io
	K(12,23) = K(23,12)
	K(23,13) = 1.60184115439382d-17*sqrt(temperature)	! 3Ii4Io + 3Ii2Io
	K(13,23) = K(23,13)
	K(23,14) = 1.61479274245206d-17*sqrt(temperature)	! 3Ii4Io + 4Ii2Io
	K(14,23) = K(23,14)
	K(23,15) = 1.62298653857220d-17*sqrt(temperature)	! 3Ii4Io + 3Io
	K(15,23) = K(23,15)
	K(23,16) = 1.60044655677337d-17*sqrt(temperature)	! 3Ii4Io + 1Ii3Io
	K(16,23) = K(23,16)
	K(23,17) = 1.60110270712207d-17*sqrt(temperature)	! 3Ii4Io + 2Ii3Io
	K(17,23) = K(23,17)
	K(23,18) = 1.61331294823566d-17*sqrt(temperature)	! 3Ii4Io + 3Ii3Io
	K(18,23) = K(23,18)
	K(23,19) = 1.63192391445304d-17*sqrt(temperature)	! 3Ii4Io + 4Ii3Io
	K(19,23) = K(23,19)
	K(23,20) = 1.60127292628362d-17*sqrt(temperature)	! 3Ii4Io + 4Io
	K(20,23) = K(23,20)
	K(23,21) = 1.60045607920488d-17*sqrt(temperature)	! 3Ii4Io + 1Ii4Io
	K(21,23) = K(23,21)
	K(23,22) = 1.61188456001933d-17*sqrt(temperature)	! 3Ii4Io + 2Ii4Io
	K(22,23) = K(23,22)
	K(23,23) = 1.63004641984493d-17*sqrt(temperature)	! 3Ii4Io + 3Ii4Io
	K(24,1) = 1.97080300953997d-17*sqrt(temperature)	! 4Ii4Io + 1Ii
	K(1,24) = K(24,1)
	K(24,2) = 1.74104240643991d-17*sqrt(temperature)	! 4Ii4Io + 2Ii
	K(2,24) = K(24,2)
	K(24,3) = 1.66691413824831d-17*sqrt(temperature)	! 4Ii4Io + 3Ii
	K(3,24) = K(24,3)
	K(24,4) = 1.64115096823312d-17*sqrt(temperature)	! 4Ii4Io + 4Ii
	K(4,24) = K(24,4)
	K(24,5) = 2.01303111286078d-17*sqrt(temperature)	! 4Ii4Io + 1Io
	K(5,24) = K(24,5)
	K(24,6) = 1.75236200650969d-17*sqrt(temperature)	! 4Ii4Io + 1Ii1Io
	K(6,24) = K(24,6)
	K(24,7) = 1.67099981364142d-17*sqrt(temperature)	! 4Ii4Io + 2Ii1Io
	K(7,24) = K(24,7)
	K(24,8) = 1.64237820032555d-17*sqrt(temperature)	! 4Ii4Io + 3Ii1Io
	K(8,24) = K(24,8)
	K(24,9) = 1.63673513658085d-17*sqrt(temperature)	! 4Ii4Io + 4Ii1Io
	K(9,24) = K(24,9)
	K(24,10) = 1.76485822671332d-17*sqrt(temperature)	! 4Ii4Io + 2Io
	K(10,24) = K(24,10)
	K(24,11) = 1.67548212627028d-17*sqrt(temperature)	! 4Ii4Io + 1Ii2Io
	K(11,24) = K(24,11)
	K(24,12) = 1.64378552393355d-17*sqrt(temperature)	! 4Ii4Io + 2Ii2Io
	K(12,24) = K(24,12)
	K(24,13) = 1.63664767387678d-17*sqrt(temperature)	! 4Ii4Io + 3Ii2Io
	K(13,24) = K(24,13)
	K(24,14) = 1.64214304912100d-17*sqrt(temperature)	! 4Ii4Io + 4Ii2Io
	K(14,24) = K(24,14)
	K(24,15) = 1.68039552109373d-17*sqrt(temperature)	! 4Ii4Io + 3Io
	K(15,24) = K(24,15)
	K(24,16) = 1.64538485938063d-17*sqrt(temperature)	! 4Ii4Io + 1Ii3Io
	K(16,24) = K(24,16)
	K(24,17) = 1.63666091347174d-17*sqrt(temperature)	! 4Ii4Io + 2Ii3Io
	K(17,24) = K(24,17)
	K(24,18) = 1.64128808149347d-17*sqrt(temperature)	! 4Ii4Io + 3Ii3Io
	K(18,24) = K(24,18)
	K(24,19) = 1.65349176354163d-17*sqrt(temperature)	! 4Ii4Io + 4Ii3Io
	K(19,24) = K(24,19)
	K(24,20) = 1.64718918272318d-17*sqrt(temperature)	! 4Ii4Io + 4Io
	K(20,24) = K(24,20)
	K(24,21) = 1.63678044426543d-17*sqrt(temperature)	! 4Ii4Io + 1Ii4Io
	K(21,24) = K(24,21)
	K(24,22) = 1.64049407655197d-17*sqrt(temperature)	! 4Ii4Io + 2Ii4Io
	K(22,24) = K(24,22)
	K(24,23) = 1.65215875329922d-17*sqrt(temperature)	! 4Ii4Io + 3Ii4Io
	K(23,24) = K(24,23)
	K(24,24) = 1.66862844350482d-17*sqrt(temperature)	! 4Ii4Io + 4Ii4Io

end subroutine get_coll_3

!-----------------------------------------------------------

subroutine get_evap_3(E,K,temperature)
	implicit none
	real(kind(1.d0)) :: E(24,24), K(24,24), temperature

	! evaporation coefficients

	E = 0.d0
	E(1,1) = 0.5d0*7.33893243358348d+27/temperature*exp(&
			 &((-23.9523596d0/temperature-(-47.68781345d0)/1.d3)&
			 &-0.d0&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(1,1)	! 2Ii -> 1Ii + 1Ii
	E(2,1) = 7.33893243358348d+27/temperature*exp(&
			 &((-46.43243564d0/temperature-(-94.70212967d0)/1.d3)&
			 &-(-23.9523596d0/temperature-(-47.68781345d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(2,1)	! 3Ii -> 2Ii + 1Ii
	E(1,2) = E(2,1)
	E(3,1) = 7.33893243358348d+27/temperature*exp(&
			 &((-77.77838269d0/temperature-(-144.4714705d0)/1.d3)&
			 &-(-46.43243564d0/temperature-(-94.70212967d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(3,1)	! 4Ii -> 3Ii + 1Ii
	E(1,3) = E(3,1)
	E(2,2) = 0.5d0*7.33893243358348d+27/temperature*exp(&
			 &((-77.77838269d0/temperature-(-144.4714705d0)/1.d3)&
			 &-(-23.9523596d0/temperature-(-47.68781345d0)/1.d3)&
			 &-(-23.9523596d0/temperature-(-47.68781345d0)/1.d3))&
			 &*5.03218937158374d+02)*K(2,2)	! 4Ii -> 2Ii + 2Ii
	E(5,1) = 7.33893243358348d+27/temperature*exp(&
			 &((-30.21939422d0/temperature-(-45.50526148d0)/1.d3)&
			 &-0.d0&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(5,1)	! 1Ii1Io -> 1Io + 1Ii
	E(1,5) = E(5,1)
	E(6,1) = 7.33893243358348d+27/temperature*exp(&
			 &((-58.58714661d0/temperature-(-95.81971414d0)/1.d3)&
			 &-(-30.21939422d0/temperature-(-45.50526148d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(6,1)	! 2Ii1Io -> 1Ii1Io + 1Ii
	E(1,6) = E(6,1)
	E(5,2) = 7.33893243358348d+27/temperature*exp(&
			 &((-58.58714661d0/temperature-(-95.81971414d0)/1.d3)&
			 &-0.d0&
			 &-(-23.9523596d0/temperature-(-47.68781345d0)/1.d3))&
			 &*5.03218937158374d+02)*K(5,2)	! 2Ii1Io -> 1Io + 2Ii
	E(2,5) = E(5,2)
	E(7,1) = 7.33893243358348d+27/temperature*exp(&
			 &((-83.27407492d0/temperature-(-142.0468805d0)/1.d3)&
			 &-(-58.58714661d0/temperature-(-95.81971414d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(7,1)	! 3Ii1Io -> 2Ii1Io + 1Ii
	E(1,7) = E(7,1)
	E(6,2) = 7.33893243358348d+27/temperature*exp(&
			 &((-83.27407492d0/temperature-(-142.0468805d0)/1.d3)&
			 &-(-30.21939422d0/temperature-(-45.50526148d0)/1.d3)&
			 &-(-23.9523596d0/temperature-(-47.68781345d0)/1.d3))&
			 &*5.03218937158374d+02)*K(6,2)	! 3Ii1Io -> 1Ii1Io + 2Ii
	E(2,6) = E(6,2)
	E(5,3) = 7.33893243358348d+27/temperature*exp(&
			 &((-83.27407492d0/temperature-(-142.0468805d0)/1.d3)&
			 &-0.d0&
			 &-(-46.43243564d0/temperature-(-94.70212967d0)/1.d3))&
			 &*5.03218937158374d+02)*K(5,3)	! 3Ii1Io -> 1Io + 3Ii
	E(3,5) = E(5,3)
	E(8,1) = 7.33893243358348d+27/temperature*exp(&
			 &((-111.5351121d0/temperature-(-194.8132841d0)/1.d3)&
			 &-(-83.27407492d0/temperature-(-142.0468805d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(8,1)	! 4Ii1Io -> 3Ii1Io + 1Ii
	E(1,8) = E(8,1)
	E(7,2) = 7.33893243358348d+27/temperature*exp(&
			 &((-111.5351121d0/temperature-(-194.8132841d0)/1.d3)&
			 &-(-58.58714661d0/temperature-(-95.81971414d0)/1.d3)&
			 &-(-23.9523596d0/temperature-(-47.68781345d0)/1.d3))&
			 &*5.03218937158374d+02)*K(7,2)	! 4Ii1Io -> 2Ii1Io + 2Ii
	E(2,7) = E(7,2)
	E(6,3) = 7.33893243358348d+27/temperature*exp(&
			 &((-111.5351121d0/temperature-(-194.8132841d0)/1.d3)&
			 &-(-30.21939422d0/temperature-(-45.50526148d0)/1.d3)&
			 &-(-46.43243564d0/temperature-(-94.70212967d0)/1.d3))&
			 &*5.03218937158374d+02)*K(6,3)	! 4Ii1Io -> 1Ii1Io + 3Ii
	E(3,6) = E(6,3)
	E(5,4) = 7.33893243358348d+27/temperature*exp(&
			 &((-111.5351121d0/temperature-(-194.8132841d0)/1.d3)&
			 &-0.d0&
			 &-(-77.77838269d0/temperature-(-144.4714705d0)/1.d3))&
			 &*5.03218937158374d+02)*K(5,4)	! 4Ii1Io -> 1Io + 4Ii
	E(4,5) = E(5,4)
	E(5,5) = 0.5d0*7.33893243358348d+27/temperature*exp(&
			 &((-30.41407318d0/temperature-(-42.69341052d0)/1.d3)&
			 &-0.d0&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(5,5)	! 2Io -> 1Io + 1Io
	E(10,1) = 7.33893243358348d+27/temperature*exp(&
			 &((-62.03219296d0/temperature-(-92.71110348d0)/1.d3)&
			 &-(-30.41407318d0/temperature-(-42.69341052d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(10,1)	! 1Ii2Io -> 2Io + 1Ii
	E(1,10) = E(10,1)
	E(6,5) = 7.33893243358348d+27/temperature*exp(&
			 &((-62.03219296d0/temperature-(-92.71110348d0)/1.d3)&
			 &-(-30.21939422d0/temperature-(-45.50526148d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(6,5)	! 1Ii2Io -> 1Ii1Io + 1Io
	E(5,6) = E(6,5)
	E(11,1) = 7.33893243358348d+27/temperature*exp(&
			 &((-88.1832749d0/temperature-(-139.0540271d0)/1.d3)&
			 &-(-62.03219296d0/temperature-(-92.71110348d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(11,1)	! 2Ii2Io -> 1Ii2Io + 1Ii
	E(1,11) = E(11,1)
	E(10,2) = 7.33893243358348d+27/temperature*exp(&
			 &((-88.1832749d0/temperature-(-139.0540271d0)/1.d3)&
			 &-(-30.41407318d0/temperature-(-42.69341052d0)/1.d3)&
			 &-(-23.9523596d0/temperature-(-47.68781345d0)/1.d3))&
			 &*5.03218937158374d+02)*K(10,2)	! 2Ii2Io -> 2Io + 2Ii
	E(2,10) = E(10,2)
	E(7,5) = 7.33893243358348d+27/temperature*exp(&
			 &((-88.1832749d0/temperature-(-139.0540271d0)/1.d3)&
			 &-(-58.58714661d0/temperature-(-95.81971414d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(7,5)	! 2Ii2Io -> 2Ii1Io + 1Io
	E(5,7) = E(7,5)
	E(6,6) = 0.5d0*7.33893243358348d+27/temperature*exp(&
			 &((-88.1832749d0/temperature-(-139.0540271d0)/1.d3)&
			 &-(-30.21939422d0/temperature-(-45.50526148d0)/1.d3)&
			 &-(-30.21939422d0/temperature-(-45.50526148d0)/1.d3))&
			 &*5.03218937158374d+02)*K(6,6)	! 2Ii2Io -> 1Ii1Io + 1Ii1Io
	E(12,1) = 7.33893243358348d+27/temperature*exp(&
			 &((-119.4487242d0/temperature-(-191.2458533d0)/1.d3)&
			 &-(-88.1832749d0/temperature-(-139.0540271d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(12,1)	! 3Ii2Io -> 2Ii2Io + 1Ii
	E(1,12) = E(12,1)
	E(11,2) = 7.33893243358348d+27/temperature*exp(&
			 &((-119.4487242d0/temperature-(-191.2458533d0)/1.d3)&
			 &-(-62.03219296d0/temperature-(-92.71110348d0)/1.d3)&
			 &-(-23.9523596d0/temperature-(-47.68781345d0)/1.d3))&
			 &*5.03218937158374d+02)*K(11,2)	! 3Ii2Io -> 1Ii2Io + 2Ii
	E(2,11) = E(11,2)
	E(10,3) = 7.33893243358348d+27/temperature*exp(&
			 &((-119.4487242d0/temperature-(-191.2458533d0)/1.d3)&
			 &-(-30.41407318d0/temperature-(-42.69341052d0)/1.d3)&
			 &-(-46.43243564d0/temperature-(-94.70212967d0)/1.d3))&
			 &*5.03218937158374d+02)*K(10,3)	! 3Ii2Io -> 2Io + 3Ii
	E(3,10) = E(10,3)
	E(8,5) = 7.33893243358348d+27/temperature*exp(&
			 &((-119.4487242d0/temperature-(-191.2458533d0)/1.d3)&
			 &-(-83.27407492d0/temperature-(-142.0468805d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(8,5)	! 3Ii2Io -> 3Ii1Io + 1Io
	E(5,8) = E(8,5)
	E(7,6) = 7.33893243358348d+27/temperature*exp(&
			 &((-119.4487242d0/temperature-(-191.2458533d0)/1.d3)&
			 &-(-58.58714661d0/temperature-(-95.81971414d0)/1.d3)&
			 &-(-30.21939422d0/temperature-(-45.50526148d0)/1.d3))&
			 &*5.03218937158374d+02)*K(7,6)	! 3Ii2Io -> 2Ii1Io + 1Ii1Io
	E(6,7) = E(7,6)
	E(13,1) = 7.33893243358348d+27/temperature*exp(&
			 &((-148.7720382d0/temperature-(-241.7960301d0)/1.d3)&
			 &-(-119.4487242d0/temperature-(-191.2458533d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(13,1)	! 4Ii2Io -> 3Ii2Io + 1Ii
	E(1,13) = E(13,1)
	E(12,2) = 7.33893243358348d+27/temperature*exp(&
			 &((-148.7720382d0/temperature-(-241.7960301d0)/1.d3)&
			 &-(-88.1832749d0/temperature-(-139.0540271d0)/1.d3)&
			 &-(-23.9523596d0/temperature-(-47.68781345d0)/1.d3))&
			 &*5.03218937158374d+02)*K(12,2)	! 4Ii2Io -> 2Ii2Io + 2Ii
	E(2,12) = E(12,2)
	E(11,3) = 7.33893243358348d+27/temperature*exp(&
			 &((-148.7720382d0/temperature-(-241.7960301d0)/1.d3)&
			 &-(-62.03219296d0/temperature-(-92.71110348d0)/1.d3)&
			 &-(-46.43243564d0/temperature-(-94.70212967d0)/1.d3))&
			 &*5.03218937158374d+02)*K(11,3)	! 4Ii2Io -> 1Ii2Io + 3Ii
	E(3,11) = E(11,3)
	E(10,4) = 7.33893243358348d+27/temperature*exp(&
			 &((-148.7720382d0/temperature-(-241.7960301d0)/1.d3)&
			 &-(-30.41407318d0/temperature-(-42.69341052d0)/1.d3)&
			 &-(-77.77838269d0/temperature-(-144.4714705d0)/1.d3))&
			 &*5.03218937158374d+02)*K(10,4)	! 4Ii2Io -> 2Io + 4Ii
	E(4,10) = E(10,4)
	E(9,5) = 7.33893243358348d+27/temperature*exp(&
			 &((-148.7720382d0/temperature-(-241.7960301d0)/1.d3)&
			 &-(-111.5351121d0/temperature-(-194.8132841d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(9,5)	! 4Ii2Io -> 4Ii1Io + 1Io
	E(5,9) = E(9,5)
	E(8,6) = 7.33893243358348d+27/temperature*exp(&
			 &((-148.7720382d0/temperature-(-241.7960301d0)/1.d3)&
			 &-(-83.27407492d0/temperature-(-142.0468805d0)/1.d3)&
			 &-(-30.21939422d0/temperature-(-45.50526148d0)/1.d3))&
			 &*5.03218937158374d+02)*K(8,6)	! 4Ii2Io -> 3Ii1Io + 1Ii1Io
	E(6,8) = E(8,6)
	E(7,7) = 0.5d0*7.33893243358348d+27/temperature*exp(&
			 &((-148.7720382d0/temperature-(-241.7960301d0)/1.d3)&
			 &-(-58.58714661d0/temperature-(-95.81971414d0)/1.d3)&
			 &-(-58.58714661d0/temperature-(-95.81971414d0)/1.d3))&
			 &*5.03218937158374d+02)*K(7,7)	! 4Ii2Io -> 2Ii1Io + 2Ii1Io
	E(10,5) = 7.33893243358348d+27/temperature*exp(&
			 &((-60.03371557d0/temperature-(-88.35020888d0)/1.d3)&
			 &-(-30.41407318d0/temperature-(-42.69341052d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(10,5)	! 3Io -> 2Io + 1Io
	E(5,10) = E(10,5)
	E(15,1) = 7.33893243358348d+27/temperature*exp(&
			 &((-89.23828631d0/temperature-(-139.4518114d0)/1.d3)&
			 &-(-60.03371557d0/temperature-(-88.35020888d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(15,1)	! 1Ii3Io -> 3Io + 1Ii
	E(1,15) = E(15,1)
	E(11,5) = 7.33893243358348d+27/temperature*exp(&
			 &((-89.23828631d0/temperature-(-139.4518114d0)/1.d3)&
			 &-(-62.03219296d0/temperature-(-92.71110348d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(11,5)	! 1Ii3Io -> 1Ii2Io + 1Io
	E(5,11) = E(11,5)
	E(10,6) = 7.33893243358348d+27/temperature*exp(&
			 &((-89.23828631d0/temperature-(-139.4518114d0)/1.d3)&
			 &-(-30.41407318d0/temperature-(-42.69341052d0)/1.d3)&
			 &-(-30.21939422d0/temperature-(-45.50526148d0)/1.d3))&
			 &*5.03218937158374d+02)*K(10,6)	! 1Ii3Io -> 2Io + 1Ii1Io
	E(6,10) = E(10,6)
	E(16,1) = 7.33893243358348d+27/temperature*exp(&
			 &((-120.2575982d0/temperature-(-186.1020181d0)/1.d3)&
			 &-(-89.23828631d0/temperature-(-139.4518114d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(16,1)	! 2Ii3Io -> 1Ii3Io + 1Ii
	E(1,16) = E(16,1)
	E(15,2) = 7.33893243358348d+27/temperature*exp(&
			 &((-120.2575982d0/temperature-(-186.1020181d0)/1.d3)&
			 &-(-60.03371557d0/temperature-(-88.35020888d0)/1.d3)&
			 &-(-23.9523596d0/temperature-(-47.68781345d0)/1.d3))&
			 &*5.03218937158374d+02)*K(15,2)	! 2Ii3Io -> 3Io + 2Ii
	E(2,15) = E(15,2)
	E(12,5) = 7.33893243358348d+27/temperature*exp(&
			 &((-120.2575982d0/temperature-(-186.1020181d0)/1.d3)&
			 &-(-88.1832749d0/temperature-(-139.0540271d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(12,5)	! 2Ii3Io -> 2Ii2Io + 1Io
	E(5,12) = E(12,5)
	E(11,6) = 7.33893243358348d+27/temperature*exp(&
			 &((-120.2575982d0/temperature-(-186.1020181d0)/1.d3)&
			 &-(-62.03219296d0/temperature-(-92.71110348d0)/1.d3)&
			 &-(-30.21939422d0/temperature-(-45.50526148d0)/1.d3))&
			 &*5.03218937158374d+02)*K(11,6)	! 2Ii3Io -> 1Ii2Io + 1Ii1Io
	E(6,11) = E(11,6)
	E(10,7) = 7.33893243358348d+27/temperature*exp(&
			 &((-120.2575982d0/temperature-(-186.1020181d0)/1.d3)&
			 &-(-30.41407318d0/temperature-(-42.69341052d0)/1.d3)&
			 &-(-58.58714661d0/temperature-(-95.81971414d0)/1.d3))&
			 &*5.03218937158374d+02)*K(10,7)	! 2Ii3Io -> 2Io + 2Ii1Io
	E(7,10) = E(10,7)
	E(17,1) = 7.33893243358348d+27/temperature*exp(&
			 &((-158.4769381d0/temperature-(-240.2133116d0)/1.d3)&
			 &-(-120.2575982d0/temperature-(-186.1020181d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(17,1)	! 3Ii3Io -> 2Ii3Io + 1Ii
	E(1,17) = E(17,1)
	E(16,2) = 7.33893243358348d+27/temperature*exp(&
			 &((-158.4769381d0/temperature-(-240.2133116d0)/1.d3)&
			 &-(-89.23828631d0/temperature-(-139.4518114d0)/1.d3)&
			 &-(-23.9523596d0/temperature-(-47.68781345d0)/1.d3))&
			 &*5.03218937158374d+02)*K(16,2)	! 3Ii3Io -> 1Ii3Io + 2Ii
	E(2,16) = E(16,2)
	E(15,3) = 7.33893243358348d+27/temperature*exp(&
			 &((-158.4769381d0/temperature-(-240.2133116d0)/1.d3)&
			 &-(-60.03371557d0/temperature-(-88.35020888d0)/1.d3)&
			 &-(-46.43243564d0/temperature-(-94.70212967d0)/1.d3))&
			 &*5.03218937158374d+02)*K(15,3)	! 3Ii3Io -> 3Io + 3Ii
	E(3,15) = E(15,3)
	E(13,5) = 7.33893243358348d+27/temperature*exp(&
			 &((-158.4769381d0/temperature-(-240.2133116d0)/1.d3)&
			 &-(-119.4487242d0/temperature-(-191.2458533d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(13,5)	! 3Ii3Io -> 3Ii2Io + 1Io
	E(5,13) = E(13,5)
	E(12,6) = 7.33893243358348d+27/temperature*exp(&
			 &((-158.4769381d0/temperature-(-240.2133116d0)/1.d3)&
			 &-(-88.1832749d0/temperature-(-139.0540271d0)/1.d3)&
			 &-(-30.21939422d0/temperature-(-45.50526148d0)/1.d3))&
			 &*5.03218937158374d+02)*K(12,6)	! 3Ii3Io -> 2Ii2Io + 1Ii1Io
	E(6,12) = E(12,6)
	E(11,7) = 7.33893243358348d+27/temperature*exp(&
			 &((-158.4769381d0/temperature-(-240.2133116d0)/1.d3)&
			 &-(-62.03219296d0/temperature-(-92.71110348d0)/1.d3)&
			 &-(-58.58714661d0/temperature-(-95.81971414d0)/1.d3))&
			 &*5.03218937158374d+02)*K(11,7)	! 3Ii3Io -> 1Ii2Io + 2Ii1Io
	E(7,11) = E(11,7)
	E(10,8) = 7.33893243358348d+27/temperature*exp(&
			 &((-158.4769381d0/temperature-(-240.2133116d0)/1.d3)&
			 &-(-30.41407318d0/temperature-(-42.69341052d0)/1.d3)&
			 &-(-83.27407492d0/temperature-(-142.0468805d0)/1.d3))&
			 &*5.03218937158374d+02)*K(10,8)	! 3Ii3Io -> 2Io + 3Ii1Io
	E(8,10) = E(10,8)
	E(18,1) = 7.33893243358348d+27/temperature*exp(&
			 &((-187.511044d0/temperature-(-287.5307016d0)/1.d3)&
			 &-(-158.4769381d0/temperature-(-240.2133116d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(18,1)	! 4Ii3Io -> 3Ii3Io + 1Ii
	E(1,18) = E(18,1)
	E(17,2) = 7.33893243358348d+27/temperature*exp(&
			 &((-187.511044d0/temperature-(-287.5307016d0)/1.d3)&
			 &-(-120.2575982d0/temperature-(-186.1020181d0)/1.d3)&
			 &-(-23.9523596d0/temperature-(-47.68781345d0)/1.d3))&
			 &*5.03218937158374d+02)*K(17,2)	! 4Ii3Io -> 2Ii3Io + 2Ii
	E(2,17) = E(17,2)
	E(16,3) = 7.33893243358348d+27/temperature*exp(&
			 &((-187.511044d0/temperature-(-287.5307016d0)/1.d3)&
			 &-(-89.23828631d0/temperature-(-139.4518114d0)/1.d3)&
			 &-(-46.43243564d0/temperature-(-94.70212967d0)/1.d3))&
			 &*5.03218937158374d+02)*K(16,3)	! 4Ii3Io -> 1Ii3Io + 3Ii
	E(3,16) = E(16,3)
	E(15,4) = 7.33893243358348d+27/temperature*exp(&
			 &((-187.511044d0/temperature-(-287.5307016d0)/1.d3)&
			 &-(-60.03371557d0/temperature-(-88.35020888d0)/1.d3)&
			 &-(-77.77838269d0/temperature-(-144.4714705d0)/1.d3))&
			 &*5.03218937158374d+02)*K(15,4)	! 4Ii3Io -> 3Io + 4Ii
	E(4,15) = E(15,4)
	E(14,5) = 7.33893243358348d+27/temperature*exp(&
			 &((-187.511044d0/temperature-(-287.5307016d0)/1.d3)&
			 &-(-148.7720382d0/temperature-(-241.7960301d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(14,5)	! 4Ii3Io -> 4Ii2Io + 1Io
	E(5,14) = E(14,5)
	E(13,6) = 7.33893243358348d+27/temperature*exp(&
			 &((-187.511044d0/temperature-(-287.5307016d0)/1.d3)&
			 &-(-119.4487242d0/temperature-(-191.2458533d0)/1.d3)&
			 &-(-30.21939422d0/temperature-(-45.50526148d0)/1.d3))&
			 &*5.03218937158374d+02)*K(13,6)	! 4Ii3Io -> 3Ii2Io + 1Ii1Io
	E(6,13) = E(13,6)
	E(12,7) = 7.33893243358348d+27/temperature*exp(&
			 &((-187.511044d0/temperature-(-287.5307016d0)/1.d3)&
			 &-(-88.1832749d0/temperature-(-139.0540271d0)/1.d3)&
			 &-(-58.58714661d0/temperature-(-95.81971414d0)/1.d3))&
			 &*5.03218937158374d+02)*K(12,7)	! 4Ii3Io -> 2Ii2Io + 2Ii1Io
	E(7,12) = E(12,7)
	E(11,8) = 7.33893243358348d+27/temperature*exp(&
			 &((-187.511044d0/temperature-(-287.5307016d0)/1.d3)&
			 &-(-62.03219296d0/temperature-(-92.71110348d0)/1.d3)&
			 &-(-83.27407492d0/temperature-(-142.0468805d0)/1.d3))&
			 &*5.03218937158374d+02)*K(11,8)	! 4Ii3Io -> 1Ii2Io + 3Ii1Io
	E(8,11) = E(11,8)
	E(10,9) = 7.33893243358348d+27/temperature*exp(&
			 &((-187.511044d0/temperature-(-287.5307016d0)/1.d3)&
			 &-(-30.41407318d0/temperature-(-42.69341052d0)/1.d3)&
			 &-(-111.5351121d0/temperature-(-194.8132841d0)/1.d3))&
			 &*5.03218937158374d+02)*K(10,9)	! 4Ii3Io -> 2Io + 4Ii1Io
	E(9,10) = E(10,9)
	E(15,5) = 7.33893243358348d+27/temperature*exp(&
			 &((-88.27882945d0/temperature-(-137.4586806d0)/1.d3)&
			 &-(-60.03371557d0/temperature-(-88.35020888d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(15,5)	! 4Io -> 3Io + 1Io
	E(5,15) = E(15,5)
	E(10,10) = 0.5d0*7.33893243358348d+27/temperature*exp(&
			 &((-88.27882945d0/temperature-(-137.4586806d0)/1.d3)&
			 &-(-30.41407318d0/temperature-(-42.69341052d0)/1.d3)&
			 &-(-30.41407318d0/temperature-(-42.69341052d0)/1.d3))&
			 &*5.03218937158374d+02)*K(10,10)	! 4Io -> 2Io + 2Io
	E(20,1) = 7.33893243358348d+27/temperature*exp(&
			 &((-119.6589794d0/temperature-(-187.4027099d0)/1.d3)&
			 &-(-88.27882945d0/temperature-(-137.4586806d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(20,1)	! 1Ii4Io -> 4Io + 1Ii
	E(1,20) = E(20,1)
	E(16,5) = 7.33893243358348d+27/temperature*exp(&
			 &((-119.6589794d0/temperature-(-187.4027099d0)/1.d3)&
			 &-(-89.23828631d0/temperature-(-139.4518114d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(16,5)	! 1Ii4Io -> 1Ii3Io + 1Io
	E(5,16) = E(16,5)
	E(15,6) = 7.33893243358348d+27/temperature*exp(&
			 &((-119.6589794d0/temperature-(-187.4027099d0)/1.d3)&
			 &-(-60.03371557d0/temperature-(-88.35020888d0)/1.d3)&
			 &-(-30.21939422d0/temperature-(-45.50526148d0)/1.d3))&
			 &*5.03218937158374d+02)*K(15,6)	! 1Ii4Io -> 3Io + 1Ii1Io
	E(6,15) = E(15,6)
	E(11,10) = 7.33893243358348d+27/temperature*exp(&
			 &((-119.6589794d0/temperature-(-187.4027099d0)/1.d3)&
			 &-(-62.03219296d0/temperature-(-92.71110348d0)/1.d3)&
			 &-(-30.41407318d0/temperature-(-42.69341052d0)/1.d3))&
			 &*5.03218937158374d+02)*K(11,10)	! 1Ii4Io -> 1Ii2Io + 2Io
	E(10,11) = E(11,10)
	E(21,1) = 7.33893243358348d+27/temperature*exp(&
			 &((-151.5999442d0/temperature-(-236.1596999d0)/1.d3)&
			 &-(-119.6589794d0/temperature-(-187.4027099d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(21,1)	! 2Ii4Io -> 1Ii4Io + 1Ii
	E(1,21) = E(21,1)
	E(20,2) = 7.33893243358348d+27/temperature*exp(&
			 &((-151.5999442d0/temperature-(-236.1596999d0)/1.d3)&
			 &-(-88.27882945d0/temperature-(-137.4586806d0)/1.d3)&
			 &-(-23.9523596d0/temperature-(-47.68781345d0)/1.d3))&
			 &*5.03218937158374d+02)*K(20,2)	! 2Ii4Io -> 4Io + 2Ii
	E(2,20) = E(20,2)
	E(17,5) = 7.33893243358348d+27/temperature*exp(&
			 &((-151.5999442d0/temperature-(-236.1596999d0)/1.d3)&
			 &-(-120.2575982d0/temperature-(-186.1020181d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(17,5)	! 2Ii4Io -> 2Ii3Io + 1Io
	E(5,17) = E(17,5)
	E(16,6) = 7.33893243358348d+27/temperature*exp(&
			 &((-151.5999442d0/temperature-(-236.1596999d0)/1.d3)&
			 &-(-89.23828631d0/temperature-(-139.4518114d0)/1.d3)&
			 &-(-30.21939422d0/temperature-(-45.50526148d0)/1.d3))&
			 &*5.03218937158374d+02)*K(16,6)	! 2Ii4Io -> 1Ii3Io + 1Ii1Io
	E(6,16) = E(16,6)
	E(15,7) = 7.33893243358348d+27/temperature*exp(&
			 &((-151.5999442d0/temperature-(-236.1596999d0)/1.d3)&
			 &-(-60.03371557d0/temperature-(-88.35020888d0)/1.d3)&
			 &-(-58.58714661d0/temperature-(-95.81971414d0)/1.d3))&
			 &*5.03218937158374d+02)*K(15,7)	! 2Ii4Io -> 3Io + 2Ii1Io
	E(7,15) = E(15,7)
	E(12,10) = 7.33893243358348d+27/temperature*exp(&
			 &((-151.5999442d0/temperature-(-236.1596999d0)/1.d3)&
			 &-(-88.1832749d0/temperature-(-139.0540271d0)/1.d3)&
			 &-(-30.41407318d0/temperature-(-42.69341052d0)/1.d3))&
			 &*5.03218937158374d+02)*K(12,10)	! 2Ii4Io -> 2Ii2Io + 2Io
	E(10,12) = E(12,10)
	E(11,11) = 0.5d0*7.33893243358348d+27/temperature*exp(&
			 &((-151.5999442d0/temperature-(-236.1596999d0)/1.d3)&
			 &-(-62.03219296d0/temperature-(-92.71110348d0)/1.d3)&
			 &-(-62.03219296d0/temperature-(-92.71110348d0)/1.d3))&
			 &*5.03218937158374d+02)*K(11,11)	! 2Ii4Io -> 1Ii2Io + 1Ii2Io
	E(22,1) = 7.33893243358348d+27/temperature*exp(&
			 &((-190.5820764d0/temperature-(-287.4570378d0)/1.d3)&
			 &-(-151.5999442d0/temperature-(-236.1596999d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(22,1)	! 3Ii4Io -> 2Ii4Io + 1Ii
	E(1,22) = E(22,1)
	E(21,2) = 7.33893243358348d+27/temperature*exp(&
			 &((-190.5820764d0/temperature-(-287.4570378d0)/1.d3)&
			 &-(-119.6589794d0/temperature-(-187.4027099d0)/1.d3)&
			 &-(-23.9523596d0/temperature-(-47.68781345d0)/1.d3))&
			 &*5.03218937158374d+02)*K(21,2)	! 3Ii4Io -> 1Ii4Io + 2Ii
	E(2,21) = E(21,2)
	E(20,3) = 7.33893243358348d+27/temperature*exp(&
			 &((-190.5820764d0/temperature-(-287.4570378d0)/1.d3)&
			 &-(-88.27882945d0/temperature-(-137.4586806d0)/1.d3)&
			 &-(-46.43243564d0/temperature-(-94.70212967d0)/1.d3))&
			 &*5.03218937158374d+02)*K(20,3)	! 3Ii4Io -> 4Io + 3Ii
	E(3,20) = E(20,3)
	E(18,5) = 7.33893243358348d+27/temperature*exp(&
			 &((-190.5820764d0/temperature-(-287.4570378d0)/1.d3)&
			 &-(-158.4769381d0/temperature-(-240.2133116d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(18,5)	! 3Ii4Io -> 3Ii3Io + 1Io
	E(5,18) = E(18,5)
	E(17,6) = 7.33893243358348d+27/temperature*exp(&
			 &((-190.5820764d0/temperature-(-287.4570378d0)/1.d3)&
			 &-(-120.2575982d0/temperature-(-186.1020181d0)/1.d3)&
			 &-(-30.21939422d0/temperature-(-45.50526148d0)/1.d3))&
			 &*5.03218937158374d+02)*K(17,6)	! 3Ii4Io -> 2Ii3Io + 1Ii1Io
	E(6,17) = E(17,6)
	E(16,7) = 7.33893243358348d+27/temperature*exp(&
			 &((-190.5820764d0/temperature-(-287.4570378d0)/1.d3)&
			 &-(-89.23828631d0/temperature-(-139.4518114d0)/1.d3)&
			 &-(-58.58714661d0/temperature-(-95.81971414d0)/1.d3))&
			 &*5.03218937158374d+02)*K(16,7)	! 3Ii4Io -> 1Ii3Io + 2Ii1Io
	E(7,16) = E(16,7)
	E(15,8) = 7.33893243358348d+27/temperature*exp(&
			 &((-190.5820764d0/temperature-(-287.4570378d0)/1.d3)&
			 &-(-60.03371557d0/temperature-(-88.35020888d0)/1.d3)&
			 &-(-83.27407492d0/temperature-(-142.0468805d0)/1.d3))&
			 &*5.03218937158374d+02)*K(15,8)	! 3Ii4Io -> 3Io + 3Ii1Io
	E(8,15) = E(15,8)
	E(13,10) = 7.33893243358348d+27/temperature*exp(&
			 &((-190.5820764d0/temperature-(-287.4570378d0)/1.d3)&
			 &-(-119.4487242d0/temperature-(-191.2458533d0)/1.d3)&
			 &-(-30.41407318d0/temperature-(-42.69341052d0)/1.d3))&
			 &*5.03218937158374d+02)*K(13,10)	! 3Ii4Io -> 3Ii2Io + 2Io
	E(10,13) = E(13,10)
	E(12,11) = 7.33893243358348d+27/temperature*exp(&
			 &((-190.5820764d0/temperature-(-287.4570378d0)/1.d3)&
			 &-(-88.1832749d0/temperature-(-139.0540271d0)/1.d3)&
			 &-(-62.03219296d0/temperature-(-92.71110348d0)/1.d3))&
			 &*5.03218937158374d+02)*K(12,11)	! 3Ii4Io -> 2Ii2Io + 1Ii2Io
	E(11,12) = E(12,11)
	E(23,1) = 7.33893243358348d+27/temperature*exp(&
			 &((-224.1158831d0/temperature-(-337.0895745d0)/1.d3)&
			 &-(-190.5820764d0/temperature-(-287.4570378d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(23,1)	! 4Ii4Io -> 3Ii4Io + 1Ii
	E(1,23) = E(23,1)
	E(22,2) = 7.33893243358348d+27/temperature*exp(&
			 &((-224.1158831d0/temperature-(-337.0895745d0)/1.d3)&
			 &-(-151.5999442d0/temperature-(-236.1596999d0)/1.d3)&
			 &-(-23.9523596d0/temperature-(-47.68781345d0)/1.d3))&
			 &*5.03218937158374d+02)*K(22,2)	! 4Ii4Io -> 2Ii4Io + 2Ii
	E(2,22) = E(22,2)
	E(21,3) = 7.33893243358348d+27/temperature*exp(&
			 &((-224.1158831d0/temperature-(-337.0895745d0)/1.d3)&
			 &-(-119.6589794d0/temperature-(-187.4027099d0)/1.d3)&
			 &-(-46.43243564d0/temperature-(-94.70212967d0)/1.d3))&
			 &*5.03218937158374d+02)*K(21,3)	! 4Ii4Io -> 1Ii4Io + 3Ii
	E(3,21) = E(21,3)
	E(20,4) = 7.33893243358348d+27/temperature*exp(&
			 &((-224.1158831d0/temperature-(-337.0895745d0)/1.d3)&
			 &-(-88.27882945d0/temperature-(-137.4586806d0)/1.d3)&
			 &-(-77.77838269d0/temperature-(-144.4714705d0)/1.d3))&
			 &*5.03218937158374d+02)*K(20,4)	! 4Ii4Io -> 4Io + 4Ii
	E(4,20) = E(20,4)
	E(19,5) = 7.33893243358348d+27/temperature*exp(&
			 &((-224.1158831d0/temperature-(-337.0895745d0)/1.d3)&
			 &-(-187.511044d0/temperature-(-287.5307016d0)/1.d3)&
			 &-0.d0)&
			 &*5.03218937158374d+02)*K(19,5)	! 4Ii4Io -> 4Ii3Io + 1Io
	E(5,19) = E(19,5)
	E(18,6) = 7.33893243358348d+27/temperature*exp(&
			 &((-224.1158831d0/temperature-(-337.0895745d0)/1.d3)&
			 &-(-158.4769381d0/temperature-(-240.2133116d0)/1.d3)&
			 &-(-30.21939422d0/temperature-(-45.50526148d0)/1.d3))&
			 &*5.03218937158374d+02)*K(18,6)	! 4Ii4Io -> 3Ii3Io + 1Ii1Io
	E(6,18) = E(18,6)
	E(17,7) = 7.33893243358348d+27/temperature*exp(&
			 &((-224.1158831d0/temperature-(-337.0895745d0)/1.d3)&
			 &-(-120.2575982d0/temperature-(-186.1020181d0)/1.d3)&
			 &-(-58.58714661d0/temperature-(-95.81971414d0)/1.d3))&
			 &*5.03218937158374d+02)*K(17,7)	! 4Ii4Io -> 2Ii3Io + 2Ii1Io
	E(7,17) = E(17,7)
	E(16,8) = 7.33893243358348d+27/temperature*exp(&
			 &((-224.1158831d0/temperature-(-337.0895745d0)/1.d3)&
			 &-(-89.23828631d0/temperature-(-139.4518114d0)/1.d3)&
			 &-(-83.27407492d0/temperature-(-142.0468805d0)/1.d3))&
			 &*5.03218937158374d+02)*K(16,8)	! 4Ii4Io -> 1Ii3Io + 3Ii1Io
	E(8,16) = E(16,8)
	E(15,9) = 7.33893243358348d+27/temperature*exp(&
			 &((-224.1158831d0/temperature-(-337.0895745d0)/1.d3)&
			 &-(-60.03371557d0/temperature-(-88.35020888d0)/1.d3)&
			 &-(-111.5351121d0/temperature-(-194.8132841d0)/1.d3))&
			 &*5.03218937158374d+02)*K(15,9)	! 4Ii4Io -> 3Io + 4Ii1Io
	E(9,15) = E(15,9)
	E(14,10) = 7.33893243358348d+27/temperature*exp(&
			 &((-224.1158831d0/temperature-(-337.0895745d0)/1.d3)&
			 &-(-148.7720382d0/temperature-(-241.7960301d0)/1.d3)&
			 &-(-30.41407318d0/temperature-(-42.69341052d0)/1.d3))&
			 &*5.03218937158374d+02)*K(14,10)	! 4Ii4Io -> 4Ii2Io + 2Io
	E(10,14) = E(14,10)
	E(13,11) = 7.33893243358348d+27/temperature*exp(&
			 &((-224.1158831d0/temperature-(-337.0895745d0)/1.d3)&
			 &-(-119.4487242d0/temperature-(-191.2458533d0)/1.d3)&
			 &-(-62.03219296d0/temperature-(-92.71110348d0)/1.d3))&
			 &*5.03218937158374d+02)*K(13,11)	! 4Ii4Io -> 3Ii2Io + 1Ii2Io
	E(11,13) = E(13,11)
	E(12,12) = 0.5d0*7.33893243358348d+27/temperature*exp(&
			 &((-224.1158831d0/temperature-(-337.0895745d0)/1.d3)&
			 &-(-88.1832749d0/temperature-(-139.0540271d0)/1.d3)&
			 &-(-88.1832749d0/temperature-(-139.0540271d0)/1.d3))&
			 &*5.03218937158374d+02)*K(12,12)	! 4Ii4Io -> 2Ii2Io + 2Ii2Io

end subroutine get_evap_3
