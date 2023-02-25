! This code is the fussion of a lot of codes that will be used to carry out
! an SCF of a HHe^+ diatomic.
!

program SCF
  implicit none
  ! +++++++++++++++++++++++++++++ Interface +++++++++++++++++++++++++++++++++
  interface

    subroutine jacobi_diagonalization(incoming_matrix, m1, diagonalized_mat, &
    total_jcbi)
      real(kind=8), dimension(:,:), intent(in) :: incoming_matrix
      real(kind=8), dimension(:,:), intent(out) :: diagonalized_mat, &
      total_jcbi
      integer, intent(in) :: m1
    end subroutine jacobi_diagonalization

    subroutine inverse_sqrt(diagonal_overlap, m1, inverse_diagonal_overlap)
      real(kind=8), dimension(:,:), intent(in) :: diagonal_overlap
      real(kind=8), dimension(:,:), intent(out) :: inverse_diagonal_overlap
      integer, intent(in) :: m1
    end subroutine inverse_sqrt

    subroutine similarity_transformation(incoming_matrix, transforming_mat, &
    m1, transformed_mat)
      real(kind=8), dimension(:,:), intent(in) :: incoming_matrix, &
      transforming_mat
      real(kind=8), dimension(:,:), intent(out) :: transformed_mat
      integer, intent(in) :: m1
    end subroutine similarity_transformation

    subroutine set_up_calculation(full_basis, system_coord, nuclear_repulsion)
      real(kind=8), dimension(3,2,2), intent(out) :: full_basis
      real(kind=8), dimension(3,2), intent(out) :: system_coord
      real(kind=8), intent(out) :: nuclear_repulsion
    end subroutine set_up_calculation

    subroutine calculate_overlap_matrix(basis_info, coord_array, ovrlpmat)
      real(kind=8), dimension(:,:,:), intent(in) :: basis_info
      real(kind=8), dimension(:,:), intent(in) :: coord_array
      real(kind=8), dimension(:,:), intent(out) :: ovrlpmat
    end subroutine calculate_overlap_matrix

    subroutine overlap_integral(expA, expB, coeffA, coeffB, basis_distance, &
      intgrl)
      real(kind=8), intent(in) :: expA, expB, coeffA, coeffB, basis_distance
      real(kind=8), intent(out) :: intgrl
    end subroutine overlap_integral

    subroutine calculate_T_matrix(basis_info, coord_array, kin_engy_mat)
      real(kind=8), dimension(:,:,:), intent(in) :: basis_info
      real(kind=8), dimension(:,:), intent(in) :: coord_array
      real(kind=8), dimension(:,:), intent(out) :: kin_engy_mat
    end subroutine calculate_T_matrix

    subroutine kinetic_integral(expA, expB, coeffA, coeffB, basis_distance, &
    intgrl)
      real(kind=8), intent(in) :: expA, expB, coeffA, coeffB, basis_distance
      real(kind=8), intent(out) :: intgrl
    end subroutine kinetic_integral

    subroutine attr_integral(expA, expB, coeffA, coeffB, basis_distance, &
    nuclear_coord, pointP_coord, nuclear_charge, intgrl)
      real(kind=8), dimension(3), intent(in) :: nuclear_coord, pointP_coord
      real(kind=8), intent(in) :: expA, expB, coeffA, coeffB, basis_distance, &
      nuclear_charge
      real(kind=8), intent(out) :: intgrl
    end subroutine attr_integral

    subroutine calculate_V_matrices(basis_info, coord_array, V1_matrix, &
    V2_matrix, Vtot_matrix)
      real(kind=8), dimension(:,:,:), intent(in) :: basis_info
      real(kind=8), dimension(:,:), intent(in) :: coord_array
      real(kind=8), dimension(:,:), intent(out) :: V1_matrix, V2_matrix, &
      Vtot_matrix
    end subroutine calculate_V_matrices

    subroutine calculate_two_elec_ints(basis_info, coord_array, two_e_ints)
      real(kind=8), dimension(:,:,:), intent(in) :: basis_info
      real(kind=8), dimension(:,:), intent(in) :: coord_array
      real(kind=8), dimension(:,:,:,:), intent(out) :: two_e_ints
    end subroutine calculate_two_elec_ints

    subroutine repulsion_integral(expA, expB, expC, expD, coeffA, coeffB, &
    coeffC, coeffD, ab_dist, cd_dist, pq_dist, intgrl)
      real(kind=8), intent(in) :: expA, expB, expC, expD, coeffA, coeffB, &
      coeffC, coeffD, ab_dist, cd_dist, pq_dist
      real(kind=8), intent(out) :: intgrl
    end subroutine repulsion_integral

    subroutine generate_G_matrix(two_e_ints, density_mat, gmat)
      real(kind=8), dimension(:,:,:,:), intent(in) :: two_e_ints
      real(kind=8), dimension(:,:), intent(in) :: density_mat
      real(kind=8), dimension(:,:), intent(out) :: gmat
    end subroutine generate_G_matrix

    subroutine generate_Hcore_matrix(kin_engy_mat, Vtot_matrix, H_core_matrix)
      real(kind=8), dimension(:,:), intent(in) :: kin_engy_mat, Vtot_matrix
      real(kind=8), dimension(:,:), intent(out) :: H_core_matrix
    end subroutine generate_Hcore_matrix

    subroutine generate_Fock_matrix(H_core_matrix, gmat, Fock_mat)
      real(kind=8), dimension(:,:), intent(in) :: H_core_matrix, gmat
      real(kind=8), dimension(:,:), intent(out) :: Fock_mat
    end subroutine generate_Fock_matrix

    subroutine eigenvalue_vector_ordering(diagonal_mat, eigenvector_mat)
      real(kind=8), dimension(:,:), intent(inout) :: diagonal_mat, &
      eigenvector_mat
    end subroutine eigenvalue_vector_ordering

    subroutine matrix_multiplication(matrixA, matrixB, resulting_matrix)
      real(kind=8), dimension(:,:), intent(in) :: matrixA, matrixB
      real(kind=8), dimension(:,:), intent(out) :: resulting_matrix
    end subroutine  matrix_multiplication

    subroutine generate_new_density_mat(coeff_mat, new_density_mat)
      real(kind=8), dimension(:,:), intent(in) :: coeff_mat
      real(kind=8), dimension(:,:), intent(out) :: new_density_mat
    end subroutine generate_new_density_mat

    subroutine determine_convergence(old_density_mat, new_density_mat, &
    convergence)
      real(kind=8), dimension(:,:), intent(in) :: old_density_mat, &
      new_density_mat
      logical, intent(out) :: convergence
    end subroutine determine_convergence

    subroutine calculate_electronic_energy(density_mat, Hcore_mat, fock_mat, &
    total_E)
      real(kind=8), dimension(:,:), intent(in) :: density_mat, Hcore_mat, &
      fock_mat
      real(kind=8), intent(out) :: total_E
    end subroutine calculate_electronic_energy

  end interface
  ! +++++++++++++++++++++++++ End of Interface ++++++++++++++++++++++++++++++
  !
  !
  !
  ! +++++++++++++++++++++++++++ Actual Program ++++++++++++++++++++++++++++++
  !
  ! --------------------- Declare stuff used in the main program ------------
  real(kind=8), dimension(:,:,:,:), allocatable :: two_e_intgrls
  real(kind=8), dimension(:,:,:), allocatable :: basis_set
  real(kind=8), dimension(:,:), allocatable :: s, diagonal_s, inv_sqrt_diag_s, &
  u_diagonalizer, final_transformed_s, coordinates, T_matrix, V_1_Matrix, &
  V_2_Matrix, V_tot_Matrix, G_matrix, density_matrix, Hcore_matrix, &
  Fock_matrix, transformed_Fock_matrix, diagonalized_Fock_mat, &
  diagonal_trans_Fock_matrix, c_prime_mat, c_matrix, old_density_matrix
  real(kind=8) :: total_energy, nuclear_rep
  integer :: m, x, y, i, j, iteration_counter
  logical :: is_converged = .false.


  ! ---------- Allocate arrays containing basis set and coordinates ---------

  allocate(basis_set(3,2,2))
  allocate(coordinates(3,2))

  ! ----------- Call subroutine to get HHe^+ system's information -----------

  call set_up_calculation(basis_set, coordinates, nuclear_rep)

    ! ********************** A tool for debugging ***********************
    ! write(*,*)"The H basis information is: "
    ! do i = lbound(basis_set,1), ubound(basis_set,1)
    !   write(*,*) (basis_set(i, j, 1), j=lbound(basis_set,2), ubound(basis_set,2))
    ! end do
    ! write(*,*)" "
    ! write(*,*)"The He basis information is: "
    ! do i = lbound(basis_set,1), ubound(basis_set,1)
    !   write(*,*) (basis_set(i, j, 2), j=lbound(basis_set,2), ubound(basis_set,2))
    ! end do
    ! write(*,*)" "
    ! write(*,*)"The systems's coordinates are:"
    ! write(*,*)"H coordinates:"
    ! do i = 1, 3
    !   write(*,*) coordinates(i,1)
    ! end do
    ! write(*,*)"He coordinates"
    ! do i = 1, 3
    !   write(*,*) coordinates(i, 2)
    ! end do
    ! *******************************************************************

  ! -------------- Allocate a bunch of matrices (S, T, V, etc.) -------------

    allocate(s(2,2))
    allocate(diagonal_s(2,2))
    allocate(inv_sqrt_diag_s(2,2))
    allocate(u_diagonalizer(2,2))
    allocate(final_transformed_s(2,2))
    allocate(T_matrix(2,2))
    allocate(V_1_Matrix(2,2))
    allocate(V_2_Matrix(2,2))
    allocate(V_tot_Matrix(2,2))
    allocate(two_e_intgrls(2,2,2,2))
    allocate(density_matrix(2,2))
    allocate(G_matrix(2,2))
    allocate(Hcore_matrix(2,2))
    allocate(Fock_matrix(2,2))
    allocate(transformed_Fock_matrix(2,2))
    allocate(diagonalized_Fock_mat(2,2))
    allocate(diagonal_trans_Fock_matrix(2,2))
    allocate(c_prime_mat(2,2))
    allocate(c_matrix(2,2))
    allocate(old_density_matrix(2,2))

    m = 2 !<--- Useful to allocate intermediates in matrix manipulation
    !                                                       subroutines

  ! -------------------- Initialize iteration counter --------------------- 

  iteration_counter = 0   !<--- counter to keep track of iterations

  ! -- Call subroutine to calculate kinetic energy integrals, get T matrix --

  call calculate_T_matrix(basis_set, coordinates, T_matrix)

    ! ********************** A tool for debugging ***********************
    ! write(*,*)"The T matrix is"
    ! do i = lbound(T_matrix, 1), ubound(T_matrix,1)
    !  write(*,*) (T_matrix(i, j), j=lbound(T_matrix,2), ubound(T_matrix,2))
    ! end do
    ! write(*,*)" "
    ! *******************************************************************


  ! --- Call subroutine calculate electron-nuclear integrals, get V matrix --

  call calculate_V_matrices(basis_set, coordinates, V_1_Matrix, V_2_Matrix, &
  V_tot_Matrix)

    ! ********************** A tool for debugging ***********************
    ! write(*,*) "The V1 Matrix is:"
    ! write(*,*) " "
    ! do i = lbound(V_1_Matrix,1), ubound(V_1_Matrix,1)
    !   write(*,*) (V_1_Matrix(i, j), j=lbound(V_1_Matrix,2), &
    !   ubound(V_1_Matrix,2))
    ! end do
    ! write(*,*) " "
    ! write(*,*) "The V2 Matrix is:"
    ! write(*,*) " "
    ! do i = lbound(V_2_Matrix,1), ubound(V_2_Matrix,1)
    !   write(*,*) (V_2_Matrix(i, j), j=lbound(V_2_Matrix,2), &
    !   ubound(V_2_Matrix,2))
    ! end do
    ! write(*,*) " "
    ! write(*,*) "The Vtotal Matrix is:"
    ! write(*,*) " "
    ! do i = lbound(V_tot_Matrix,1), ubound(V_tot_Matrix,1)
    !   write(*,*) (V_tot_Matrix(i, j), j=lbound(V_tot_Matrix,2), &
    !   ubound(V_tot_Matrix,2))
    ! end do
    ! *******************************************************************


  ! --------------- Call subroutine to generate H_core matrix ---------------

  call generate_Hcore_matrix(T_matrix, V_tot_Matrix, Hcore_matrix)

    ! ********************** A tool for debugging ***********************
    ! write(*,*) "The H_core Matrix is:"
    ! write(*,*) " "
    ! do i = lbound(Hcore_matrix,1), ubound(Hcore_matrix,1)
    !   write(*,*) (Hcore_matrix(i, j), j=lbound(Hcore_matrix,2), &
    !   ubound(Hcore_matrix,2))
    ! end do
    ! *******************************************************************


  ! -------- Call subroutine to calculate all 2 electron integrals ----------

  call calculate_two_elec_ints(basis_set, coordinates, two_e_intgrls)

    ! ********************** A tool for debugging ***********************
    ! write(*,*)"The following integrals should be te same: "
    ! write(*,*)" "
    ! write(*,*) two_e_intgrls(1,2,1,2) ! 1,2,3,4
    ! write(*,*)" "
    ! write(*,*) two_e_intgrls(1,2,1,2) ! 3,2,1,4
    ! write(*,*)" "
    ! write(*,*) two_e_intgrls(1,2,1,2) ! 1,4,3,2
    ! write(*,*)" "
    ! write(*,*) two_e_intgrls(1,2,1,2) ! 3,4,1,2
    ! write(*,*)" "
    ! write(*,*) two_e_intgrls(2,1,2,1) ! 2,1,4,3
    ! write(*,*)" "
    ! write(*,*) two_e_intgrls(2,1,2,1) ! 4,1,2,3
    ! write(*,*)" "
    ! write(*,*) two_e_intgrls(2,1,2,1) ! 2,3,4,1
    ! write(*,*)" "
    ! write(*,*) two_e_intgrls(2,1,2,1) ! 4,3,2,1
    ! *******************************************************************


  ! ------------- Call subroutine to calculate overlap matrix ---------------

  call calculate_overlap_matrix(basis_set, coordinates, s)

    ! ********************** A tool for debugging ***********************
    ! write(*,*)"The overlap matrix is"
    ! do i = lbound(s, 1), ubound(s,1)
    !  write(*,*) (s(i, j), j=lbound(s,2), ubound(s,2))
    ! end do
    ! write(*,*)" "
    ! *******************************************************************


  ! ------------------ Call subroutine for diagonalization ------------------

  call jacobi_diagonalization(s, m, diagonal_s, u_diagonalizer)

    ! ********************** A tool for debugging ***********************
    ! WRITE(*,*)"The diagonalization gives:"
    ! do i = lbound(diagonal_s,1), ubound(diagonal_s,1)
    !  WRITE(*,*) (diagonal_s(i,y), y = lbound(diagonal_s,1), ubound(diagonal_s,1))
    ! end do
    ! WRITE(*,*)" "
    ! *******************************************************************

    ! ********************** A tool for debugging ***********************
    ! WRITE(*,*)"Final eigenvector matrix" ! Print out final Jtot matrix
    ! do i = lbound(u_diagonalizer,1), ubound(u_diagonalizer,1)
    !  WRITE(*,*) (u_diagonalizer(i,y), y = lbound(u_diagonalizer,1), &
    !  &ubound(u_diagonalizer,1))
    ! end do
    ! WRITE(*,*)" "
    ! *******************************************************************


  ! --- Call subroutine to take inverse square roots of diagonal elements ---

  call inverse_sqrt(diagonal_s, m, inv_sqrt_diag_s)

    ! ********************** A tool for debugging ***********************
    ! WRITE(*,*)"Inverse square root values:"
    ! do i = lbound(inv_sqrt_diag_s,1), ubound(inv_sqrt_diag_s,1)
    !  WRITE(*,*) (inv_sqrt_diag_s(i,y), y = lbound(inv_sqrt_diag_s,1), &
    !  &ubound(diagonal_s,1))
    ! end do
    ! WRITE(*,*)" "
    ! *******************************************************************


  ! ----------- Call subroutine to perform U·s^(-1/2)·U^t -------------------

  call similarity_transformation(inv_sqrt_diag_s, u_diagonalizer, m, &
  final_transformed_s)

    ! ********************** A tool for debugging ***********************
    ! WRITE(*,*)"Final S to the -1/2:"
    ! do i = lbound(final_transformed_s,1), ubound(final_transformed_s,1)
    !   WRITE(*,*) (final_transformed_s(i,y), y = lbound(final_transformed_s,1), &
    !   &ubound(final_transformed_s,1))
    ! end do
    ! WRITE(*,*)" "
    ! *******************************************************************


  ! -- Set initial density matrix to the null matrix (so that F = H_core) ---

  density_matrix = 0.0

  ! ------------------------- Perform iterations ----------------------------

  do while (.not.is_converged)

    old_density_matrix = density_matrix


    ! --------------- Call subroutine to generate the G matrix ----------------

    call generate_G_matrix(two_e_intgrls, density_matrix, G_matrix)

      ! ********************** A tool for debugging ***********************
      ! write(*,*)"The G matrix is: "
      ! write(*,*)" "
      ! do i = lbound(G_matrix, 1), ubound(G_matrix, 1)
      !   write(*,*) (G_matrix(i, j), j=lbound(G_matrix,2), &
      !   ubound(G_matrix,2))
      ! end do
      ! *******************************************************************


    ! ---------------- Call subroutine to generate Fock matrix ----------------

    call generate_Fock_matrix(Hcore_matrix, G_matrix, Fock_matrix)

      ! ********************** A tool for debugging ***********************
      ! write(*,*)"The Fock matrix is: "
      ! write(*,*)" "
      ! do i = lbound(Fock_matrix, 1), ubound(Fock_matrix, 1)
      !   write(*,*) (Fock_matrix(i, j), j=lbound(Fock_matrix,2), &
      !   ubound(Fock_matrix,2))
      ! end do
      ! write(*,*)" "
      ! *******************************************************************


    ! ------------ Call subroutine to calculate electronic energy -------------

    call calculate_electronic_energy(density_matrix, Hcore_matrix, &
    Fock_matrix, total_energy)

    write(*,*)"---------- Electronic energy:", total_energy, "------------------"
    write(*,*)" "
    write(*,*)" "
    write(*,*)"********* Iteration:", iteration_counter, "*********"

    ! ---------------- Call subroutine to transform Fock matrix ---------------

    call similarity_transformation(Fock_matrix, final_transformed_s, m, &
    transformed_Fock_matrix)

      ! ********************** A tool for debugging ***********************
      ! write(*,*)"The transformed Fock matrix is: "
      ! write(*,*)" "
      ! do i = lbound(transformed_Fock_matrix, 1), ubound(transformed_Fock_matrix, 1)
      !   write(*,*) (transformed_Fock_matrix(i, j), &
      !   j=lbound(transformed_Fock_matrix,2), &
      !   ubound(transformed_Fock_matrix,2))
      ! end do
      ! write(*,*)" "
      ! *******************************************************************


    ! ------ Call subroutine to diagonalize the transformed Fock matrix -------

    call jacobi_diagonalization(transformed_Fock_matrix, m, &
    diagonal_trans_Fock_matrix, c_prime_mat)

      ! ********************** A tool for debugging ***********************
      ! write(*,*)"The diagonalized transformed Fock matrix is: "
      ! write(*,*)" "
      ! do i = lbound(diagonal_trans_Fock_matrix, 1), &
      !   ubound(diagonal_trans_Fock_matrix, 1)
      !   write(*,*) (diagonal_trans_Fock_matrix(i, j), &
      !   j=lbound(diagonal_trans_Fock_matrix,2), &
      !   ubound(diagonal_trans_Fock_matrix,2))
      ! end do
      ! write(*,*)" "
      ! write(*,*)"The unordered ordered C' matrix is: "
      ! write(*,*)" "
      ! do i = lbound(c_prime_mat,1), ubound(c_prime_mat,1)
      !   write(*,*) (c_prime_mat(i,j), j=lbound(c_prime_mat,2), &
      !   ubound(c_prime_mat,2))
      ! end do
      ! write(*,*)" "
      ! *******************************************************************


    ! - Call subroutine to order eigenvalues/vectors in diagonal Fock Matrix --

    call eigenvalue_vector_ordering(diagonal_trans_Fock_matrix, c_prime_mat)

      ! ********************** A tool for debugging ***********************
      ! write(*,*)"The ordered transformed Fock matrix is: "
      ! write(*,*)" "
      ! do i = lbound(diagonal_trans_Fock_matrix, 1), &
      !   ubound(diagonal_trans_Fock_matrix, 1)
      !   write(*,*) (diagonal_trans_Fock_matrix(i, j), &
      !   j=lbound(diagonal_trans_Fock_matrix,2), &
      !   ubound(diagonal_trans_Fock_matrix,2))
      ! end do
      ! write(*,*)" "
      ! write(*,*)"The ordered C' matrix is: "
      ! write(*,*)" "
      ! do i = lbound(c_prime_mat,1), ubound(c_prime_mat,1)
      !   write(*,*) (c_prime_mat(i,j), j=lbound(c_prime_mat,2), &
      !   ubound(c_prime_mat,2))
      ! end do
      ! write(*,*)" "
      ! *******************************************************************


    ! -------------- Call subroutine to obtain C from S^-1/2*C'----------------

    call matrix_multiplication(final_transformed_s, c_prime_mat, c_matrix)

      ! ********************** A tool for debugging ***********************
      ! write(*,*)"The C matrix is: "
      ! write(*,*)" "
      ! do i = lbound(c_matrix,1), ubound(c_matrix,1)
      !   write(*,*) (c_matrix(i,j), j=lbound(c_matrix,2), &
      !   ubound(c_matrix,2))
      ! end do
      ! write(*,*)" "
      ! *******************************************************************


    ! -------- Call subroutine to obtain new density matrix from C ------------

    call generate_new_density_mat(c_matrix, density_matrix)

      ! ********************** A tool for debugging ***********************
      write(*,*)"The new density matrix is: "
      write(*,*)" "
      do i = lbound(density_matrix,1), ubound(density_matrix,1)
        write(*,*) (density_matrix(i,j), j=lbound(density_matrix,2), &
        ubound(density_matrix,2))
      end do
      write(*,*)" "
      ! *******************************************************************

    ! -------------- Call subroutine to evaluate convergence ------------------

    call determine_convergence(old_density_matrix, density_matrix, is_converged)

      ! ********************** A tool for debugging ***********************
      ! write(*,*)"Is converged: "
      ! write(*,*) is_converged
      ! *******************************************************************

    iteration_counter = iteration_counter + 1

  end do

  write(*,*) " "
  write(*,*) "************** Density converged **************"
  write(*,*) " "
  write(*,*) "+++++++++++++++++++++++++++++++++"
  write(*,*) "+          Final data           +"
  write(*,*) "+++++++++++++++++++++++++++++++++"
  write(*,*) " "
  write(*,*)"Orbitals:"
  write(*,*)" "
  do i = lbound(c_matrix,1), ubound(c_matrix,1)
    write(*,*) (c_matrix(i,j), j=lbound(c_matrix,2), &
    ubound(c_matrix,2))
  end do
  write(*,*)" "
  write(*,*) "Orbital energies:"
  write(*,*)" "
  do i = lbound(diagonal_trans_Fock_matrix, 1), &
    ubound(diagonal_trans_Fock_matrix, 1)
    write(*,*) (diagonal_trans_Fock_matrix(i, j), &
    j=lbound(diagonal_trans_Fock_matrix,2), &
    ubound(diagonal_trans_Fock_matrix,2))
  end do
  write(*,*)" "
  write(*,*)"Final electronic energy:"
  write(*,*) total_energy
  write(*,*)" "
  write(*,*)"Final system energy:"
  write(*,*) total_energy + nuclear_rep

end program SCF

!
!
!            It's not a phase mom, this is who I really am
!
!

subroutine jacobi_diagonalization(incoming_matrix, m1, diagonalized_mat, &
total_jcbi)
  ! ------------------- Declare stuff coming in and out ---------------------
  real(kind=8), dimension(:,:), intent(in) :: incoming_matrix
  real(kind=8), dimension(:,:), intent(out) :: diagonalized_mat, &
  total_jcbi
  integer, intent(in) :: m1
  !
  !
  ! ------------------ Declare stuff used in the subroutine -----------------
  real(kind=8), dimension(:,:), allocatable :: J, Jt, C, E, Jtot, Temp, A
  real(kind=8), dimension(:), allocatable :: VA, VE, DV, EIG, Tem1
  integer :: i, j2, k, l, x, y, z, deltacount
  real(kind=8) :: pp, pq, qp, qq, phi, jpp, jpq, jqp, product, delta, aux1


  ! ---------------- Allocate all the intermediates to use ------------------
  allocate (A(m1,m1))
  allocate (J(m1,m1))
  allocate (Jt(m1,m1))
  allocate (C(m1,m1))
  allocate (E(m1,m1))
  allocate (Jtot(m1,m1))
  allocate (Temp(m1,m1))
  allocate (VA(m1))
  allocate (VE(m1))
  allocate (DV(m1))
  allocate (EIG(m1))
  allocate (Tem1(m1))

  A = incoming_matrix         !<- Put incoming matrix into A

  ! -------------------- Begin diagonalization process -----------------------

  !                   Initialize Jtot (Eigenvector matrix) as an identity matrix
  !                   so that the first time it is multiplied it does not alter
  !                   the incoming matrix

  do x = lbound(Jtot,1), ubound(Jtot,1)   !<- Row counter
    do y = lbound(Jtot,2), ubound(Jtot,2) !<- Column counter
      if (x == y) then                    !<- Elements in the diagonal are 1
        Jtot(x,y) = 1
      else                                !<- Elements anywhere else are 0
        Jtot(x,y) = 0
      end if
    end do
  end do

  deltacount = 0 ! Counter for convergence

  !                When there is virutally no change in the elemtns of
  !                the rotated matrix, and the previous one in 5 consecutive
  !                loops, the code will stop and give the resulting matrix back

  do while (deltacount < 5)

  ! ------------------ Begin sweep thorugh the matrix ------------------------
  !
  ! Since we only want to zero the off-diagonal elements of matrix A, counter k
  ! will go from i+1 to the end, ensuring that no element in the diagonal will
  ! be zeroed.

    do i = lbound(A,1), ubound(A,1) - 1   !<- Row counter - p
      do k = i+1, ubound(A,2)             !<- Column counter - q

        ! ********************** A tool for debugging ***********************
        ! WRITE(*,*)"Matrix A to use now"
        ! do x = lbound(A,1), ubound(A,1)
        !   WRITE(*,*) (A(x,y), y = lbound(A,2), ubound(A,2))
        ! end do
        ! WRITE(*,*)" "
        ! *******************************************************************

        pp = A(i, i) !<- Element pp
        pq = A(i, k) !<- Element pq
        qp = A(k, i) !<- Element qp
        qq = A(k, k) !<- Element qq

        ! ********************** A tool for debugging ***********************
        ! WRITE(*,*)"Coefficients to use"
        ! WRITE(*,*) pp, pq
        ! WRITE(*,*) qp, qq
        ! WRITE(*,*)" "
        ! *******************************************************************

        !                                     We have extracted pp, pq, qp, qq
        !                                     Now calculate phi
        phi = 0                              !<- Clean phi
        phi = (0.5) * (atan((2*pq)/(qq-pp))) ! Compute phi

        ! ********************** A tool for debugging ***********************
        ! WRITE(*,*)"Phi to use", phi
        ! WRITE(*,*)" "
        ! *******************************************************************

        !                                      Now create J
        !                                      J is an indentity matrix first

        do x = lbound(J,1), ubound(J,1)      !<- Row counter
          do y = lbound(J,2), ubound(J,2)    !<- Column counter
            if (x == y) then                 !<- Elements in the diagonal are 1
              J(x,y) = 1
            else                             ! Elements anywhere else are 0
              J(x,y) = 0
            end if
          end do
        end do

        !                                Now plug in the sines and cosines in J

        J(i,i) = cos(phi)       !<- Same as pp
        J(i,k) = sin(phi)       !<- pq
        J(k,i) = -1 * sin(phi)  !<- qp
        J(k,k) = cos(phi)       !<- and qq

        ! ********************** A tool for debugging ***********************
        ! WRITE(*,*)"Matrix J"
        ! do x = lbound(J,1), ubound(J,1)
        !    WRITE(*,*) (J(x,y), y = lbound(J,2), ubound(J,2))
        ! end do
        !WRITE(*,*)"One do"
        ! *******************************************************************


        !                                         Begin Jtot * J multiplication
        product = 0
        do z = lbound(Jtot,1), ubound(Jtot,1)     !<- This counter is the "slowest", it will allow us to move the rows in Jtot and Temp
         do x = lbound(J,2), ubound(J,2)          !<- This counter sets the column in matrix J
           do y = lbound(Jtot,2), ubound(Jtot,2)  !<- This counter is the "fastest", it sets the column in Jtot and row in J
             product = Jtot(z,y)*J(y,x) + product !<- Multiply each element and add
           end do
           Temp(z,x) = product                    !<- Set the result in C
           product = 0                            !<- Clear variable value
         end do
        end do

        !                              We now have matrix Jtot(old)*J in Temp,
        !                              now set Jtot(new) = temp
        !                              for next iteration

        Jtot = Temp

        !                          We now have matrix Jtot(old)*J in Temp,
        !                          now set Jtot(new) = temp, for next iteration
        !                          Jtot = Temp

        !                          End Jtot * J multiplication

        !                          Generate J (transpose)

        do x = lbound(J,1), ubound(J,1)
          do y = lbound(J,1), ubound(J,1)
            Jt(y,x) = J(x,y)
          end do
        end do

        ! ********************** A tool for debugging ***********************
        ! WRITE(*,*)"Matrix Jt"
        ! do x = lbound(Jt,1), ubound(Jt,1)
        !    WRITE(*,*) (Jt(x,y), y = lbound(Jt,2), ubound(Jt,2))
        ! end do
        ! *******************************************************************

        product = 0                   !<- Clear variable

        !                         Obtain matrix C, which is the product of A·J

        do z = lbound(A,1), ubound(A,1)    !<- This counter is the "slowest", it will allow us to move the rows in A and C
         do x = lbound(J,2), ubound(J,2)   !<- This counter sets the column in matrix B
           do y = lbound(A,2), ubound(A,2) !<- This counter is the "fastest", it sets the column in A and row in B
             product = A(z,y)*J(y,x) + product !<- Multiply each element and add
           end do
           C(z,x) = product                !<- Set the result in C
           product = 0                     !<- Clear variable value
         end do
        end do
        !                            Once we have C, multiply it by J (tranpose)

        ! ********************** A tool for debugging ***********************
        ! WRITE(*,*)"A·J"
        ! do x = lbound(C,1), ubound(C,1)
        !   WRITE(*,*) (C(x,y), y = lbound(C,2), ubound(C,2))
        ! end do
        ! WRITE(*,*)" "
        ! *******************************************************************

        product = 0         !<- Clear variable

        !        Perform J·C to complete the rotation/similarity transformation

        do z = lbound(Jt,1), ubound(Jt,1)       !<- This counter is the "slowest", it will allow us to move the rows in Pt and B
         do x = lbound(C,2), ubound(C,2)        !<- This counter sets the column in matrix C
           do y = lbound(Jt,2), ubound(Jt,2)    !<- This counter is the "fastest", it sets the column in Pt and row in C
             product = Jt(z,y)*C(y,x) + product !<- Multiply each element and add
           end do
           E(z,x) = product                     !<- Set the result in C
           product = 0                          !<- Clear variable value
         end do
        end do

        ! ********************** A tool for debugging ***********************
        ! WRITE(*,*)"Jt·A·J"
        ! do x = lbound(E,1), ubound(E,1)
        !   WRITE(*,*) (E(x,y), y = lbound(E,2), ubound(E,2))
        ! end do
        ! WRITE(*,*)" "
        ! *******************************************************************

        do x = lbound(A,1), ubound(A,1) !<- Extract diagonal values of A
          VA(x) = A(x,x)
        end do

        do x = lbound(E,1), ubound(E,1) !<- Extract diagonal values of E
          VE(x) = E(x,x)
        end do

        DV = ABS(VE)-ABS(VA)    ! Take the difference in the norm of the diagonal elements of A and B
        delta = SUM(DV)         ! Sum all elements in DV. If all differences were small, delta should be small

        ! ********************** A tool for debugging ***********************
        ! WRITE(*,*) delta
        ! *******************************************************************

        if (delta < 1E-10) then        !<- When delta is less than E-10, add 1 to the counter
          deltacount = deltacount + 1
        else if (delta > 1E-10) then   !<- When delta is bigger than E-10, clean counter
          deltacount = 0               !<- This ensures convergence through multiple consecutive identical results
        end if

        A = E                          !<- Take the rotated matrix and redo

      end do
    end do                             !<- End sweep
  end do                               !<- End do while

  diagonalized_mat = A
  total_jcbi = Jtot

end subroutine jacobi_diagonalization

!
!
!                 Dad, why do I look like the mail man?
!
!

subroutine inverse_sqrt(diagonal_overlap, m1, inverse_diagonal_overlap)

  ! ------------------- Declare stuff coming in and out ---------------------
  real(kind=8), dimension(:,:), intent(in) :: diagonal_overlap
  real(kind=8), dimension(:,:), intent(out) :: inverse_diagonal_overlap
  integer, intent(in) :: m1
  !
  !
  !
  ! ------------------ Declare stuff used in the subroutine -----------------
  integer :: i, j

  ! Take inverse square root of diagonal elements in diagonalized overlap matrix

  do i = lbound(diagonal_overlap, 1), ubound(diagonal_overlap, 1)
    do j = lbound(diagonal_overlap, 2), ubound(diagonal_overlap, 2)
      if (i == j) then
        inverse_diagonal_overlap(i, j) = 1/sqrt(diagonal_overlap(i,j))
      else
        inverse_diagonal_overlap(i,j) = 0.0
      end if
    end do
  end do

end subroutine inverse_sqrt

!
!
!                       Does this look infected?
!
!

subroutine similarity_transformation(incoming_matrix, transforming_mat, m1, &
transformed_mat)
  ! ------------------- Declare stuff coming in and out ---------------------
  real(kind=8), dimension(:,:), intent(in) :: incoming_matrix, &
  transforming_mat
  real(kind=8), dimension(:,:), intent(out) :: transformed_mat
  integer, intent(in) :: m1
  ! ------------------ Declare stuff used in the subroutine -----------------
  real(kind=8), dimension(:,:), allocatable :: A, P, C, Pt, B
  real(kind=8) :: product
  integer :: x, y, i, j

  ! ---------------- Allocate all the intermediates to use ------------------
  allocate (A(m1,m1))
  allocate (P(m1,m1))
  allocate (C(m1,m1))
  allocate (Pt(m1,m1))
  allocate (B(m1,m1))

  ! ---------- Set variables so that I don't have to rewrite this code ------
  !                                           This code performs P(trans)*A*P
  A = incoming_matrix
  P = transforming_mat

  ! --------------- Begin similarity transformation process -----------------

  !                                                       Generate P transpose

  do x = lbound(P,1), ubound(P,1)
    do y = lbound(P,2), ubound(P,2)
      Pt(y,x) = P(x,y)
    end do
  end do

  product = 0

  !                               Obtain matrix C, which is the product of A·Pt

  do i = lbound(A,1), ubound(A,1)           !<- This counter is the "slowest", it will allow us to move the rows in A and C
    do x = lbound(P,2), ubound(P,2)         !<- This counter sets the column in matrix B
      do y = lbound(A,2), ubound(A,2)       !<- This counter is the "fastest", it sets the column in A and row in B
        product = A(i,y)*Pt(y,x) + product  !<- Multiply each element and add
      end do
      C(i,x) = product                      !<- Set the result in C
      product = 0                           !<- Clear variable value
    end do
  end do

  !                                     Once we have matrix C, multiply P by it

  do i = lbound(P,1), ubound(P,1)           !<- This counter is the "slowest", it will allow us to move the rows in P and C
    do x = lbound(C,2), ubound(C,2)         !<- This counter sets the column in matrix C
      do y = lbound(P,2), ubound(P,2)       !<- This counter is the "fastest", it sets the column in P and row in C
        product = P(i,y)*C(y,x) + product   !<- Multiply each element and add
      end do
      B(i,x) = product                      !<- Set the result in C
      product = 0                           !<- Clear variable value
    end do
  end do

  transformed_mat = B

end subroutine similarity_transformation

!
!
!          You should respect my opinion! *their opinion: 2+3 = 7*
!
!

subroutine set_up_calculation(full_basis, system_coord, nuclear_repulsion)
  ! ------------------- Declare stuff coming in and out ---------------------
  real(kind=8), dimension(3,2,2), intent(out) :: full_basis
  real(kind=8), dimension(3,2), intent(out) :: system_coord
  real(kind=8), intent(out) :: nuclear_repulsion
  ! ------------------ Declare stuff used in the subroutine -----------------
  ! The ranks of these arrays are specific for STO-3G
  real(kind=8), DIMENSION(3,2) :: basis_H, basis_He
  REAL(kind=8), DIMENSION(3) :: coordinates_H, coordinates_He
  real(kind=8) :: internuclear_dist
  INTEGER :: i, j, k, contraction_H, contraction_He

  ! ------------------ Setting up the system: Basis set ---------------------
  ! Basis function info is stored in a 2D array. Column 1 contains exponents,
  ! column 2 contains the corresponding contraction coefficients.

  ! ------------------------ BSE STO-3G basis set ---------------------------
  ! basis_H = reshape((/0.3425250914E+01, 0.1543289673E+00, 0.6239137298E+00, &
  ! 0.5353281423E+00, 0.1688554040E+00, 0.4446345422E+00/), (/3,2/), ORDER=(/2,1/))
  
  ! basis_He = reshape((/0.6362421394E+01, 0.1543289673E+00, 0.1158922999E+01, &
  ! 0.5353281423E+00, 0.3136497915E+00, 0.4446345422E+00/), (/3,2/), ORDER=(/2,1/))

  ! ---------------------- Szabo & Ostlund basis set ------------------------
  basis_H = reshape((/0.168856157, 0.444635, 0.62391349, 0.535328, &
  3.425250016, 0.154329/), (/3,2/), ORDER=(/2,1/))

  basis_He = reshape((/0.48084429, 0.444635, 1.776691148, 0.535328, &
  9.753934616, 0.154329/), (/3,2/), ORDER=(/2,1/))



  !                             Set up the array that packs the basis set info

  do i = lbound(full_basis,3), ubound(full_basis,3)
    do j = lbound(full_basis,1), ubound(full_basis,1)
      do k = lbound(full_basis,2), ubound(full_basis,2)
        if (i==1) then
          full_basis(j,k,i) = basis_H(j,k)
        elseif (i==2) then
          full_basis(j,k,i) = basis_He(j,k)
        end if
      end do
    end do
  end do

  ! ----------------- Setting up the system: Coordinates --------------------
  !
  ! Array contains coordinates at x, y, z in positions 1, 2, 3 respectively.

  coordinates_H(1) = 0.7316
  coordinates_H(2) = 0.0
  coordinates_H(3) = 0.0

  coordinates_He(1) = -0.7316
  coordinates_He(2) = 0.0
  coordinates_He(3) = 0.0

  !                         Setting the coordinates of the diatomic in an array

  do i = 1, 2                     !<--- Two atoms, two columns
    do j = 1, 3                   !<--- Three coordinates
      if (i == 1) then
        system_coord(j,i) = coordinates_H(j)
      elseif (i==2) then
        system_coord(j,i) = coordinates_He(j)
      end if
    end do
  end do

  internuclear_dist = 0.0

  do i = 1, 3                            !<--- Calculate internuclear distance
    internuclear_dist = internuclear_dist + &
    ((coordinates_H(i)-coordinates_He(i))*&
    (coordinates_H(i)-coordinates_He(i)))
  end do

  internuclear_dist = sqrt(internuclear_dist)

  nuclear_repulsion = 2/internuclear_dist

end subroutine set_up_calculation

!
!
!        When life gives you lemons... you start doing tequila shots
!
!

subroutine calculate_overlap_matrix(basis_info, coord_array, ovrlpmat)
  ! ------------------- Declare stuff coming in and out ---------------------
  real(kind=8), dimension(:,:,:), intent(in) :: basis_info
  real(kind=8), dimension(:,:), intent(in) :: coord_array
  real(kind=8), dimension(:,:), intent(out) :: ovrlpmat

  !
  !
  ! ------------------ Declare stuff used in the subroutine -----------------
  real(kind=8) :: coefA, coefB, exponentA, exponentB, integral, &
  total_integral, dist_ab2, function_distance
  integer :: i, j, k, l

  ! --------- Calculate distance H-He^2 for Gaussian Product Rule -----------
  dist_ab2 = 0.0

  do i = 1, 3
    dist_ab2 = dist_ab2 + ((coord_array(i,1) - coord_array(i,2)) * &
    (coord_array(i,1) - coord_array(i,2)))
  end do

  ! ---------------------- Generate overlap matrix --------------------------
  !
  ! The idea is the following: have do loops whose counters will correspond
  ! to shells in the system. This is will allow to fill the overlap matrix
  ! with those indices. Example, i=1 and j=2; in this case, the element to be
  ! calculated will be S_12. Inside this structure, there will be do loops
  ! that will carry out the integration of the primitives and add the results
  ! to the corresponding overlap matrix element.


   do i = 1, 2
     do j = 1, 2
       ovrlpmat(i,j) = 0.0       !<--- Set the element of the S matrix to 0
       total_integral = 0.0      !<--- Set the integral to 0 before computing it

       if (i == j) then                 !<--- If we are to evaluate an element
         function_distance = 0.0        !     of the diagonal, the distance
       else                             !     between the basis functions is
         function_distance = dist_ab2   !     0, otherwise it is the calculated
       end if                           !     squared distance.

       !    Loop around number of primitives to calculate individual integrals
       !                                                    between primitives

       do k = lbound(basis_info,1), ubound(basis_info,1)
         do l = lbound(basis_info,1), ubound(basis_info,1)
           exponentA = basis_info(k, 1, i)!<-- Pick exponent of primitive k, shell i
           exponentB = basis_info(l, 1, j)
           coefA = basis_info(k, 2, i)    !<-- Pick coefficient of primitive k, shell i
           coefB = basis_info(l, 2, j)

           ! ******************* A tool for debugging ************************
           ! write(*,*)"For element: ", i, j
           ! write(*,*)"Exponents to use, A: ", exponentA, "Exponent B: ", exponentB
           ! write(*,*)"Coefficients to use, A: ", coefA, "Coefficient B: ", coefB
           ! *****************************************************************



           ! Keep this code readable, and calculate the integral somwhere else
           call overlap_integral(exponentA, exponentB, coefA, coefB, &
           function_distance, integral)

           !     Add the value of the integral between the two primitves to the
           !                                      total integral between shells
           total_integral = total_integral + integral


         end do
       end do

       !      After looping over primitives, put the total integral between the
       !                                   shells in the overlap matrix element
       ovrlpmat(i,j) = total_integral

     end do
   end do

end subroutine calculate_overlap_matrix

!
!
!           When life turns its back on you... you grab its butt
!
!

subroutine overlap_integral(expA, expB, coeffA, coeffB, basis_distance, intgrl)
  ! ------------------- Declare stuff coming in and out ---------------------
  real(kind=8), intent(in) :: expA, expB, coeffA, coeffB, basis_distance
  real(kind=8), intent(out) :: intgrl
  ! ------------------ Declare stuff used in the subroutine -----------------
  real(kind=8) :: kp, gamma, pi, dummy, norm_coefA, norm_coefB
  integer :: i
  !
  ! The idea is to have this subroutine calculate the overlap integral when
  ! given the information of two basis functions (i.e., coefficients,
  ! exponents and distance between functions)
  !
  ! We will calculate the integral d_a·d_b·∫Xa·XbdV = da*dBKp·
  ! (pi/(alfa_a + alfa_b)), where Kp is obtained from the Gaussian
  ! Product Rule, d_a and alfa_a represent the contraction coefficient and
  ! exponent of primitive a, and respectively for b.

  pi = 3.141592653589793238462643383279
  gamma = expA + expB       !<--- Place the sum of exponents in gamma variable

  kp = exp(-((expA*expB*basis_distance)/gamma))

  ! Dummy represents the term to elevate to the 3/2

  dummy = ((pi*pi*pi)/(gamma*gamma*gamma))

  dummy = sqrt(dummy)

  norm_coefA = coeffA*((2*expA/pi)**(0.75))  ! Add normalization to contraction
  !                                                                 coefficient
  norm_coefB = coeffB*((2*expB/pi)**(0.75))

  intgrl = kp*dummy*norm_coefA*norm_coefB    ! <--- Put it all toghether

end subroutine overlap_integral

!
!
!                          Mommy says I am special :)
!
!

subroutine calculate_T_matrix(basis_info, coord_array, kin_engy_mat)
  ! ------------------- Declare stuff coming in and out ---------------------
  real(kind=8), dimension(:,:,:), intent(in) :: basis_info
  real(kind=8), dimension(:,:), intent(in) :: coord_array
  real(kind=8), dimension(:,:), intent(out) :: kin_engy_mat
  !
  !
  !
  ! ------------------ Declare stuff used in the subroutine -----------------
  real(kind=8) :: coefA, coefB, exponentA, exponentB, integral, &
  total_integral, dist_ab2, function_distance
  integer :: i, j, k, l

  ! --------- Calculate distance H-He^2 for Gaussian Product Rule -----------
  dist_ab2 = 0.0

  do i = 1, 3
    dist_ab2 = dist_ab2 + ((coord_array(i,1) - coord_array(i,2)) * &
    (coord_array(i,1) - coord_array(i,2)))
  end do

  ! ------------------- Generate kinetic energy matrix ----------------------
  !
  ! The idea is the following: have do loops whose counters will correspond
  ! to shells in the system. This is will allow to fill the kinetic energy
  ! matrix with those indices. Example: i=1 and j=2; in this case, the element
  ! to be calculated will be T_12. Inside this structure, there will be do
  ! loops that will carry out the integration of the primitives and add the
  ! results to the corresponding kinetic energy matrix element.

  do i = 1, 2
    do j = 1, 2
      kin_engy_mat(i,j) = 0.0   !<--- Set the element of the T matrix to 0
      total_integral = 0.0      !<--- Set the integral to 0 before computing it

      if (i == j) then                 !<--- If we are to evaluate an element
        function_distance = 0.0        !     of the diagonal, the distance
      else                             !     between the basis functions is
        function_distance = dist_ab2   !     0, otherwise it is the calculated
      end if                           !     squared distance.

      !    Loop around number of primitives to calculate individual integrals
      !                                                    between primitives

      do k = lbound(basis_info,1), ubound(basis_info,1)
        do l = lbound(basis_info,1), ubound(basis_info,1)
          exponentA = basis_info(k, 1, i)!<-- Pick exponent of primitive k, shell i
          exponentB = basis_info(l, 1, j)
          coefA = basis_info(k, 2, i)    !<-- Pick coefficient of primitive k, shell i
          coefB = basis_info(l, 2, j)

          ! ******************* A tool for debugging ************************
          ! write(*,*)"For element: ", i, j
          ! write(*,*)"Exponents to use, A: ", exponentA, "Exponent B: ", exponentB
          ! write(*,*)"Coefficients to use, A: ", coefA, "Coefficient B: ", coefB
          ! *****************************************************************



          ! Keep this code readable, and calculate the integral somewhere else
          call kinetic_integral(exponentA, exponentB, coefA, coefB, &
          function_distance, integral)

          !     Add the value of the integral between the two primitves to the
          !                                      total integral between shells
          total_integral = total_integral + integral


        end do
      end do

      !      After looping over primitives, put the total integral between the
      !                                   shells in the overlap matrix element
      kin_engy_mat(i,j) = total_integral

    end do
  end do

end subroutine calculate_T_matrix

!
!
!                 Silence is golden... duct tape is silver
!
!

subroutine kinetic_integral(expA, expB, coeffA, coeffB, basis_distance, &
intgrl)
  ! ------------------- Declare stuff coming in and out ---------------------
  real(kind=8), intent(in) :: expA, expB, coeffA, coeffB, basis_distance
  real(kind=8), intent(out) :: intgrl
  ! ------------------ Declare stuff used in the subroutine -----------------
  real(kind=8) :: kp, gamma, pi, dummy, norm_coefA, norm_coefB, dummy2, dummy3
  !
  ! This subroutine calculates the kinetic energy integral for two given
  ! primitive gaussian functions (i.e., contraction coefficients, exponents,
  ! distance between functions).
  !
  ! The integral is calculated according to the expression provided in the
  ! book Modern Quantum Chemistry by A. Szabo & N. S. Ostlund (First edition,
  ! Appendix 1, equation A.11):
  !
  ! ∫Xa·[-(1/2)·/nabla^(1/2)]·Xb·dV = [alfa_a·alfa_b/(alfa_a+alfa_b)]·[3 -
  !                                    2·alfa_a·alfa_b/(alfa_a+alfa_b)·
  !                                    |Ra-Rb|^2]·[(pi/(alfa_a+alfa_b))^(3/2]·
  !                                    [exp(-((alfa_a·alfa_b/(alfa_a+alfa_b))·
  !                                    |Ra-Rb|^2))]
  !

  pi = 3.141592653589793238462643383279
  gamma = expA + expB

  kp = exp(-((expA*expB*basis_distance)/gamma)) !<-- From gaussian product rule,
  !                                                  last [] in the expression

  ! Dummy represents the term to elevate to the 3/2

  dummy = ((pi*pi*pi)/(gamma*gamma*gamma))

  dummy = sqrt(dummy)

  ! Dummy2 represents the other term with |Ra - Rb|^2

  dummy2 = 3 - ((2*expA*expB*basis_distance)/gamma)

  ! Dummy3 is just alfa_a·alfa_b/(alfa_a + alfa_b)

  dummy3 = (expA*expB)/gamma

  ! Include normalization in coefficients

  norm_coefA = coeffA*((2*expA/pi)**(0.75))

  norm_coefB = coeffB*((2*expB/pi)**(0.75))

  ! Bring everything toghether

  intgrl = norm_coefA*norm_coefB*dummy3*dummy2*dummy*kp

end subroutine kinetic_integral

!
!
!                   At least the PhD means I get a job, right?
!
!

subroutine calculate_V_matrices(basis_info, coord_array, V1_matrix, V2_matrix, &
Vtot_matrix)
  ! ------------------- Declare stuff coming in and out ---------------------
  real(kind=8), dimension(:,:,:), intent(in) :: basis_info
  real(kind=8), dimension(:,:), intent(in) :: coord_array
  real(kind=8), dimension(:,:), intent(out) :: V1_matrix, V2_matrix, Vtot_matrix
  !
  !
  !
  ! ------------------ Declare stuff used in the subroutine -----------------
  real(kind=8) :: coefA, coefB, exponentA, exponentB, integral, &
  total_integral, dist_ab2, function_distance, H_nuc_charge, He_nuc_charge
  real(kind=8), dimension(3) :: H_atom_coordinates, He_atom_coordinates, &
  P_coordinates
  integer :: i, j, k, l, m

  ! --------- Calculate distance H-He^2 for Gaussian Product Rule -----------
  dist_ab2 = 0.0

  do i = 1, 3
    dist_ab2 = dist_ab2 + ((coord_array(i,1) - coord_array(i,2)) * &
    (coord_array(i,1) - coord_array(i,2)))
  end do

  ! ------------------ Prepare data specific to each atom -------------------

  H_nuc_charge = 1.0    !<--- Set atomic numbers
  He_nuc_charge = 2.0   !     of H and He

  do i = 1, 2           !<--- Two atoms, two columns
    do j = 1, 3         !<--- Three coordinates
      if (i==1) then    !<--- Column 1 = Hydrogen, Column 2 = Helium
        H_atom_coordinates(j) = coord_array(j, i)
      elseif (i==2) then
        He_atom_coordinates(j) = coord_array(j, i)
      end if
    end do
  end do


  ! ------------- Generate attraction potential energy matrix ---------------
  !
  ! The idea is the following: have do loops whose counters will correspond
  ! to shells in the system. This is will allow to fill each V matrix with
  ! those indices. Example; i=1 and j=2; in this case, the element to be
  ! calculated will be S_12. Inside this structure, there will be do loops that
  ! will carry out the integration over the primitives and add the results to
  ! the corresponding V matrix element. This process will be carried out once
  ! for every atom in the system.



  ! ------------------- Generate V1 (Hydrogen) matrix -----------------------

  do i = 1, 2
    do j = 1, 2
      V1_matrix(i,j) = 0.0      !<--- Set the element of the V1 matrix to 0
      total_integral = 0.0      !<--- Set the integral to 0 before computing it

      if (i == j) then                 !<--- If we are to evaluate an element
        function_distance = 0.0        !     of the diagonal, the distance
      else                             !     between the basis functions is
        function_distance = dist_ab2   !     0, otherwise it is the calculated
      end if                           !     squared distance.

      !    Loop around number of primitives to calculate individual integrals
      !                                                    between primitives

      do k = lbound(basis_info,1), ubound(basis_info,1)
        do l = lbound(basis_info,1), ubound(basis_info,1)
          exponentA = basis_info(k, 1, i)!<-- Pick exponent of primitive k, shell i
          exponentB = basis_info(l, 1, j)
          coefA = basis_info(k, 2, i)    !<-- Pick coefficient of primitive k, shell i
          coefB = basis_info(l, 2, j)

          ! Compute coordinates of point P between gaussian functions A and B.
          ! Choose coordinates of atoms on which functions are centered
          ! according to indices i, j (e.g., i = 1, j = 2, function A is
          ! centered on H and function B is centered on He)
          P_coordinates = 0.0
          do m = 1, 3
            P_coordinates(m) = ((exponentA*coord_array(m, i)) + &
            (exponentB*coord_array(m, j)))/(exponentA+exponentB)
          end do

          ! We have our data, keep this readable and compute the integral
          !                                                somewhere else
          call attr_integral(exponentA, exponentB, coefA, coefB, &
          function_distance, H_atom_coordinates, P_coordinates, H_nuc_charge, &
          integral)

          !     Add the value of the integral between the two primitves to the
          !                                      total integral between shells
          total_integral = total_integral + integral


        end do
      end do

      !      After looping over primitives, put the total integral between the
      !                                   shells in the overlap matrix element
      V1_matrix(i,j) = total_integral

    end do
  end do




  ! ------------------- Generate V2 (Helium) matrix -----------------------




  do i = 1, 2
    do j = 1, 2
      V2_matrix(i,j) = 0.0      !<--- Set the element of the V1 matrix to 0
      total_integral = 0.0      !<--- Set the integral to 0 before computing it

      if (i == j) then                 !<--- If we are to evaluate an element
        function_distance = 0.0        !     of the diagonal, the distance
      else                             !     between the basis functions is
        function_distance = dist_ab2   !     0, otherwise it is the calculated
      end if                           !     squared distance.

      !    Loop around number of primitives to calculate individual integrals
      !                                                    between primitives

      do k = lbound(basis_info,1), ubound(basis_info,1)
        do l = lbound(basis_info,1), ubound(basis_info,1)
          exponentA = basis_info(k, 1, i)!<-- Pick exponent of primitive k, shell i
          exponentB = basis_info(l, 1, j)
          coefA = basis_info(k, 2, i)    !<-- Pick coefficient of primitive k, shell i
          coefB = basis_info(l, 2, j)

          ! Compute coordinates of point P between gaussian functions A and B.
          ! Choose coordinates of atoms on which functions are centered
          ! according to indices i, j (e.g., i = 1, j = 2, function A is
          ! centered on H and function B is centered on He)
          P_coordinates = 0.0
          do m = 1, 3
            P_coordinates(m) = ((exponentA*coord_array(m, i)) + &
            (exponentB*coord_array(m, j)))/(exponentA+exponentB)
          end do

          ! We have our data, keep this readable and compute the integral
          !                                                somewhere else
          call attr_integral(exponentA, exponentB, coefA, coefB, &
          function_distance, He_atom_coordinates, P_coordinates, &
          He_nuc_charge, integral)

          !     Add the value of the integral between the two primitves to the
          !                                      total integral between shells
          total_integral = total_integral + integral


        end do
      end do

      !      After looping over primitives, put the total integral between the
      !                                   shells in the overlap matrix element
      V2_matrix(i,j) = total_integral

    end do
  end do

  Vtot_matrix = V1_matrix + V2_matrix

end subroutine calculate_V_matrices

!
!
!                       Don't feed me after midnight
!
!

subroutine attr_integral(expA, expB, coeffA, coeffB, basis_distance, &
nuclear_coord, pointP_coord, nuclear_charge, intgrl)
  ! ------------------- Declare stuff coming in and out ---------------------
  real(kind=8), dimension(3), intent(in) :: nuclear_coord, pointP_coord
  real(kind=8), intent(in) :: expA, expB, coeffA, coeffB, basis_distance, &
  nuclear_charge
  real(kind=8), intent(out) :: intgrl
  !
  !
  !
  ! ------------------ Declare stuff used in the subroutine -----------------
  real(kind=8) :: kp, gamma, pi, norm_coefA, norm_coefB, dummy1, dummy2, &
  dist_pc2, argument
  integer :: i
  !
  ! This subroutine calculates electron - nuclear potential energy integrals
  ! for two given primitive gaussian functions (i.e., contraction coefficients,
  ! exponents, distance between functions, coordinates of point P between
  ! functions), and a given nucleus (i.e., nuclear charge and coordinates of
  ! the nucleus).
  !
  ! The integral is calculated according to the expression provided in the book
  ! Modern Quantum Chemistry by A. Szabo & N. S. Ostlund (First edition,
  ! Appendix 1, equation A.33).
  !
  !

  pi = 3.141592653589793238462643383279
  gamma = expA + expB

  !                                      Compute first elements of the integral

  dummy1 = (2*pi*nuclear_charge)/gamma

  kp = exp(-((expA*expB*basis_distance)/gamma))

  !          Now comes the fun part, the F_0 function and |Rp-Rc|^2

  !                                                               Get |Rp-Rc|^2
  dist_pc2 = 0.0

  do i = 1, 3
    dist_pc2 = dist_pc2 + ((pointP_coord(i) - nuclear_coord(i)) * &
    (pointP_coord(i) - nuclear_coord(i)))
  end do

  !                                   Compute the argument of the F_0 function
  argument = gamma * dist_pc2

  !                                                Dummy2 will be F_0(argument)
  ! Turns out that if the argument is very small, you will get a NaN. To
  ! account for this, when that is the case, you evaluate for the asymptotic
  ! value.
  if (argument < 1E-6) then
    dummy2 = 1.0 - (argument/3.0)
  else
    dummy2 = (0.5*sqrt((pi/argument))) * erf((sqrt(argument)))
  end if
  !                                       Include normalization in coefficients

  norm_coefA = coeffA*((2*expA/pi)**(0.75))

  norm_coefB = coeffB*((2*expB/pi)**(0.75))

  !                                                   Put everything toghether
  intgrl = norm_coefA*norm_coefB*(-1)*dummy1*kp*dummy2

end subroutine attr_integral

!
!
!                         So, what exactly are we?
!
!

subroutine calculate_two_elec_ints(basis_info, coord_array, two_e_ints)
  ! ------------------- Declare stuff coming in and out ---------------------
  real(kind=8), dimension(:,:,:), intent(in) :: basis_info
  real(kind=8), dimension(:,:), intent(in) :: coord_array
  real(kind=8), dimension(:,:,:,:), intent(out) :: two_e_ints
  !
  !
  !
  ! ------------------ Declare stuff used in the subroutine -----------------
  real(kind=8) :: coefA, coefB, coefC, coefD, exponentA, exponentB, exponentC, &
  exponentD, integral, total_integral, interatomic_dist, dist_ab2, dist_cd2, &
  dist_pq2
  real(kind=8), dimension(3) :: P_coordinates, Q_coordinates, &
  H_atom_coordinates, He_atom_coordinates
  integer :: i, j, k, l, m, n, o, p, r

  ! --------- Calculate distance H-He^2 for Gaussian Product Rule -----------
  interatomic_dist = 0.0

  do i = 1, 3
    interatomic_dist = interatomic_dist + &
    ((coord_array(i,1) - coord_array(i,2)) * &
    (coord_array(i,1) - coord_array(i,2)))
  end do

  ! ------------- Set atom coordinates for gaussian product rule ------------

  do i = 1, 2           !<--- Two atoms, two columns
    do j = 1, 3         !<--- Three coordinates
      if (i==1) then    !<--- Column 1 = Hydrogen, Column 2 = Helium
        H_atom_coordinates(j) = coord_array(j, i)
      elseif (i==2) then
        He_atom_coordinates(j) = coord_array(j, i)
      end if
    end do
  end do

  ! -------- Generate 4 dimensional tensor for 2e integral storage ----------
  !
  ! The idea is the following: have do loops whose counters (i, j, k, l)
  ! correspond to the shells in the system. This will allow to compute
  ! every possible combination of functions (i.e., 1111, 1112, etc.). Inside
  ! this structure there will be do loops that will carry out integration over
  ! primitives and add the results to the corresponding (ij|kl) element.
  !

  do i = 1, 2
    do j = 1, 2
      do k = 1, 2
        do l = 1, 2

         two_e_ints(i,j,k,l) = 0.0 !<--- Set the element to 0
         total_integral = 0.0      !<--- Set the integral initially to 0

         if (i == j) then          !<--- If i, j share center, the distance
           dist_ab2 = 0.0          !     between them is 0. If i, j have
         else if (i /= j) then     !     different centers, then they are
           dist_ab2 = interatomic_dist ! separated by the distance H-He
         end if

         if (k == l) then          !<--- If k, l share center, the distance
           dist_cd2 = 0.0          !     between them is 0. If k, l have
         else if (k /= l) then     !     different centers, then they are
           dist_cd2 = interatomic_dist ! separated by the distance H-He
         end if

         !      Now we loop over primitives to calculate individual integrals
         !                  between primitives. Four primitives, four indices
         do m = lbound(basis_info, 1), ubound(basis_info, 1)
           do n = lbound(basis_info, 1), ubound(basis_info, 1)
             do o = lbound(basis_info, 1), ubound(basis_info, 1)
               do p = lbound(basis_info, 1), ubound(basis_info, 1)
                 exponentA = basis_info(m, 1, i) !<--- Pick exponent of
                 exponentB = basis_info(n, 1, j) !     primitive m, taken
                 exponentC = basis_info(o, 1, k) !     from shell i.
                 exponentD = basis_info(p, 1, l)

                 coefA = basis_info(m, 2, i)     !<--- Pick coefficient of
                 coefB = basis_info(n, 2, j)     !     primitive m, taken
                 coefC = basis_info(o, 2, k)     !     from shell i.
                 coefD = basis_info(p, 2, l)

                 ! We will perform the gaussian product rule between
                 ! primitives m and n, and between primitives o and p.
                 ! Results are termed P for m-n and Q for o-p. Choose
                 ! coordinates of atoms on which functions are centered
                 ! according to indices i, j, k, l (e.g., i=1, j=2, k=2,
                 ! l=2, function m is centered on H, function n is centered
                 ! on atom He, function o is centered on atom He, function p
                 ! is centered on atom He)
                 P_coordinates = 0.0
                 Q_coordinates = 0.0
                 do r = 1, 3
                   P_coordinates(r) = ((exponentA*coord_array(r, i)) + &
                   (exponentB*coord_array(r, j)))/(exponentA+exponentB)

                   Q_coordinates(r) = ((exponentC*coord_array(r, k)) + &
                   (exponentD*coord_array(r, l)))/(exponentC+exponentD)
                 end do

                 ! We have the coordinates of P and Q, we now compute the
                 ! distance^2 between P and Q
                 dist_pq2 = 0.0
                 do r = 1, 3
                   dist_pq2 = dist_pq2 + &
                   ((P_coordinates(r) - Q_coordinates(r)) * &
                   (P_coordinates(r) - Q_coordinates(r)))
                 end do

                 ! *********************** Debugging ************************
                 ! write(*,*)"P-Q distance: ", dist_pq2
                 ! **********************************************************


                 ! We have everything we need for our integral, keep this
                 !       readable and compute the integral somewhere else
                 call repulsion_integral(exponentA, exponentB, exponentC, &
                 exponentD, coefA, coefB, coefC, coefD, dist_ab2, dist_cd2, &
                 dist_pq2, integral)

                 ! Add the value of the integral between the four primitives
                 !                      to the total integral between shells
                 total_integral = total_integral + integral

               end do
             end do
           end do
         end do !<---- Finish looping over primitives

         ! After looping over primitives, put the total integral between the
         !                              four shells in the 4D tensor element

          ! ******************** A tool for debugging ************************
          ! write(*,*)"4D tensor element: ", i, j, k, l
          ! write(*,*)total_integral
          ! ******************************************************************

         two_e_ints(i, j, k, l) = total_integral

       end do
     end do
   end do
 end do

end subroutine calculate_two_elec_ints

!
!
!                       At least I don't have diarrhea...
!
!

subroutine repulsion_integral(expA, expB, expC, expD, coeffA, coeffB, coeffC, &
coeffD, ab_dist, cd_dist, pq_dist, intgrl)
  ! ------------------- Declare stuff coming in and out ---------------------
  real(kind=8), intent(in) :: expA, expB, expC, expD, coeffA, coeffB, &
  coeffC, coeffD, ab_dist, cd_dist, pq_dist
  real(kind=8), intent(out) :: intgrl
  ! ------------------ Declare stuff used in the subroutine -----------------
  real(kind=8) :: pi, kpq, dummy1, dummy2, argument, gamma_ab, gamma_cd, &
  gamma_abcd, norm_coefA, norm_coefB, norm_coefC, norm_coefD
  !
  ! This subroutine calculates electron - electron repulsion integrals for four
  ! given primitive gaussians (i.e., contraction coefficients, exponents,
  ! distances, etc.)
  !
  ! The integral is calculated according to the expression provided in the book
  ! Modern Quantum Chemistry by A. Szabo & N. S. Ostlund (First edition,
  ! Appendix 1, equation A.41).
  !
  !

  pi = 3.141592653589793238462643383279
  gamma_ab = expA + expB
  gamma_cd = expC + expD
  gamma_abcd = expA + expB + expC + expD

  !                                                              Compute K_pq

  kpq = exp((-1)*((expA*expB*ab_dist/gamma_ab)+(expC*expD*cd_dist/gamma_cd)))

  !        Dummy1 represents the term including 2pi^(5/2)/[bunch_of_exponents]

  dummy1 = (2*sqrt(pi*pi*pi*pi*pi))/((gamma_ab)*(gamma_cd)*(sqrt(gamma_abcd)))

  !                                  Compute the argument of the F_0 function

  argument = (gamma_ab*gamma_cd*pq_dist)/(gamma_abcd)

  !                                              Dummy2 will be F_0(argument)
  !                Account for very small arguments with asymptotic F_0 value
  if (argument < 1E-6) then
    dummy2 = 1.0 - (argument/3.0)
  else
    dummy2 = (0.5*sqrt((pi/argument))) * erf((sqrt(argument)))
  end if

  !                                     Include normalization in coefficients

  norm_coefA = coeffA*((2*expA/pi)**(0.75))

  norm_coefB = coeffB*((2*expB/pi)**(0.75))

  norm_coefC = coeffC*((2*expC/pi)**(0.75))

  norm_coefD = coeffD*((2*expD/pi)**(0.75))

  !                                                  Put everything toghether
  intgrl = norm_coefA*norm_coefB*norm_coefC*norm_coefD*kpq*dummy1*dummy2

end subroutine repulsion_integral

!
!
!                        Never google your symptoms
!
!

subroutine generate_G_matrix(two_e_ints, density_mat, gmat)
  ! ------------------- Declare stuff coming in and out ---------------------
  real(kind=8), dimension(:,:,:,:), intent(in) :: two_e_ints
  real(kind=8), dimension(:,:), intent(in) :: density_mat
  real(kind=8), dimension(:,:), intent(out) :: gmat
  !
  !
  !
  ! ------------------ Declare stuff used in the subroutine -----------------
  real(kind=8) :: dummy1, totalG_element
  integer :: i, j, k, l
  !
  !
  ! This subroutine generates the G matrix (two electron contribution to
  ! Fock matrix), given the density matrix and the 4D tensor that stores
  ! all possible 2-electron integrals.
  !
  ! The idea is the following: For the G matrix element (i,j), use those
  ! same indices to chose two shells. These shells remain fixed while looping
  ! over k,l. k,l set the density matrix element to choose, and also complete
  ! the shell quartet to perform the operations between ijkl elements (that is:
  ! P_kl[(ij|lk)-1/2(ik|lj)])
  !
  !
  do i = 1, 2
    do j =1, 2
      totalG_element = 0.0    ! <--- Set variable that stores G_ij to zero

      do k = 1, 2
        do l = 1, 2
          dummy1 = density_mat(k,l) * (two_e_ints(i,j,l,k) - &
          (0.5*two_e_ints(i,k,l,j)))  ! <--- Calculate P_kl[(ij|lk)-1/2(ik|lj)])

          totalG_element = totalG_element + dummy1 !<--- Add k,l element to
          !                                              total i,j element
        end do
      end do

      gmat(i,j) = totalG_element   !<--- Put the element in the G matrix
    end do
  end do

end subroutine generate_G_matrix

!
!
!                           Scotty doesn't know
!
!

subroutine generate_Hcore_matrix(kin_engy_mat, Vtot_matrix, H_core_matrix)
  ! ------------------- Declare stuff coming in and out ---------------------
  real(kind=8), dimension(:,:), intent(in) :: kin_engy_mat, Vtot_matrix
  real(kind=8), dimension(:,:), intent(out) :: H_core_matrix
  !

  H_core_matrix = kin_engy_mat + Vtot_matrix!<--- Hcore is just T+V (1 electron)

end subroutine generate_Hcore_matrix

!
!
!                                   42
!
!

subroutine generate_Fock_matrix(H_core_matrix, gmat, Fock_mat)
  ! ------------------- Declare stuff coming in and out ---------------------
  real(kind=8), dimension(:,:), intent(in) :: H_core_matrix, gmat
  real(kind=8), dimension(:,:), intent(out) :: Fock_mat
  !

  Fock_mat = H_core_matrix + gmat

end subroutine generate_Fock_matrix

!
!
! All men are created equal... then they are immediately sorted out by society
!
!

subroutine eigenvalue_vector_ordering(diagonal_mat, eigenvector_mat)
  ! ------------------- Declare stuff coming in and out ---------------------
  real(kind=8), dimension(:,:), intent(inout) :: diagonal_mat, &
  eigenvector_mat
  !
  !
  ! ------------------ Declare stuff used in the subroutine -----------------
  real(kind=8), dimension(:), allocatable :: eigenvalues, dummy2
  real(kind=8) :: dummy1
  integer i, j, k

  allocate(eigenvalues(ubound(diagonal_mat,1)))
  allocate(dummy2(ubound(eigenvector_mat,2)))

  ! This subroutine orders the eigenvectors in a previously-diagonalized matrix.
  ! Eigenvalues are extracted from the diagonal of 'diagonal_mat', stored in
  ! eigenvector_array, ordered from smallest to largest. This heirarchy is
  ! used to order the eigenvector matrix too.
  !
  do i = lbound(diagonal_mat,1), ubound(diagonal_mat,1) !<--- Extract diagonal
    eigenvalues(i) = diagonal_mat(i,i)
  end do

  !                               Begin sorting of eigenvalues and eigenvectors

  do i = 1, ubound(eigenvalues,1)-1
    do j = 1, ubound(eigenvalues,1)-1
      if (eigenvalues(j) > eigenvalues(j+1)) then

        !                                                 Eigenvalue ordering

        dummy1 = eigenvalues(j+1)           ! The idea is this: if a given
        eigenvalues(j+1) = eigenvalues(j)   ! element is larger than the element
        eigenvalues(j) = dummy1             ! to its right, their positions are
        !                                     swapped

        !                                                Eigenvector ordering
        !      Eigenvectors will be moved in the same way as eigenvalues were

        do k = 1, ubound(eigenvector_mat,1)     ! Save column j+1 to dummy
          dummy2(k) = eigenvector_mat(k, j+1)   ! variable.
        end do

        do k = 1, ubound(eigenvector_mat,1)     ! Move column j to the right
          eigenvector_mat(k, j+1) = eigenvector_mat(k, j)
        end do

        do k = 1, ubound(eigenvector_mat,1)     ! Put column j+1 in column j
          eigenvector_mat(k, j) = dummy2(k)
        end do

      end if
    end do
  end do

  diagonal_mat = 0.0 !<--- Clear diagonalized matrix
  !                                      Put back eigenvalues in the diagonal
  do i = lbound(diagonal_mat,1), ubound(diagonal_mat,1)
    diagonal_mat(i,i) = eigenvalues(i)
  end do

end subroutine eigenvalue_vector_ordering

!
!
!         Every pizza is individual size if you believe in yourself
!
!

subroutine matrix_multiplication(matrixA, matrixB, resulting_matrix)
  ! ------------------- Declare stuff coming in and out ---------------------
  real(kind=8), dimension(:,:), intent(in) :: matrixA, matrixB
  real(kind=8), dimension(:,:), intent(out) :: resulting_matrix
  !
  !
  ! ------------------ Declare stuff used in the subroutine -----------------
  real(kind=8) :: product
  integer :: x, y, i

  !                                        This code performs matrixA*matrixB

  do i = lbound(matrixA,1), ubound(matrixA,1)     !<- This counter is the "slowest", it will allow us to move the rows in A and C
    do x = lbound(matrixB,2), ubound(matrixB,2)   !<- This counter sets the column in matrix B
      do y = lbound(matrixA,2), ubound(matrixA,2)      !<- This counter is the "fastest", it sets the column in A and row in B
        product = matrixA(i,y)*matrixB(y,x) + product  !<- Multiply each element and add
      end do
      resulting_matrix(i,x) = product                      !<- Set the result in C
      product = 0                           !<- Clear variable value
    end do
  end do

end subroutine matrix_multiplication

!
!
!                      One must imagine Sisyphus happy
!
!

subroutine generate_new_density_mat(coeff_mat, new_density_mat)
  ! ------------------- Declare stuff coming in and out ---------------------
  real(kind=8), dimension(:,:), intent(in) :: coeff_mat
  real(kind=8), dimension(:,:), intent(out) :: new_density_mat
  !
  !
  ! ------------------ Declare stuff used in the subroutine -----------------
  real(kind=8) :: product
  integer :: i, j, k

  do i = lbound(coeff_mat,1), ubound(coeff_mat,1)
    do j = lbound(coeff_mat,2), ubound(coeff_mat,2)
      product = 0.0
      do k = lbound(coeff_mat,1), (ubound(coeff_mat,1)/2)
        product = product + (2*coeff_mat(i,k)*coeff_mat(j,k))
      end do

      new_density_mat(i,j) = product
    end do
  end do

end subroutine generate_new_density_mat

!
!
!  Started from the bottom, now we are here (still pretty much at the bottom)
!
!

subroutine determine_convergence(old_density_mat, new_density_mat, &
convergence)
  ! ------------------- Declare stuff coming in and out ---------------------
  real(kind=8), dimension(:,:), intent(in) :: old_density_mat, &
  new_density_mat
  logical, intent(out) :: convergence
  !
  !
  ! ------------------ Declare stuff used in the subroutine -----------------
  real(kind=8) :: total_difference
  integer :: i, j

  total_difference = 0.0

  ! Sum the differences in all matrix elements. If total is less than a certain
  ! threshold, density is converged

  do i = lbound(new_density_mat,1), ubound(new_density_mat,1)
    do j = lbound(new_density_mat,2), ubound(new_density_mat,2)
      total_difference = total_difference + &
      abs(old_density_mat(i,j)-new_density_mat(i,j))
    end do
  end do

  if (total_difference < 1E-5) then !<--- 1E-5 is the threshold, may be changed
    convergence = .true.
  else
    convergence = .false.
  end if

end subroutine determine_convergence

!
!
!       Started from the top, now you are at the bottom (of this file)
!
!

subroutine calculate_electronic_energy(density_mat, Hcore_mat, fock_mat, &
total_E)
  ! ------------------- Declare stuff coming in and out ---------------------
  real(kind=8), dimension(:,:), intent(in) :: density_mat, Hcore_mat, &
  fock_mat
  real(kind=8), intent(out) :: total_E
  !
  !
  ! ------------------ Declare stuff used in the subroutine -----------------
  integer :: i, j

  !                Compute total electronic energy as in Szabo equation 3.274

  total_E = 0.0

  do i = lbound(Hcore_mat,1), ubound(Hcore_mat,1)
    do j = lbound(Hcore_mat,2), ubound(Hcore_mat,2)
      total_E = total_E + density_mat(j,i)*(Hcore_mat(i,j)+fock_mat(i,j))
    end do
  end do

  total_E = 0.5*total_E

end subroutine
