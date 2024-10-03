module mod_problem
  use Iterative_Solver_Problem
  private
  type, extends(Problem), public :: our_problem
    double precision, dimension(:, :), pointer, public :: m_rhs
  contains
    procedure, pass :: diagonals
    procedure, pass :: action
    procedure, pass :: RHS
  end type our_problem
contains
  !> @brief Optionally provide the diagonal elements of the underlying kernel. If
  !> implemented and returning true, the provided diagonals will be used by
  !> the base precondition() function for preconditioning (and therefore the precondition() function does
  !> not need to be reimplemented), and, in the case of linear problems, for selection of
  !> the P space. Otherwise, preconditioning will be done with precondition(), and any P
  !> space has to be provided manually.
  !> @param d On exit, contains the diagonal elements
  !> @return Whether diagonals have been provided.
  logical function diagonals(this, d)
    class(our_problem), intent(in) :: this
    double precision, intent(inout), dimension(:) :: d
    call A_diagonal(d, size(d))
    diagonals = .true.
  end function diagonals

  !> @brief Calculate the action of the kernel matrix on a set of parameters. Used by
  !> linear solvers, but not by the non-linear solvers (NonLinearEquations, Optimize).
  !> @param parameters The trial solutions for which the action is to be calculated
  !> @param actions The action vectors
  !> @param range The range of the space for which actions should be computed. It's OK to provide also the values outside this range (which will happen in a multiprocessing context), but they will be ignored by the solver.
  subroutine action(this, parameters, actions, range)
    class(our_problem), intent(in) :: this
    double precision, intent(in), dimension(:, :) :: parameters
    double precision, intent(inout), dimension(:, :) :: actions
    integer, dimension(2), intent(in) :: range
    call A_times_x(actions, parameters, ubound(actions, 1), ubound(actions, 2))
    ! note that only actions(range(1)+1:range(2),:) is required. Following code proves that in parallel case
    actions(:range(1),:) = 999
    actions(range(2)+1:,:) = 999
  end subroutine action
  !> @brief Provide the inhomogeneous part of one of the sets of linear equations. Implementation required only for linear equation solver.
  !> @param vector Will contain the requested RHS on exit
  !> @param instance Which equation set for which the RHS should be provided
  !> @param range The range of the space for which actions should be computed. It's OK to provide also the values outside this range (which will happen in a multiprocessing context), but they will be ignored by the solver.
  !> @return Whether the requested instance exists
  logical function RHS(this, vector, instance, range)
    class(our_problem), intent(in) :: this
    double precision, intent(inout), dimension(:) :: vector
    integer, dimension(2), intent(in) :: range
    integer, intent(in) :: instance
    RHS = .false.
    if (instance .ge. lbound(this%m_rhs, 2) .and. instance .le. ubound(this%m_rhs, 2)) then
      RHS = .true.
      vector(range(1)+1:range(2)) = this%m_rhs(range(1)+1:range(2), instance)
    end if
  end function RHS

end module mod_problem
program test_solver
  use iso_fortran_env, only : dp => real64
  use mod_problem
  use Iterative_Solver
  implicit double precision (a-h, o-z)
  real(dp), allocatable, target :: rhs(:, :)
  real(dp), allocatable, target :: diag(:)
  real(dp), allocatable :: A(:, :)
  real(dp), allocatable :: unit(:, :)
  integer, allocatable :: ipiv(:)
  integer :: i, nvec, ivec
  type(our_problem) :: problem
  real(dp), allocatable, dimension(:, :) :: c, g

  call mpi_init
  if (mpi_rank_global() .gt. 0) close(6)

  !size of matrix
  n = 10

  !rhs vectors
  nvec = 4
  allocate(rhs(n, nvec))
  call rhs_vectors(rhs, n, nvec)

  !matrix diagonal
  allocate(diag(n))
  call A_diagonal(diag, n)
  call outsqr(diag, n, n, 1, 'diagonal')

  !iterative solver stuff to be inserted here
  !dimension:       n
  !number of rhs's: nvec
  !rhs:             corresponding variable
  !preconditioner:  diag holds the matrix diagonal
  !matrix*vector:   A_times_x subroutine
  problem%m_rhs => rhs
  n_buffers = nvec ! could be more (not useful) or fewer
  n_buffers = 2
  allocate(c(n, n_buffers), g(n, n_buffers))
  call Solve_Linear_Equations(c, g, problem, verbosity = 2)
  call outsqr(c, n, n, min(nvec, n_buffers), 'solution to A.x=rhs')
  call outsqr(g, n, n, min(nvec, n_buffers), 'residual')
  call Iterative_Solver_Solution([nvec], c, g)
  write (6, *) 'last solution', c(:, 1)
  print *, 'Errors', Iterative_Solver_Errors()
  call Iterative_Solver_Finalize
  deallocate(c, g)

  call mpi_finalize


  !compile with -D_DEBUG -llapack to get the Lapack solution for comparison
#ifdef _DEBUG
  allocate(A(n,n))
  allocate(unit(n,n))
  unit(:,:)=0d0
  do i=1,n
     unit(i,i)=1d0
  enddo
  call A_times_x(A,unit,n,n)
  call outsqr(A,n,n,n,'A')
  call outsqr(rhs,n,n,nvec,'rhs vectors')

  !solve with lapack
  allocate(ipiv(n))
  info=0
  call dgesv(n,nvec,A,n,ipiv,rhs,n,info)
  if(info.ne.0) stop 'problem in dgesv'

  !print solution
  call outsqr(rhs,n,n,nvec,'solution to A.x=rhs')

  !restore matrix and test
  call A_times_x(A,unit,n,n)
  call outsqr(matmul(A,rhs),n,n,nvec,'test A.x')

  deallocate(A,unit,ipiv)
#endif

  deallocate(rhs, diag)
end program test_solver

!---------------------------------------------------------------------
subroutine A_times_x(Ax, x, n, nvec)
  !---------------------------------------------------------------------
  !apply A.x on nvec many vectors
  use iso_fortran_env, only : dp => real64
  implicit none
  integer, intent(in) :: n, nvec
  real(dp), intent(in) :: x(n, nvec)
  real(dp), intent(out) :: Ax(n, nvec)
  integer :: i, j, ivec
  real(dp) :: aij
  Ax(:, :) = 0d0
  do ivec = 1, nvec
    do i = 1, n
      do j = max(i - 2, 1), min(i + 2, n)
        aij = 10d0 * dble(i) + dble(j)
        Ax(i, ivec) = Ax(i, ivec) + aij * x(j, ivec)
      enddo
    enddo
  enddo
  return
end subroutine A_times_x

!---------------------------------------------------------------------
subroutine A_diagonal(diag, n)
  !---------------------------------------------------------------------
  !get matrix diagonal
  use iso_fortran_env, only : dp => real64
  implicit none
  integer, intent(in) :: n
  real(dp), intent(out) :: diag(n)
  real(dp), allocatable :: A(:, :), unit(:, :)
  integer :: i

  !just for convenience, we compute the whole matrix and then extract the
  !diagonal
  allocate(A(n, n))
  allocate(unit(n, n))
  unit(:, :) = 0d0
  do i = 1, n
    unit(i, i) = 1d0
  enddo
  call A_times_x(A, unit, n, n)
  do i = 1, n
    diag(i) = A(i, i)
  enddo
  deallocate(A, unit)
  return
end subroutine A_diagonal

!---------------------------------------------------------------------
subroutine rhs_vectors(rhs, n, nvec)
  !---------------------------------------------------------------------
  use iso_fortran_env, only : dp => real64
  implicit none
  integer, intent(in) :: n, nvec
  real(dp), intent(out) :: rhs(n, nvec)
  integer :: ivec, i
  rhs(:, :) = 0d0
  do ivec = 1, nvec
    !    do i = max(i - 3, 1), min(i + 3, n) ! I don't understand this at all !!!
    do i = max(ivec - 3, 1), min(ivec + 3, n) ! my guess at what was intended
      rhs(i, ivec) = dble(i)
    enddo
  enddo
  return
end subroutine rhs_vectors

!---------------------------------------------------------------------
subroutine outsqr(a, n, ia, ib, name)
  !---------------------------------------------------------------------
  implicit none
  integer :: i, j, k, l
  integer, intent(in) :: n, ia, ib
  character(*), intent(in) :: name
  real(8), intent(in) :: a(n, n)
  write(6, '(/,1x,a)') name
  k = 1
  l = 8
  do while(k<=ib)
    l = min0(l, ib)
    write(6, '(6x,8i16)') (i, i = k, l)
    do j = 1, ia
      write(6, '(6x,i4,8f16.8)') j, (a(j, i), i = k, l)
    enddo
    write(6, *)
    k = k + 8
    l = l + 8
  enddo
  return
end subroutine outsqr
