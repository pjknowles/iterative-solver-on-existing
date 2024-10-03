program test_solver
  use iso_fortran_env, only: dp => real64
  implicit double precision (a-h,o-z)
  real(dp), allocatable :: rhs(:,:)
  real(dp), allocatable :: diag(:)
  real(dp), allocatable :: A(:,:)
  real(dp), allocatable :: unit(:,:)
  integer, allocatable :: ipiv(:)
  integer :: i,nvec,ivec


  !size of matrix
  n=10

  !rhs vectors
  nvec=4
  allocate(rhs(n,nvec))
  call rhs_vectors(rhs,n,nvec)

  !matrix diagonal
  allocate(diag(n))
  call A_diagonal(diag,n)
  call outsqr(diag,n,n,1,'diagonal')

  !iterative solver stuff to be inserted here
  !dimension:       n
  !number of rhs's: nvec  
  !rhs:             corresponding variable
  !preconditioner:  diag holds the matrix diagonal
  !matrix*vector:   A_times_x subroutine
  



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

  deallocate(rhs,diag)
end program test_solver



!---------------------------------------------------------------------
subroutine A_times_x(Ax,x,n,nvec)
!---------------------------------------------------------------------
  !apply A.x on nvec many vectors
  use iso_fortran_env, only: dp => real64
  implicit none
  integer, intent(in) :: n,nvec
  real(dp), intent(in) :: x(n,nvec)
  real(dp), intent(out) :: Ax(n,nvec)
  integer :: i,j,ivec
  real(dp) :: aij

  Ax(:,:)=0d0
  do ivec=1,nvec
     do i=1,n
        do j=max(i-2,1),min(i+2,n)
           aij=10d0*dble(i)+dble(j)
           Ax(i,ivec)=Ax(i,ivec)+aij*x(j,ivec)
        enddo
     enddo
  enddo

  return
end subroutine A_times_x


!---------------------------------------------------------------------
subroutine A_diagonal(diag,n)
!---------------------------------------------------------------------
  !get matrix diagonal
  use iso_fortran_env, only: dp => real64
  implicit none
  integer, intent(in) :: n
  real(dp), intent(out) :: diag(n)
  real(dp), allocatable :: A(:,:),unit(:,:)
  integer :: i

  !just for convenience, we compute the whole matrix and then extract the
  !diagonal
  allocate(A(n,n))
  allocate(unit(n,n))
  unit(:,:)=0d0
  do i=1,n
     unit(i,i)=1d0
  enddo
  call A_times_x(A,unit,n,n)
  
  do i=1,n
     diag(i)=A(i,i)
  enddo
  
  deallocate(A,unit)

  return
end subroutine A_diagonal



!---------------------------------------------------------------------
subroutine rhs_vectors(rhs,n,nvec)
!---------------------------------------------------------------------
  use iso_fortran_env, only: dp => real64
  implicit none
  integer, intent(in) :: n,nvec
  real(dp), intent(out) :: rhs(n,nvec)
  integer :: ivec,i

  rhs(:,:)=0d0
  do ivec=1,nvec
!     do i=max(i-3,1),min(i+3,n) ! I don't understand this
     do i=max(ivec-3,1),min(ivec+3,n)
        rhs(i,ivec)=dble(i)
     enddo
  enddo

  return
end subroutine rhs_vectors


!---------------------------------------------------------------------
subroutine outsqr(a,n,ia,ib,name)
!---------------------------------------------------------------------
  implicit none
  integer :: i,j,k,l
  integer, intent(in)      :: n,ia,ib
  character(*), intent(in) :: name
  real(8), intent(in) :: a(n,n)

  write(6,'(/,1x,a)') name
  k=1
  l=8
  do while(k<=ib)
     l=min0(l,ib)
     write(6,'(6x,8i16)') (i,i=k,l)
     do j=1,ia
        write(6,'(6x,i4,8f16.8)') j,(a(j,i),i=k,l)
     enddo
     write(6,*)
     k=k+8
     l=l+8
  enddo

  return
end subroutine outsqr
