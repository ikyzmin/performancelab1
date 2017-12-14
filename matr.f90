module matr

	interface Populate
		module procedure rand_int
		module procedure rand_int_vec
		module procedure rand_real
		module procedure rand_real_vec
	end interface

	interface APrint
		module procedure print_matrix
		module procedure print_real_matrix
		module procedure print_vec
		module procedure print_real_vec	
	end interface

	interface PrintVecAsMat
		module procedure print_vec_as_mat_real
		module procedure print_vec_as_mat_int
	end interface 
	
	contains
		
		subroutine rand_int(a)
			integer,dimension(:,:),allocatable,intent(inout) :: a
			real :: u
			integer :: n,m
			n = size(a,1)
			do i = 1,n
				do j=1,n
					call RANDOM_NUMBER(u)
					a(i,j) = FLOOR((100)*u)+1
				enddo
			enddo
		end subroutine rand_int

		subroutine rand_int_vec(veca)
			integer,dimension(:),allocatable,intent(inout) :: veca
 			real :: u
			integer :: n
			n = size(veca)
			do i = 1,n
				call RANDOM_NUMBER(u)
				veca(i) = FLOOR((100)*u)+1
			enddo
		end subroutine rand_int_vec


		subroutine rand_real_vec(veca)
			real,dimension(:),allocatable,intent(inout) :: veca
 			real :: u
			integer :: n
			n = size(veca)
			do i = 1,n
				call RANDOM_NUMBER(veca(i))
			enddo
		end subroutine rand_real_vec

		subroutine rand_real(a)
			real,dimension(:,:),allocatable,intent(inout) :: a
			real :: u
			integer :: n,m
			n = size(a,1)
			do i = 1,n
				do j=1,n
					call RANDOM_NUMBER(a(i,j))
				enddo
			enddo
		end subroutine rand_real




		subroutine print_matrix(a)
			integer,dimension(:,:),allocatable,intent(in) :: a
 			integer :: length
			length = size(a,1)
			do i=1,length
				print *,(a(i,j), j=1,length)
			enddo
			print *,""
		end subroutine print_matrix

		subroutine print_real_matrix(r_a)
			real,dimension(:,:),allocatable,intent(in) :: r_a
 			integer :: length
			length = size(r_a,1)
			do i=1,length
				print *,(r_a(i,j), j=1,length)
			enddo
			print *,""
		end subroutine print_real_matrix
		
		subroutine print_real_vec(r_veca)
			real,dimension(:),allocatable,intent(in) :: r_veca
 			integer :: length
			length = size(r_veca)
			do i=1,length
				print *,r_veca(i)
			enddo
			print *,""	
		end subroutine print_real_vec

		subroutine print_vec(veca)
			integer,dimension(:),allocatable,intent(in) :: veca
 			integer :: length
			length = size(veca)
			do i=1,length
				print *,veca(i)
			enddo
			print *,""
		end subroutine print_vec

		subroutine print_vec_as_mat_real(matasvec)
			real,dimension(:),allocatable,intent(in) :: matasvec
 			integer :: length
			length = size(matasvec)
			do i=1,length
				do j=1,length
					print *,matasvec(i*length+j)
				enddo
				print *,""
			enddo
			print *,""
		end subroutine print_vec_as_mat_real

		subroutine print_vec_as_mat_int(matasvec)
			integer,dimension(:),allocatable,intent(in) :: matasvec
 			integer :: length
			length = size(matasvec)
			do i=1,length
				do j=1,length
					print *,matasvec(i*length+j)
				enddo
				print *,""
			enddo
			print *,""
		end subroutine print_vec_as_mat_int

		subroutine print_results_formatted(n,time,isvector)
			integer,intent(in) :: n
			real,intent(in)::time
			logical,optional,intent(in):: isvector
			if (present(isvector)) then
				print '("vector(",i8,") multiplication takes ",f9.4," seconds")',n*n,time
			else
				print '(i4,"x",i4" matrix multiplication takes ",f9.4," seconds")',n,n,time
			endif
		end subroutine print_results_formatted

		function matr_mul(a,b)
			integer,dimension(:,:),allocatable,intent(in) :: a,b
			integer,dimension(:,:),allocatable :: matr_mul
			allocate(matr_mul(n,n))
			matr_mul = matmul(a,b)
		end function matr_mul

		function matr_mul_dot(a,b)
			integer,dimension(:,:),allocatable,intent(in) :: a,b
			integer,dimension(:,:),allocatable :: matr_mul_dot
			n  = size (a,1)
			allocate(matr_mul_dot(n,n))
			do i=1,n
				do j=1,n
					matr_mul_dot(i,j) = dot_product(a(i,:),b(:,j))
				enddo
			enddo
		end function matr_mul_dot

		function stdmatmul(a,b)
			real,dimension(:,:),allocatable,intent(in) :: a,b
			real,dimension(:,:),allocatable :: stdmatmul
			n  = size (a,1)
			allocate(stdmatmul(n,n))
			do i=1,n
				do j=1,n
					do k=1,n
 						stdmatmul(i,j) = stdmatmul(i,j) + a(i,k)*b(k,j)
					enddo
				enddo
			enddo
		end function stdmatmul

		function stdmatmulnonallocatable(a,b)
			real,dimension(:,:),intent(in) :: a,b
			real,dimension(:,:),allocatable :: stdmatmulnonallocatable
			n  = size (a,1)
			allocate(stdmatmulnonallocatable(n,n))
			do i=1,n
				do j=1,n
					do k=1,n
 						stdmatmulnonallocatable(i,j) = stdmatmulnonallocatable(i,j) + a(i,k)*b(k,j)
					enddo
				enddo
			enddo
		end function stdmatmulnonallocatable

		function matasvecmul(a,b)
			real,dimension(:),allocatable,intent(in) :: a,b
			real,dimension(:),allocatable :: matasvecmul
			n  = size (a)
			allocate(matasvecmul(n))
			do i=1,length
				do j=1,length
					do k=1,length
 						matasvecmul(i*n+j) = matasvecmul(i*n+j) + a(i*n+k)*b(k*n+j)
					enddo
				enddo
			enddo
		end function matasvecmul

		function matmultranspose(a,b)
			real,dimension(:,:),allocatable,intent(in) :: a,b
			real,dimension(:,:),allocatable :: matmultranspose,at
			n  = size (a,1)
			allocate(matmultranspose(n,n),at(n,n))
			at = transpose(a)
			do i=1,length
				do j=1,length
					do k=1,length
 						matmultranspose(i,j) = matmultranspose(i,j) + at(k,i)*b(k,j)
					enddo
				enddo
			enddo
		end function matmultranspose

		function matmulblock(a,b,block_size)
			real,dimension(:,:),allocatable :: a,b
			real,dimension(:,:),allocatable :: matmulblock
			integer :: block_size , n,block_count
			n  = size (a,1)
			allocate(matmulblock(n,n))
			blocks_count = n / block_size
			do i=1,n-block_size+1,block_size 
				do j =1,n-block_size+1,block_size  
					do k=1,n-block_size+1,block_size 
						matmulblock(i:i+block_size-1,j:j+block_size-1)=matmulblock(i:i+block_size-1,j:j+block_size-1)+stdmatmulnonallocatable(a(j:j+block_size-1,k:k+block_size-1),b(k:k+block_size-1,j:j+block_size-1))
					enddo
				enddo		
			enddo
		end function matmulblock

		function blasmul(a,b)
			real,dimension(:,:),allocatable :: blasmul
			real,dimension(:,:),allocatable :: a,b
			character :: TRANSA,TRANSB
			integer :: M,N,K,LDA,LDB
			real :: ALPHA,BETA
			allocate(blasmul(n,n))
			ALPHA =1.0
			BETA =0.0
			TRANSA = 'N'
			TRANSB = 'N'
			M=size(a,1)
			N=M
			K=M
			call SGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,N,B,N,BETA,blasmul,N)
			
		end function blasmul

end module matr

