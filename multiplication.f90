	program multiplication
	use matr
	real,allocatable,dimension(:,:) :: a,b,res
	real,allocatable,dimension(:) :: veca,vecb,vecres
	integer:: lengths(8),block_lengths(8),lengths_for_blocked(8)
	lengths = (/4,8,10,100,200,500,1000,2000/)
	lengths_for_blocked = (/4,8,16,128,256,512,1024,2048/)
	block_lengths = (/2,4,8,32,64,128,256,512/)
	n = size (lengths)
	n_block = size(block_lengths)


	do i=1,n
		allocate(a(lengths(i),lengths(i)),b(lengths(i),lengths(i)),res(lengths(i),lengths(i)))
		call Populate(a)
		call Populate(b)
		call cpu_time(startTime)
		res = stdmatmul(a,b)
		call cpu_time(endTime)
		call print_results_formatted(lengths(i),endTime-startTime)
		deallocate(a,b,res)
	enddo 
	do i=1,n
		allocate(a(lengths(i),lengths(i)),b(lengths(i),lengths(i)),res(lengths(i),lengths(i)))
		call Populate(a)
		call Populate(b)
		call cpu_time(startTime)
		res = matmultranspose(a,b)
		call cpu_time(endTime)
		call print_results_formatted(lengths(i),endTime-startTime)
		deallocate(a,b,res)
	enddo
	do i=1,n
		allocate(veca(lengths(i)*lengths(i)),vecb(lengths(i)*lengths(i)),vecres(lengths(i)*lengths(i)))
		call Populate(veca)
		call Populate(vecb)
		call cpu_time(startTime)
		vecres = matasvecmul(veca,vecb)
		call cpu_time(endTime)
		call print_results_formatted(lengths(i),endTime-startTime,.TRUE.)
		deallocate(veca,vecb,vecres)
	enddo
	do i=1,n_block
		allocate(a(lengths_for_blocked(i),lengths_for_blocked(i)),b(lengths_for_blocked(i),lengths_for_blocked(i)),res(lengths_for_blocked(i),lengths_for_blocked(i)))
		call Populate(a)
		call Populate(b)
		
		call cpu_time(startTime)
		res = matmulblock(a,b,block_lengths(i))
		call cpu_time(endTime)
		call print_results_formatted(lengths_for_blocked(i),endTime-startTime)
		deallocate(a,b,res)
	enddo

	do i=1,n
		allocate(a(lengths(i),lengths(i)),b(lengths(i),lengths(i)),res(lengths(i),lengths(i)))
		call Populate(a)
		call Populate(b)
		call cpu_time(startTime)
		res = blasmul(a,b)
		call cpu_time(endTime)
		call print_results_formatted(lengths(i),endTime-startTime)
		deallocate(a,b,res)
	enddo
	end program
