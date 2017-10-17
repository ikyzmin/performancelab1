	program multiplication
	use matr
	real,allocatable,dimension(:,:) :: a,b,res
	real,allocatable,dimension(:) :: veca,vecb,vecres
	integer:: lengths(8)
	lengths = (/2,8,10,100,200,500,1000,2000/)
	n = size (lengths)
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
	end program
