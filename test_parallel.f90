program parallel_test
    use omp_lib
    integer :: omp_rank,i
    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14,300)
    real(kind=dp):: start, finish
    INTEGER :: partial_Sum, total_Sum

    ! call cpu_time(start)

    start=omp_get_wtime()
    partial_Sum = 0;
    total_Sum = 0;
    DO i=1,1000000
        partial_Sum=partial_Sum+i
    end do    
    total_Sum = total_Sum + partial_Sum
    finish=omp_get_wtime()     
    write(*,*) 'total sum:', total_Sum
    write(*,*) 'Execution time in seconds for single= ', (finish-start)
    
    start=omp_get_wtime()
    !$omp parallel private (i,partial_sum) shared(total_sum)
        partial_Sum = 0;
        total_Sum = 0;
        omp_rank=omp_get_max_threads()
        ! write(*,*) omp_rank,omp_get_thread_num()
        !$omp  do
            DO i=1,1000000
                ! write(*,*) 'i*2=', i*2, ' on thread ', omp_rank
                partial_Sum=partial_Sum+i
                ! if (omp_get_thread_num()==1) then
                !     write(*,*) 'loop:', omp_get_thread_num(), i
                ! end if
            end do    
        !$omp end do  
    
        !!$OMP CRITICAL
            total_Sum = total_Sum + partial_Sum
        !!$OMP END CRITICAL
    
    !$omp end parallel 
    finish=omp_get_wtime()     
            ! call cpu_time(finish)  
    write(*,*) 'total sum:', total_Sum

    write(*,*) 'Execution time in seconds for parallel= ', (finish-start)
end program parallel_test