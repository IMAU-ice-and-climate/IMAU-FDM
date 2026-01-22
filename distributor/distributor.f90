program main
    use mpi
    implicit none

    double precision, dimension(60000) :: lat, lon
    double precision cur_lat, cur_lon, new_lat, new_lon
    integer io, i_row, cur_point, i_point, ierror, rank, job_request, n_proc_done, size, cur_i_point
    integer, dimension(60000) :: pointlist
    integer, dimension(MPI_Status_size) :: status
    double precision :: jobs_done = -99999999.0

    job_request = 42

    call MPI_Init(ierror)
    call MPI_Comm_size(MPI_Comm_World, size, ierror)
    call MPI_Comm_Rank(MPI_Comm_World, rank, ierror)

    ! TKTKTK read/set paths
        ! read in domain from toml
        ! domain = 
        ! pad_ref = "../reference/{domain}/IN_ll_{domain}.txt"
        ! pad_pointlist = "../pointlist/....txt"
        

    if (rank == 0) then

        ! Read coordinate file
        open (unit=20, file="../reference/FGRN055/IN_ll_FGRN055.txt", action="read")

        i_row = 1
        do
            read (20, *, iostat=io) cur_lat, cur_lon
            if (io/=0) exit
            lat(i_row) = cur_lat
            lon(i_row) = cur_lon
            i_row = i_row + 1
        end do
        close(20)

        ! Read point list
        open (unit=30, file="rundir/pointlists/pointlist_run-1957-2023-FGRN055-era055-2.txt", action="read")
        i_point = 1
        do
            read (30, *, iostat=io) pointlist(i_point)
            if (io/=0) exit
            i_point = i_point + 1
        end do
        close(30)
        i_point = i_point-1  ! The last point is 0, not sure, but should be skipped
    end if

    ! Worker threads
    if (rank /= 0) then
        do
            ! Send job request
            call MPI_send(job_request, 1, MPI_Integer, 0, 1, MPI_Comm_World, ierror)
            ! Receive new job
            call MPI_recv(cur_lat, 1, MPI_Double_Precision, 0, 2, MPI_Comm_World, status, ierror)
            call MPI_recv(cur_lon, 1, MPI_Double_Precision, 0, 2, MPI_Comm_World, status, ierror)
            print*, rank, ": doing job with lat: ", cur_lat, "lon: ", cur_lon
            ! call sleep(1)
            if (cur_lat == jobs_done) then  ! Exit the loop when no more jobs are available
                exit
            end if
        end do
    ! Distributor thread
    else
        n_proc_done = 0
        cur_i_point = 1
        do
            ! Check for any job requests and receive them.
            call MPI_Probe(MPI_Any_Source, 1, MPI_Comm_World, status, ierror)
            ! print*, "Received job from PID", status(MPI_Source)
            call MPI_recv(job_request, 1, MPI_Integer, status(MPI_Source), 1, MPI_Comm_World, status, ierror)
            if (cur_i_point > i_point) then  ! Send job_done message
                new_lat = jobs_done
                new_lon = jobs_done
                n_proc_done = n_proc_done + 1
            else  ! Send new point on the list
                new_lat = lat(pointlist(cur_i_point))
                new_lon = lon(pointlist(cur_i_point))
                cur_i_point = cur_i_point + 1
            end if
            call MPI_send(new_lat, 1, MPI_Double_Precision, status(MPI_Source), 2, MPI_Comm_World, ierror)
            call MPI_send(new_lon, 1, MPI_Double_Precision, status(MPI_Source), 2, MPI_Comm_World, ierror)
            if (n_proc_done == size-1) exit  ! Exit the loop when all job_done messages are sent
        end do
    end if

    call MPI_Finalize(ierror)
end program main