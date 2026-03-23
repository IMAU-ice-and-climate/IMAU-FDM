program main
    
    use mpi
    use model_main
    implicit none

    double precision, dimension(60000) :: lat, lon, ind_lats, ind_lons
    double precision :: cur_lat, cur_lon, new_lat, new_lon
    integer :: ind_lon, ind_lat, new_ind_lat, new_ind_lon
    integer :: io, i_row, cur_point, i_point, ierror, rank, job_request, n_proc_done, size, cur_i_point, cur_int_lat, cur_int_lon
    integer, dimension(60000) :: pointlist
    integer, dimension(MPI_Status_size) :: status
    double precision :: jobs_done = -99999999.0
    double precision :: ignore_ref1, ignore_ref2, ignore_ref3 !ignored
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
            read (20, *, iostat=io) cur_lat, cur_lon, ignore_ref1, ignore_ref2, ignore_ref3, cur_int_lat, cur_int_lon
            if (io/=0) exit
            lat(i_row) = cur_lat
            lon(i_row) = cur_lon
            ind_lats(i_row) = cur_int_lat
            ind_lons(i_row) = cur_int_lon
            i_row = i_row + 1
        end do
        close(20)

        ! Read point list
        open (unit=30, file="/ec/res4/scratch/nld4814/test-new-distributor/pointlist_1.txt", action="read")
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
            ! TKTK: can we also send ind_lat, ind_lon 
            !-> eventually we can remove function in FDM to use ind_lat, ind_lon to find true lat/lon but at the moment, need integer coordinates
            call MPI_recv(cur_lat, 1, MPI_Double_Precision, 0, 2, MPI_Comm_World, status, ierror)
            call MPI_recv(cur_lon, 1, MPI_Double_Precision, 0, 2, MPI_Comm_World, status, ierror)
            call MPI_recv(cur_int_lat, 1, MPI_Integer, 0, 2, MPI_Comm_World, status, ierror)
            call MPI_recv(cur_int_lon, 1, MPI_Integer, 0, 2, MPI_Comm_World, status, ierror)
            call MPI_recv(point_numb, 1, MPI_Integer, 0, 2, MPI_Comm_World, status, ierror)

            print*, rank, ": doing job with lat: ", cur_lat, "lon: ", cur_lon
            ! call sleep(1)
            if (cur_lat == jobs_done) then  ! Exit the loop when no more jobs are available
                exit
            end if

            write(point_numb, '(I0)') pointlist(cur_i_point)   ! converts integer to string
            call Run_Model(cur_lat, cur_lon, cur_int_lat, cur_int_lon)

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
                new_ind_lat = ind_lats(pointlist(cur_i_point))
                new_ind_lon = ind_lons(pointlist(cur_i_point))
                cur_i_point = cur_i_point + 1
            end if

            ! distributor sends lat, lon, point numb to worker
            call MPI_send(new_lat, 1, MPI_Double_Precision, status(MPI_Source), 2, MPI_Comm_World, ierror)
            call MPI_send(new_lon, 1, MPI_Double_Precision, status(MPI_Source), 2, MPI_Comm_World, ierror)
            call MPI_send(new_ind_lat, 1, MPI_Double_Precision, status(MPI_Source), 2, MPI_Comm_World, ierror)
            call MPI_send(new_ind_lon, 1, MPI_Double_Precision, status(MPI_Source), 2, MPI_Comm_World, ierror)
            call MPI_send(pointlist(cur_i_point), 1, MPI_Integer, status(MPI_Source), 2, MPI_Comm_World, ierror)

            if (n_proc_done == size-1) exit  ! Exit the loop when all job_done messages are sent
        end do
    end if

    call MPI_Finalize(ierror)
end program main