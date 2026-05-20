program main

    use mpi
    use model_main
    implicit none

    double precision, dimension(60000) :: lat, lon, ind_lats, ind_lons
    double precision :: cur_lat, cur_lon, new_lat, new_lon
    integer :: new_ind_lat, new_ind_lon
    integer :: io, i_row, i_point, ierror, rank, job_request, n_proc_done, size, cur_i_point, cur_int_lat, cur_int_lon
    integer :: recv_point_numb_int
    integer, dimension(60000) :: pointlist
    integer, dimension(MPI_Status_size) :: status
    double precision :: jobs_done = -99999999.0
    double precision :: ignore_ref1, ignore_ref2, ignore_ref3 !ignored
    character(len=512) :: ref_pointlist_path, work_pointlist_path
    character(len=512) :: point_log_dir, dist_log_dir
    character(len=512) :: log_fname, dist_log_fname
    integer :: dist_log_unit, point_log_unit, devnull_unit
    job_request = 42

    call MPI_Init(ierror)
    call MPI_Comm_size(MPI_Comm_World, size, ierror)
    call MPI_Comm_Rank(MPI_Comm_World, rank, ierror)

    call get_command_argument(1, work_pointlist_path)
    call get_command_argument(2, settings_dir)
    call get_command_argument(3, point_log_dir)
    call get_command_argument(4, dist_log_dir)

    ! All ranks suppress Define_Job startup output (paths, averages) from going to SLURM log
    open(newunit=devnull_unit, file='/dev/null', action='write', status='old')
    log_unit = devnull_unit

    call Define_Job()
    ref_pointlist_path = trim(reference_dir)//"IN_ll_"//trim(domain)//".txt"

    if (rank == 0) then

        ! Open distributor log and redirect rank 0 output there
        write(dist_log_fname, '(4A)') trim(dist_log_dir), "/distributor_", trim(project_name), ".log"
        open(newunit=dist_log_unit, file=trim(dist_log_fname), action='write', status='replace', iostat=io)
        if (io /= 0) then
            write(*, '(2A)') "ERROR: could not open distributor log: ", trim(dist_log_fname)
            call MPI_Abort(MPI_Comm_World, 1, ierror)
        end if
        log_unit = dist_log_unit
        write(dist_log_unit, '(2A)')   "Distributor started: ", trim(project_name)
        write(dist_log_unit, '(A,I0)') "MPI size: ", size
        flush(dist_log_unit)

        write(dist_log_unit, '(2A)') "Reading coordinate file: ", trim(REF_POINTLIST_PATH)
        flush(dist_log_unit)

        ! Read coordinate file
        open(unit=20, file=trim(REF_POINTLIST_PATH), action="read", iostat=io)
        if (io /= 0) then
            write(dist_log_unit, '(2A)') "ERROR: cannot open ref pointlist: ", trim(REF_POINTLIST_PATH)
            flush(dist_log_unit)
            call MPI_Abort(MPI_Comm_World, 1, ierror)
        end if

        i_row = 1
        do
            read (20, *, iostat=io) cur_lon, cur_lat, ignore_ref1, ignore_ref2, ignore_ref3, cur_int_lat, cur_int_lon
            if (io/=0) exit
            lat(i_row) = cur_lat
            lon(i_row) = cur_lon
            ind_lats(i_row) = cur_int_lat
            ind_lons(i_row) = cur_int_lon
            i_row = i_row + 1
        end do
        close(20)

        write(dist_log_unit, '(2A)') "Reading work pointlist: ", trim(work_pointlist_path)
        flush(dist_log_unit)

        ! Read point list
        open(unit=30, file=trim(work_pointlist_path), action="read", iostat=io)
        if (io /= 0) then
            write(dist_log_unit, '(2A)') "ERROR: cannot open work pointlist: ", trim(work_pointlist_path)
            flush(dist_log_unit)
            call MPI_Abort(MPI_Comm_World, 1, ierror)
        end if
        i_point = 1
        do
            read (30, *, iostat=io) pointlist(i_point)
            if (io/=0) exit
            i_point = i_point + 1
        end do
        close(30)
        i_point = i_point-1
        write(dist_log_unit, '(A,I0)') "Total points to run: ", i_point
        flush(dist_log_unit)

    end if

    ! Worker ranks
    if (rank /= 0) then
        do
            ! Send job request
            call MPI_send(job_request, 1, MPI_Integer, 0, 1, MPI_Comm_World, ierror)

            ! Receive new job
            call MPI_recv(cur_lat, 1, MPI_Double_Precision, 0, 2, MPI_Comm_World, status, ierror)
            call MPI_recv(cur_lon, 1, MPI_Double_Precision, 0, 2, MPI_Comm_World, status, ierror)
            call MPI_recv(cur_int_lat, 1, MPI_Integer, 0, 2, MPI_Comm_World, status, ierror)
            call MPI_recv(cur_int_lon, 1, MPI_Integer, 0, 2, MPI_Comm_World, status, ierror)
            call MPI_recv(recv_point_numb_int, 1, MPI_Integer, 0, 2, MPI_Comm_World, status, ierror)
            write(point_numb, '(I0)') recv_point_numb_int

            if (cur_lat == jobs_done) exit

            ! Open per-point log and redirect model output into it
            write(log_fname, '(4A,I0,A)') trim(point_log_dir), "/log_IMAU-FDM_", &
                trim(project_name), "_", recv_point_numb_int, ".out"
            open(newunit=point_log_unit, file=trim(log_fname), action='write', status='replace', iostat=io)
            if (io /= 0) then
                print *, "WARNING rank ", rank, ": could not open log file: ", trim(log_fname)
                log_unit = 6  ! fallback to stdout
            else
                log_unit = point_log_unit
            end if

            call Run_Model(cur_lat, cur_lon, cur_int_lat, cur_int_lon)

            if (io == 0) then
                flush(point_log_unit)
                close(point_log_unit)
            end if
            log_unit = devnull_unit  ! suppress output between points

        end do

        close(devnull_unit)

    ! Distributor rank
    else
        n_proc_done = 0
        cur_i_point = 1
        do
            ! Check for any job requests and receive them.
            call MPI_Probe(MPI_Any_Source, 1, MPI_Comm_World, status, ierror)
            call MPI_recv(job_request, 1, MPI_Integer, status(MPI_Source), 1, MPI_Comm_World, status, ierror)

            if (cur_i_point > i_point) then  ! No more jobs — send done signal
                new_lat = jobs_done
                new_lon = jobs_done
                n_proc_done = n_proc_done + 1
                write(dist_log_unit, '(A,I0,A)') "Rank ", status(MPI_Source), " <- done"
                flush(dist_log_unit)
            else  ! Send new point
                new_lat = lat(pointlist(cur_i_point))
                new_lon = lon(pointlist(cur_i_point))
                new_ind_lat = ind_lats(pointlist(cur_i_point))
                new_ind_lon = ind_lons(pointlist(cur_i_point))
                write(dist_log_unit, '(A,I0,A,I0,A,F8.4,A,F8.4)') &
                    "Rank ", status(MPI_Source), " <- point ", pointlist(cur_i_point), &
                    "  lat=", new_lat, "  lon=", new_lon
                flush(dist_log_unit)
                cur_i_point = cur_i_point + 1
            end if

            call MPI_send(new_lat, 1, MPI_Double_Precision, status(MPI_Source), 2, MPI_Comm_World, ierror)
            call MPI_send(new_lon, 1, MPI_Double_Precision, status(MPI_Source), 2, MPI_Comm_World, ierror)
            call MPI_send(new_ind_lat, 1, MPI_Integer, status(MPI_Source), 2, MPI_Comm_World, ierror)
            call MPI_send(new_ind_lon, 1, MPI_Integer, status(MPI_Source), 2, MPI_Comm_World, ierror)
            call MPI_send(pointlist(cur_i_point-1), 1, MPI_Integer, status(MPI_Source), 2, MPI_Comm_World, ierror)

            if (n_proc_done == size-1) exit
        end do

        write(dist_log_unit, '(A)') "All points dispatched. Distributor exiting."
        close(dist_log_unit)
    end if

    call MPI_Finalize(ierror)
end program main
