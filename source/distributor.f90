program main
    ! MPI distributor for IMAU-FDM.
    !
    ! Layout: rank 0 is the distributor; all other ranks are workers.
    !   - Distributor reads the reference pointlist (all grid points) and the work
    !     pointlist (subset to run), then sends one point at a time to whichever
    !     worker requests a job.
    !   - Workers run Run_Model() for each received point, writing output and a
    !     per-point log file, then immediately request the next job.
    !
    ! MPI message tags:
    !   tag 1 — worker → distributor: job request (integer, value=job_request)
    !   tag 2 — distributor → worker: point data (lat, lon, ind_lat, ind_lon, point number)
    !
    ! Sentinel value: distributor sends cur_lat=jobs_done to signal no more work.

    use mpi
    use model_main
    implicit none

    ! Reference pointlist arrays (rank 0 only): full grid coordinates and 1-based indices
    double precision, allocatable :: lat(:), lon(:)
    integer, allocatable :: ind_lats(:), ind_lons(:)

    ! Per-point coordinates sent from distributor to worker
    double precision :: cur_lat, cur_lon, new_lat, new_lon
    integer :: new_ind_lat, new_ind_lon

    ! Loop / MPI bookkeeping
    integer :: io, i_row, i_point, ierror, rank, job_request, n_proc_done, size, cur_i_point, cur_int_lat, cur_int_lon
    integer :: recv_point_numb_int, n_ref
    integer, allocatable :: pointlist(:)           ! 1-based row numbers into the reference pointlist
    integer, dimension(MPI_Status_size) :: status

    double precision :: jobs_done = -99999999.0    ! sentinel lat value signalling no more work
    double precision :: ignore_ref1, ignore_ref2, ignore_ref3  ! unused columns in reference pointlist
    character(len=512) :: ref_pointlist_path, work_pointlist_path
    character(len=512) :: point_log_dir, dist_log_dir
    character(len=512) :: log_fname, dist_log_fname
    integer :: dist_log_unit, point_log_unit, devnull_unit
    character(8)  :: ts_date
    character(10) :: ts_time
    job_request = 42

    call MPI_Init(ierror)
    call MPI_Comm_size(MPI_Comm_World, size, ierror)
    call MPI_Comm_Rank(MPI_Comm_World, rank, ierror)

    ! Command-line arguments (set by launch_job.sh):
    !   1: path to work pointlist
    !   2: path to settings directory (contains run.toml, model.toml, constants.toml)
    !   3: directory for per-point log files
    !   4: directory for the distributor log file
    call get_command_argument(1, work_pointlist_path)
    call get_command_argument(2, settings_dir)
    call get_command_argument(3, point_log_dir)
    call get_command_argument(4, dist_log_dir)

    ! Suppress Read_Job startup output on all ranks so it doesn't clutter the SLURM log
    open(newunit=devnull_unit, file='/dev/null', action='write', status='old')
    log_unit = devnull_unit

    call Read_Job()
    ref_pointlist_path = trim(reference_dir)//"IN_ll_"//trim(domain)//".txt"

    ! -----------------------------------------------------------------------
    ! Rank 0: read pointlists and open distributor log
    ! -----------------------------------------------------------------------
    if (rank == 0) then

        write(dist_log_fname, '(4A)') trim(dist_log_dir), "/distributor_", trim(project_name), ".log"
        open(newunit=dist_log_unit, file=trim(dist_log_fname), action='write', status='replace', iostat=io)
        if (io /= 0) then
            write(*, '(2A)') "ERROR: could not open distributor log: ", trim(dist_log_fname)
            call MPI_Abort(MPI_Comm_World, 1, ierror)
        end if
        log_unit = dist_log_unit
        call date_and_time(ts_date, ts_time)
        write(dist_log_unit, '(A)') "Started: "//ts_date(1:4)//"-"//ts_date(5:6)//"-"//ts_date(7:8)// &
                                    " "//ts_time(1:2)//":"//ts_time(3:4)//":"//ts_time(5:6)
        write(dist_log_unit, '(2A)')   "Distributor started: ", trim(project_name)
        write(dist_log_unit, '(A,I0)') "MPI size: ", size
        flush(dist_log_unit)

        ! --- Read reference pointlist (all grid points) ---
        ! Format: lon, lat, <3 ignored cols>, rlat_idx, rlon_idx
        ! rlat_idx/rlon_idx are 0-based (Python convention); +1 converts to 1-based Fortran indices.
        write(dist_log_unit, '(2A)') "Reading coordinate file: ", trim(REF_POINTLIST_PATH)
        flush(dist_log_unit)

        open(unit=20, file=trim(REF_POINTLIST_PATH), action="read", iostat=io)
        if (io /= 0) then
            write(dist_log_unit, '(2A)') "ERROR: cannot open ref pointlist: ", trim(REF_POINTLIST_PATH)
            flush(dist_log_unit)
            call MPI_Abort(MPI_Comm_World, 1, ierror)
        end if

        ! First pass: count rows so we can allocate
        n_ref = 0
        do
            read(20, *, iostat=io) cur_lon, cur_lat, ignore_ref1, ignore_ref2, ignore_ref3, cur_int_lat, cur_int_lon
            if (io /= 0) exit
            n_ref = n_ref + 1
        end do
        rewind(20)
        allocate(lat(n_ref), lon(n_ref), ind_lats(n_ref), ind_lons(n_ref))

        ! Second pass: store coordinates and convert indices to 1-based
        do i_row = 1, n_ref
            read(20, *) cur_lon, cur_lat, ignore_ref1, ignore_ref2, ignore_ref3, cur_int_lat, cur_int_lon
            lat(i_row)      = cur_lat
            lon(i_row)      = cur_lon
            ind_lats(i_row) = cur_int_lat + 1   ! 0-based → 1-based
            ind_lons(i_row) = cur_int_lon + 1
        end do
        close(20)

        ! --- Read work pointlist (subset of points to run) ---
        ! Each line is a 1-based row number into the reference pointlist.
        write(dist_log_unit, '(2A)') "Reading work pointlist: ", trim(work_pointlist_path)
        flush(dist_log_unit)

        open(unit=30, file=trim(work_pointlist_path), action="read", iostat=io)
        if (io /= 0) then
            write(dist_log_unit, '(2A)') "ERROR: cannot open work pointlist: ", trim(work_pointlist_path)
            flush(dist_log_unit)
            call MPI_Abort(MPI_Comm_World, 1, ierror)
        end if

        ! First pass: count points
        i_point = 0
        do
            read(30, *, iostat=io) cur_int_lat  ! dummy variable; only used to count lines
            if (io /= 0) exit
            i_point = i_point + 1
        end do
        rewind(30)
        allocate(pointlist(i_point))

        ! Second pass: store point numbers
        do i_row = 1, i_point
            read(30, *) pointlist(i_row)
        end do
        close(30)
        write(dist_log_unit, '(A,I0)') "Total points to run: ", i_point
        flush(dist_log_unit)

    end if

    ! -----------------------------------------------------------------------
    ! Worker ranks (rank /= 0): request and run points until done signal
    ! -----------------------------------------------------------------------
    if (rank /= 0) then
        do
            ! Request next job from distributor
            call MPI_send(job_request, 1, MPI_Integer, 0, 1, MPI_Comm_World, ierror)

            ! Receive point data: lat, lon, grid indices, point number
            call MPI_recv(cur_lat, 1, MPI_Double_Precision, 0, 2, MPI_Comm_World, status, ierror)
            call MPI_recv(cur_lon, 1, MPI_Double_Precision, 0, 2, MPI_Comm_World, status, ierror)
            call MPI_recv(cur_int_lat, 1, MPI_Integer, 0, 2, MPI_Comm_World, status, ierror)
            call MPI_recv(cur_int_lon, 1, MPI_Integer, 0, 2, MPI_Comm_World, status, ierror)
            call MPI_recv(recv_point_numb_int, 1, MPI_Integer, 0, 2, MPI_Comm_World, status, ierror)
            write(point_numb, '(I0)') recv_point_numb_int

            ! Check for done sentinel
            if (cur_lat == jobs_done) exit

            ! Open per-point log file
            write(log_fname, '(4A,I0,A)') trim(point_log_dir), "/log_IMAU-FDM_", &
                trim(project_name), "_", recv_point_numb_int, ".out"
            open(newunit=point_log_unit, file=trim(log_fname), action='write', status='replace', iostat=io)
            if (io /= 0) then
                print *, "WARNING rank ", rank, ": could not open log file: ", trim(log_fname)
                log_unit = 6  ! fallback to stdout
            else
                log_unit = point_log_unit
                call date_and_time(ts_date, ts_time)
                write(log_unit, '(A)') "Started: "//ts_date(1:4)//"-"//ts_date(5:6)//"-"//ts_date(7:8)// &
                                       " "//ts_time(1:2)//":"//ts_time(3:4)//":"//ts_time(5:6)
            end if

            call Run_Model(cur_lat, cur_lon, cur_int_lat, cur_int_lon)

            if (io == 0) then
                call date_and_time(ts_date, ts_time)
                write(log_unit, '(A)') "Finished: "//ts_date(1:4)//"-"//ts_date(5:6)//"-"//ts_date(7:8)// &
                                       " "//ts_time(1:2)//":"//ts_time(3:4)//":"//ts_time(5:6)
                flush(point_log_unit)
                close(point_log_unit)
            end if
            log_unit = devnull_unit  ! suppress output between points

        end do

        close(devnull_unit)

    ! -----------------------------------------------------------------------
    ! Distributor rank (rank == 0): dispatch points to workers as they become free
    ! -----------------------------------------------------------------------
    else
        n_proc_done = 0
        cur_i_point = 1
        do
            ! Wait for any worker to request a job
            call MPI_Probe(MPI_Any_Source, 1, MPI_Comm_World, status, ierror)
            call MPI_recv(job_request, 1, MPI_Integer, status(MPI_Source), 1, MPI_Comm_World, status, ierror)

            if (cur_i_point > i_point) then  ! All points dispatched — send done sentinel
                new_lat = jobs_done
                new_lon = jobs_done
                n_proc_done = n_proc_done + 1
                write(dist_log_unit, '(A,I0,A)') "Rank ", status(MPI_Source), " <- done"
                flush(dist_log_unit)
            else  ! Send next point
                new_lat     = lat(pointlist(cur_i_point))
                new_lon     = lon(pointlist(cur_i_point))
                new_ind_lat = ind_lats(pointlist(cur_i_point))
                new_ind_lon = ind_lons(pointlist(cur_i_point))
                write(dist_log_unit, '(A,I0,A,I0,A,F8.4,A,F8.4)') &
                    "Rank ", status(MPI_Source), " <- point ", pointlist(cur_i_point), &
                    "  lat=", new_lat, "  lon=", new_lon
                flush(dist_log_unit)
                cur_i_point = cur_i_point + 1
            end if

            ! Send point data (workers ignore indices when lat==jobs_done)
            call MPI_send(new_lat,     1, MPI_Double_Precision, status(MPI_Source), 2, MPI_Comm_World, ierror)
            call MPI_send(new_lon,     1, MPI_Double_Precision, status(MPI_Source), 2, MPI_Comm_World, ierror)
            call MPI_send(new_ind_lat, 1, MPI_Integer,          status(MPI_Source), 2, MPI_Comm_World, ierror)
            call MPI_send(new_ind_lon, 1, MPI_Integer,          status(MPI_Source), 2, MPI_Comm_World, ierror)
            call MPI_send(pointlist(cur_i_point-1), 1, MPI_Integer, status(MPI_Source), 2, MPI_Comm_World, ierror)

            if (n_proc_done == size-1) exit  ! All workers have received the done sentinel
        end do

        write(dist_log_unit, '(A)') "All points dispatched. Distributor exiting."
        call date_and_time(ts_date, ts_time)
        write(dist_log_unit, '(A)') "Finished: "//ts_date(1:4)//"-"//ts_date(5:6)//"-"//ts_date(7:8)// &
                                    " "//ts_time(1:2)//":"//ts_time(3:4)//":"//ts_time(5:6)
        close(dist_log_unit)
    end if

    call MPI_Finalize(ierror)
end program main
