program distribute_points
implicit none
integer       :: nargs, iargc, zargs
character*255 :: pointlistfile
character*255 :: argv
character*255 :: path2request
logical       :: usetimeguess
character*255 :: settingsfile


integer, parameter :: ngridpointsmax = 60000
integer       :: pointlist(ngridpointsmax)
logical       :: ptodo(ngridpointsmax)
integer       :: expruntime(ngridpointsmax)
real          :: exprunhours
real          :: totalhours
integer       :: it, nthreads, mintreads
integer, allocatable:: p4thread(:)
logical, allocatable:: threadok(:)
integer       :: io_open, io_read, io_inline
integer       :: io_write, io_close
integer       :: npoints, ipoint
character*8   :: cpoint
character*5   :: cit
character*255 :: reqfile
character*10  :: request
integer       :: nsecleft
character*255 :: settingsline
integer       :: nabort, nover
integer       :: ic, is, ie, nc
integer       :: ntimesort, itact
integer, allocatable:: sortedlist(:)
integer       :: itodo


logical       :: lgoon
logical       :: lfatal
integer*8     :: nanosleep
integer       :: sleeperror

! the idea is that this program scans requests from the different threads
! commands are: in: "provide" -> provide a number to this thread
!               in: "nonew"   -> This thread stops working, last point is ready
!               in: "abort"   -> This thread stops working, current point is not ready    
!              out: "wait"    -> request calling script to stop 
!              out: "stop"    -> stop script 
!              out: <number>  -> the grid point to be done   

! threads here run from 1 to nthreads 
! threads outside (system, logfiles) run from 0 to nthreads-1

nanosleep    = 500000000
usetimeguess = .false.
lfatal       = .false.

nargs   = iargc()         ! determine number of command line args

if ( nargs < 3 ) then
  write(6,'(A)') 'distribute_points.x has as input:'
  write(6,'(A)') '1) file with list of gridpoints'
  write(6,'(A)') '2) number of threads running'
  write(6,'(A)') '3) path to requests'
  write(6,'(A)') 'Options: -t <settingsfile> (to read the expected duration)'
  write(0,'(A)') '           We assume a layout of lon, lat, <int>, <int>, duration (hr), y, x'
  stop
endif  

call getarg(2,argv)
read(argv(1:8),'(I5)') nthreads
write(6,'(A,I5,A)') 'DP will feed ',nthreads,' threads.'
allocate(p4thread(nthreads))
allocate(threadok(nthreads))
threadok = .true.
mintreads = nthreads/11
write(6,'(A,I5,A)') 'DP will terminate once less than ',mintreads,' are active.'

call getarg(3,path2request)

zargs = 4
do while ( zargs <= nargs )
  call getarg(zargs, argv)
  select case ( argv )
  case ('-t')
    usetimeguess = .true.
    zargs = zargs + 1
    call getarg(zargs, settingsfile)
  case default  
    write(0,'(2A)') 'DP: Unknown action provided: ',trim(argv)
  end select
  zargs = zargs + 1
enddo
    
!-------- read file with list point numbers ------ (of the settingsfile which is not nessairly provided here)
call getarg(1,pointlistfile)
open(111,file=pointlistfile,action='read',status='old',iostat=io_open)
if (io_open /= 0 ) then
  write(0,'(A,I8,A)') "DP: Error ",io_open," occurred during opening pointlistfile."
  write(0,'(2A)') "pointlistfile = ",trim(pointlistfile)
  stop
endif

io_read = 0
npoints = 0
do while (io_read == 0 .and. npoints <= ngridpointsmax )
 read(111, '(A8)', iostat = io_read) cpoint
 if ( io_read == 0 ) then
   npoints = npoints + 1
   read(cpoint,'(I8)') pointlist(npoints)
 endif  
enddo
close(111)

if ( io_read == 0 ) then
  write(0,'(A,I8,A)') 'DP: More than ',ngridpointsmax,&
  & ' points in file, adjust ngridpointsmax'
  stop
endif

!-------- read settingsfile if provided
if ( usetimeguess ) then
 open(111,file=settingsfile,action='read',status='old',iostat=io_open)
 if (io_open /= 0 ) then
  write(0,'(A,I8,A)') "DP: Error ",io_open," occurred during opening the settingsfile."
  write(0,'(2A)') "settingsfile = ",trim(settingsfile)
  stop
 endif
 
 ipoint = 0
 read(111,'(A)', iostat = io_read) settingsline
 do while (io_read == 0)
  ipoint = ipoint + 1
  nc     = 0
  ic     = 1
  do while ( nc < 5)
   if ( settingsline(ic:ic) .eq. "," ) then
    nc = nc + 1
    is = ie
    ie = ic
   endif
   ic = ic + 1
  enddo
  read(settingsline(is+2:ie-1),*, iostat=io_inline) exprunhours
  if ( io_inline /= 0 ) &
&   write(0,'(A,I6,3A)') "DP: Error while reading duration from line ",ipoint,", segment:",&
&     settingsline(is+2:ie-1),", line:",trim(settingsline)
  expruntime(ipoint)= int(exprunhours*3600.)
  
  read(111,'(A)', iostat = io_read) settingsline  
 enddo
 
 close (111)
 
! some info output
 totalhours = sum(expruntime(pointlist(1:npoints))/3600.)
 write(6,'(A,I6,A)')  'DP: ', ipoint, ' durations readed from settingsfile.'
 write(6,'(A,F8.1,A)')'DP: Total estimated calculation time is ',totalhours,' hours,'
 write(6,'(A,F8.1,A,I4,A)') 'DP: and ',totalhours/nthreads,' hours using ',&
 & nthreads,' threads.'
 
 expruntime(ipoint+1:) = 0
 
 ntimesort = 5*nthreads
 allocate(sortedlist(ntimesort))
 sortedlist = -1
 itact      = 0
else
 expruntime = -60 
endif


! ----------- start to be distributor
ptodo(:npoints)   = .true.
ptodo(npoints+1:) = .false.

ipoint = 1
do while ( count(threadok).ge.mintreads )
 lgoon = .true.
 do it=1,nthreads
! note that the ranks go from 0 to nthread-1  
    write(cit,'(I5.5)') it-1
    reqfile = trim(path2request) // "/" // cit
    
    open(112, file=reqfile, action='read', status='old', iostat=io_open)
    if ( io_open == 0 ) then
      read(112,'(A)', iostat = io_read) request
      if ( io_read /= 0 ) request = ""
      select case (trim(request))
      case ( 'provide')
       if ( usetimeguess) then
        read(112,'(I8)', iostat = io_inline) nsecleft
	close(112)
	if ( io_inline /= 0 ) write(0,'(A,I6,A)') 'DP: Error ',io_inline,' while reading nsecleft'
!	write(6,'(A,I4,A,I6,A)') 'Get a point for ',it,' with ', nsecleft, ' seconds left.'

	if ( ipoint <= npoints ) call update_list(&
&         ipoint, npoints, ntimesort, itact, sortedlist, &
&         ngridpointsmax, pointlist, expruntime)
        if ( itact >= 1 ) then
	  call get_point_list(ntimesort, itact, sortedlist, &
&           ngridpointsmax, pointlist, expruntime, nsecleft, itodo)
        else 
	  itodo = -1
	endif

	if ( itodo > 0 ) then
	  open(112, file=reqfile, action='write')
	  write(112,'(I8)') pointlist(itodo)
	  write(112,'(I8)') expruntime(pointlist(itodo))/60
	  close(112)
	  ptodo(itodo) = .false.
	  p4thread(it)  = pointlist(itodo)
          write(6,'(A,I4,A,I6,A,I4,A,I4,A)') &
	  & 'DP: Give thread ',it-1,' # ',pointlist(itodo),&
	  & ' , exprt ',expruntime(pointlist(itodo))/60, &
	  & ' m; ',nsecleft/60, ' m left.'
	  
	else
	  close(112)
	  open(112, file=reqfile, action='write')
	  write(112,'(A)') 'wait'
	  close(112)
	  threadok(it) = .false.
	  p4thread(it) = -1
	  write(6,'(A,I4,A,I4,A,I4,A)') &
	  & 'DP: No points for thread ',it-1,'; no points available suitable with ',&
	  & nsecleft/60, ' minutes left. ', count(threadok), ' threads still active.'
	
	endif  
       else
        close(112)
        if ( ipoint <= npoints ) then
	  open(112, file=reqfile, action='write')
	  write(112,'(I8)') pointlist(ipoint)
	  close(112)
	  ptodo(ipoint) = .false.
	  p4thread(it)  = pointlist(ipoint)
	  ipoint        = ipoint + 1
	   
	else
	  open(112, file=reqfile, action='write')
	  write(112,'(A)') 'wait'
	  close(112)
	  threadok(it) = .false.
	  p4thread(it) = -1
	endif
	lgoon = .false. 
       endif	 
      case ( 'abort' )
        close(112, status='delete')
        ! terminate
	threadok(it) = .false.
	lgoon        = .false.
	write(6,'(A,I4,A,I4)') 'DP: Thread ',it-1,' aborts, ', count(threadok), &
	& ' threads still active.'
	
      case ( 'nonew' )
	p4thread(it) = -1
        close(112, status='delete')
        ! terminate
	threadok(it) = .false.
	lgoon        = .false.

	write(6,'(A,I4,A,I4)') 'DP: Thread ',it-1,' do not nead a new point, ', count(threadok), &
	& ' threads still active.'

      case ( 'fatal' )
        close(112, status='delete')
        write(6,'(A,I4,A)') 'DP: Thread ',it-1,' sends fatal error - abort without restart.'
	threadok = .false. ! so simply stop
	lgoon    = .false.
	lfatal   = .true.
	
      case default ! = wait, stop, a number
        ! don't do anything, leave file there
	close(112)

      end select	
    ! else no file thus no request
    endif    
  enddo
  
  if ( lgoon ) then
    call ftnsleep(nanosleep, sleeperror)
    if ( sleeperror /= 0 ) write(0,'(A,I8)') 'DP: ftnsleep gave error ',sleeperror 
  endif  
enddo

! ----- end of normal operations
! notifiy the scripts they need to stop
! as the threads may write as well, some error-safety is needed
do it=1,nthreads
! note that the ranks go from 0 to nthread-1  
  write(cit,'(I5.5)') it-1
  reqfile = trim(path2request) // "/" // cit
  
  is = 0
  lgoon = .true.
  do while ( lgoon )
    open(112, file=reqfile, action='write', iostat=io_open)
    if ( io_open == 0 ) then
      write(112,'(A)', iostat=io_write) 'stop'
      close(112, iostat=io_close)
    endif
    if ( io_open==0 .and. io_write==0 .and. io_close==0 ) then
      lgoon = .false.
    elseif ( is == 2 ) then
      lgoon = .false.
    else  
      is = is + 1
    endif
  enddo       
enddo

! make the list of points yet to do
if ( any(ptodo) .or. any(p4thread > 1) ) then
  write(6,'(A)') 'distribute_points has still points to do.'
  nabort = 0
  nover  = 0
  open(111,file=pointlistfile,action='write')
  do it=1,nthreads
    if ( p4thread(it)>0 ) then
      write(111,'(I8)') p4thread(it)
      nabort = nabort + 1
    endif
  enddo    
  do it=1,npoints
    if ( ptodo(it) ) then
      write(111,'(I8)') pointlist(it)
      nover = nover + 1
    endif  
  enddo
  close(111)

  write(6,'(I6,A,I6,A,I6,A)') nabort,' points will be aborted, ', &
  & nover, ' points are not yet started. In total ',nabort+nover, &
  & ' points.'

  request = "continue"

else
  write(6,'(A)') 'distribute_points sees that everyting is completed!'
  request = "done"  

endif

if ( lfatal ) then
  write(6,'(A)') 'DP: Fatal error - no continuation'
  request="killed"
endif  

! write I'm ready
reqfile = trim(path2request) // "/DP"
open(113, file=reqfile, action="write")
write(113,'(A)') trim(request)
close(113)

write(6,'(A)') 'distribute_points is ready and teminates.'
end
  

!------------------------------------------------------------------
subroutine update_list(ipoint, npoint, nlist, ilist, sortedlist, nexpt, pointlist, expruntime)
implicit none
integer, intent(inout):: ipoint
integer, intent(in)   :: npoint
integer, intent(in)   :: nlist
integer, intent(inout):: ilist
integer, intent(inout):: sortedlist(nlist)
integer, intent(in)   :: nexpt
integer, intent(in)   :: pointlist(nexpt)
integer, intent(in)   :: expruntime(nexpt)

integer :: j, im
! sort from longest to shortest
do while (ilist < nlist .and. ipoint <= npoint )
 ! add
 ilist = ilist + 1
 sortedlist(ilist) = ipoint
 ipoint = ipoint + 1
 
 ! sort
 if ( ilist > 1 ) then
  j = ilist
  do while ( j>1 )
   if (expruntime(pointlist(sortedlist(j-1))) < &
   &   expruntime(pointlist(sortedlist(j))) ) then
   ! swap numbers
    im              = sortedlist(j-1)
    sortedlist(j-1) = sortedlist(j)
    sortedlist(j)   = im
    j               = j - 1
   else
    j               = -1
   endif  
  enddo
 endif
enddo

return
end


!-------------------------------------
subroutine get_point_list(nlist, ilist, sortedlist, &
& nexpt, pointlist, expruntime, nsecleft, itodo)
implicit none
integer, intent(in)   :: nlist
integer, intent(inout):: ilist
integer, intent(inout):: sortedlist(nlist)
integer, intent(in)   :: nexpt
integer, intent(in)   :: pointlist(nexpt)
integer, intent(in)   :: expruntime(nexpt)
integer, intent(in)   :: nsecleft
integer, intent(out)  :: itodo

integer :: jtodo
real    :: upfc, lwfc

upfc = 1.3
lwfc = 1.0

if ( nsecleft > expruntime(pointlist(sortedlist(1)))*upfc ) then
  itodo = sortedlist(1)
  sortedlist(1:ilist-1) = sortedlist(2:ilist)
  ilist = ilist - 1
elseif ( nsecleft < expruntime(pointlist(sortedlist(ilist)))*lwfc ) then
  ! not enough time
  itodo = -1
elseif ( nsecleft < expruntime(pointlist(sortedlist(ilist)))*upfc ) then
  ! find the largest one that is upfc of the remaining time
  jtodo = ilist
  do while ( nsecleft > expruntime(pointlist(sortedlist(jtodo-1)))*lwfc ) 
   jtodo = jtodo - 1
  enddo       
  itodo = sortedlist(jtodo)
  sortedlist(jtodo:ilist-1) = sortedlist(jtodo+1:ilist) 
  ilist = ilist -1
  
else
  ! find the largest one that is no more than half of the remaining time
  jtodo = ilist
  do while ( nsecleft > expruntime(pointlist(sortedlist(jtodo-1)))*upfc ) 
   jtodo = jtodo - 1
  enddo       
  itodo = sortedlist(jtodo)
  sortedlist(jtodo:ilist-1) = sortedlist(jtodo+1:ilist) 
  ilist = ilist -1

endif

return
end
