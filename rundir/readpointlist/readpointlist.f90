program readpointlist
implicit none
integer       :: nargs, zargs
integer       :: iaction
character*255 :: pointlistfile
character*255 :: argv
integer, parameter :: ngridpointsmax = 60000
integer       :: pointlist(ngridpointsmax)
integer       :: io_open, io_read
integer       :: iargc
integer       :: npoints, ipoint
character*8   :: cpoint

iaction = 0
nargs   = iargc()         ! determine number of command line args

if ( nargs == 0 ) then
  write(0,'(A)') 'readpointlist.x has as input:'
  write(0,'(A)') '1) file with list of gridpoints'
  write(0,'(A)') '2 optional) action'
  write(0,'(A)') '  default action is read list, prompt first and rewrite list without first.'
  write(0,'(A)') ' -count: read list and prompt number of entries.'
  stop
endif  

call getarg(1,pointlistfile)
open(111,file=pointlistfile,action='read',status='old',iostat=io_open)
if (io_open /= 0 ) then
  write(0,'(A,I8,A)') "Error ",io_open," occurred during opening pointlistfile."
  write(0,'(2A)') "pointlistfile = ",trim(pointlistfile)
  stop
endif

zargs = 2
do while ( zargs <= nargs )
  call getarg(zargs, argv)
  select case ( argv )
  case ('-count')
    iaction = 1
  case default  
    write(0,'(2A)') 'Unknown action provided: ',trim(argv)
  end select
  zargs = zargs + 1
enddo
    
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
  write(0,'(A,I8,A)') 'More than ',ngridpointsmax,&
  & ' points in file, adjust ngridpointsmax'
  stop
endif

if ( iaction == 0 ) then 
  ! write point-to-do to screen
  print *, pointlist(1)

  ! rewrite pointlist file with one less
  open(111,file=pointlistfile,action='write')
  if ( npoints > 1 ) then
    do ipoint=2,npoints
      write(111,'(I8)') pointlist(ipoint)
    enddo
  else
    write(111,'(A)') 'Ready'
  endif
  close(111)  

elseif ( iaction == 1 ) then
  ! write number of points to screen
  print *, npoints

else
  write(0,'(A,I6)') 'Unknown action: ',iaction
endif

end
  

