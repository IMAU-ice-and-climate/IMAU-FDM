program toml_main
  use, intrinsic :: iso_fortran_env, only: error_unit
  use model_settings, only: read_constants, Constants, read_settings, Settings
  use tomlf

  implicit none

  integer :: fu, rc
  logical :: file_exists
  character(len=*), parameter :: constants_filename = 'settings/model_variables.toml'
  character(len=*), parameter :: settings_filename = 'settings/model_settings.toml'
  type(toml_table), allocatable :: table
  type(Constants), allocatable :: const
  type(Settings), allocatable :: config

  inquire(file=constants_filename, exist=file_exists)
  if (.not. file_exists) then
    write(error_unit, '("Error: TOML file ",a," not found")') trim(constants_filename)
    error stop 1
  end if

  open(action='read', file=constants_filename, iostat=rc, newunit=fu)
  if (rc /= 0) then
    write(error_unit, '("Error: opening TOML file ",a," failed; iostat=",i0)') trim(constants_filename), rc
    error stop 1
  end if

  call toml_parse(table, fu)
  close(fu)

  if (.not. allocated(table)) then
    write(error_unit, '("Error: TOML parsing failed for ",a)') trim(constants_filename)
    error stop 1
  end if

  call read_constants(table, const)
  deallocate(table)

  print *, const%g

  inquire(file=settings_filename, exist=file_exists)
  if (.not. file_exists) then
    write(error_unit, '("Error: TOML file ",a," not found")') trim(settings_filename)
    error stop 1
  end if

  open(action='read', file=settings_filename, iostat=rc, newunit=fu)
  if (rc /= 0) then
    write(error_unit, '("Error: opening TOML file ",a," failed; iostat=",i0)') trim(settings_filename), rc
    error stop 1
  end if

  call toml_parse(table, fu)
  close(fu)

  if (.not. allocated(table)) then
    write(error_unit, '("Error: TOML parsing failed for ",a)') trim(settings_filename)
    error stop 1
  end if

  call read_settings(table, config)

  print *, config%minimum_values%ts_minimum

end program toml_main