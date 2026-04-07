program toml_main
    use model_settings, only: read_constants, Constants
    use tomlf

    implicit none

    integer :: fu, rc
    type(toml_table), allocatable :: table
    type(Constants), allocatable :: const

    open (action="read", file="settings/model_variables.toml", iostat=rc, newunit=fu)
    call toml_parse(table, fu)
    call read_constants(table, const)
    print *, const%g
end program toml_main