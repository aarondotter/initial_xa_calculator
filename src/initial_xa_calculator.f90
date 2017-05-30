program initial_xa_calculator

  use const_def
  use const_lib, only: const_init
  use chem_def
  use chem_lib
  use net_def
  use net_lib
  use rates_lib, only: rates_init

  implicit none

  logical, parameter :: verbose = .true.
  logical, parameter :: do_OPAL = .false.
  logical, parameter :: do_MESA = .true.
  logical, parameter :: use_ATLAS_atm = .false.
  logical, parameter :: do_XYZ = .true.
  logical :: scale_Y
  !for the mass fractions, etc.
  real(dp) :: log_num_frac(num_chem_elements)
  real(dp) :: mass_frac(num_chem_elements), num_frac(num_chem_elements)
  real(dp) :: alpha_div_Fe = 0.0d0 
  real(dp) :: Fe_div_H = 0.0d0
  real(dp) :: X_div_H(num_chem_elements)
  real(dp) :: dY_div_dZ = 0d0, Z_div_X, X, Y, Z
  real(dp), parameter :: Y_BBN = 0.249d0
  integer :: el, i, id, ierr, j, k, handle, species
  integer, pointer :: chem_id(:), net_iso(:)

  !for the net
  integer :: current, num_net_elements

  !file names
  character (len=128) :: net_file, atm_file, arg

  !for the net
  type(Net_General_info), pointer :: g      

  !new element type
  type element
     integer :: id, num_isos
     integer, allocatable :: isos(:)
     character(len=iso_name_length), allocatable :: iso_names(:)
     real(dp), allocatable :: iso_fracs(:)
  end type element
  type(element), allocatable :: e(:)


  if(command_argument_count() < 2) then
     write(0,*) ' usage:'
     write(0,*) '   ./initial_xa_calculator [net_name] [Fe/H] [alpha/Fe] {dY/dZ}'
     write(0,*) '    net_name = mesa/net network '
     write(0,*) '    [Fe/H] = 0 for solar'
     write(0,*) '    [alpha/Fe] = 0 for solar'
     write(0,*) '    dY/dZ = 1.5 by default, OPTIONAL'
     stop
  else
     call get_command_argument(1,net_file)
     if(use_ATLAS_atm) then
        call get_command_argument(2,atm_file)
     else
        call get_command_argument(2,arg)
        read(arg,*) Fe_div_H
        call get_command_argument(3,arg)
        read(arg,*) alpha_div_Fe
        if(command_argument_count()==4)then
           call get_command_argument(4,arg)
           read(arg,*) dY_div_dZ
        else
           dY_div_dZ = 1.5d0
        endif
     endif
  endif

  !initialize MESA modules
  ierr = 0
  call initialize_mesa('', ierr)
  if (ierr /= 0) stop '   oh no!'         
  call setup_net(net_file, handle, species, chem_id, net_iso, ierr)
  if (ierr /= 0) stop '   uh oh!'
  call get_net_ptr(handle, g, ierr)

  if(use_ATLAS_atm)then
     !read A09 file, adjust ratios, compute number and mass fractions
     call init_atlas_abunds
  else     
     call set_abunds(Fe_div_H, alpha_div_Fe)
  endif

  !determine number and names of elements in net g
  call net_elements_and_isotopes


  if(do_OPAL) call output_for_OPAL

  if(do_MESA) call output_for_MESA

  if(do_XYZ) call output_XYZ


contains

  subroutine set_abunds(Fe_div_H, alpha_div_Fe)
    real(dp), intent(in) :: Fe_div_H, alpha_div_Fe
    real(dp) :: Zold, Xold, Yold
    log_num_frac = -20d0
    log_num_frac(e_H ) = 12.0
    log_num_frac(e_He) = 10.93
    ! first store log abundances from the paper
    log_num_frac(e_li) = 3.26
    log_num_frac(e_be) = 1.38
    log_num_frac(e_b ) = 2.70
    log_num_frac(e_c ) = 8.43
    log_num_frac(e_n ) = 7.83
    log_num_frac(e_o ) = 8.69
    log_num_frac(e_f ) = 4.56
    log_num_frac(e_ne) = 7.93
    log_num_frac(e_na) = 6.24
    log_num_frac(e_mg) = 7.60
    log_num_frac(e_al) = 6.45
    log_num_frac(e_si) = 7.51
    log_num_frac(e_p ) = 5.41
    log_num_frac(e_s ) = 7.12
    log_num_frac(e_cl) = 5.50
    log_num_frac(e_ar) = 6.40
    log_num_frac(e_k ) = 5.03
    log_num_frac(e_ca) = 6.34
    log_num_frac(e_sc) = 3.15
    log_num_frac(e_ti) = 4.95
    log_num_frac(e_v ) = 3.93
    log_num_frac(e_cr) = 5.64
    log_num_frac(e_mn) = 5.43
    log_num_frac(e_fe) = 7.50
    log_num_frac(e_co) = 4.99
    log_num_frac(e_ni) = 6.22
    log_num_frac(e_cu) = 4.19
    log_num_frac(e_zn) = 4.56
    log_num_frac(e_ga) = 3.04
    log_num_frac(e_ge) = 3.65
    log_num_frac(e_as) = 2.30 !meteor
    log_num_frac(e_se) = 3.34 !meteor
    log_num_frac(e_br) = 2.54 !meteor
    log_num_frac(e_kr) = 3.25 !indirect
    log_num_frac(e_rb) = 2.52
    log_num_frac(e_sr) = 2.87
    log_num_frac(e_y ) = 2.21
    log_num_frac(e_zr) = 2.58
    log_num_frac(e_nb) = 1.46
    log_num_frac(e_mo) = 1.88
    log_num_frac(e_Ru) = 1.75
    log_num_frac(e_Rh) = 0.91
    log_num_frac(e_Pd) = 1.57
    log_num_frac(e_Ag) = 0.94
    log_num_frac(e_Cd) = 1.71
    log_num_frac(e_In) = 0.80
    log_num_frac(e_Sn) = 2.04
    log_num_frac(e_Sb) = 1.01
    log_num_frac(e_Te) = 2.18
    log_num_frac(e_I ) = 1.55
    log_num_frac(e_Xe) = 2.24
    log_num_frac(e_Cs) = 1.08
    log_num_frac(e_Ba) = 2.18
    log_num_frac(e_La) = 1.10
    log_num_frac(e_Ce) = 1.58
    log_num_frac(e_Pr) = 0.72
    log_num_frac(e_Nd) = 1.42
    log_num_frac(e_Sm) = 0.96
    log_num_frac(e_Eu) = 0.52
    log_num_frac(e_Gd) = 1.07
    log_num_frac(e_Tb) = 0.30
    log_num_frac(e_Dy) = 1.10
    log_num_frac(e_Ho) = 0.48
    log_num_frac(e_Er) = 0.92
    log_num_frac(e_Tm) = 0.10
    log_num_frac(e_Yb) = 0.84
    log_num_frac(e_Lu) = 0.10
    log_num_frac(e_Hf) = 0.85
    log_num_frac(e_Ta) =-0.12
    log_num_frac(e_W ) = 0.85
    log_num_frac(e_Re) = 0.26
    log_num_frac(e_Os) = 1.40
    log_num_frac(e_Ir) = 1.38
    log_num_frac(e_Pt) = 1.62
    log_num_frac(e_Au) = 0.92
    log_num_frac(e_Hg) = 1.17
    log_num_frac(e_Tl) = 0.90
    log_num_frac(e_Pb) = 1.75
    log_num_frac(e_Bi) = 0.65
    log_num_frac(e_Th) = 0.02
    log_num_frac(e_U) = -0.54

    call set_X_div_H(alpha_div_Fe,Fe_div_H, X_div_H)
    
    log_num_frac = log_num_frac + X_div_H
    
    !set number fractions and normalize
    num_frac = exp(ln10*log_num_Frac)
    num_frac = num_frac/sum(num_frac)

    !set mass fractions and normalize
    mass_frac = num_frac*element_atomic_weight
    mass_frac = mass_frac/sum(mass_frac)

    X = mass_frac(1)
    Y = mass_frac(2)
    Z = 1.0d0-(X+Y)
    
    if( dY_div_dZ /=0.0d0 )then

       Xold= X
       Yold= Y
       Zold= Z

       !reset X, Y, and Z starting from the initial values
       !and now accounting for the desired Y value
       Z_div_X = Z/X

       !new X
       X = (1.0d0 - Y_BBN)/(1.0d0+ Z_div_X*(1.0d0 + dY_div_dZ))

       !new Z
       Z = Z_div_X * X

       !new Y
       Y = 1.0d0 - (X + Z)

       !now reset mass_frac
       mass_frac(1) = X
       mass_frac(2) = Y
       mass_frac(3:num_chem_elements) = (Z/Zold)*mass_frac(3:num_chem_elements)

       num_frac = mass_frac/element_atomic_weight
       num_frac = num_frac/sum(num_frac)
       
    endif
  end subroutine set_abunds

  subroutine set_X_div_H(alpha_div_Fe,Fe_div_H,X_div_H)
    real(dp), intent(in) :: alpha_div_Fe, Fe_div_H
    real(dp), intent(out) :: X_div_H(num_chem_elements)
    real(dp) :: X_div_Fe(num_chem_elements)
    X_div_H = 0d0

    X_div_H(e_O) = alpha_div_Fe
    X_div_H(e_Ne) = alpha_div_Fe
    X_div_H(e_Mg) = alpha_div_Fe
    X_div_H(e_Si) = alpha_div_Fe
    X_div_H(e_S) = alpha_div_Fe
    X_div_H(e_Ar) = alpha_div_Fe
    X_div_H(e_Ca) = alpha_div_Fe
    X_div_H(e_Ti) = alpha_div_Fe

    X_div_H(e_Li:num_chem_elements) = X_div_H(e_Li:num_chem_elements) + Fe_div_H
    
  end subroutine set_X_div_H
  
  subroutine init_atlas_abunds
    real(dp), parameter :: altas_offset = 12.04d0
    real(dp) :: Num_Frac_Z
    character(len=128) :: stuff
    integer :: i, ilo, ihi
    open(1,file=trim(atm_file))
    do i=1,4
       read(1,*) !skip header lines
    enddo

    num_frac = 0d0
    log_num_frac = -9.9d1
    read(1,'(a)') stuff

    read(stuff(45:52),*) num_frac(1)
    read(stuff(55:62),*) num_frac(2)

    if(verbose) write(0,*) ' H, He = ', num_frac(1:2)
       
    do i=1,16
       ilo = 6*(i-1)+3 ! start at Li
       ihi = ilo+5
       read(1,'(17x,6(3x,f7.2))') log_num_frac(ilo:ihi)
    enddo
    read(1,'(17x,3x,f7.2)') log_num_frac(ihi+1:ihi+1)
    close(1)

    num_Frac_Z = 1d0 - sum(num_frac(1:2))
    num_frac(3:num_chem_elements) = 1d1**log_num_frac(3:num_chem_elements)
    num_frac(3:num_chem_elements) = num_frac(3:num_chem_elements)*num_frac_Z/sum(num_frac(3:num_chem_elements))

    !set mass fractions
    mass_frac = num_frac*element_atomic_weight
    mass_frac = mass_frac/sum(mass_frac)

    if(verbose)then
       do i=1,num_chem_elements
          write(0,*) chem_element_name(i), log_num_frac(i), num_frac(i), mass_frac(i)
       enddo
       write(0,*) '   sum(X) = ', sum(num_frac)
    endif
  end subroutine init_atlas_abunds

  subroutine net_elements_and_isotopes
    integer :: pass
    character(len=iso_name_length) :: name
    !print out some details of the current net and count the elements
    !first pass counts, second pass assigns to net_elements array
    do pass = 1,3
       current = -1
       num_net_elements = 0 
       do i=1,g% num_isos
          id = g% chem_id(i)
          el = chem_isos% Z(id)
          name = chem_isos% name(id)
          if(pass==1.and.id > 0.and.verbose) write(0,'(2i4,a8,f10.6,2i8)') i, &
               chem_isos% chem_id(id), trim(name), &
               chem_isos% W(id), el, chem_isos% N(id)
          !count elements present:
          if(trim(name) /= 'neut' .and. trim(name) /= 'prot') then
             if(el == current)then
                if(pass==1) cycle
                if(pass>1) then
                   e(num_net_elements)% num_isos = e(num_net_elements)% num_isos +1
                   if(pass==3) e(num_net_elements)% isos(e(num_net_elements)% num_isos) = chem_isos% chem_id(id)
                endif
             else
                num_net_elements = num_net_elements + 1
                if(pass>1) then
                   e(num_net_elements)% id = el
                   e(num_net_elements)% num_isos = 1
                   if(pass==3) e(num_net_elements)% isos(1) = chem_isos% chem_id(id)
                endif
                current = el
             endif
          endif
       enddo
       if (pass==1) then
          if(verbose) write(0,*) ' num_elements = ', num_net_elements
          allocate(e(num_net_elements))
          e(:)% num_isos = 0
          e(:)% id = 0
       else if (pass==2) then
          do i=1,num_net_elements
             allocate(e(i)% isos(e(i)% num_isos))
             e(i)% isos = 0
          enddo
       endif
    enddo

    !now assign names to the isotopes for each element
    do i=1,num_net_elements
       allocate(e(i)% iso_names(e(i)% num_isos))
       allocate(e(i)% iso_fracs(e(i)% num_isos))
       do j=1,e(i)% num_isos
          e(i)% iso_names(j) = chem_isos% name(e(i)% isos(j))
          e(i)% iso_fracs(j) = lodders03_element_atom_percent(e(i)% iso_names(j))
       enddo
       call fix_isotope_element_fractions(e(i)% iso_fracs)
    enddo
  end subroutine net_elements_and_isotopes

  subroutine fix_isotope_element_fractions(fracs)
    real(dp), intent(inout) :: fracs(:)
    real(dp) :: sum_fracs, diff
    real(dp), parameter :: eps = 1d-10
    integer :: i, n, count
    n=size(fracs)
    fracs=fracs*1d-2 !originally sum to 100%
    sum_fracs=sum(fracs)
    if(sum_fracs==0d0 .or. n==1 )then
       fracs=1d0
    else if(sum_fracs > 0d0 .and. sum_fracs < 1d0)then
       count=0
       diff=1d0-sum_fracs
       do i=1,n
          if(fracs(i)==0d0) count=count+1
       enddo
       where(fracs==0d0) fracs = diff/real(count,kind=dp)
    endif
    sum_fracs = sum(fracs)
    fracs = fracs/sum_fracs
    if(abs(sum(fracs) - 1d0) > eps) then
       write(0,*) '  fix_isotope_element_fractions  '
       write(0,*) fracs
       write(0,*) sum_fracs
       stop 77
    endif
  end subroutine fix_isotope_element_fractions


  subroutine initialize_mesa(mesa_dir, ierr)
    character (len=*), intent(in) :: mesa_dir
    integer, intent(out) :: ierr
    ierr = 0
    call const_init(mesa_dir,ierr)     
    if (ierr /= 0) then
       write(0,*) 'const_init failed'
       stop 1
    end if
    call chem_init('isotopes.data', ierr)
    if (ierr /= 0) then
       write(0,*) 'chem_init failed'
       return
    end if

!      subroutine rates_init( &
!           reactionlist_filename, jina_reaclib_filename, &
!           rates_table_dir_in, &
!           use_special_weak_rates, &
!           special_weak_states_file, &
!           special_weak_transitions_file, &
!           cache_dir, ierr)
!         use rates_def
!         use reaclib_input, only: do_read_reaclib
!         use load_weak, only: load_weak_data
!         use load_ecapture, only: load_ecapture_data
!         use rates_initialize, only: init_rates_info
!         
!         character (len=*), intent(in) :: reactionlist_filename, jina_reaclib_filename, rates_table_dir_in
!         logical, intent(in) :: use_special_weak_rates
!         character (len=*), intent(in) :: special_weak_states_file, special_weak_transitions_file
!         character (len=*), intent(in) :: cache_dir ! '' means use default
!         integer, intent(out) :: ierr ! 0 means AOK. 

    call rates_init('reactions.list', '', 'rate_tables', .false., '', '', '', ierr)
    if (ierr /= 0) then
       write(0,*) 'rates_init failed'
       return
    end if
    call net_init(ierr)
    if (ierr /= 0) then
       write(0,*) 'net_init failed'
       return
    end if
  end subroutine initialize_mesa


  subroutine setup_net( net_file, handle, species, chem_id, net_iso, ierr)

    character (len=*), intent(in) :: net_file
    integer, pointer :: chem_id(:), net_iso(:) ! set, but not allocated
    integer, intent(out) :: handle, species, ierr

    ierr = 0
    handle = alloc_net_handle(ierr)
    if (ierr /= 0) then
       write(0,*) 'alloc_net_handle failed'
       return
    end if

    call net_start_def(handle, ierr)
    if (ierr /= 0) then
       write(0,*) 'net_start_def failed'
       return
    end if

    if(verbose) write(0,*) 'load ' // trim(net_file)
    call read_net_file(net_file, handle, ierr)
    if (ierr /= 0) then
       write(0,*) 'read_net_file failed ', trim(net_file)
       return
    end if

    call net_finish_def(handle, ierr)
    if (ierr /= 0) then
       write(0,*) 'net_finish_def failed'
       return
    end if

    species = net_num_isos(handle, ierr)
    if (ierr /= 0) then
       write(0,*) 'failed in net_num_isos'
       return
    end if

    call get_chem_id_table_ptr(handle, chem_id, ierr)
    if (ierr /= 0) then
       write(0,*) 'failed in get_chem_id_table_ptr'
       return
    end if

    call get_net_iso_table_ptr(handle, net_iso, ierr)
    if (ierr /= 0) then
       write(0,*) 'failed in get_net_iso_table_ptr'
       return
    end if

  end subroutine setup_net

  subroutine output_XYZ
    character(len=9) :: fmt = '(a9,f8.6)'
    if(verbose)then
       write(0,fmt) ' X Y Z '
       write(0,fmt) ' X = ', mass_frac(1)
       write(0,fmt) ' Y = ', mass_frac(2)
       write(0,fmt) ' Z = ', sum(mass_frac(3:))
       write(0,fmt) ' X+Y+Z = ', sum(mass_frac)
    endif
    open(98,file='input_XYZ')
    write(98,'(1p,e20.10)') mass_frac(1)
    write(98,'(1p,e20.10)') mass_frac(2)
    write(98,'(1p,e20.10)') sum(mass_frac(3:num_chem_elements))
    close(98)
  end subroutine output_XYZ
  
  subroutine output_for_OPAL
    !for OPAL
    integer, parameter :: num_OPAL_elements = 21
    integer :: OPAL_elements(num_OPAL_elements)
    integer :: op
    real(dp) :: Z_frac
    !---------------H He C N O Ne Na Mg Al Si  P  S Cl Ar  K Ca Ti Cr Mn Fe Ni
    OPAL_elements = [1,2,6,7,8,10,11,12,13,14,15,16,17,18,19,20,22,24,25,26,28]
    Z_frac = sum(num_frac(3:num_chem_elements))
    if(verbose)then
       write(0,*)
       write(0,*) '  OPAL metals  '
       do i=3,num_OPAL_elements
          op = OPAL_elements(i)
          write(0,'(i3,1x,a4,1x,f8.6)') i, chem_element_name(op), num_frac(op)/Z_frac
       enddo
       write(0,*)
       write(0,*) ' sum = ', sum(num_frac(OPAL_elements))
       write(0,*)
    endif
  end subroutine output_for_OPAL

  subroutine output_for_MESA
    integer :: i, j, k, num_elements, pass
    real(dp) :: sum, extra
    !for MESA/star file_for_uniform_xa
    num_elements = size(e)
    sum = 0d0
    open(99,file='input_initial_xa.data')
    do pass=1,2
       do i=1,num_elements
          do j=1,num_chem_elements
             if(e(i)% id == j) then
                do k=1,e(i)% num_isos
                   if(pass==1) then
                      sum = sum + e(i)% iso_fracs(k)*mass_frac(j)
                   elseif(pass==2)then
                      if(i==num_elements) then
                         write(99,'(a8,1p,e20.10)') e(i)% iso_names(k), e(i)% iso_fracs(k) * (mass_frac(j) + extra)
                      else
                         write(99,'(a8,1p,e20.10)') e(i)% iso_names(k), e(i)% iso_fracs(k) * mass_frac(j)
                      endif
                   endif
                enddo
             endif
          enddo
       enddo
       !here, if sum/=1, add the leftovers to the heaviest element
       if(pass==1) extra= 1d0 - sum
    enddo
    close(99)

  end subroutine output_for_MESA

end program initial_xa_calculator
