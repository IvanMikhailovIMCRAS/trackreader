!> @brief	Потенциалы взаимодействий
!!              модуль для прикладных программистов
!! @date		03.06.2018
!!
!! @author	Михайлов И.В., Шавыкин О.В.
Module Potential
    Use CommonModule
    Use TestErrors
    Use Surface
    Use BoxMuller
    
    !todo:
    integer(4), parameter :: num_layers = 200
    real(8) num(num_layers)
    real(8) num_norm(num_layers)
    !!!!!!!!!!!!
    
    !> индетификаторы поверхностей:
    integer(4), parameter, public :: FORCE_CODE_HARM = 1, &
                                     FORCE_CODE_HARM_ANGL = 2, &
                                     FORCE_CODE_LJ = 3, &
                                     FORCE_CODE_FIELD = 4, &
                                     FORCE_CODE_PARABOLIC = 5, &
                                     FORCE_CODE_PARABOLIC_RANDOM = 6, &
                                     FORCE_CODE_PARABOLIC_PHI = 7, &
                                     FORCE_CODE_PARABOLIC_SELF = 8, &
                                     FORCE_CODE_PARABOLIC_FORCE = 9, &
                                     FORCE_CODE_PARABOLIC_FORCE2 = 10
                                     
    
    type force_info 
        !>>> эти параметры задает пользователь
        integer(4) type_code
        logical(1) bunbond
        integer(4) num_param
        integer(4) num_rows
        real(8), allocatable, dimension(:) :: param
        !>>> эти параметры программа определяет сама
        ! сколько типов объектов (частиц, поверхностей) охватывает такой потенциал:
        integer(4) num_type
        integer(4), allocatable, dimension(:) :: types
        logical(1) bsurface 
    end type

Contains 

!>
subroutine PotentialRecognition(n_file, num_type, type_name, num_surface, surf, force_name, force_str)
    ! input 
    integer(4) n_file
    integer(4) num_type
    character(LEN_TYPE_NAME) type_name(num_type)
    integer(4) num_surface
    type(surface_info) surf(num_surface)
    ! output
    character(255) force_name
    type(force_info) force_str
    ! local
    character(255) temp
    logical(1) bfind, bfind2, bfind3
    integer(4) ioer
    character(255), allocatable, dimension(:) :: temp_type
    integer(4) i, t
    
    read(n_file,*) force_name
    
    bfind = .false.

    ! гармонический без учета обрезки:
    if (force_name(1:len_trim(force_name)).eq.'harm') then
        force_str%type_code = FORCE_CODE_HARM
        force_str%bunbond = .false.
        force_str%num_type = 2
        force_str%num_param = 2
        force_str%num_rows = 1
        
        allocate(force_str%param(1:force_str%num_param*force_str%num_rows))
        backspace(n_file)
        read(n_file, *, iostat=ioer) temp, temp, temp, force_str%param
        if (ioer.ne.0) then
            call ERROR(ERR_CODE_READ_POTENTIAL)
        endif
        
        ! нашли:
        bfind = .true.
    endif
    
    ! леннард-джонса с учетом обрезки:
    if (force_name(1:len_trim(force_name)).eq.'lj') then
        force_str%type_code = FORCE_CODE_LJ
        force_str%bunbond = .true.
        force_str%num_type = 2
        ! rcut=2.8 sigma=0.8 eps=1.0 solvent=0.0
        force_str%num_param = 4
        force_str%num_rows = 1
        
        allocate(force_str%param(1:force_str%num_param*force_str%num_rows))
        backspace(n_file)
        read(n_file, *, iostat=ioer) temp, temp, temp, force_str%param
        if (ioer.ne.0) then
            call ERROR(ERR_CODE_READ_POTENTIAL)
        endif
        
        ! нашли:
        bfind = .true.
    endif
    
    ! поле:
    if (force_name(1:len_trim(force_name)).eq.'field') then
        force_str%type_code = FORCE_CODE_FIELD
        force_str%bunbond = .false.
        force_str%num_type = 1
        force_str%num_param = 6
        
        !>>узнаем имя входного файла
        backspace(n_file)
        read(n_file,*) force_name, temp, temp
        
        call CalcNumRows(force_str%num_param, temp, force_str%num_rows)
        allocate(force_str%param(1:force_str%num_param*force_str%num_rows))
        force_str%param(:) = 0.0
        !теперь читаем:
        open(n_field,file=temp(1:len_trim(temp)))
	do i = 1, force_str%num_rows
            read(n_field, *) force_str%param(force_str%num_param*(i-1)+1:force_str%num_param*i)
	enddo
	close(n_field)
        
        ! нашли:
        bfind = .true.
    endif
    
    ! поле:
    if (force_name(1:len_trim(force_name)).eq.'parabolic') then
        force_str%type_code = FORCE_CODE_PARABOLIC
        force_str%bunbond = .false.
        force_str%num_type = 1
        ! k, H1
        force_str%num_param = 2
        force_str%num_rows = 1
        
        allocate(force_str%param(1:force_str%num_param*force_str%num_rows))
        backspace(n_file)
        read(n_file, *, iostat=ioer) temp, temp, force_str%param
        if (ioer.ne.0) then
            call ERROR(ERR_CODE_READ_POTENTIAL)
        endif
        
        ! нашли:
        bfind = .true.
    endif
    
    ! поле:
    if (force_name(1:len_trim(force_name)).eq.'parabolic_random') then
        force_str%type_code = FORCE_CODE_PARABOLIC_RANDOM
        force_str%bunbond = .false.
        force_str%num_type = 1
        ! k, H1, dp, rcut, np_minus
        force_str%num_param = 5
        force_str%num_rows = 1
        
        allocate(force_str%param(1:force_str%num_param*force_str%num_rows))
        backspace(n_file)
        read(n_file, *, iostat=ioer) temp, temp, force_str%param(1:force_str%num_param-1)
        ! находим число точек разбиения полуинтервала (в отрицательную сторону):
        force_str%param(5) = force_str%param(4)/force_str%param(3)
        if (ioer.ne.0) then
            call ERROR(ERR_CODE_READ_POTENTIAL)
        endif
        
        ! нашли:
        bfind = .true.
    endif

    ! поле: 
    !  sense of density
    if (force_name(1:len_trim(force_name)).eq.'parabolic_phi') then
        force_str%type_code = FORCE_CODE_PARABOLIC_PHI
        force_str%bunbond = .false.
        force_str%num_type = 1
        ! k, H1, x_shift,  abs_Fx
        force_str%num_param = 4
        force_str%num_rows = 1
        
        allocate(force_str%param(1:force_str%num_param*force_str%num_rows))
        backspace(n_file)
        read(n_file, *, iostat=ioer) temp, temp, force_str%param(1:force_str%num_param)
        if (ioer.ne.0) then
            call ERROR(ERR_CODE_READ_POTENTIAL)
        endif
        
        ! нашли:
        bfind = .true.
    endif

    ! поле: 
    !  sense of density
    if (force_name(1:len_trim(force_name)).eq.'parabolic_self') then
        force_str%type_code = FORCE_CODE_PARABOLIC_SELF
        force_str%bunbond = .false.
        force_str%num_type = 1
        ! k, H1, x_shift,  abs_Fx
        force_str%num_param = 4
        force_str%num_rows = 1
        
        allocate(force_str%param(1:force_str%num_param*force_str%num_rows))
        backspace(n_file)
        read(n_file, *, iostat=ioer) temp, temp, force_str%param(1:force_str%num_param)
        if (ioer.ne.0) then
            call ERROR(ERR_CODE_READ_POTENTIAL)
        endif
        
        ! нашли:
        bfind = .true.
    endif

    !поле:
    if (force_name(1:len_trim(force_name)).eq.'parabolic_force') then
        force_str%type_code = FORCE_CODE_PARABOLIC_FORCE
        force_str%bunbond = .false.
        force_str%num_type = 1
        ! k, H1
        force_str%num_param = 2
        force_str%num_rows = 1
        
        allocate(force_str%param(1:force_str%num_param*force_str%num_rows))
        backspace(n_file)
        read(n_file, *, iostat=ioer) temp, temp, force_str%param(1:force_str%num_param)
        if (ioer.ne.0) then
            call ERROR(ERR_CODE_READ_POTENTIAL)
        endif
        
        ! нашли:
        bfind = .true.
    endif
    
    !поле:
    if (force_name(1:len_trim(force_name)).eq.'parabolic_force2') then
        force_str%type_code = FORCE_CODE_PARABOLIC_FORCE2
        force_str%bunbond = .false.
        force_str%num_type = 1
        ! k, H1
        force_str%num_param = 2
        force_str%num_rows = 1
        
        allocate(force_str%param(1:force_str%num_param*force_str%num_rows))
        backspace(n_file)
        read(n_file, *, iostat=ioer) temp, temp, force_str%param(1:force_str%num_param)
        if (ioer.ne.0) then
            call ERROR(ERR_CODE_READ_POTENTIAL)
        endif
        
        ! нашли:
        bfind = .true.
    endif
    
    !>>>автоматическое определение типов частиц и поверхностей
    if (bfind.eqv..true.) then
        allocate(force_str%types(1:force_str%num_type))
        allocate(temp_type(1:force_str%num_type))
        !считаем типы для которых работает этот потенциал:
        backspace(n_file)
        read(n_file, *, iostat=ioer) temp, temp_type(1:force_str%num_type)
        
        do t = 1, force_str%num_type
            bfind2 = .false.
            do i = 1, num_type
                temp=type_name(i)
                if (temp_type(t).eq.temp(1:len_trim(temp))) then
                    bfind2 = .true.
                    force_str%types(t) = i
                endif
            enddo
            
            if (bfind2.eqv..false.) then
                !поверхность может идти первой:
                if (force_str%num_type.eq.2.and.t.eq.1) then
                    bfind3 = .false.
                    !ищем среди поверхностей:
                    do i = 1, num_surface
                        temp = surf(i)%name
                        if (temp_type(t).eq.temp(1:len_trim(temp))) then
                            bfind3 = .true.
                            force_str%types(t) = i
                        endif
                    enddo
                    
                    if (bfind3.eqv..true.) then
                        ! потенциал с поверхностью
                        force_str%bsurface = .true.
                    else
                        call ERROR(ERR_CODE_POTENTIAL_UNKNOWN_PART_TYPE2)
                    endif
                else
                    call ERROR(ERR_CODE_POTENTIAL_UNKNOWN_PART_TYPE)
                endif
            else
                ! потенциал без поверхности
                force_str%bsurface = .false.
            endif
        enddo
    else
        call ERROR(ERR_CODE_POTENTIAL_UNKNOWN)
    endif
    
    !
    !force_str%param(:) = 0.0

    return
end subroutine PotentialRecognition

!>> подсчитаем число строк во внешнем файле
subroutine CalcNumRows(num_param, file_name, num_rows)
    ! input
    integer(4) num_param
    character(255) file_name
    ! output
    integer(4) num_rows
    ! local
    real(8) param(num_param)
    logical(1) file_exists
    integer(4) ioer
    
    num_rows = 0
    
    INQUIRE(FILE=file_name(1:len_trim(file_name)), EXIST=file_exists) 
    if (file_exists.eqv..true.) then
	open(n_field,file=file_name(1:len_trim(file_name)))
	ioer = 0
	do while (ioer.eq.0)
            read(n_file, *, iostat=ioer) param(:)
            if (ioer.eq.0) num_rows = num_rows + 1
	enddo
	
	close(n_field)
    else
        call ERROR(ERR_CODE_POTENTIAL_OPEN_FIELD)
    endif

    return
end subroutine CalcNumRows

!>
subroutine PotentialRcut(force_str, r_cut)
    ! input 
    type(force_info) force_str
    ! output
    real(8) r_cut
    
    if (force_str%bunbond.eqv..true.) then
        r_cut = force_str%param(1)
    else
        call ERROR(ERR_CODE_POTENTIAL_NO_UNBOND)
    endif

    return
end subroutine PotentialRcut


!>
subroutine PotentialForceOneParticle(x, y, z, force_str, Fx, Fy, Fz)
    ! input 
    real(8) x, y, z
    type(force_info) force_str
    ! output
    real(8) Fx, Fy, Fz
    ! 
    select case(force_str%type_code)
        case(FORCE_CODE_FIELD)
            Fx = 0.0 
            Fy = 0.0
            Fz = 0.0
        case(FORCE_CODE_PARABOLIC)
            call Parabolic_dx_Force(x, y, z, force_str%param(1), force_str%param(2), &
                                & Fx, Fy, Fz)
        case(FORCE_CODE_PARABOLIC_RANDOM)
            ! (1) k (2) H1 (3) dp (5) np_minus
            call Parabolic_dx_Force_Random(x, y, z, force_str%param(1), force_str%param(2), &
                                & force_str%param(3), int(force_str%param(5)), Fx, Fy, Fz)
        case(FORCE_CODE_PARABOLIC_PHI)
            call Parabolic_dx_Force_Phi(x, y, z, force_str%param(1), force_str%param(2), &
                                & force_str%param(3), force_str%param(4), Fx, Fy, Fz)
        case(FORCE_CODE_PARABOLIC_SELF)
            call Parabolic_dx_Force_Self(x, y, z, force_str%param(1), force_str%param(2), &
                                & force_str%param(3), force_str%param(4), Fx, Fy, Fz)
        case(FORCE_CODE_PARABOLIC_FORCE)
            call Parabolic_dx_Force_dUdr(x, y, z, force_str%param(1), force_str%param(2), &
                                & Fx, Fy, Fz)
        case(FORCE_CODE_PARABOLIC_FORCE2)
            call Parabolic_dx_Force_dUdr2(x, y, z, force_str%param(1), force_str%param(2), &
                                & Fx, Fy, Fz)
        case default
            call ERROR(ERR_CODE_UNKNOWN_1PART_POTENTIAL)
    end select 

    return
end subroutine PotentialForceOneParticle

!>
subroutine PotentialForceTwoParticle(x1, y1, z1, x2, y2, z2, force_str, Fx, Fy, Fz)
    ! input 
    real(8) x1, y1, z1, x2, y2, z2
    type(force_info) force_str
    ! output
    real(8) Fx, Fy, Fz
    ! local
    real(8) dx, dy, dz

    select case(force_str%type_code)
        case(FORCE_CODE_HARM)
            call PartToPartDistance(x1, y1, z1, x2, y2, z2, dx, dy, dz)
            ! kbond, lbond
            call Harmonic_Force(dx, dy, dz, force_str%param(1), force_str%param(2), &
                                & Fx, Fy, Fz)
        case(FORCE_CODE_LJ)
            call PartToPartDistance(x1, y1, z1, x2, y2, z2, dx, dy, dz)
            !sig, rcut, eps, solvent
            call Lennard_Jones_Force(dx, dy, dz,  force_str%param(2), force_str%param(1), &
                                    & force_str%param(3), force_str%param(4), &
                                    & Fx, Fy, Fz)
        case default
            call ERROR(ERR_CODE_UNKNOWN_2PART_POTENTIAL)
    end select 

    return
end subroutine PotentialForceTwoParticle

!>
subroutine PotentialForceThreeParticle(x1, y1, z1, x2, y2, z2, x3, y3, z3, force_str, Fx, Fy, Fz)
    ! input 
    real(8) x1, y1, z1, x2, y2, z2, x3, y3, z3
    type(force_info) force_str
    ! output
    real(8) Fx, Fy, Fz
    ! local
    real(8) dx, dy, dz

    select case(force_str%type_code)
        case(FORCE_CODE_HARM_ANGL)
            call PartToPartDistance(x1, y1, z1, x3, y3, z3, dx, dy, dz)
            call Angle_Force(force_str%param(1), force_str%param(2), &
                            & dx, dy, dz, &
                            & x1, y1, z1, & 
                            & x2, y2, z2, &
                            & x3, y3, z3, &
                            & Fx, Fy, Fz)
        case default
            call ERROR(ERR_CODE_UNKNOWN_3PART_POTENTIAL)
    end select 

    return
end subroutine PotentialForceThreeParticle

!>
subroutine PotentialForceSurfaceParticle(surf, x, y, z, force_str, Fx, Fy, Fz)
    ! input 
    type(surface_info) surf
    real(8) x, y, z
    type(force_info) force_str
    ! output
    real(8) Fx, Fy, Fz
    ! local
    real(8) dx, dy, dz

    select case(force_str%type_code)
        case(FORCE_CODE_HARM)
           ! call PartToSurfaceDistance(x, y, z, surf, dx, dy, dz)
            call DistanceToSurface(x, y, z, surf, dx, dy, dz)
            ! kbond, lbond
            call Harmonic_Force(dx, dy, dz, force_str%param(1), force_str%param(2), &
                                & Fx, Fy, Fz)
        case(FORCE_CODE_LJ)
           ! call PartToSurfaceDistance(x, y, z, surf, dx, dy, dz)
            call DistanceToSurface(x, y, z, surf, dx, dy, dz)
            !sig, rcut, eps, solvent
            call Lennard_Jones_Force(dx, dy, dz,  force_str%param(2), force_str%param(1), &
                                    & force_str%param(3), force_str%param(4), &
                                    & Fx, Fy, Fz)
        case default
            call ERROR(ERR_CODE_UNKNOWN_SURFACE_PART_POTENTIAL)
    end select 

    return
end subroutine PotentialForceSurfaceParticle

!>>>>>>>>>>>>>>>>>>>>>
subroutine Parabolic_dx_Force(x, y, z, k, H1, &
                                & Fx, Fy, Fz)
    ! input
    real(8) x, y, z      
    real(8) k, H1
    ! output
    real(8) Fx, Fy, Fz
    ! local
    
    if (x.gt.0) then
        if (abs(x).le.H1) then
            Fx = (3.0/2.0)*k**2*(H1**2 - x**2)
            !Fx = (9.0/2.0)*(k**4)*x*(H1**2 - x**2)
        else
            Fx = 0.0
        endif
    else
        Fx = 1000   !*k**2*H1**2
    endif
    Fy = 0.0
    Fz = 0.0
          
    return
end subroutine Parabolic_dx_Force

!>>>>>>>>>>>>>>>>>>>>>
!>>>
subroutine Boltzmann(k, H1, x, u, G1)
    ! input
    real(8) k, H1, x
    ! output 
    real(8) u, G1
    !
    u = (3.0/2.0)*k**2*(H1**2 - x**2) 
    G1 = exp( - u)
    
    return 
end subroutine Boltzmann

!>>>>>>>>>>>>>>>>>>>>>
subroutine Parabolic_dx_Force_Random(x, y, z, k, H1, dp, np_minus, &
                                & Fx, Fy, Fz)
    ! input
    real(8) x, y, z      
    real(8) k, H1
    real(8) dp
    integer(4) np_minus
    ! output
    real(8) Fx, Fy, Fz
    ! local
    real(8) xtemp
    real(8) G1(2*np_minus), u(2*np_minus)
    real(8) sumg, prob_rand
    integer(4) i
    
    if (x.gt.0) then
        !if (abs(x).le.H1) then
            ! находим распределение Больцмана на выбранном участке:
            do i = 1, 2*np_minus
                xtemp = x - (np_minus - i)*dp
                call Boltzmann(k, H1, xtemp, u(i), G1(i))
            enddo
            sumg = sum(G1)
            G1(:) = G1(:) / sumg
            
            prob_rand = randm()
            
            sumg = 0
            !ищем интервал:
            do i = 1, 2*np_minus
                sumg = sumg + G1(i)
                
                if (prob_rand.le.sumg) then
                   Fx = u(i)!/20.0
                   exit
                endif
            enddo
        !else
            !Fx = 0.0
        !endif
    else
        Fx = 1000   !*k**2*H1**2
    endif
    Fy = 0.0
    Fz = 0.0
          
    return
end subroutine Parabolic_dx_Force_Random

!>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>
subroutine Parabolic_dx_Force_Phi(x, y, z, k, H1, x_shift, abs_Fx, &
                                & Fx, Fy, Fz)
    ! input
    real(8) x, y, z      
    real(8) k, H1
    real(8) x_shift
    real(8) abs_Fx
    ! output
    real(8) Fx, Fy, Fz
    ! local
    real(8) phi
    real(8) rand_phi
    real(8) u_m, u_p
    real(8) G1_p, G1_m
    
    if (x.gt.0) then
        !определяем плотность "мнимых" (размазанных) частиц
        if (abs(x).le.H1) then
            phi = (3.0/2.0)*k**2*(H1**2 - x**2)
        else
            phi = 0.0
        endif
        ! ищем вероятность столкновения с облаком размазанных частиц:
        rand_phi = randm()
        
        ! если столкновение произошло, то необходимо найти понять в какую сторону будем двигаться:
        if (rand_phi.ge.phi) then
           call Boltzmann(k, H1, x + x_shift, u_p, G1_p) 
           if (x.ge.x_shift) then
                call Boltzmann(k, H1, x - x_shift, u_m, G1_m)
           else
                call Boltzmann(k, H1, x, u_m, G1_m)
           endif 
           G1_m = G1_m / (G1_m + G1_p)
           u_m = u_m / (u_m + u_p)
           
           ! на сколько в долевом отношении мы ушли от критической плотности:
           rand_phi = (rand_phi - phi) / (1.0 - phi)
           
!           write(*,*) rand_phi, u_m, G1_m
           ! узнаем, в какую сторону пойдет частица из-за столкновения:
           if (rand_phi.le.u_m) then
                Fx = abs_Fx
            else
                Fx = -abs_Fx  ! -0.243
           endif
           
        endif
    else
        Fx = 1000   !*k**2*H1**2
    endif
    Fy = 0.0
    Fz = 0.0
          
    return
end subroutine Parabolic_dx_Force_Phi


!>>>>>>>>>>>>>>>>>>>>>
subroutine Parabolic_dx_Force_Self(x, y, z, k, H1, x_shift, abs_Fx, &
                                & Fx, Fy, Fz)
    ! input
    real(8) x, y, z      
    real(8) k, H1
    real(8) x_shift
    real(8) abs_Fx
    ! output
    real(8) Fx, Fy, Fz
    ! local
    real(8) phi
    real(8) rand_phi
    real(8) sum_rad
    real(8) u_m, u_p
    real(8) G1_p, G1_m
    real(8) phi_p, phi_m
    !integer(4), parameter :: num_layer = 200
    !real(8) rad(num_layer), phi_now(num_layer)
    real(8), parameter :: dr = 0.1
    integer(4) i
    
    if (x.gt.0) then
        !определяем плотность "мнимых" (размазанных) частиц
        if (abs(x).le.H1) then
            phi = (3.0/2.0)*k**2*(H1**2 - x**2)
        else
            phi = 0.0
        endif

        
        sum_rad = 0.0
        
        i = int(x/dr) + 1
        ! если моментальная больше мнимой, то заменяем:
        if (num_norm(i).gt.phi) then
            phi = num_norm(i)
        endif
        
        ! ищем вероятность столкновения с облаком размазанных частиц:
        rand_phi = randm()
        
        
        ! если столкновение произошло, то необходимо найти понять в какую сторону будем двигаться:
        if (rand_phi.ge.phi) then
           call Boltzmann(k, H1, x + x_shift, u_p, G1_p) 
           i = int((x + x_shift)/dr) + 1
           phi_p = num_norm(i)
           if (x.ge.x_shift) then
                call Boltzmann(k, H1, x - x_shift, u_m, G1_m)
                i = int((x - x_shift)/dr) + 1
                phi_p = num_norm(i)
           else
                call Boltzmann(k, H1, x, u_m, G1_m)
                i = int(x/dr) + 1
                phi_p = num_norm(i)
           endif 
           !G1_m = G1_m / (G1_m + G1_p)
           u_m = u_m / (u_m + u_p)
           phi_m = phi_m / (phi_m + phi_p)
           
           ! на сколько в долевом отношении мы ушли от критической плотности:
           rand_phi = (rand_phi - phi) / (1.0 - phi)
           
!           write(*,*) rand_phi, u_m, G1_m
           ! узнаем, в какую сторону пойдет частица из-за столкновения:
           if (rand_phi.le.u_m) then
                Fx = abs_Fx
            else
                Fx = -abs_Fx
           endif
           
        endif
    else
        Fx = 1000   !*k**2*H1**2
    endif
    Fy = 0.0
    Fz = 0.0
          
    return
end subroutine Parabolic_dx_Force_Self

!>>>>>>>>>>>>>>>>>>>>>
subroutine Parabolic_dx_Force_dUdr(x, y, z, k, H1, &
                                & Fx, Fy, Fz)
    ! input
    real(8) x, y, z      
    real(8) k, H1
    ! output
    real(8) Fx, Fy, Fz
    ! local
    real(8) phi
    real(8) rand_phi
    !integer(4), parameter :: num_layer = 200
    !real(8) rad(num_layer), phi_now(num_layer)
    real(8), parameter :: dr = 0.1
    integer(4) i
    
    if (x.gt.0) then
        !определяем плотность "мнимых" (размазанных) частиц
        if (abs(x).le.H1) then
            phi = (3.0/2.0)*k**2*(H1**2 - x**2)
        else
            phi = 0.0
        endif
    
        i = int(x/dr) + 1
        ! если моментальная больше мнимой, то заменяем:
        if (num_norm(i).gt.phi) then
            phi = num_norm(i)
        endif
        
        ! ищем вероятность столкновения с облаком размазанных частиц:
        rand_phi = randm()
    
        ! если столкновение произошло, то необходимо найти понять в какую сторону будем двигаться:
        if (rand_phi.ge.phi) then
           ! производная от параболического потенциала:
           Fx = (3.0*k**2)*x
        endif
    else
        Fx = 1000   !*k**2*H1**2
    endif
    Fy = 0.0
    Fz = 0.0
          
    return
end subroutine Parabolic_dx_Force_dUdr

!>>>>>>>>>>>>>>>>>>>>>
subroutine Parabolic_dx_Force_dUdr2(x, y, z, k, H1, &
                                & Fx, Fy, Fz)
    ! input
    real(8) x, y, z      
    real(8) k, H1
    ! output
    real(8) Fx, Fy, Fz
    ! local
    !integer(4), parameter :: num_layer = 200
    !real(8) rad(num_layer), phi_now(num_layer)
    real(8), parameter :: dr = 0.1
    integer(4) i
    
    if (x.gt.0) then
        if (x.le.H1) then
            ! производная от параболического потенциала:
            Fx = (3.0*k**2)*x
        else
            ! выше - нет никакой растягивающей силы!
            Fx = 0.0
        endif
    else
        Fx = 1000   !*k**2*H1**2
    endif
    Fy = 0.0
    Fz = 0.0
          
    return
end subroutine Parabolic_dx_Force_dUdr2


!>
subroutine PartToPartDistance(x1, y1, z1, x2, y2, z2, dx, dy, dz)
    ! input
    real(8) x1, y1, z1, x2, y2, z2
    ! output
    real(8) dx, dy, dz
    
    dx = x2 - x1
    dy = y2 - y1
    dz = z2 - z1
    
! #ifdef BOUNDARY_CONDITION_CUBIC                        
!     ! граничные условия
!     call BoundaryConditionCubic(dx, dy, dz, box)
! #endif   
    
    return
end subroutine PartToPartDistance


!************************************************************
!>	@brief	Подсчет сил объемного взаимодействия
!!	@param	dx		[IN]	
!!	@param	dy		[IN]	расстояние между атомами
!!	@param	dz		[IN]	
!!	@param	sig 		[IN]	диаметр атомов
!!	@param	rcut		[IN]	радиус обрезки
!!	@param	eps		[IN]	энергия невалентного взаимодействия
!!	@param	solvent		[IN]	качество расворителя
!!	@param	Fx1		[OUT]	
!!	@param	Fy1		[OUT]	проекции силы действующей на первую частицу
!!	@param	Fz1		[OUT]	
  subroutine Lennard_Jones_Force(dx, dy, dz,  sig, rcut, eps, solvent, &
				  & Fx1, Fy1, Fz1)
!************************************************************
    implicit none
    ! входные
    real(8) dx, dy, dz
    real(8) sig,rcut,eps
    real(8) solvent
    ! выходные
    real(8) Fx1, Fy1, Fz1 !, Fx2, Fy2, Fz2
    ! рабочие
    real(8) r,F
    
    ! предупредать исключительные случаи:
    F = 0.0
    
    ! находим расстояние между частицами:
    r = sqrt(dx**2 + dy**2 + dz**2)
    
    !if (r.le.sig*2**(1/6)) then
    !	if (r.gt.0.78) then
    !		F = 4.0*eps*(-12*(sig**12/r**13) + 6*(sig**6/r**7))
    !	else
    !		F = 1000.0 * (r - 1.0)
    !	endif
    !endif
    if (r.le.1.0) F = eps * (r - 1.0)
    !if (r.le.sig*2**(1/6)) F = eps * (r - 1.0)
    
    ! sigma*2^(1/6) < r <= rcut
    !if (r.gt.sig*2**(1/6).and.r.le.rcut) F = solvent*sig*(1./r**2 + 2./rcut**2 - (3.*r**2)/rcut**4)
    if (r.gt.1.0.and.r.le.rcut) F = solvent*sig*(1./r**2 + 2./rcut**2 - (3.*r**2)/rcut**4)
    !if (r.gt.sig*2**(1/6).and.r.le.rcut) F = 4.0*solvent*(-12*(sig**12/r**13) + 6*(sig**6/r**7))

    F = F/r
    
    ! находим проекции силы:
    Fx1 = F*dx
    Fy1 = F*dy
    Fz1 = F*dz

    return
    end subroutine Lennard_Jones_Force
!**************************************************************

!**************************************************************
!<	@brief	Подсчет силы валентных взаимодействий
!!	@param	dx	[IN]
!!	@param	dy	[IN]	проекции расстояния между атомами
!!	@param	dz	[IN]
!!	@param	k_bond	[IN]	жесткость связи
!!	@param	lbond	[IN]	равновесная длина связи
!!	@param	Fx1	[OUT]
!!	@param	Fy1	[OUT] проекции силы на левую точку
!!	@param	Fz1	[OUT]
    subroutine Harmonic_Force(dx, dy, dz, kbond, lbond, &
			      & Fx1, Fy1, Fz1)
!**************************************************************
        implicit none
        ! входные
        real(8) dx, dy, dz, kbond, lbond
        ! выходные
        real(8) Fx1, Fy1, Fz1
        ! вспомогательные, для подсчета расстояния и силы:
        real(8) r, F 
        
        ! находим расстояние между частицами:
        r = sqrt(dx**2 + dy**2 + dz**2)

        ! считаем нормированную силу:
        F = kbond * (r - lbond)/r

        ! иситаем проекции силы:
        Fx1 = F*dx
        Fy1 = F*dy
        Fz1 = F*dz
         
        return
    end subroutine Harmonic_Force
!**************************************************************

!
!>	@brief	расчет силы угловых взаимодействий
!!	@param	k_ang	[IN]	жесткость угла
!!	@param	av_ang	[IN]	равновесный угол
!!  	@param 	x_l 	[IN]
!!  	@param	y_l 	[IN] левая точка
!!  	@param 	z_l 	[IN] 
!!  	@param 	x_c 	{IN]
!!  	@param 	y_c 	[IN] средняя точка
!!  	@param 	z_c 	[IN] 
!!  	@param 	x_r 	[IN]
!!  	@param 	y_r 	[IN] правая точка
!!  	@param 	z_r 	[IN]
!!	@param	Fx1	[OUT]
!!	@param	Fy1	[OUT] проекции силы на левую точку
!!	@param	Fz1	[OUT]
!!	@param	Fx2	[OUT]
!!	@param	Fy2	[OUT] проекции силы на правую точку
!!	@param	Fz2	[OUT]

subroutine Angle_Force(k_ang, av_ang, &
                     & dx, dy, dz, &
		     & x_l, y_l, z_l, & 
		     & x_c, y_c, z_c, &
		     & x_r, y_r, z_r, &
		     & Fx1, Fy1, Fz1)
      ! входные
      real(8) k_ang, av_ang
      real(8) x_l, y_l, z_l
      real(8) x_c, y_c, z_c
      real(8) x_r, y_r, z_r
      !real(8) box
      ! выходные
      real(8) Fx1, Fy1, Fz1
      ! локальные
      real(8) ang_val, F, r
      real(8) dx, dy, dz
      
      ! расстояние между двумя точками:
!       dx = x_l - x_r
!       dy = y_l - y_r
!       dz = z_l - z_r
!       
! #ifdef BOUNDARY_CONDITION_CUBIC                 
!       ! граничные условия
!       call BoundaryCondition(dx, dy, dz, box)      
! #endif
      
      r = sqrt(dx**2 + dy**2 + dz**2)
	
      ! получаем угол
      call GetAngle(x_l, y_l, z_l, &
		  & x_c, y_c, z_c, &
		  & x_r, y_r, z_r, ang_val)
	    
      F = k_ang*(ang_val - av_ang)/r
	    
          
      ! находим проекции сил:  
        Fx1 = F*dx
        Fy1 = F*dy
        Fz1 = F*dz
 	    
      return
end subroutine Angle_Force

!!!!!!!!!!!!!!!!!
!>  @brief получаем косинус угла по трем точкам
!!  @param x_l [IN]
!!  @param y_l [IN] левая точка
!!  @param z_l [IN] 
!!  @param x_c {IN]
!!  @param y_c [IN] средняя точка
!!  @param z_c [IN] 
!!  @param x_r [IN]
!!  @param y_r [IN] правая точка
!!  @param z_r [IN]
!!  @param cos_val	[OUT] косинус угла между двумя векторам
!!!!!!!!!!!!!!!!!
subroutine GetCos(x_l, y_l, z_l, x_c, y_c, z_c, x_r, y_r, z_r, cos_val)
	real(8)	x_l, y_l, z_l, x_c, y_c, z_c, x_r, y_r, z_r
	real(8) cos_val
	! work
	real(8)	r_x, r_x2, r_y, r_y2, r_z, r_z2
	real(8) d1, d2
	
	! находим координаты векторов:
	r_x = x_l - x_c 
	r_y = y_l - y_c
	r_z = z_l - z_c
							
	r_x2 = x_r - x_c
	r_y2 = y_r - y_c
	r_z2 = z_r - z_c

	! находим длины векторов:
	d1 = sqrt(r_x*r_x + r_y*r_y + r_z*r_z)
	d2 = sqrt(r_x2*r_x2 + r_y2*r_y2 + r_z2*r_z2)
	
	! находим  
	cos_val = (r_x*r_x2 + r_y*r_y2 + r_z*r_z2)/(d1*d2)
end subroutine GetCos


!!!!!!!!!!!!!!!!!
!>  @brief получаем угол по трем точкам
!!  @param x_l [IN]
!!  @param y_l [IN] левая точка
!!  @param z_l [IN] 
!!  @param x_c {IN]
!!  @param y_c [IN] средняя точка
!!  @param z_c [IN] 
!!  @param x_r [IN]
!!  @param y_r [IN] правая точка
!!  @param z_r [IN]
!!  @param ang_val	[OUT] угол между двумя векторам
!!!!!!!!!!111
subroutine GetAngle(x_l, y_l, z_l, x_c, y_c, z_c, x_r, y_r, z_r, ang_val)
	real(8)	x_l, y_l, z_l, x_c, y_c, z_c, x_r, y_r, z_r
	real(8) ang_val
	! work
	real(8)	cos_val
	
	! находим координаты векторов:
	call GetCos(x_l, y_l, z_l, x_c, y_c, z_c, x_r, y_r, z_r, cos_val)
	
	! находим угол между векторами 
	ang_val = acos(cos_val)*180.0/3.1415926535897934
	
end subroutine GetAngle

end module Potential
