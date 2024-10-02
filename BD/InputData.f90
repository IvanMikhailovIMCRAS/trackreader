!> @brief	Работа с входными файлами и данными
!! @date		09.03.2016
!!
!! @author	Михайлов И.В., Шавыкин О.В.
! Расскоментировать если нужно убирать вращение как целого на каждом шаге:
!#define ANTI_GYRATION
Module InputData
    Use CommonModule
    Use TestErrors 
    Use Surface
    Use Potential
Contains 
#ifdef HI_CHOLESKY
#define HYDRODYNAMIC_INTERACTION
#endif
#ifdef HI_FIXMAN
#define HYDRODYNAMIC_INTERACTION
#endif
#ifdef HI_TEA
#define HYDRODYNAMIC_INTERACTION
#endif

!> @brief получаем информацию о входных файлах
!!        а также сколько нам нужно памяти для чтения данных
subroutine ScanInputFiles(num_atom, num_bond, num_angles, &
                        & num_group, max_group_elem, num_taboo, num_surface, num_potential)
    integer(4) num_atom, num_bond, num_angles
    integer(4) num_group, num_taboo
    integer(4) num_surface, num_potential
    integer(8) max_group_elem
    
    ! локальные:
    logical(1) file_exists
    character(5) temp_word
    integer(4) temp_int, temp_int2
    
    ! Открываем выходной файл
    open(n_infor,file='INFOR')
    write(n_infor,'(A,A,A,A,A)') '<<< ', PROGRAM_NAME, ' version (', &
                & VERSION_STR(1:len_trim(VERSION_STR)), ') /I. Mikhailov, O. Shavykin/ (c) >>>'				      
    close(n_infor)
    
    ! зануляем значения всех переменных:
    num_atom = 0
    num_bond = 0
    num_angles = 0
    num_group = 0
    max_group_elem = 0
    num_taboo = 0
    
    !  1. считываем файл COORD
    INQUIRE(FILE="COORD", EXIST=file_exists) 
    if (file_exists.eqv..true.) then
	open(n_coord,file='COORD')
	read(n_coord,*) temp_word, num_atom
        close(n_coord)
    end if

    ! 2. считываем файл BONDS
    INQUIRE(FILE="BONDS", EXIST=file_exists) 
    if (file_exists.eqv..true.) then
	open(n_bonds,file='BONDS')
	read(n_bonds,*) temp_word, num_bond
	close(n_bonds)
    end if

    ! 3. считываем ANGL
    INQUIRE(FILE="ANGL", EXIST=file_exists) 
    if (file_exists.eqv..true.) then
      open(n_angl, file='ANGL')
      read(n_angl,*) temp_word, num_angles
      close(n_angl)
    end if
    
    ! 4. считываем GROUP    
    ! читаем группы атомов:
    INQUIRE(FILE="GROUP", EXIST=file_exists) 
    if (file_exists.eqv..true.) then
 	open(n_group, file='GROUP')
 	ioer = 0
 	! сканируем файл и считываем число групп
        do while (ioer.eq.0)
            ! ищем ключевое слово
            read(n_group,*,iostat=ioer) temp_word, temp_int
            
            if (temp_word(1:len_trim(temp_word)).eq.'GROUP'.and.ioer.eq.0) then
                ! нашли группу атомов
                num_group = num_group + 1
                if (temp_int.gt.max_group_elem) max_group_elem = temp_int
            
                ! пролистаем элементы
                do i=1, temp_int
                    read(n_group,*,iostat=ioer) temp_int2
                    if (ioer.ne.0) call ERROR(ERR_CODE_READ_GROUP)
                end do
            end if
        end do
        
        close(n_group)
    end if
    
    ! 5. считываем TABOO
    INQUIRE(FILE="TABOO", EXIST=file_exists) 
    if (file_exists.eqv..true.) then
	open(n_taboo, file='TABOO')
	read(n_taboo,*) temp_word, num_taboo
	close(n_taboo)
    end if
    
    ! 6. считываем SURFACE
    INQUIRE(FILE="SURFACE", EXIST=file_exists) 
    if (file_exists.eqv..true.) then
      open(n_surface, file='SURFACE')
      read(n_surface,*) temp_word, num_surface
      close(n_surface)
    end if
    
    ! 7. считываем FORCE
    INQUIRE(FILE="FORCE", EXIST=file_exists) 
    if (file_exists.eqv..true.) then
      open(n_potential, file='FORCE')
      read(n_potential,*) temp_word, num_potential
      close(n_potential)
    end if
    
    return
end subroutine ScanInputFiles


! считываем координаты частиц и их типы
subroutine ScanCoordType(num_atom, num_type, type_name)
    ! входной:
    integer(4) num_atom
    
    ! выходные:
    integer(4) num_type
    character(2), dimension(:), allocatable ::  type_name
    ! локальные:
    real(8) x, y, z
    logical(1) file_exists
    character(5) temp_word
    character(2) temp_type
    character(2) temp_type_arr(1:num_atom)
    integer(4) i, j
    logical(1) bfind
    integer(4) ioer
    !
    ioer = 0
    
    num_type = 0
    ! считываем координаты частиц
    INQUIRE(FILE="COORD", EXIST=file_exists) 
    if (file_exists.eqv..true.) then
	open(n_coord,file='COORD')
	
	read(n_coord,*) temp_word, num_atom
	
	do i=1,num_atom
	    read(n_coord,*,iostat=ioer) x, y, z, temp_type
	    
	    if (ioer.ne.0) then
                call ERROR(ERR_CODE_READ_COORD)
	    endif
	    
	    bfind = .false.
	    do j = 1, num_type
                if (temp_type(1:len_trim(temp_type)).eq.temp_type_arr(j)) then
                    bfind = .true.
                endif
	    enddo
	    
	    ! если нет такого, то добавляем!
	    if (bfind.eqv..false.) then
                num_type = num_type + 1
                temp_type_arr(num_type) = temp_type
	    endif
	enddo
	
	close(n_coord)
    else	
	call ERROR(ERR_CODE_OPEN_COORD)
    end if
    
    allocate(type_name(1:num_type))
    do j = 1, num_type
         type_name(j) = temp_type_arr(j)
    enddo
    
    open(n_infor,file='INFOR', access='append')
    write(n_infor, *) 
    write(n_infor,'(a,1x,I0)') 'число типов:',num_type
    write(n_infor, *) 
    close(n_infor)

    return
end subroutine ScanCoordType


! считываем координаты частиц и их типы
subroutine ReadCoord(num_atom, x, y, z, num_type, type_name, atom_type)
    ! входной:
    integer(4) num_atom
    integer(4) num_type
    character(2) type_name(num_type)
    ! выходные:
    real(8), dimension(:), allocatable :: x, y, z !x(num_atom), y(num_atom), z(num_atom)
    integer(4), dimension(:), allocatable :: atom_type
    ! локальные:
    logical(1) file_exists
    character(5) temp_word
    character(2) temp_type
    integer(4) i, j

    ! выделяем память:
    allocate(x(1:num_atom))
    allocate(y(1:num_atom))
    allocate(z(1:num_atom))
    allocate(atom_type(1:num_atom))
    
 ! считываем координаты частиц
    INQUIRE(FILE="COORD", EXIST=file_exists) 
    if (file_exists.eqv..true.) then
	open(n_coord,file='COORD')
	open(n_infor,file='INFOR', access='append')
	
	read(n_coord,*) temp_word, num_atom
	write(n_infor,'(a,1x,I0)') 'число атомов:',num_atom
	
	do i=1,num_atom
	    read(n_coord,*) x(i),y(i),z(i), temp_type
	    write(n_infor,*) x(i),y(i),z(i), temp_type
	    
	    do j = 1, num_type
                if (temp_type(1:len_trim(temp_type)).eq.type_name(j)) then
                    atom_type(i) = j
                endif
	    enddo
	enddo
	write(n_infor,'(a)') 'атомы считаны'
	write(n_infor, *) 
	close(n_coord)
	close(n_infor)
    else	
	call ERROR(ERR_CODE_OPEN_COORD)
    end if

    return
end subroutine ReadCoord

! считываем координаты частиц и их типы
subroutine ReadSurface(num_surface, surf)
    ! входной:
    integer(4) num_surface
    ! выходные:
    type(surface_info), dimension(:), allocatable :: surf 
    ! локальные:
    logical(1) file_exists
    character(255) temp_surf
    integer(4) ioer
    integer(4) i, j

    ! выделяем память:
    allocate(surf(1:num_surface))
    
    open(n_infor,file='INFOR', access='append')
 ! считываем координаты частиц
    INQUIRE(FILE="SURFACE", EXIST=file_exists) 
    if (file_exists.eqv..true.) then
	open(n_surface,file='SURFACE')
	
	read(n_surface,*) temp_surf, num_surface
	write(n_infor,'(a,1x,I0)') 'число поверхностей:',num_surface
	
	do i=1,num_surface	    
	    call SurfaceRecognition(n_surface, temp_surf, surf(i))
	    
	    if (ioer.ne.0) call ERROR(ERR_CODE_READ_SURFACE)
	    
	    write(n_infor,*) temp_surf(1:len_trim(temp_surf)), surf(i)%type_code, surf(i)%name(1:len_trim(surf(i)%name)), surf(i)%param
	    
	    do j = 1, i-1
                if (surf(j)%name(1:len_trim(surf(j)%name)).eq.surf(i)%name(1:len_trim(surf(i)%name))) then
                    call ERROR(ERR_CODE_SAME_SURFACE)
                endif
	    enddo
	enddo
	write(n_infor,'(a)') 'поверхности считаны'
	write(n_infor, *) 
	close(n_surface)
	
    else	
	num_surface = 0
	write(n_infor, '(a)') 'файл SURFACE отсутствует'
    end if
    
    close(n_infor)

    return
end subroutine ReadSurface

! считываем координаты частиц и их типы
subroutine ReadPotential(num_potential, num_type, type_name, num_surface, surf, forces)
    ! входной:
    integer(4) num_potential
    integer(4) num_type
    character(LEN_TYPE_NAME) type_name(num_type)
    integer(4) num_surface
    type(surface_info) surf(num_surface)
    ! выходные:
    type(force_info), dimension(:), allocatable :: forces
    ! локальные:
    logical(1) file_exists
    character(255) temp_force
    integer(4) ioer
    integer(4) i, j

    ! выделяем память:
    allocate(forces(1:num_potential))
    
    open(n_infor,file='INFOR', access='append')
 ! считываем координаты частиц
    INQUIRE(FILE="FORCE", EXIST=file_exists) 
    if (file_exists.eqv..true.) then
	open(n_force,file='FORCE')
	
	read(n_force,*) temp_force, num_force
	write(n_infor,'(a,1x,I0)') 'число потенциалов:',num_force
	
	do i=1,num_force	    
	    call PotentialRecognition(n_force, num_type, type_name, num_surface, surf, temp_force, forces(i))
	    
	    if (ioer.ne.0) then
                close(n_infor)
                call ERROR(ERR_CODE_READ_POTENTIAL)
            endif
	    
	    write(n_infor,'(I0,A)',advance='no') i, ')'
	    write(n_infor,*) temp_force(1:len_trim(temp_force)), forces(i)%type_code 
	    write(n_infor,'(A)',advance='no') 'типы:'
	    do j=1, forces(i)%num_type
                write(n_infor,'(3x,A)',advance='no') type_name(forces(i)%types(j))
            enddo
            write(n_infor,*)
	    write(n_infor,'(A,3x)',advance='no') 'параметры:'
	    write(n_infor,*) forces(i)%param
	    
! 	    do j = 1, i-1
!                 if (surf(j)%name(1:len_trim(surf(j)%name)).eq.surf(i)%name(1:len_trim(surf(i)%name))) then
!                     close(n_infor)
!                     call ERROR(ERR_CODE_SAME_POTENTIAL)
!                 endif
! 	    enddo
	enddo
	write(n_infor,'(a)') 'потенциалы считаны'
	write(n_infor, *) 
	close(n_force)
	
    else	
	num_force = 0
	write(n_infor, '(a)') 'файл FORCE отсутствует'
    end if
    
    close(n_infor)

    return
end subroutine ReadPotential

! считываем координаты частиц и их типы
subroutine ReadFriction(num_atom, atom_type, num_type, type_name, friction)
    ! входной:
    integer(4) num_atom
    integer(4) atom_type(num_atom)
    integer(4) num_type
    character(LEN_TYPE_NAME) type_name(num_type)
    ! выходные:
    real(8), dimension(:), allocatable :: friction
    ! локальные:
    logical(1) file_exists
    character(255) temp_word, temp_type
    real(8) temp_friction
    integer(4) ioer
    integer(4) i, num
    logical(1) bfind2
    real(8) friction_type(num_type)

    ! выделяем память:
    allocate(friction(1:num_atom))
    
    friction(:) = 0.0
    
    open(n_infor,file='INFOR', access='append')
 ! считываем координаты частиц
    INQUIRE(FILE="FORCE", EXIST=file_exists) 
    if (file_exists.eqv..true.) then
	open(n_force,file='FORCE')
	
	write(n_infor,'(a,1x,I0)') 'коэффициенты трения:'
	
	ioer = 0
	do while (ioer.eq.0)	    
	    read(n_force, *, iostat=ioer) temp_word
	    
	    if (temp_word(1:len_trim(temp_word)).eq.'friction') then
                backspace(n_force)
                read(n_force, *, iostat=ioer) temp_word, temp_type, temp_friction
        
                bfind2 = .false.
                do i = 1, num_type
                    temp_word=type_name(i)
                    if (temp_type(1:len_trim(temp_type)).eq.temp_word(1:len_trim(temp_word))) then
                        bfind2 = .true.
                        num = i ! запомним
                        friction_type(num) = temp_friction
                        
                        write(n_infor,'(A,3x, F20.10)') temp_type(1:len_trim(temp_type)), friction_type(num)
                    endif
                enddo
            
                if (bfind2.eqv..false.) then
                    call ERROR(ERR_CODE_FRICTION_UNKNOWN_PART_TYPE)
                endif
	    endif
	enddo
	write(n_infor,'(a)') 'коэффициенты трения считаны'
	write(n_infor, *) 
	close(n_force)
	
    else	
	num_force = 0
	write(n_infor, '(a)') 'файл FORCE отсутствует'
    end if
    
    close(n_infor)
    
    do i = 1, num_atom
        friction(i) = friction_type(atom_type(i))
    enddo

    return
end subroutine ReadFriction


! запоминаем пары связанных частиц
subroutine ReadBonds(num_bond, bond1, bond2)
    ! входной:
    integer(4) num_bond
    ! выходные:
    integer(4), dimension(:), allocatable :: bond1, bond2!bond1(num_bond), bond2(num_bond)
    ! локальные:
    logical(1) file_exists
    character(5) temp_word
    integer(4) i
    
    ! выделяем память:
    allocate(bond1(1:num_bond))
    allocate(bond2(1:num_bond))
    ! считываем связи
    INQUIRE(FILE="BONDS", EXIST=file_exists) 
    if (file_exists.eqv..true.) then
	open(n_bonds,file='BONDS')
	open(n_infor,file='INFOR', access='append')
	
	read(n_bonds,*) temp_word, num_bond
	write(n_infor,'(a,1x,I0)') 'число связей:',num_bond
	do i=1,num_bond
	    read(n_bonds,*) bond1(i), bond2(i)
	    write(n_infor,*) bond1(i),bond2(i)
	enddo
	write(n_infor,'(a)') 'связи считаны'
	write(n_infor, *) 
	close(n_bonds)
	close(n_infor)
    else
	call ERROR(ERR_CODE_OPEN_BONDS)
    end if

    return
end subroutine ReadBonds

! считываем углы
subroutine ReadAngl(num_angles, angles_l, angles_c, angles_r)
    ! входной:
    integer(4) num_angles
    ! выходные:
    integer(4), dimension(:), allocatable :: angles_l, angles_c, angles_r
    !angles_l(num_angles), angles_c(num_angles), angles_r(num_angles)
    ! локальные:
    logical(1) file_exists
    character(5) temp_word
    integer(4) i
    
    ! выделяем память:
    allocate(angles_r(1:num_angles))
    allocate(angles_c(1:num_angles))
    allocate(angles_l(1:num_angles))
    
    open(n_infor,file='INFOR', access='append')
    ! считываем углы
    INQUIRE(FILE="ANGL", EXIST=file_exists) 
    if (file_exists.eqv..true.) then
      open(n_angl, file='ANGL')
      
      read(n_angl,*) temp_word, num_angles
      write(n_infor, '(a,1x,I0)') 'число углов:', num_angles
      
      do i=1, num_angles
	  read(n_angl,*) angles_r(i), angles_c(i), angles_l(i)
	  write(n_infor,*) angles_r(i), angles_c(i), angles_l(i)
      end do
    
      close(n_angl)
    
      write(n_infor,'(a)') 'углы считаны'
      write(n_infor, *) 
    else
	num_angles = 0
	write(n_infor, '(a)') 'файл ANGL отсутствует'
    end if
    close(n_infor)
    
    return
end subroutine ReadAngl


subroutine ReadGroup(num_atom, num_group, max_group_elem, natom_group, &
                    & group_x, group_y, group_z, &
                    &  num_atom_group, atom_group, bgroup_atom)
    integer(4) num_atom
    integer(4) num_group
    integer(8) max_group_elem
    ! выходные
    ! центры масс:
    real(8), dimension(:), allocatable :: group_x, group_y, group_z !group_x(num_group), group_y(num_group), group_z(num_group)
    ! число атомов в группе:
    integer(4), dimension(:), allocatable :: num_atom_group!(num_group)
    ! списки номеров атомов в группах
    integer(4), dimension(:,:), allocatable :: atom_group !(num_group, max_group_elem)
    ! свойство "быть в группе":
    logical(4), dimension(:), allocatable :: bgroup_atom!(num_atom)
    integer(4) natom_group 
    ! локальные
    logical(1) file_exists
    character(5) temp_word
    integer(4) ioer
    integer(4) j, i
    integer(4) temp_int
    
    
    ! свойство "быть в группе"
    allocate(bgroup_atom(1:num_atom))
    bgroup_atom(:) = .false.
    ! центры масс  групп атомов
    allocate(group_x(1:num_group))
    allocate(group_y(1:num_group))
    allocate(group_z(1:num_group))
    ! число атомов в группах
    allocate(num_atom_group(1:num_group))
    ! списки номеров атомов в группах
    allocate(atom_group(1:num_group,1:max_group_elem))
    
    natom_group = 0
    open(n_infor,file='INFOR', access='append')
    ! читаем группы атомов:
    INQUIRE(FILE="GROUP", EXIST=file_exists) 
    if (file_exists.eqv..true.) then
        write(n_infor, '(a,1x,I0)') 'число групп с заданным центром масс:', num_group
 	open(n_group, file='GROUP')
 	        
        j = 0
        ioer = 0

        ! сканируем файл
        do while (ioer.eq.0)
            ! ищем ключевое слово
            read(n_group,*,iostat=ioer) temp_word
            if (temp_word.eq.'GROUP'.and.ioer.eq.0) then
                ! нашли группу атомов
                j = j + 1
                ! возвращаемся, чтобы прочитать информацию о группе
                backspace(n_group)
                ! считываем число элементов и центр масс группы
                read(n_group,*,iostat=ioer) temp_word, temp_int, group_x(j), group_y(j), group_z(j)
            
                ! запоминаем число элементов в группе:
                num_atom_group(j) = temp_int
                ! читаем номера атомов
                do i=1, num_atom_group(j)
                    read(n_group,*,iostat=ioer) atom_group(j,i)
                    
                    if (ioer.ne.0) call ERROR(ERR_CODE_READ_GROUP_ELEM)
                    ! запоминаем, что атом состоит в группе:
                    bgroup_atom(atom_group(j,i)) = .true.
                end do
            end if
        end do
        close(n_group)
        
        ! подсчитываем число атомов входящих в группы
        do j=1, num_group 
            natom_group = natom_group + num_atom_group(j)
        end do
        write(n_infor,'(a)') 'группы с заданным центром масс считаны'
	write(n_infor, *)
    else
	! нет файла
	write(n_infor, '(a)') 'файл GROUP отсутствует'	
    end if
    close(n_infor)

    return
end subroutine ReadGroup




! читаем запрещенные объемные взаимодействия
subroutine ReadTaboo(num_taboo, taboo1, taboo2)
    ! входной:
    integer(4) num_taboo
    ! выходные:
    integer(4), dimension(:), allocatable :: taboo1, taboo2!taboo1(num_taboo), taboo2(num_taboo)
    ! локальные:
    logical(1) file_exists
    character(5) temp_word
    integer(4) i
    
    ! выделяем память:
    allocate(taboo1(1:num_taboo))
    allocate(taboo2(1:num_taboo))
    
    !======== Список запрещенных взаимодействий ===========
    open(n_infor,file='INFOR', access='append')
    ! читаем запрещенные объемные взаимодействия:
    INQUIRE(FILE="TABOO", EXIST=file_exists) 
    if (file_exists.eqv..true.) then
	open(n_taboo, file='TABOO')
	read(n_taboo,*) temp_word, num_taboo
	write(n_infor, '(a,1x,I0)') 'число запрещенных объемных взаимодействий:', num_taboo
	
	do i=1,num_taboo
	    read(n_taboo,*) taboo1(i),taboo2(i)
	    write(n_infor,*) taboo1(i),taboo2(i)
	enddo    
	close(n_taboo)
	
        write(n_infor,'(a)') 'запрещенные объёмные взаимодействия считаны'
	write(n_infor, *)
    else
	num_taboo = 0
	write(n_infor, '(a)') 'файл TABOO отсутствует'	
	write(n_infor, *)
    end if
    
    close(n_infor)
    
    return
end subroutine ReadTaboo

!> @brief читаем параметры системы:
subroutine ReadCONTR(dt, step, kT, &
		    & delta_r_sphere, max_unbond_atoms, &
		    & bondary_condition, h_hi, hradius, &
		    & recalc_diffusion_tensor, hi_method, n_anti_gyration,&
		    & print_track, print_stats, print_picture, print_tremor)
    ! выходные
    real(8)  dt, kT
    integer(8) step
    integer(4) max_unbond_atoms
    real(8)  delta_r_sphere 
    real(8)	bondary_condition, h_hi, hradius
    integer(4) n_anti_gyration
    integer(4) recalc_diffusion_tensor
    integer(4) hi_method
    integer(8)	print_track, print_tremor
    integer(8) print_stats, print_picture
    ! рабочий
    integer(4) ioer, ioer2
    real(8)  r_sphere
    logical(1) file_exists
    character(255) temp
    
    dt = 0.0  ! /шаг моделирования/
    step = 0  ! /число шагов моделирования/
    kT = 0.0  ! /энергия в дисперсии случайной силы/ 
    delta_r_sphere = 0.0  ! /параметр метода двойных списков/
    print_track = 0 ! /периодичность вывода координат частиц/
    print_picture = 0 ! /периодичность вывода картинок/
    print_stats = 0 ! /периодичность вывода энергии системы/
    print_tremor = 0  
    bondary_condition = 0.0   !  /размер коробки/
    h = 0.0 ! /гидродинамический радиус/
    hi_tensor_period =  1 ! /период перерасчета тензора инерции (только для метода cholesky)/
    hi_method  = HI_METHOD_CHOLESKY ! /метод учета гидродинамических взаимодействий |cholesky | chebyshev | tea | krylov |/
    n_anti_gyr = 0 ! /относительно какого атома будет происходить удаление вращения/
    
    ! читаем параметры системы:
    INQUIRE(FILE="CONTR", EXIST=file_exists) 
    
    if (file_exists.eqv..true.) then
        
	ioer = 0
	ioer2 = 0
	open(n_infor, file='INFOR', access='append')
	write(n_infor,'(a)') 'параметры моделирования:'
	open(n_contr, file='CONTR')
	
	! считываем параметры моделирования
	do while(ioer.eq.0)
                write(temp,'(A)') ' '
		read(n_contr,*,iostat=ioer) temp
		
		if (temp(1:len_trim(temp)).eq.'dt') then
			backspace(n_contr)
			read(n_contr,*,iostat=ioer2) temp, dt
                        write(n_infor,'(F15.7, a)') dt, ' - шаг моделирования'
                endif
                if (temp(1:len_trim(temp)).eq.'step') then
			backspace(n_contr)
                        read(n_contr,*,iostat=ioer2) temp, step
                        write(n_infor,'(I15, a)') step, ' - число шагов моделирования'
                endif
                if (temp(1:len_trim(temp)).eq.'kT') then
			backspace(n_contr)
                        read(n_contr,*,iostat=ioer2) temp, kT
                        write(n_infor,'(F15.3, a)') kT, ' - энергия в дисперсии случайной силы'
                endif
        
                if (temp(1:len_trim(temp)).eq.'delta_r_sphere') then
			backspace(n_contr)
                        read(n_contr,*,iostat=ioer2) temp, delta_r_sphere
                        write(n_infor,'(F15.3, a)') delta_r_sphere, &
                                &  ' - параметр метода двойных списков'
                endif
                
                if (temp(1:len_trim(temp)).eq.'limit_packing') then
                        ! максимальное число шариков в сфере:
			backspace(n_contr)
                        read(n_contr,*,iostat=ioer2) temp, max_unbond_atoms
                        !if (max_unbond_atoms.gt.num_atom-1) max_unbond_atoms = num_atom - 1
                        !обычно: max_unbond_atoms=857
                        write(n_infor,'(I15, a)') max_unbond_atoms, &
                                & ' - верхняя граница числа шариков в сфере радиуса r_sphere'
                endif
                    
                if (temp(1:len_trim(temp)).eq.'print_track') then
			backspace(n_contr)
                        read(n_contr,*,iostat=ioer2) temp, print_track
                        write(n_infor,'(I15, a)') print_track, &
                                        &  ' шагов - интервал между печатью траектории'
                endif
                if (temp(1:len_trim(temp)).eq.'print_picture') then
			backspace(n_contr)
                        read(n_contr,*,iostat=ioer2) temp, print_picture
                        write(n_infor,'(I15, a)') print_picture, &
                                        &  ' шагов - интервал между печатью снимков'
                endif
                if (temp(1:len_trim(temp)).eq.'print_stats') then
			backspace(n_contr)
                        read(n_contr,*,iostat=ioer2) temp, print_stats
                        write(n_infor,'(I15, a)') print_stats, &
                                        &   ' шагов - интервал между печатью статистики'
                endif
                if (temp(1:len_trim(temp)).eq.'print_tremor') then
			backspace(n_contr)
                        read(n_contr,*,iostat=ioer2) temp, print_tremor
                        write(n_infor,'(I15, a)') print_tremor, & 
                            & ' шагов - интервал между печатью сил, &
                            & действующих на фиксированные'
                endif
#ifdef BOUNDARY_CONDITION_CUBIC
                if (temp(1:len_trim(temp)).eq.'box') then
			backspace(n_contr)
                        read(n_contr,*,iostat=ioer2) temp, bondary_condition

                        write(n_infor, '(F15.3, a)') bondary_condition, ' - размер коробки'
                endif
#endif         
 
        

		
#ifdef HYDRODYNAMIC_INTERACTION
		hi_method = HI_METHOD_CHOLESKY
		if (temp(1:len_trim(temp)).eq.'hi_method') then
			backspace(n_contr)
			read(n_contr,*,iostat=ioer2) temp, temp
			if (temp(1:len_trim(temp)).eq.'cholesky') hi_method = HI_METHOD_CHOLESKY
			if (temp(1:len_trim(temp)).eq.'chebyshev') hi_method = HI_METHOD_CHEBYSHEV
			if (temp(1:len_trim(temp)).eq.'fixman') hi_method = HI_METHOD_FIXMAN
			if (temp(1:len_trim(temp)).eq.'tea') hi_method = HI_METHOD_TEA
			if (temp(1:len_trim(temp)).eq.'krylov') hi_method = HI_METHOD_KRYLOV
			if (temp(1:len_trim(temp)).eq.'rpy_test') hi_method = HI_METHOD_RPY_TEST
			
			select case(hi_method)
                            case(HI_METHOD_CHOLESKY)
                                write(n_infor, *) 'алгоритм: разложение Холецкого'
                            case(HI_METHOD_CHEBYSHEV)
                                write(n_infor, *) 'алгоритм: метод Чебышева'
                            case(HI_METHOD_FIXMAN)
                                write(n_infor, *) 'алгоритм: метод Чебышева-Фиксмана'
                            case(HI_METHOD_TEA)
                                write(n_infor, *) 'алгоритм: метод Гейера (TEA)'
                            case(HI_METHOD_KRYLOV)
                                write(n_infor, *) 'алгоритм: метод Крылова'
                            case(HI_METHOD_RPY_TEST)
                                write(n_infor, *) 'алгоритм: тестирование тензора Ротне-Прагер-Ямакавы'
                            case default
                                write(n_infor, *) 'алгоритм: неизвестный'
                                stop
                        end select
		endif
                if (temp(1:len_trim(temp)).eq.'h') then
                        backspace(n_contr)
			read(n_contr,*,iostat=ioer2) temp, h_hi
                        write(n_infor, '(F15.3, a)') h_hi, ' - гидродинамический параметр'
                        hradius = h_hi*sqrt(3.1415926535897934/3)*lbond
                        write(n_infor, '(F15.3, a)') hradius, ' - гидродинамический радиус'
                endif
                if (temp(1:len_trim(temp)).eq.'hi_tensor_period') then
                    backspace(n_contr)
                    read(n_contr,*,iostat=ioer2) temp, recalc_diffusion_tensor
                    write(n_infor, '(I15, a)') recalc_diffusion_tensor, &
                    & ' шагов - интервал между пересчётом тензора Ротне-Прагер-Ямакавы'
                endif
#endif
                 
#ifdef ANTI_GYRATION
                if (temp(1:len_trim(temp)).eq.'n_anti_gyr') then
                    backspace(n_contr)
                    read(n_contr,*,iostat=ioer2) temp, n_anti_gyration
                    if (n_anti_gyration.eq.0) then
                        write(n_infor, '(a)') 'удалять вращение относительно центра масс' 
                    else
                        write(n_infor, '(a, I15, a)') 'удалять вращение относительно ', & 
                                     & n_anti_gyration,'-го атома'
                    end if
                endif
#endif
                ! если файл CONTR содержит неправильный параметр:
                if (ioer2.lt.0) then
                    call ERROR(ERR_CODE_BAD_CONTR)
                end if
        enddo
        
	close(n_contr)
    
        write(n_infor,'(a)') 'параметры моделирования считаны'
        write(n_infor, *) 
        close(n_infor)
      
    else
      call ERROR(ERR_CODE_OPEN_CONTR)
    end if
    
    return
end subroutine ReadCONTR

!> @brief читаем параметры системы:
subroutine ReadMoveCONTR(num_group, group_x, group_y, group_z, &
                    &   num_motion, motion_group1, motion_group2, &
                    &   motion_start_x1, motion_start_y1, motion_start_z1,  &
                    &   motion_start_x2, motion_start_y2, motion_start_z2, &
                    &   motion_dr1, motion_dr2, motion_min_dr, motion_max_dr, &
                    &   motion_dx, motion_dy, motion_dz)
    ! входные
    integer(4) num_group
    real(8) group_x(num_group), group_y(num_group), group_z(num_group)
    ! выходные
    integer(4) num_motion
    integer(4), dimension(:), allocatable ::  motion_group1, motion_group2
    real(8), dimension(:), allocatable ::  motion_start_x1, motion_start_y1, motion_start_z1
    real(8), dimension(:), allocatable ::  motion_start_x2, motion_start_y2, motion_start_z2
    real(8), dimension(:), allocatable ::  motion_dr1, motion_dr2
    real(8), dimension(:), allocatable ::  motion_min_dr, motion_max_dr
    real(8), dimension(:), allocatable ::  motion_dx, motion_dy, motion_dz
    
    ! рабочий
    integer(4) ioer
    logical(1) file_exists
    character(255) temp
    
    integer(4) n1, n2
    real(8) dr1, dr2, min_dr, max_dr
    real(8) dx, dy, dz, dr
    
    integer(4) counter
    
    ! читаем параметры системы:
    INQUIRE(FILE="CONTR", EXIST=file_exists) 
    
    if (file_exists.eqv..true.) then
        
	ioer = 0
	open(n_contr, file='CONTR')
	
	num_motion = 0
	! считываем параметры моделирования
	do while(ioer.eq.0)
                write(temp, '(A)')  ' '
		read(n_contr,*,iostat=ioer) temp
		
		if (temp(1:len_trim(temp)).eq.'speed') then
                    backspace(n_contr)
                    read(n_contr,*,iostat=ioer) temp, n1, n2, dr1, dr2, min_dr, max_dr  
			
                    num_motion = num_motion + 1
                endif
        enddo
        close(n_contr)
        
        !> память
        allocate(motion_group1(num_motion))
        allocate(motion_group2(num_motion))
                        
        allocate(motion_start_x1(num_motion))
        allocate(motion_start_y1(num_motion))
        allocate(motion_start_z1(num_motion))
        allocate(motion_start_x2(num_motion))
        allocate(motion_start_y2(num_motion))
        allocate(motion_start_z2(num_motion))
                        
        allocate(motion_dr1(num_motion))
        allocate(motion_dr2(num_motion))
        
        allocate(motion_min_dr(num_motion))
        allocate(motion_max_dr(num_motion))
        
        allocate(motion_dx(num_motion))
        allocate(motion_dy(num_motion))
        allocate(motion_dz(num_motion))
        !>>>>>>>>>>>>>>
        
	open(n_infor, file='INFOR', access='append')
	write(n_infor, *)
	write(n_infor,'(a, 2x, I0)') 'число пар двигающихся групп:', num_motion
	
	open(n_contr, file='CONTR')
	
	counter = 0
	ioer = 0
	! считываем параметры моделирования
	do while(ioer.eq.0)
                write(temp, '(A)')  ' '
		read(n_contr,*,iostat=ioer) temp
		
		if (temp(1:len_trim(temp)).eq.'speed') then
                    backspace(n_contr)
                    read(n_contr,*,iostat=ioer) temp, n1, n2, dr1, dr2, min_dr, max_dr
                    
                    if  (min_dr.gt.max_dr) then
                        call ERROR(ERR_CODE_SPEED_MIN_MAX)
                    endif
			
                    counter = counter + 1
                    if (n1.ne.n2.and.ioer.eq.0) then
                        if (n1.gt.0.and.n1.le.num_group.and.n2.gt.0.and.n2.le.num_group) then
                            motion_group1(counter) = n1
                            motion_group2(counter) = n2
                        
                            motion_start_x1(counter) = group_x(n1)
                            motion_start_y1(counter) = group_y(n1)
                            motion_start_z1(counter) = group_z(n1)
                            motion_start_x2(counter) = group_x(n2)
                            motion_start_y2(counter) = group_y(n2)
                            motion_start_z2(counter) = group_z(n2)
                        
                            motion_dr1(counter) = dr1
                            motion_dr2(counter) = dr2
                            
                            motion_min_dr(counter) = min_dr
                            motion_max_dr(counter) = max_dr
                            
                           ! узнаем расстояние между группами
                            dx = group_x(n2) - group_x(n1)
                            dy = group_y(n2) - group_y(n1)
                            dz = group_z(n2) - group_z(n1)
            
                            dr = sqrt( dx**2 + dy**2 + dz**2 )
                            
                            ! начальное положение центров масс групп на недопустимом расстоянии:
                            if  (dr.lt.min_dr.or.dr.gt.max_dr) then
                                call ERROR(ERR_CODE_SPEED_ABROAD)
                            endif
            
                            ! это определяет также вклады направлений:
                            motion_dx(counter) = dx/dr
                            motion_dy(counter) = dy/dr
                            motion_dz(counter) = dz/dr 
                            
                            write(n_infor,'(A, I0, 2x, I0, 2x, A, F15.7, 2x, F15.7)') &
                                    & 'groups:', n1, n2, 'dr:', dr1, dr2
                            write(n_infor,'(10x, A, F15.7, 2x, F15.7, 2x, F15.7, 2x, A, F15.7, 2x, F15.7)') &
                                    & 'dx-dy-dz', motion_dx(counter), motion_dy(counter), motion_dz(counter),  &
                                    & 'min-max', min_dr, max_dr
                        else
                            call ERROR(ERR_CODE_SPEED_WRONG_GROUP)
                        endif
                    else
                        call ERROR(ERR_CODE_READ_SPEED)
                    endif
                        
                endif
        enddo
        
	close(n_contr)
    
        write(n_infor,'(a)') 'параметры пар двигающихся групп считаны'
        write(n_infor, *) 
        close(n_infor)
        
    else
      call ERROR(ERR_CODE_OPEN_CONTR)
    end if
    
    return
end subroutine ReadMoveCONTR

! заполняем метки <не иметь объемных взаимодействий>
subroutine LabeledFantomPairs(num_atom, num_bond, bond1, bond2, &
                            & num_angles, angles_l, angles_c, angles_r, &
                            & num_taboo, taboo1, taboo2, &
                            & atom_type, num_type, type_name, P_aa_matrix, label_taboo)
    ! input
    integer(4) num_atom, num_bond
    integer(4) bond1(num_bond), bond2(num_bond)
    integer(4) num_angles
    integer(4) angles_l(num_angles), angles_c(num_angles), angles_r (num_angles)
    integer(4) num_taboo
    integer(4), dimension(:), allocatable ::  taboo1, taboo2!taboo1(num_taboo), taboo2(num_taboo)
    integer(4) atom_type(num_atom)
    integer(4) num_type
    character(LEN_TYPE_NAME) type_name(num_type)
    integer(4) P_aa_matrix(num_type, num_type)
    ! output
    logical(4), dimension(:, :), allocatable :: label_taboo!(num_atom, num_atom)
    ! локальная
    integer(4) i, j
    
    
    ! метка <не иметь объемных взаимодействий>
    allocate(label_taboo(1:num_atom,1:num_atom))
    
    label_taboo(:,:) = .false. 
    ! запрещаем связям
    do i=1,num_bond
        label_taboo(bond1(i),bond2(i)) = .true.
        label_taboo(bond2(i),bond1(i)) = .true.
    enddo
    ! запрещаем углам:
    do i=1, num_angles
	! атомы являются концами одного угла:
	label_taboo(angles_l(i),angles_r(i)) = .true.
	label_taboo(angles_r(i), angles_l(i)) = .true.
    end do
    ! запрещаем "запрещенные пары"
    do i=1,num_taboo
        label_taboo(taboo1(i),taboo2(i)) = .true.
        label_taboo(taboo2(i),taboo1(i)) = .true.
    enddo  
    ! если между этими парами нет взаимодействий:
    do i=1,num_atom
        do j=i,num_atom
            if (P_aa_matrix(atom_type(i), atom_type(j)).eq.0) then
                label_taboo(i, j) = .true.
                label_taboo(j, i) = .true.
            endif
        enddo
    enddo 
    !>>>
    do i=1,num_type
        do j=i,num_type
            if (P_aa_matrix(i, j).eq.0) then
                write(*, '(A,A,A,A,A)') &
                        & 'WARNING: unbond potential between ', type_name(i), ' and ', type_name(j), 'does not exist'
            endif
        enddo
    enddo
    
    ! здесь можно освободить массивы taboo1, taboo2
    if (num_taboo.gt.0) then
        deallocate(taboo1)
        deallocate(taboo2)
    end if
    
    return
end subroutine LabeledFantomPairs    


End Module InputData
