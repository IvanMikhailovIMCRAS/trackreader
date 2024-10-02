!> @brief BrownianDynamic 2014 version (1.8)
!! 
!! Программа по моделированию крупномасштабных моделей одиночных полимерных молекул
!! методом Броуновской динамики.
!!
!! @n Для запуска программы необходимы следующие входные файлы:
!! @n COORD - содержит исходные координаты атомов
!! @n BONDS - содержит информацию о связях в молекуле
!! @n CONTR - задает параметры моделирования
!! @n Необязательные (!) входные файлы: 
!! @n GROUP - список групп атомов с постоянным центром масс
!! @n ANGL - список углов
!! @n TABOO - список запрещенных объёмных взаимодействий
!! @n Программа генерирует следующие выходные файлы:
!! @n INFOR - информация о запуске и выполнении программы
!! @n TRACK - траектории движения атомов
!! @n STATS - энергия системы, процент принятых конфигураций
!! @n TREMOR - силы действующие на группы частиц с постоянным центром масс
!!
!!
!! @author I. Mikhailov, O. Shavykin
!! @date 04.02.2015 -
!! @todo продумать возможные ошибки (функция ERROR(n))
!! @todo добавить тестирование параметров системы
!! @todo вывод STATS = 3 энергии + общая (вывод температуры, проверка импульса)
!! @todo проверка все ли частицы в коробке
!! @todo проверить существование циклов
!
!
! изменения в версиях:
! version (1) 
! version (1.1) заданы валентные, невалентные, угловые потенциалы
! version (1.2) силы для всех атомов считаются в одной функции
!	      UnBond_List теперь заносит всё в один общий список 
!
! version (1.3) добавлены граничные условия, разбит на несколько модулей
!		- BoxMuller.f90 
!    		- CommonModule.f90 
!		- TestErrors.f90 
!		- ForceEnergy.f90
!		- BrownianDynamic.f90
!		добавлен вывод в STATIS энергии и температуры
!		@date 18.12.2014
! 		
!		добавлены команды препроцесору на граничные условия
!		restart траектории
!		@date 03.02.2015
!
! version (1.4) добавлен файл TABOO - список запрещенных объёмных взаимодействий
!		label_bond заменен на label_taboo и теперь в нём 
!			отмечвются все пары, которые не имеют объёмных взаимодействий
!		@date	01.03.2015
!
! version (1.5) добавлены гидродинамические взаимодействия
!		@date	19.04.2015
!
!		добавлены фиксированные атомы, 
!		силы действующие на фиксированные выводятся в файл TREMOR,
!		добавлена функция для вывода снимков Picture_Ex
!		@date	05.05.2015
!
! version (1.6) добавлен учёт групп атомов с постоянным центром массив
!               Поскольку фиксированные атомы можно моделировать с помощью 
!               атома с постоянным центром массы, прежний механизм фиксированных
!               атомов был убран из программы
!               Файл FIX заменен на файл GROUP
!               Выходной файл TREMOR теперь будет выводить силы, действующие на
!               группу атомов
!               Исправлена ошибка с граничными условиями (добавлен #define в ForceEnergy.f90)
!               @date   21.10.2015
! version (1.7) добавлена возможность убирать вращение как целого на каждом шаге
!               модуль AntiGyration.f90
!               добавлены и затем удалены двойные списки для гидродинамических взаимодействий 
!               (оказались неэффективными)
!               @date   27.11.2015
!
!               Для гидродинамических взаимодействий добавлена частота пересчета 
!               тензора Ротне-Прагера-Ямакавы (считывается из файла CONTR)
!               Цикл для проверки необходимости пересчета двойных списков изменён 
!               (сделано ускорение)
!               @date   17.12.2015               
!               
!               Частота пересчёта тензора диффузии заменена алгоритмом пересчета тензора 
!               в случаях аналогичным пересчёту двойных списков (поздее отклонено, т.к. нельзя так делать)
!               Добавлен Метод Фиксмана расчета случайных перемещений при гидродинамических взаимодействиях
!               @date   19.12.2015
!               
!               Добавлен Метод Гейера (TEA) расчета случайных перемещений при гидродинамических взаимодействиях
!               @date   21.12.2015
!
!               замена директивы кубических граничных условий на BOUNDARY_CONDITION_CUBIC
!               добавлены цилиндрические граничные условия (директива BOUNDARY_CONDITION_CYLINDRICAL)
!               @date   03.03.2016
!
!               найдены ошибки при использование граничных условий: построение файлов .ENT и 
!               добавлен учёт цилиндрических условий при составлении двойных списков
!               для функций реализующих периодические граничные условия создан файл BoundaryCondition.f90
!               @TODO: для угловой энергии нет периодических условий
!               @date   04.03.2016
!
!               добавлены три новых модуля: 
!               - InputData - перенесены функции по работе с входными данными
!               - OutputData - перенесены функции по работе с выводом данных
!               - BoundaryCondition - перенесены функции для граничных условий
!               @date 09.03.2016
!
!               в файле InputData.f90 не хватало строки: !#define ANTI_GYRATION (для чтения из файла CONTR)
!               обновлены гидродинамические взаимодействия:
!               1. обновлен Makefile
!               2. для всех алгоритмов выделены общие части:
!               - Выделение памяти
!               - Подготовка перед стартом итераций 
!               - Перемещение на каждом шаге итераций
!               @date 09.09.2016
!
!               устранены проблемы:
!               1. в файле ForceEnergy.f90 неправильно рассчитывалось расстояние между частицами при расчете энергии
!               2. при отсутствии объемных взаимодействий не устанавливался num_bond = 0 (на некоторых ПК были ошибки)
!               обновлены гидродинамические взаимодействия:
!               в методе Фиксмана поиск собственных чисел заменен только поиском максимального и минимального 
!               собственных значений, что алгоритмически более простая задача (т.е. увеличилась скорость вычислений)
!               @date 30.05.2017
!
! version (1.8) Изменения коснулись учета гидродинамических взаимодействий:
!		1. Теперь вместо разных макросов для каждого метода остался один макрос HYDRODYNAMIC_INTERACTION
!		2. Добавлен числовой индетификатор для распознавания типа метода для учета гидродинамических взаимодействий
!			(числовые коды для него см. в CommonModule.f90)
!		3. Обновлены процедуры для работы со всеми методами учета HI
!			+ добавлен метод основанный на подстранствах Крылова
!			Доступные методы:
!			- с разложением Холецкого (cholesky)
!			- с разложением Чебышева (chebyshev)
!				метод Фиксмана разновидность метода Чебышева 
!				при условии оценки собственных чисел по способу предложенному Фиксманом
!			- разложение TEA (tea)
!			- подпространства Крылова (krylov)
!		4. По результатам тестов выявлено:
!		4.1.	Фантомные цепи
!			- метод TEA некорректен для фантомных взаимодействий
!			- метод Cholesky не всегда корректно работает
!			- метод Krylov совсем плох с QR разложением, иногда с Cholesky
!			- метод Фиксмана работает!
!               @date   08.03.2018-10.03.2018
!
! version (2.0) Новая эпоха Броуновской динамики! 
!               Анонс: будут типы атомов, возможность задавать поверхности и разные потенциалы, внешние поля
!
!               beta-версия
!               Типы атомов в файле COORD
!               @date   03.06.2018
!
!               beta-версия
!               - Удалено: 
!               Цилиндрические условия извлечены временно из кода для избежания путаницы
!               - Добавлены: 
!               Surface.f90
!               Potential.f90
!               - Изменения: 
!               BrownianDynamic.f90
!               ForceEnergy.f90
!               InputData.f90
!               @date   04.06.2018-07.06.2018
!
!               перемещение пар групп
!               новый CONTR - счет параметров по ключевым словам, убраны параметры потенциалов
!               TODO
!               потенциалы,углы
!               surface
!               @date   18.06.2018-20.06.2018
!
!               Добавлено:
!               - силы для углов
!               - обработка ошибок при чтении файла GROUP
!               - потенциалы в Potential
!               - в CONTR параметр limit_packing - максимальное число шариков в сфере
!               предупреждение: возможно до этого введения могли быть проблемы с объемными взаимодействиями
!               Исправлено:
!               - ошибка в CONTR при чтении speed
!               - ошибка в случае существования 'лишних' строк в CONTR
!               @date   25.09.2018
!               
!
! Расскоментировать если нужно убирать вращение как целого на каждом шаге:
!#define ANTI_GYRATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program BrownianDynamic
    ! список используемых модулей:
    use CommonModule	! Константы
    use BoxMuller ! генератор случайных чисел
    use TestErrors	! Обработка ошибок
    use Surface     ! библиотека поверхностей
    use Potential   ! библиотека потенциалов взаимодействий
    use InputData   ! обработка входнных данных из файлов
    use OutputData  ! работа с данными на выход
    use PeriodicCondition ! граничные условия
    use ForceEnergy	! Подсчет Энергии и сил

#ifdef HYDRODYNAMIC_INTERACTION
    use HydroDynamic
#endif

#ifdef ANTI_GYRATION
    use AntiGyration
#endif

    implicit none
               
    integer(4) i, j
 
    integer(4) num_atom, num_bond, num_group, num_angles, num_unbond, num_taboo
    integer(4) natom_group ! общее число "групповых атомов" 
    
    real(8), dimension(:), allocatable :: x, y, z
    integer(4), dimension(:), allocatable :: bond1, bond2
    integer(4), dimension(:), allocatable :: unbond1, unbond2
    integer(4), dimension(:), allocatable :: taboo1, taboo2
    integer(4), dimension(:,:), allocatable :: atom_group   
    integer(4), dimension(:), allocatable :: angles_r, angles_c, angles_l
    
    !
    integer(4), dimension(:), allocatable :: angles_p
    integer(4), dimension(:), allocatable :: bond_p
    integer(4), dimension(:,:), allocatable :: P_aa_matrix  ! номера потенциалов невалентных взаимодействий
    ! одначастичные взаимодействия:
    integer(4) num_onepart
    integer(4), dimension(:), allocatable :: onepart
    integer(4), dimension(:), allocatable :: onepart_p
    
    ! для поверхность-частица
    !integer(4) num_bond_sa
    integer(4) num_unbond_sa
    !
    !integer(4), dimension(:), allocatable :: bond_sa1, bond_sa2
    integer(4), dimension(:), allocatable :: unbond_sa1, unbond_sa2
    ! потенциалы
    !integer(4), dimension(:), allocatable :: bond_sa_p
    integer(4), dimension(:,:), allocatable :: P_sa_matrix  ! номера потенциалов невалентных взаимодействий
    
    logical, dimension(:), allocatable :: bgroup_atom
    logical, dimension(:,:), allocatable :: label_taboo
 
    real(8), dimension(:), allocatable :: Fx, Fy, Fz
    real(8), dimension(:), allocatable :: x_mem, y_mem, z_mem
    
    real(8), dimension(:), allocatable :: group_x, group_y, group_z
    integer(4), dimension(:), allocatable :: num_atom_group
    !
    integer(4) num_type
    integer(4), dimension(:), allocatable :: atom_type
    character(LEN_TYPE_NAME), dimension(:), allocatable :: type_name
    !
    integer(4) num_surface
    type(surface_info), dimension(:), allocatable :: surf
    
    integer(4) num_potential
    type(force_info), dimension(:), allocatable :: forces
    !!!, 
    integer(4) num_motion
    integer(4), dimension(:), allocatable ::  motion_group1, motion_group2
    real(8), dimension(:), allocatable ::  motion_start_x1, motion_start_y1, motion_start_z1
    real(8), dimension(:), allocatable ::  motion_start_x2, motion_start_y2, motion_start_z2
    real(8), dimension(:), allocatable ::  motion_dr1, motion_dr2
    real(8), dimension(:), allocatable ::  motion_min_dr, motion_max_dr
    real(8), dimension(:), allocatable ::  motion_dx, motion_dy, motion_dz
    !!!!!!!!!!!!!!!
    
    
    ! параметры системы из CONTR
    real(8) dt,kt,lbond,kbond,sig,eps,rcut
    real(8)	av_ang, k_ang
    real(8) solvent
    real(8) friction_contr
    real(8), dimension(:), allocatable :: friction
    integer(8) step,nstep
    integer(4) n_anti_gyration ! номер атома относительно, которого убираем
    integer(8) print_track, print_picture, print_stats, print_tremor
    integer(4) recalc_diffusion_tensor ! периодичность обновления тензора
    real(8) delta_r_sphere, r_sphere_contr
    !new
    real(8), dimension(:,:), allocatable :: r_sphere2
    !>>
    real(8)	box
    ! периодический угол:
    real(8) alpha_box
    real(8)	hradius, h_hi
    integer(4) hi_method
    
    ! =======================
    
    integer(4)	max_unbond_atoms, max_unbond_num
    integer(8) max_group_elem
    ! центр масс группы атомов:
    real(8) center_x, center_y, center_z 

    ! для чтения файлов
    character(5) temp_word
    integer(8)	temp_int
    integer(4) ioer
    
    ! перемещение
    real(8), dimension(:), allocatable ::  delta_coord_x, delta_coord_y, delta_coord_z
    real(8) sum_delta_x, sum_delta_y, sum_delta_z, sum_delta

#ifdef HYDRODYNAMIC_INTERACTION
    ! тензор диффузии (используется при наличие гидродинамических взаимодействий):
    real(8), dimension(:,:), allocatable :: D

    real(8), dimension(:,:), allocatable :: B

    integer(4) L_fixman
    real(8), dimension(:,:), allocatable :: fixman_matrix
    real(8), dimension(:), allocatable :: cheb_coeff

    real(8) eps_tea, beta_tea
    real(8), dimension(:), allocatable :: C_tea

    integer(4) m_lanczos
#endif    
    !отклонение от фиксированной конфигурации
    real(8) d_max
    ! подсчет пересчета двойных списков:
    integer(8) count_unbond_refresh
    !множитель случайной силы:
    real(8)   lambda_contr
    real(8), dimension(:), allocatable :: lambda
    
    ! степени свободы:
    integer(4) nfree
    
    ! энергия
    real(8) E_v, E_lj, E_ang, sumE
    
    ! действие случайных сил:
    real(8), dimension(:), allocatable :: grand_x, grand_y, grand_z
    
    ! флаги:
    !logical(1)  brefresh ! пересчет списка несвязанных атомов
    logical(1)	file_exists  
!    logical(1)  bcalc_diffusion_tensor ! пересчет тензор
     
    ! для подсчета времени:
    real(8) time_0, time_start, time_finish, time_cal    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 0. получаем начальную информацию:
    call ScanInputFiles(num_atom, num_bond, num_angles, &
                        & num_group, max_group_elem, num_taboo, num_surface, num_potential)
    
    ! 1. считываем координаты атомов
    call ScanCoordType(num_atom, num_type, type_name)
    call ReadCoord(num_atom, x, y, z, num_type, type_name, atom_type)

    ! 1+++. Считываем поверхности 
    call ReadSurface(num_surface, surf)
    ! 1+++++++. Считываем потенциалы
    call ReadPotential(num_potential, num_type, type_name, num_surface, surf, forces)
    ! Считываем: friction
    call ReadFriction(num_atom, atom_type, num_type, type_name, friction)
    
    ! 2. считываем связи:
    call ReadBonds(num_bond, bond1, bond2)
    ! проводим тест на корректность связей:
    call TestBonds(num_atom, num_bond, bond1, bond2)
    
    ! 3. считываем углы:
    call ReadAngl(num_angles, angles_l, angles_c, angles_r)
    ! тестируем правильность углов:
    call TestAngles(num_atom, num_angles, angles_l, angles_c, angles_r)
    
    
    ! 4. считываем группы частиц с фиксированным центром масс    
    call ReadGroup(num_atom, num_group, max_group_elem, natom_group, &
                    & group_x, group_y, group_z, &
                    &  num_atom_group, atom_group, bgroup_atom)
                    
    call ReadMoveCONTR(num_group, group_x, group_y, group_z, &
                    &   num_motion, motion_group1, motion_group2, &
                    &   motion_start_x1, motion_start_y1, motion_start_z1,  &
                    &   motion_start_x2, motion_start_y2, motion_start_z2, &
                    &   motion_dr1, motion_dr2, motion_min_dr, motion_max_dr, &
                    &   motion_dx, motion_dy, motion_dz)
    
    ! 5. считываем запрещённые взаимодействия
    call ReadTaboo(num_taboo, taboo1, taboo2)
    
    ! 4.---. Задать потенциалы по связям, углам и матрицы невалентных взаимодействий
    call OneParticleInteractionList(num_atom, atom_type, num_potential, forces, num_onepart, onepart, onepart_p)
    call ParticleParticleBondInteractionList(num_bond, bond1, bond2, num_atom, atom_type, &
                                        & num_potential, forces, bond_p)
    call ParticleParticleAnglInteractionList(num_angles, angles_r, angles_c, angles_l, num_atom, atom_type, &
                                & num_potential, forces, angles_p)
    call ParticleParticleUnBondInteractionMatrix(num_atom, atom_type, &
                                & num_potential, forces, num_type, P_aa_matrix)
         
    !TODO:
    !call SurfaceParticleBondInteractionList()
    call SurfaceParticleUnBondInteractionMatrix(num_atom, atom_type, &
                                & num_potential, forces, num_surface, num_type, P_sa_matrix)
  
    ! =========================================
    ! 6. читаем параметры системы:
    call ReadCONTR(dt, step, kT, &
		      & delta_r_sphere, max_unbond_atoms, &
                      & box, &
                      & h_hi, hradius, recalc_diffusion_tensor, hi_method, n_anti_gyration, &
		      & print_track, print_stats, print_picture, print_tremor)
    call TestCONTRParam(dt, step, kT, &
		    & delta_r_sphere, max_unbond_atoms, &
		    & box, h_hi, hradius, &
		    & recalc_diffusion_tensor, hi_method, n_anti_gyration,&
		    & print_track, print_stats, print_picture, print_tremor)

#ifdef BOUNDARY_CONDITION_CUBIC    
    ! тестирование размерности коробки:
    !call TestBox(num_atom, num_bond, bond1, bond2, box, rcut + delta_r_sphere, lbond)
#endif    
    
    ! 7. рассчитываем nfree и lambda
    ! считаем число степеней свободы
    nfree = 3*num_atom - num_bond - DIMENSION_SYSTEM 
		
    !TODO:
    friction_contr =  friction(1)
    ! получаем множитель случайной силы:
    lambda_contr = sqrt(2.0*kT*dt/friction_contr)
    allocate(lambda(1:num_atom))
    lambda(:) = sqrt(2.0*kT*dt/friction(:))
    
    ! 8. TRACK, STATS, TREMOR
    ! проверяем наличие файлов TRACK, STATS, TREMOR
    ! если нет, то создаём. если есть файл TRACK, то считать последний снимок
    call CheckTRACK(num_atom, x, y, z, nstep, &
                    & box, alpha_box)
    call CheckSTATS()
    call CheckTREMOR(num_group, natom_group, step, print_tremor)
        
    ! после того как мы прочитали файл TRACK, то у нас есть начальные координаты системы:
    ! - либо прочитанные из COORD
    ! - либо последний снимок из TRACK
    ! теперь можно делать подготовительные действия:
    ! 1) для цилиндрических граничных условий
    ! 2) для гидродинамических
    ! 3) двойные списки

    ! ============================================
    !9. Метка "не иметь объёмных взаимодействий", см. InputData.f90
    call LabeledFantomPairs(num_atom, num_bond, bond1, bond2, &
                            & num_angles, angles_l, angles_c, angles_r, &
                            & num_taboo, taboo1, taboo2, &
                            & atom_type, num_type, type_name, P_aa_matrix, label_taboo)

    !9.++ Узнаем массив парных Rsphere:
    call GetRsphere(delta_r_sphere, num_potential, forces, num_type, P_aa_matrix, r_sphere2)
    
    
    !10. Двойные списки
    !todo:
    eps = 1.0
    ! если есть объемные взаимодействия 
    if ((eps.eq.0.0).eqv..false.) then
	! для систем с маленьким числом атомов:
	if (max_unbond_atoms.gt.num_atom - 1) max_unbond_atoms = num_atom - 1    
	! максимльное число несвязанных взаимодействий:
	max_unbond_num = max_unbond_atoms * num_atom/2
	allocate(unbond1(1:max_unbond_num))
	allocate(unbond2(1:max_unbond_num))

    	! находим список контактов для несвязанных атомов
	call UnBond_List(num_atom,x,y,z, max_unbond_num, &
			  & num_unbond, unbond1, unbond2, box, &
			  & atom_type, num_type, r_sphere2,label_taboo)

	! запомним начальную конфигурацию
	allocate(x_mem(1:num_atom))
	allocate(y_mem(1:num_atom))
	allocate(z_mem(1:num_atom))
	
	do i=1, num_atom
	  x_mem(i) = x(i)
	  y_mem(i) = y(i)
	  z_mem(i) = z(i)
	end do
    else
      max_unbond_num = 0
      num_unbond = 0
    end if 
    

    !TODO:
    num_unbond_sa = 0
    allocate(unbond_sa1(1:num_unbond_sa))
    allocate(unbond_sa2(1:num_unbond_sa))
    
    !todo:
    num(:) = 0.0
    
!    ============================================== 
     !11. Выделяем память
#ifdef HYDRODYNAMIC_INTERACTION
    call HI_AllocateMemory(num_atom, D, hi_method, B)
#endif  

    allocate(Fx(1:num_atom))
    allocate(Fy(1:num_atom))
    allocate(Fz(1:num_atom))
    Fx(:) = 0.0
    Fy(:) = 0.0
    Fz(:) = 0.0
   
    allocate(grand_x(1:num_atom))
    allocate(grand_y(1:num_atom))
    allocate(grand_z(1:num_atom))

    allocate(delta_coord_x(1:num_atom))
    allocate(delta_coord_y(1:num_atom))
    allocate(delta_coord_z(1:num_atom))


    !12. Подготовка к HI
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! если есть гидродинамические взаимодействия:
    ! подготовка к расчёту случайной силы:
#ifdef HYDRODYNAMIC_INTERACTION
    ! считаем тензор диффузии:
    call RotnePragerYamakawa(num_atom, X, Y, Z, kT, friction_contr, hradius, D)
    ! в зависимости от метода:
    select case(hi_method)
        case(HI_METHOD_CHOLESKY)
		call CholeskyDecomposition(3*num_atom, D, B)
        case(HI_METHOD_CHEBYSHEV)
        	!call ChebyshevKrugerPrepare(num_atom, D, L_fixman, fixman_matrix, cheb_coeff)
	case(HI_METHOD_FIXMAN)
        case(HI_METHOD_TEA)
        	!call HI_TEA_Prepare(num_atom, D, eps_tea, beta_tea, C_tea)
	case(HI_METHOD_KRYLOV)
        	m_lanczos = int(3*num_atom/2)
    end select
#endif  
    
    ! 13. строим молекулу
    call PictureTypeAtom(num_atom, X, Y, Z, num_bond, bond1, bond2, &
              & num_type, type_name, atom_type, nstep)

    !^^^^^^^^^^^^^^^^ место для основного цикла программы

    ! счетчик обновления списка несвязанных атомов
    count_unbond_refresh = 0
    ! определяем время перед началом запуска алгоритма
    ! дважды запускаем секундомер
    call timer(time_0)
    call timer(time_start)

    do while (nstep.lt.step)
	  nstep = nstep + 1
	  ! рассчитываем случайные силы:
	  call GetRandomForce(num_atom, grand_x, grand_y, grand_z)
	  
	  ! рассчитываем силы действующие на частицы:
!           call Calculation_Force(num_atom, X, Y, Z, &
!                                 & num_bond, bond1, bond2, &
!                                 & num_unbond, unbond1, unbond2, &
!                                 & num_angles, angles_r, angles_c, angles_l,&
!                                 & kbond, lbond, eps, sig, rcut,solvent, &
!                                 & box, &
!                                 & k_ang, av_ang, &
!                                 & Fx, Fy, Fz)
        
        !TODO:
        !call CalculationPhi_RightNow(num_atom, X, Y, Z, num_layers, num, num_norm, nstep)
        !
        call Calculation_Force(num_atom, X, Y, Z, atom_type, num_type, &
                                & num_onepart, onepart, onepart_p, &
                                & num_bond, bond1, bond2, bond_p, &
                                & num_unbond, unbond1, unbond2, P_aa_matrix, &
                                & num_angles, angles_r, angles_c, angles_l, angles_p, &
                                & num_unbond_sa, unbond_sa1, unbond_sa2, P_sa_matrix, &
                                & num_potential, forces, &
                                & num_surface, surf, &
                                & box, &
                                & Fx, Fy, Fz)

	  ! если нужно учитывать гидродинамические взаимодействия:
#ifdef HYDRODYNAMIC_INTERACTION
	select case(hi_method)
        	case(HI_METHOD_CHOLESKY)
			call HI_Displacement_Cholesky(nstep, recalc_diffusion_tensor, num_atom, X, Y, Z, &
    				& grand_x, grand_y, grand_z, Fx, Fy, Fz, dt, kT, friction_contr, hradius, D, &   
    				& B, delta_coord_x, delta_coord_y, delta_coord_z)
        	case(HI_METHOD_CHEBYSHEV)
        		call HI_Displacement_Chebyshev(nstep, recalc_diffusion_tensor, num_atom, X, Y, Z, &
    				& grand_x, grand_y, grand_z, Fx, Fy, Fz, dt, kT, friction_contr, hradius, D, &   
    				& delta_coord_x, delta_coord_y, delta_coord_z)
		case(HI_METHOD_FIXMAN)
        		call HI_Displacement_Fixman(nstep, recalc_diffusion_tensor, num_atom, X, Y, Z, &
    				& grand_x, grand_y, grand_z, Fx, Fy, Fz, dt, kT, friction_contr, hradius, D, &   
    				& delta_coord_x, delta_coord_y, delta_coord_z)
        	case(HI_METHOD_TEA)
        		call HI_Displacement_TEA(nstep, recalc_diffusion_tensor, num_atom, X, Y, Z, &
    				& grand_x, grand_y, grand_z, Fx, Fy, Fz, dt, kT, friction_contr, hradius, D, &   
    				& delta_coord_x, delta_coord_y, delta_coord_z)
		case(HI_METHOD_KRYLOV)
        		call HI_Displacement_Krylov(nstep, recalc_diffusion_tensor, num_atom, X, Y, Z, &
    				& grand_x, grand_y, grand_z, Fx, Fy, Fz, dt, kT, friction_contr, hradius, D, m_lanczos, &   
    				& delta_coord_x, delta_coord_y, delta_coord_z)
		case(HI_METHOD_RPY_TEST)
			call HI_Displacement_Test(nstep, recalc_diffusion_tensor, num_atom, X, Y, Z, &
    					& grand_x, grand_y, grand_z, Fx, Fy, Fz, dt, kT, friction_contr, hradius, lambda_contr, &   
    					& delta_coord_x, delta_coord_y, delta_coord_z)
    	end select
#else 
! если нет гидродинамических взаимодействий:
	  do i=1, num_atom
    	    ! находим перемещение
 	    delta_coord_x(i) = Fx(i)*dt/friction(i) + lambda(i)*grand_x(i)
  	    delta_coord_y(i) = Fy(i)*dt/friction(i) + lambda(i)*grand_y(i)
   	    delta_coord_z(i) = Fz(i)*dt/friction(i) + lambda(i)*grand_z(i)
 	  enddo 
 	  
#endif 	  

         ! из итоговых приращений убираем вращение
#ifdef ANTI_GYRATION
         if (n_anti_gyration.gt.0) then
            call AntiRotation(num_atom, x, y, z, delta_coord_x, delta_coord_y, delta_coord_z, &
                             & x(n_anti_gyration), y(n_anti_gyration), z(n_anti_gyration))
         else
            ! относительно центра масс
            call CenterMass(num_atom, x, y, z, center_x, center_y, center_z)
            call AntiRotation(num_atom, x, y, z, delta_coord_x, delta_coord_y, delta_coord_z, &
                             & center_x, center_y, center_z)
         end if
#endif
    
 	  do i=1, num_atom
 	    ! и наконец, перемещаем частичку!
 	   
            x(i) = x(i) + delta_coord_x(i)
 	    y(i) = y(i) + delta_coord_y(i)
 	    z(i) = z(i) + delta_coord_z(i)
 	  ! кубические граничнык условия:
#ifdef BOUNDARY_CONDITION_CUBIC 	  
	    call BoundaryConditionCubic(x(i), y(i), z(i), box)
#endif
	  end do
 	
        
 	 !! >> В случае если есть группы атомов с фиксированным центром масс, ТО
 	 ! когда закончились перемещения, мы возвращаем группу атомов к 
 	 ! заданному центру масс 
 	 call DisplacementGroupToCenterMass(num_group, num_atom_group, max_group_elem, &
                                    & atom_group, group_x, group_y, group_z, num_atom, x, y, z)
 	 
        call MoveGroupCenterMass(num_motion, motion_group1, motion_group2, &
                            & motion_start_x1, motion_start_y1, motion_start_z1, &
                            & motion_start_x2, motion_start_y2, motion_start_z2, &
                            & motion_dr1, motion_dr2, motion_min_dr, motion_max_dr, &
                            & motion_dx, motion_dy, motion_dz, &
                            & num_group, group_x, group_y, group_z)
    
        !=========== ПЕРЕСЧЕТ ДВОЙНЫХ СПИСКОВ  ======>
          ! если есть объемные взаимодействия 
	  if ((eps.eq.0.0).eqv..false.) then	    
	      i = 1
	      ! цикл по всем атомам:
	      do while(i.le.num_atom)
		! смотрим отклонение от фиксированной конфигурации
		d_max = sqrt((x_mem(i) - x(i) )**2 + &
			    &(y_mem(i) - y(i) )**2 + &
			    &(z_mem(i) - z(i) )**2)
	
		if (d_max.ge.delta_r_sphere*0.5) then
                  count_unbond_refresh = count_unbond_refresh + 1
		  ! находим список контактов для несвязанных атомов

                    ! находим список контактов для несвязанных атомов
                    call UnBond_List(num_atom,x,y,z, max_unbond_num, &
			  & num_unbond, unbond1, unbond2, box, &
			  & atom_type, num_type, r_sphere2,label_taboo)
			  
                    ! запоминаем текущую конфигурацию
                    x_mem(:) = x(:)
                    y_mem(:) = y(:)
                    z_mem(:) = z(:)
                
                    ! на этом шаге больше не надо обновлять
                    i = num_atom + 1
		end if
		i = i + 1
            end do ! while
	    
          end if
   
          !=======================================================
          ! каждые print_track печатаем траекторию
            if (modulo(nstep,print_track).eq.0) then
                open(n_track,file='TRACK',ACCESS='append')
                write(n_track,'(a,I0)') 'step ',nstep
                do i=1,num_atom
                    write(n_track,*) i,x(i),y(i),z(i)
                enddo
                close(n_track)                 
            endif
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            ! каждые print_track печатаем снимки симпатичных молекулок
            if (modulo(nstep, print_picture).eq.0) then
                call PictureTypeAtom(num_atom, X, Y, Z, num_bond, bond1, bond2, &
                        &  num_type, type_name, atom_type, nstep)
            endif
              
            !TODO:
             ! выводим энергию, температуру
!             if (modulo(nstep, print_stats).eq.0) then
!                 E_v = 0.0
! ! 		if (num_bond.gt.0) then
! ! 		  ! считаем энергию:
! ! 		  call Calculation_Harmonic_Energy(num_atom, X, Y, Z, &
! !                                 & num_bond, bond1, bond2, &
! !                                 & kbond, lbond, box, alpha_box, E_v)
! ! 		  ! усредняем:
! ! 		  E_v = E_v / num_bond
! !                 end if
! 		
! 		 E_lj = 0.0
! ! 		call Calculation_Lennard_Jones_Energy(num_atom, X, Y, Z, &
! !                                     & num_unbond, unbond1, unbond2, &
! !                                     & eps, sig, rcut,solvent, box, alpha_box, E_lj)
! !                 ! усредняем
! !                 E_lj = E_lj / num_atom
! 		
! 		E_ang = 0.0
! ! 		if (num_angles.gt.0) then
! ! 		  call Calculation_Angle_Energy(num_atom, X, Y, Z, &
! !                                 & num_angles, angles_r, angles_c, angles_l,&
! !                                 & k_ang, av_ang, box, alpha_box, E_ang)
! !                 
! ! 		  E_ang = E_ang / num_angles
! !                 end if
! 		! потенциальная энергия:
! 		sumE = E_v + E_lj + E_ang		
! 		
! 		open(n_stats,file='STATS',ACCESS='append')
! 		write(n_stats,'(I10,F20.10,F20.10,F20.10,F20.10, F20.10, F20.10)') nstep, E_v, E_lj, E_ang, sumE
! 		
! 		close(n_stats)
!             end if
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            ! каждые print_track печатаем траекторию
            if (modulo(nstep,print_tremor).eq.0) then
                open(n_tremor,file='TREMOR',ACCESS='append')
                write(n_tremor,'(a,I0)') 'step ',nstep
                do i=1,num_atom
		    if (bgroup_atom(i).eqv..true.) then
		      write(n_tremor,*) i, Fx(i),Fy(i),Fz(i)
		    end if
                enddo
                close(n_tremor)                 
            endif
        
    enddo
    !*****************************
    call timer(time_finish) !останавливаем секудомер
    time_cal=time_finish-2.0_8*time_start+time_0

    open(n_infor,file='INFOR',ACCESS='append')
    write(n_infor,'(A,F20.3)') 'время работы алгоритма:', time_cal
    write(n_infor, '(A,I10)') 'число обновлений двойных списков:',count_unbond_refresh
    write(n_infor,'(A,I19)') 'общее число перемещений:', nstep
    
    call PictureTypeAtom(num_atom, X, Y, Z, num_bond, bond1, bond2, &
              & num_type, type_name, atom_type, step+1)
    !::::::::::::::::::::::::::::::::::::::::::::::::::::
    write(n_infor,'(a,a,a)') 'Программа ', PROGRAM_NAME ,' выполнена успешно!'
    close(n_infor)

    ! создаем файл с координатами конечной конфигурации
    call PrintCoordFinal(num_atom, X, Y, Z, num_type, type_name, atom_type)
!     open(n_coordf, file='COORDF')
!     write(n_coordf,'(A, I0)') 'num_atom ', num_atom
! 
!     do i=1, num_atom
!         write(n_coordf,*) x(i), y(i), z(i), type_name(atom_type(i))
!     enddo
! 
!     close(n_coordf)

    stop
end program BrownianDynamic
!*******************************************************************


subroutine CalculationPhi_RightNow(num_atom, X, Y, Z, num_layers, num, num_norm, nstep)
    integer(4) num_atom
    real(8) x(num_atom), y(num_atom), z(num_atom)
    integer(4) num_layers
    real(8) num(num_layers)
    real(8) num_norm(num_layers)
    integer(8) nstep
    !
    real(8) dist
    integer(4) i, j
    real(8), parameter :: dr = 0.1
    
    
    do i = 1, num_atom
        dist = abs(x(i))
        
        do j = 1, num_layers
            if (dist.ge.(j-1)*dr.and.dist.lt.j*dr) then
                num(j) = num(j) + 1
            endif
        enddo
    enddo
    
    do j = 1, num_layers
       num_norm(j) = num(j)/dble(nstep)
    enddo
    
    
    return
end subroutine CalculationPhi_RightNow

subroutine MoveGroupCenterMass(num_motion, motion_group1, motion_group2, &
                            & motion_start_x1, motion_start_y1, motion_start_z1, &
                            & motion_start_x2, motion_start_y2, motion_start_z2, &
                            & motion_dr1, motion_dr2, motion_min_dr, motion_max_dr, &
                            & motion_dx, motion_dy, motion_dz, &
                            & num_group, group_x, group_y, group_z)
    implicit none
    ! INPUT
    integer(4) num_motion
    integer(4) motion_group1(num_motion)
    integer(4) motion_group2(num_motion)
    real(8) motion_start_x1(num_motion), motion_start_y1(num_motion), motion_start_z1(num_motion)
    real(8) motion_start_x2(num_motion), motion_start_y2(num_motion), motion_start_z2(num_motion)
    real(8) motion_dr1(num_motion), motion_dr2(num_motion)
    real(8) motion_min_dr(num_motion), motion_max_dr(num_motion)
    real(8) motion_dx(num_motion), motion_dy(num_motion), motion_dz(num_motion)
    integer(4) num_group
    !OUTPUT
    real(8) group_x(num_group), group_y(num_group), group_z(num_group)
    
    
    !local
      save     bback ! делаем переменную статической
      logical  bback
      data     bback  / .false. /
      
    real(8) dx,dy,dz,dr
!     real(8) dx1,dy1,dz1,dr1
!     real(8) dx2,dy2,dz2,dr2
    
    integer(4) i

 	 !! >> В случае если есть движение пар групп атомов с фиксированным центром масс, ТО
 	 do i=1, num_motion
            ! узнаем расстояние между группами
            dx = group_x(motion_group2(i)) - group_x(motion_group1(i))
            dy = group_y(motion_group2(i)) - group_y(motion_group1(i))
            dz = group_z(motion_group2(i)) - group_z(motion_group1(i))
            
            dr = sqrt( dx**2 + dy**2 + dz**2 )
            
            ! узнаем расстояние первой группы до своей начальной позиции
!             dx1 = group_x(motion_group1(i)) - motion_start_x1(i)
!             dy1 = group_y(motion_group1(i)) - motion_start_y1(i)
!             dz1 = group_z(motion_group1(i)) - motion_start_z1(i)
!             
!             dr1 = sqrt( dx1**2 + dy1**2 + dz1**2 )
!             
!             ! узнаем расстояние второй группы до начальной позиции первой группы
!             dx2 = group_x(motion_group2(i)) - motion_start_x1(i)
!             dy2 = group_y(motion_group2(i)) - motion_start_y1(i)
!             dz2 = group_z(motion_group2(i)) - motion_start_z1(i)
!             
!             dr2 = sqrt( dx2**2 + dy2**2 + dz2**2 )
            
            ! сравниваем - это определяет куда двигаться дальше:
            ! сблизились достаточно-пора назад:
            if (dr.le.motion_min_dr(i)) then
                bback = .true.
            else
                !вернулись-значит опять сближаем:
                if (dr.ge.motion_max_dr(i).and.bback.eqv..true.) then
                    bback = .false.
                endif
            endif
            if (bback.eqv..false.) then
                ! первая группа ближе к своему начальному положению 
                ! продолжаем сближение
                !>
                group_x(motion_group1(i)) = group_x(motion_group1(i)) + motion_dx(i)*motion_dr1(i)
                group_y(motion_group1(i)) = group_y(motion_group1(i)) + motion_dy(i)*motion_dr1(i)
                group_z(motion_group1(i)) = group_z(motion_group1(i)) + motion_dz(i)*motion_dr1(i)
                !>
                group_x(motion_group2(i)) = group_x(motion_group2(i)) - motion_dx(i)*motion_dr2(i)
                group_y(motion_group2(i)) = group_y(motion_group2(i)) - motion_dy(i)*motion_dr2(i)
                group_z(motion_group2(i)) = group_z(motion_group2(i)) - motion_dz(i)*motion_dr2(i)
            else
                ! наоборот
                ! отдаляемся
                !>
                group_x(motion_group1(i)) = group_x(motion_group1(i)) - motion_dx(i)*motion_dr1(i)
                group_y(motion_group1(i)) = group_y(motion_group1(i)) - motion_dy(i)*motion_dr1(i)
                group_z(motion_group1(i)) = group_z(motion_group1(i)) - motion_dz(i)*motion_dr1(i)
                !>
                group_x(motion_group2(i)) = group_x(motion_group2(i)) + motion_dx(i)*motion_dr2(i)
                group_y(motion_group2(i)) = group_y(motion_group2(i)) + motion_dy(i)*motion_dr2(i)
                group_z(motion_group2(i)) = group_z(motion_group2(i)) + motion_dz(i)*motion_dr2(i)
            endif
            
            
 	 enddo
 	
    return
end subroutine MoveGroupCenterMass

!! >> В случае если есть группы атомов с фиксированным центром масс, ТО
 	 ! когда закончились перемещения, мы возвращаем группу атомов к 
 	 ! заданному центру масс 
subroutine DisplacementGroupToCenterMass(num_group, num_atom_group, max_group_elem, &
                    & atom_group, group_x, group_y, group_z, num_atom, x, y, z)
    ! INPUT
    integer(4) num_group
    integer(4) num_atom_group(num_group)
    integer(8) max_group_elem
    integer(4) atom_group(1:num_group,1:max_group_elem)
    real(8) group_x(num_group), group_y(num_group), group_z(num_group)
    integer(4) num_atom
    ! OUTPUT
    real(8) x(num_atom), y(num_atom), z(num_atom)
    ! local
    integer(4) i, j
    real(8) center_x, center_y, center_z
    
 	 do i=1, num_group

            call CenterGroup(num_atom, x, y, z, num_atom_group(i), max_group_elem, &
                            &  atom_group(i,:), center_x, center_y, center_z)

            ! перетаскиваем группу в заданный центр масс
            do j=1, num_atom_group(i)
                x(atom_group(i,j)) = x(atom_group(i,j)) - center_x + group_x(i)
                y(atom_group(i,j)) = y(atom_group(i,j)) - center_y + group_y(i)
                z(atom_group(i,j)) = z(atom_group(i,j)) - center_z + group_z(i)
            end do
         
 	 end do
 	 
    return
end subroutine DisplacementGroupToCenterMass

!>	@brief	Находим центр масс
!!	@param	num_atom     [IN] число атомов	
!!	@param  x            [IN] 	
!!	@param	y            [IN] координаты атомов	
!!	@param	z            [IN]	
!!	@param	center_x     [OUT]	
!!	@param	center_y     [OUT] координаты центра масс	
!!	@param	center_z     [OUT]	
!!	@date	28.11.2015
subroutine CenterMass(num_atom, x, y, z, center_x, center_y, center_z)
    ! INPUT
    integer(4) num_atom
    real(8) x(num_atom), y(num_atom), z(num_atom)
    ! OUTPUT
    real(8) center_x, center_y, center_z
    ! work
    integer(4) i

    center_x = 0.0
    center_y = 0.0
    center_z = 0.0
    ! идем по атомам в группе
    do i=1, num_atom
        center_x = center_x + x(i)
        center_y = center_y + y(i)
        center_z = center_z + z(i)
    end do

    center_x = center_x / num_atom
    center_y = center_y / num_atom
    center_z = center_z / num_atom

    return
end subroutine CenterMass


!>	@brief	Находим центр масс
!!	@param	num_atom     [IN] число атомов	
!!	@param  x            [IN] 	
!!	@param	y            [IN] координаты атомов	
!!	@param	z            [IN]	
!!	@param	natom_group  [IN] число групп с заданным центром масс	
!!	@param	max_elem     [IN] максимальное число атомов в таких группах	
!!	@param	atom_group   [IN] списки атомов, состоящих в группах	
!!	@param	center_x     [OUT]	
!!	@param	center_y     [OUT] координаты центра масс	
!!	@param	center_z     [OUT]	
!!	@date	21.10.2015
subroutine CenterGroup(num_atom, x, y, z, natom_group, max_elem, &
                    &  atom_group, center_x, center_y, center_z)
    ! INPUT
    integer(4) num_atom
    real(8) x(num_atom), y(num_atom), z(num_atom)
    integer(4) natom_group ! число атомов в этой группе
    integer(8) max_elem ! максимальное число элементов среди всех групп
    integer(4) atom_group(max_elem) ! список
    ! OUTPUT
    real(8) center_x, center_y, center_z
    ! work
    integer(4) i

    center_x = 0.0
    center_y = 0.0
    center_z = 0.0
    ! идем по атомам в группе
    do i=1, natom_group
        center_x = center_x + x(atom_group(i))
        center_y = center_y + y(atom_group(i))
        center_z = center_z + z(atom_group(i))
    end do

    center_x = center_x / natom_group
    center_y = center_y / natom_group
    center_z = center_z / natom_group

    return
end subroutine CenterGroup


!********************************************************************
!> 	@brief Определение процессорного времени timer
!!	@param	time	[OUT]	настоящее время в секундах с точностью до миллисекунд
subroutine timer(time) 
    real(8) time ! в секундах с точностью до миллисекунд
    integer(4) ival(8)
    call date_and_time(values = ival)

    time = dble(ival(8))*0.001_8 + dble(ival(7)) + &
          & dble(ival(6))*60.0_8 + dble(ival(5))*3600.0_8 + dble(ival(3))*86400.0_8

    select case (ival(2))
        case(1,3,5,7,8,10,12)
            time=time+dble(ival(2))*2678400.0_8
        case(4,6,9,11)
            time=time+dble(ival(2))*2592000.0_8
        case(2)
            if (modulo(ival(1),4).eq.0) then
                time=time+dble(ival(2))*2505600.0_8
            else
                time=time+dble(ival(2))*2419200.0_8
            endif
    end select

    time = time+dble(ival(1)-1)*31536000.0_8+dble((ival(1)-1)/4)*86400.0_8

    return
end subroutine timer
