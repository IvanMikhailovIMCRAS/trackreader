!
!>	@brief	список общих параметров
!!
!!	@n содержит:
!!	@n 1. Дескрипторы файлов
!!	@n 2. Коды ошибок
!!

module CommonModule
    Implicit None
    character(20), parameter,	public :: PROGRAM_NAME = 'BrownianDynamic 2014'
    character(3), parameter,	public :: VERSION_STR = '2.0'
    
    integer, parameter, public :: LEN_TYPE_NAME = 2
    !
    integer, parameter, public :: DIMENSION_SYSTEM = 3
    
    ! дексрипторы входных файлов:
    integer, parameter,	public :: n_coord = 11, &
				  n_bonds = 12, &
				  n_contr = 13, &
				  n_fix = 17, &
				  n_angl = 21, &
				  n_taboo = 22, &
				  n_group = 23, &
				  n_surface = 24, &
				  n_force = 25
    ! дескрипторы выходных файлов:
    integer, parameter,	public :: n_infor = 34, &
				  n_track = 35, &
				  n_stats = 36, &
				  n_coordf = 38, &
				  n_tremor = 39
    
    
    ! коды ошибок:
    integer, parameter,	public :: ERR_CODE_DISTANCE = 1, &
				  ERR_CODE_LENGTH = 2, &
				  ERR_CODE_PERCENT = 3, &
				  ERR_CODE_OPEN_COORD = 4, &
				  ERR_CODE_OPEN_BONDS = 5, &
				  ERR_CODE_OPEN_CONTR = 6, &
				  ERR_CODE_UNCORRECT_BONDS = 7, &
				  ERR_CODE_IDENTIC_BONDS = 8, &
				  ERR_CODE_CICLE = 9, &
				  ERR_CODE_BIG_NUM = 10, &
				  ERR_CODE_SMALL_NUM = 11, &
				  ERR_CODE_UNCORRECT_ANGLES = 12, &
				  ERR_CODE_IDENTIC_ANGLES = 13, &
				  ERR_CODE_SMALL_BOX = 14, &
				  ERR_CODE_NONE_ATOMS = 15, &
				  ERR_CODE_IDENTIC_ATOMS = 16, &
				  ERR_CODE_BAD_CONTR = 17, &
				  ERR_CODE_BAD_PARAM = 18, &
				  ERR_CODE_CONTR_PARAM = 19, &
				  ERR_CODE_READ_COORD = 20, &
				  ERR_CODE_READ_SURFACE = 21, &
                                  ERR_CODE_SAME_SURFACE = 22, &
                                  ERR_SURFACE_UNKNOWN = 23, &
                                  ERR_CODE_READ_POTENTIAL = 24, &
                                  ERR_CODE_SAME_POTENTIAL = 25, &
                                  ERR_CODE_POTENTIAL_UNKNOWN = 26, &
                                  ERR_CODE_POTENTIAL_OPEN_FIELD = 27, &
                                  ERR_CODE_POTENTIAL_UNKNOWN_PART_TYPE = 28, &
                                  ERR_CODE_POTENTIAL_UNKNOWN_PART_TYPE2 = 29, &
                                  ERR_CODE_NO_BOND_POTENTIAL = 30, &
                                  ERR_CODE_NO_ANGL_POTENTIAL = 31, &
                                  ERR_CODE_FRICTION_UNKNOWN_PART_TYPE = 32, &
                                  ERR_CODE_POTENTIAL_NO_UNBOND = 33, &
                                  ERR_CODE_UNKNOWN_1PART_POTENTIAL = 34, &
                                  ERR_CODE_UNKNOWN_2PART_POTENTIAL = 35, &
                                  ERR_CODE_UNKNOWN_3PART_POTENTIAL = 36, &
                                  ERR_CODE_READ_SPEED = 37, &
                                  ERR_CODE_SPEED_WRONG_GROUP = 38, &
                                  ERR_CODE_SPEED_ABROAD = 39, &
                                  ERR_CODE_SPEED_MIN_MAX = 40, &
                                  ERR_CODE_UNKNOWN_SURFACE_PART_POTENTIAL = 41, &
                                  ERR_CODE_READ_GROUP = 42, &
                                  ERR_CODE_READ_GROUP_ELEM = 43

    ! методы учета гидродинамических взаимодействий:
    integer, parameter,	public :: HI_METHOD_CHOLESKY = 1, &
				& HI_METHOD_CHEBYSHEV = 2, &
				& HI_METHOD_FIXMAN = 3, &
				& HI_METHOD_TEA = 4, &
				& HI_METHOD_KRYLOV = 5 , &
				& HI_METHOD_RPY_TEST = 6 !,&
				!& HI_METHOD_NONE = 0
				  
    real(8), parameter, public :: BD_PI = 3.1415926535897932
Contains
    

end module CommonModule
