!> @brief       Поверхность, расстояние до поверхностей
!!
!!              модуль для прикладных программистов
!!              
!!              Для ввода новой поверхности 'X' необходимо:
!!              1. Код поверхности 'X' (SURF_CODE_X)
!!              2. SurfaceRecognition
!!              3. Функция расчета расстояния до поверхности 'X'
!!              4. DistanceToSurface
!!
!! @date		03.06.2018
!!
!! @author	Михайлов И.В., Шавыкин О.В.
Module Surface
    Use CommonModule
    Use TestErrors
    
!TODO: dx+dy+dz

    !> индетификаторы поверхностей:
    integer(4), parameter, public :: SURF_CODE_PLANE = 1, &
                                     SURF_CODE_CYLINDR = 2, &
                                     SURF_CODE_SPHERE = 3, &
                                     SURF_CODE_ELLIPSOID = 4, &
                                     SURF_CODE_TORUS = 5
                                     
    !character, parameter, public :: SURFACE_NAME = 'fdfdf'
    
    type surface_info 
        character(3) name
        integer(4) type_code
        integer(4) num_param
        real(8), allocatable, dimension(:) :: param
    end type

Contains 


!>
subroutine SurfaceRecognition(n_file, surf_name, surf)
    ! input 
    integer(4) n_file
    ! output
    character(255) surf_name
    type(surface_info) surf
    ! local
    character(255) temp_word
    logical(1) bfind
    integer(4) ioer
    
    bfind = .false.

    read(n_file,*) surf_name
    ! 
    if (surf_name(1:len_trim(surf_name)).eq.'plane') then
        surf%type_code = SURF_CODE_PLANE
       
        surf%num_param = 4
        allocate(surf%param(1:surf%num_param))
        surf%param(:) = 0.0
        ! нашли:
        bfind = .true.
    endif
    
    if (surf_name(1:len_trim(surf_name)).eq.'sphere') then
        surf%type_code = SURF_CODE_SPHERE
        
        surf%num_param = 4
        allocate(surf%param(1:surf%num_param))
        surf%param(:) = 0.0
        ! нашли:
        bfind = .true.
    endif
    
    if (surf_name(1:len_trim(surf_name)).eq.'ellipsoid') then
        surf%type_code = SURF_CODE_ELLIPSOID
       
        surf%num_param = 6
        allocate(surf%param(1:surf%num_param))
        surf%param(:) = 0.0
        ! нашли:
        bfind = .true.
    endif
    
    
    ! если нашли такую поверхность, то считаем её параметры:
    if (bfind.eqv..true.) then
        backspace(n_file)
        read(n_file, *, iostat=ioer) temp_word, surf%name, surf%param
    else
         call ERROR(ERR_SURFACE_UNKNOWN)
    endif

    return
end subroutine SurfaceRecognition

!>
subroutine DistanceToSurface(x, y, z, surf, dx, dy, dz)
    ! input 
    real(8) x, y, z
    type(surface_info) surf
    ! output
    real(8) dx, dy, dz

    select case(surf%type_code)
        case(SURF_CODE_PLANE)
            call SurfaceDistancePlane(x, y, z, surf%param(1), &
                            & surf%param(2), surf%param(3), surf%param(4), dx, dy, dz)
        case(SURF_CODE_SPHERE)
            call SurfaceDistanceSphere(x, y, z, surf%param(1), &
                            & surf%param(2), surf%param(3), surf%param(4), dx, dy, dz)
        case(SURF_CODE_ELLIPSOID)
            call SurfaceDistanceEllipsoid(x, y, z, surf%param(1), surf%param(2), surf%param(3), &
                            & surf%param(4), surf%param(5), surf%param(6), dx, dy, dz)
        case default
            !call 
    end select 

    return
end subroutine DistanceToSurface

!>>>>>>>>>>>>>>>>>>>>>>> Список функций рассчитывающих расстояние от точки до поверхностей
!>
subroutine SurfaceDistancePlane(x, y, z, A, B, C, D, dx, dy, dz)
    ! input
    real(8) x, y, z
    real(8) A, B, C, D
    ! output
    real(8) dx, dy, dz
    real(8) dist
    
    ! находим расстояние от точки до сферы:
    dist = abs(A*x + B*y + C*z + D)/sqrt(A**2 + B**2 + C**2)
    dx = 0
    dy = 0
    dz = 0

    return
end subroutine SurfaceDistancePlane

!>
subroutine SurfaceDistanceSphere(x, y, z, x0, y0, z0, r_sph, dx, dy, dz)
    ! input
    real(8) x, y, z
    real(8) x0, y0, z0
    real(8) r_sph
    ! output
    real(8) dx, dy, dz
    !local
    real(8) x_surf, y_surf, z_surf
    real(8) x_point, y_point, z_point
    real(8) t
    
    ! находим расстояние от точки до сферы:
   ! dist = abs(sqrt((x-x0)**2 - (y-y0)**2 - (z-z0)**2) - r_sph)
    
    x_point = x - x0
    y_point = y - y0
    z_point = z - z0
    
    ! находим расстояние от точки до сферы:
    t = (x_point/r_sph)**2 + (y_point/r_sph)**2 + (z_point/r_sph)**2 
    
    x_surf = x_point/t
    y_surf = y_point/t
    z_surf = z_point/t
    
    dx = x_surf - x_point
    dy = y_surf - y_point
    dz = z_surf - z_point

    return
end subroutine SurfaceDistanceSphere

!>
!http://rsdn.org/forum/alg/2773935.flat
subroutine SurfaceDistanceEllipsoid(x, y, z, x0, y0, z0, r_x, r_y, r_z, dx, dy, dz)
    ! input
    real(8) x, y, z
    real(8) x0, y0, z0
    real(8) r_x, r_y, r_z
    ! output
    real(8) dx, dy, dz
    ! local
    real(8) x_surf, y_surf, z_surf
    real(8) x_point, y_point, z_point
    real(8) t
    
    ! приводим координаты к системе координат с центром в (x0, y0, z0)
    x_point = x - x0
    y_point = y - y0
    z_point = z - z0
    
    ! находим расстояние от точки до сферы:
    t = (x_point/r_x)**2 + (y_point/r_y)**2 + (z_point/r_z)**2 
    
    x_surf = x_point/t
    y_surf = y_point/t
    z_surf = z_point/t
    
    dx = x_surf - x_point
    dy = y_surf - y_point
    dz = z_surf - z_point

    return
end subroutine SurfaceDistanceEllipsoid


!>>
subroutine SurfaceDistanceTorus(x, y, z, x0, y0, z0, r_major, r_minor, dx, dy, dz)
    ! input
    real(8) x, y, z
    real(8) x0, y0, z0
    real(8) r_major, r_minor
    ! output
    real(8) dx, dy, dz
    real(8) dist
    ! local
    real(8) t
    
    ! приводим координаты к системе координат с центром в (x0, y0, z0)
    x_point = x - x0
    y_point = y - y0
    z_point = z - z0
    
    dx = 0
    dy = 0 
    dz = 0
    dist = 0

    return
end subroutine SurfaceDistanceTorus

end module Surface