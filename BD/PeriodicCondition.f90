Module PeriodicCondition
  Use CommonModule
  Use TestErrors
  
  implicit none

  
Contains
#ifdef BOUNDARY_CONDITION_CUBIC
!******************************************************************************
!         КУБИЧЕСКИЕ ГРАНИЧНЫЕ УСЛОВИЯ      
!******************************************************************************

!>	@brief	Округление до ближайшего целого
!!	@detailed собственная реализация функции nint()
!!	@param	d	[IN]	число которое нужно округлить
!!	@return	ближайшее целое
!!	@date	07.01.2015
real(8) function bint(d)
    implicit none
    ! INPUT
    real(8)	d
    ! WORK
    real(8)	b

    b = int(d); 

    ! если нужно было округлить  в большую сторону:
    if ((10 * d - 10 * real(b)).ge.5) b = b + 1

    bint = real(b)
    return
end function bint

!
!>	@brief	проверка граничных условий
!
!!	@param	x		[IN/OUT}	координаты атомов
!!	@param	y		[IN/OUT}	координаты
!!	@param	z		[IN/OUT}	координаты
!!	@param	box		[IN}	размер коробки
!
subroutine BoundaryConditionCubic(x, y, z, box)
    ! входные
    real(8)	x,y,z
    real(8)	box
    ! функция
    real(8) :: bint
   
    ! если |x| вышел за границу, то x = x - box * sign(x)
    if (abs(x).gt.box/2.0)  x = x - sign(box, x)
    if (abs(y).gt.box/2.0)  y = y - sign(box, y)
    if (abs(z).gt.box/2.0)  z = z - sign(box, z)
!     x = x - box*bint(x/box)
!     y = y - box*bint(y/box)
!     z = z - box*bint(z/box)  
!       
    return
end subroutine BoundaryConditionCubic
#endif

#ifdef BOUNDARY_CONDITION_CYLINDRICAL
!******************************************************************************
!         ЦИЛИНДРИЧЕСКИЕ ГРАНИЧНЫЕ УСЛОВИЯ      
!******************************************************************************
subroutine CylindrToCartesian(num_atom, rho, alpha, x, y)
    ! на вход:
    integer(4) num_atom
    real(8) rho(num_atom)
    real(8) alpha(num_atom)
    ! на выход:
    real(8) x(num_atom), y(num_atom)
    ! рабочие:
    integer(4) i
    
    ! переводим координаты из цилиндрической системы в декартову:
    do i = 1, num_atom
        x(i) = rho(i)*cos(alpha(i))
        y(i) = rho(i)*sin(alpha(i))
    end do
    
    return

end subroutine CylindrToCartesian


subroutine CartesianToCylindr(num_atom, x, y, rho, alpha)
    ! на вход:
    integer(4) num_atom
    real(8) x(num_atom), y(num_atom)
    ! на выход:
    real(8) rho(num_atom)
    real(8) alpha(num_atom)
    ! рабочие:
    integer(4) i
    
    ! переводим координаты из декартовой в цилиндрическую:
    do i = 1, num_atom
        rho(i) = sqrt(x(i)**2 + y(i)**2)
        ! TODO: учитывать деление на ноль rho(i) <> 0
        alpha(i) = sign(acos(x(i)/rho(i)), y(i))
    end do
    
    return

end subroutine CartesianToCylindr


!******************************************************************************
      subroutine distance_cylindrical_box(alfa_box, x1, y1, z1, x2, y2, z2, dist)
      !******************************************************************************
      real(8) alfa_box, x1, y1, z1, x2, y2, z2
      real(8) rr, r1, r2, alfa
      real(8) dist

      r1=sqrt(x1*x1+y1*y1)
      r2=sqrt(x2*x2+y2*y2)

      alfa=acos((x1*x2+y1*y2)/r1/r2)
      if (alfa.gt.alfa_box/2.0) alfa=alfa_box-alfa

      rr=r1**2+r2**2-2.0*r1*r2*cos(alfa)

      dist=sqrt(rr+(z1-z2)**2)

      end subroutine distance_cylindrical_box

!
!>	@brief	проверка граничных условий
!
!!	@param	x		[IN/OUT}	координаты атомов
!!	@param	y		[IN/OUT}	координаты
!!	@param	alpha_box        [IN}	        периодический угол
!
subroutine BoundaryConditionCylindrical(x, y, alpha_box)
      ! входные
      real(8) alpha_box, x, y, z
      ! вспомогательные
      real(8) r, alfa     
      real(8) pi
      parameter( pi=3.1415926535897932)

      r=sqrt(x*x+y*y)

      alfa=sign(acos(x/r),y/r) !
      
      if (alfa.gt.(pi/2.0+alpha_box).or.alfa.lt.(pi/2.0-2*alpha_box)) then
        write(*,*) 'ERROR! periodic box'
        stop
      endif

      if (alfa.gt.pi/2.0) alfa=alfa-alpha_box
      if (alfa.lt.(pi/2.0-alpha_box)) alfa=alfa+alpha_box

      x=r*cos(alfa)
      y=r*sin(alfa)
      
    return
end subroutine BoundaryConditionCylindrical

!******************************************************************************
      subroutine angle_box(alpha_box, x, y, r, alpha)
      !******************************************************************************
      ! на вход:
      real(8) alpha_box, x, y
      ! на выход
      real(8) alpha
      real(8) r

      r=sqrt(x*x+y*y)
      alpha=sign(acos(x/r),y/r)

      return
      end subroutine angle_box

!
!>	@brief	проверка граничных условий
!
!!	@param	alpha		[IN/OUT}	угол
!!	@param	alpha_box        [IN}	        периодический угол
!
subroutine BoundaryConditionCylindr(alpha, alpha_box)
      ! входные
      real(8) alpha_box, alpha
      ! вспомогательные
      real(8) pi
      parameter( pi=3.1415926535897932)
      
      if (alpha.gt.(pi/2.0+alpha_box).or.alpha.lt.(pi/2.0-2*alpha_box)) then
        write(*,*) 'ERROR! periodic box'
        stop
      endif

      ! проверяем граничные условия:
      if (alpha.gt.pi/2.0) alpha=alpha-alpha_box
      if (alpha.lt.(pi/2.0-alpha_box)) alpha=alpha+alpha_box

    return
end subroutine BoundaryConditionCylindr
      
      
      !******************************************************************************
      !	@brief расчет расстояния между двумя частицами
      !        для расчёта силы действующей на первую (рассматриваемую) частицу
      !        со стороны второй (действующая)
      !
      ! @param alpha_box [IN]	граничные условия
      ! @param alpha1 	 [IN]	угол рассматриваемой частицы
      ! @param alpha2 	 [IN]	угол действующей частицы
      ! @param distance  [OUT] расстояние
      subroutine DistanceCylindr(alpha_box, rho1, alpha1, &
                                & rho2, alpha2, dx, dy)
      	 ! входные:
	 real(8) alpha_box
	 real(8) alpha1, alpha2
	 real(8) rho1, rho2
	 ! выходные:
	 real(8) dx, dy
	 ! рабочие:
	 real(8) alpha12, alpha2_s

	 ! угол между частицами:
	 alpha12 = alpha1 - alpha2
	 ! как дальше поступить?
	 if (abs(alpha12).le.alpha_box/2.0) then
                !write(*,*) rho1, rho2, alpha1, alpha2 
	 	! угол небольшой, значит ничего не меняем
		dx = rho2*cos(alpha2) - rho1*cos(alpha1)
		dy = rho2*sin(alpha2) - rho1*sin(alpha1)
		! выходим:
		return
	 end if
		
        !write(*,*) 'alpha = ', alpha12
	! пересчитываем положение второй частицы:
	alpha2_s = alpha2 + sign(alpha_box, alpha12)
	
	! пересчитываем расстояние:
	dx = rho2*cos(alpha2_s) - rho1*cos(alpha1)
        dy = rho2*sin(alpha2_s) - rho1*sin(alpha1)
	
	return
      end subroutine DistanceCylindr

      
      !******************************************************************************
      !	@brief расчет расстояния между двумя частицами
      !        для расчёта силы действующей на первую (рассматриваемую) частицу
      !        со стороны второй (действующая)
      !
      ! @param alpha_box [IN]	граничные условия
      ! @param alpha1 	 [IN]	угол рассматриваемой частицы
      ! @param alpha2 	 [IN]	угол действующей частицы
      ! @param x1 	 [IN]	координаты рассматриваемой
      ! @param y1 	 [IN]
      ! @param x2 	 [IN]	координаты действующей
      ! @param y2 	 [IN]
      ! @param r2 	 [IN]	радиусс-вектор действующей
      ! @param dx 	 [OUT]	проекция расстояния между частицами
      ! @param dy 	 [OUT]		с учётом граничных условий
      subroutine DistanceCorrect(alpha_box, alpha1, alpha2, &
			      & x1, y1, x2, y2, r2, dx, dy)
      	 ! входные:
	 real(8) alpha_box
	 real(8) alpha1, alpha2
	 real(8) x1, y1, x2, y2, r2
	 ! выходные:
	 real(8) dx, dy
	 ! рабочие:
	 real(8) alpha12, alpha_s, x2_s, y2_s

	 ! угол между частицами:
	 alpha12 = alpha1 - alpha2
	 ! как дальше поступить?
	 if (abs(alpha12).le.alpha_box/2.0) then
	 	! угол небольшой, значит ничего не меняем
		dx = x2 - x1
		dy = y2 - y1
		! выходим:
		return
	 end if
		
	! пересчитываем положение второй частицы:
	alpha_s = alpha2 + sign(alpha_box, alpha12)
	
	! радиус-вектор:
	!r2 = sqrt(x2**2 + y2**2)
	! вспомогательные координаты (с учётом граничных условий):
	x2_s = r2*cos(alpha_s)
	y2_s = r2*sin(alpha_s)

	! расстояние между частицами (с учётом граничных условий):
	dx = x2_s - x1
	dy = y2_s - y1

	return
      end subroutine DistanceCorrect
#endif


end Module PeriodicCondition