!>	@brief	Функции тестирования корректности работы и обработка кодов ошибок
Module TestErrors
    ! загружаем коды ошибок и дескрипторы файлов:
    Use CommonModule
    
    Implicit None
    
Contains
    
    !*************************************************************
    !>	@brief	Интерпретация кодов ошибок
    !!	@param	n	[IN]	код ошибки
    subroutine ERROR(n)
!*************************************************************
    integer(4) n
 
    open(n_infor, file='INFOR', access='append')
    write(n_infor,'(a)') 'ОШИБКА:'
    select case (n)
     case(ERR_CODE_DISTANCE)
 	write(n_infor,'(a)') 'Расстояние между несвязанными частицами меньше 0.5*sig'
 	write(n_infor,'(a)') 'выберите более маленький шаг dt '
    case(ERR_CODE_LENGTH)
	write(n_infor,'(a)') 'Длина связи отклоняется на больше чем 0.5*lbond'
	write(n_infor,'(a)') 'выберите более маленький шаг dt '
    case(ERR_CODE_PERCENT)
	write(n_infor,'(a)') 'Число принятых конфигураций слишком мало'
	write(n_infor,'(a)') 'выберите более маленький шаг dt '
    case(ERR_CODE_OPEN_COORD)
	write(n_infor,'(a)') 'отсутствует файл COORD'
    case(ERR_CODE_READ_COORD)
	write(n_infor,'(a)') 'проблема при чтении файла COORD'
    case(ERR_CODE_OPEN_BONDS)
	write(n_infor,'(a)') 'отсутствует файл BONDS'
    case(ERR_CODE_OPEN_CONTR)
	write(n_infor,'(a)') 'отсутствует файл CONTR'
    case(ERR_CODE_UNCORRECT_BONDS)
	write(n_infor,'(a)') 'содержит петлю (частица соединена сама с собой)'
    case(ERR_CODE_IDENTIC_BONDS)
	write(n_infor,'(a)') 'одиннаковые связи в BONDS'
    case(ERR_CODE_BIG_NUM)
	write(n_infor,'(a)') 'в BONDS номер частиц превосходит число частиц'
    case(ERR_CODE_SMALL_NUM)
 	write(n_infor,'(a)') 'в BONDS номер атома меньше 1'
    case (ERR_CODE_IDENTIC_ANGLES)
	write(n_infor,'(a)') 'в ANGL одиннаковые углы'
    case (ERR_CODE_SMALL_BOX)
	write(n_infor,'(a)') 'размер коробки слишком маленький'
    case (ERR_CODE_NONE_ATOMS)
	write(n_infor,'(a)') 'в COORD нет частиц'
    case (ERR_CODE_CICLE)
	write(n_infor,'(a)') 'в BONDS связи образует граф с циклом'
    case (ERR_CODE_BAD_CONTR)
	write(n_infor,'(a)') 'не хватает параметров в CONTR'
    case (ERR_CODE_BAD_PARAM)
      write(n_infor,'(a)') 'ошибка чтения параметра в CONTR'
    case (ERR_CODE_CONTR_PARAM)
      write(n_infor,'(a)') 'недопустимое значение параметра в CONTR'    
    case (ERR_CODE_READ_SURFACE)
      write(n_infor,'(a)') 'проблема при чтении файла SURFACE'
    case (ERR_CODE_SAME_SURFACE)
      write(n_infor,'(a)') 'одно и тоже имя у разных поверхностей (SURFACE)'
    case (ERR_SURFACE_UNKNOWN)
      write(n_infor,'(a)') 'неизвестная поверхность (SURFACE)'
    case (ERR_CODE_READ_POTENTIAL)
      write(n_infor,'(a)') 'проблема при чтении файла POTENTIAL' 
    case (ERR_CODE_POTENTIAL_OPEN_FIELD)
      write(n_infor,'(a)') 'отсутствует файл с внешним полем' 
    case (ERR_CODE_POTENTIAL_UNKNOWN)
      write(n_infor,'(a)') 'неизвестный тип потенциала' 
    case (ERR_CODE_POTENTIAL_UNKNOWN_PART_TYPE)
      write(n_infor,'(a)') 'неизвестный тип частицы для потенциала'
    case (ERR_CODE_POTENTIAL_UNKNOWN_PART_TYPE2)
      write(n_infor,'(a)') 'неизвестный тип частицы (или поверхности) &
                            & для потенциала'
    case (ERR_CODE_NO_BOND_POTENTIAL)
      write(n_infor,'(a)') 'не задан по крайней мере один потенциал связи'
    case (ERR_CODE_NO_ANGL_POTENTIAL)
      write(n_infor,'(a)') 'не задан по крайней мере один потенциал угла'
    case (ERR_CODE_FRICTION_UNKNOWN_PART_TYPE)
      write(n_infor,'(a)') 'неизвестный тип частицы для коэффициента трения'
    case (ERR_CODE_POTENTIAL_NO_UNBOND)
      write(n_infor,'(a)') 'попытка использовать неверный потенциал &
                            & для невалентных взаимодейтсвий'
    case (ERR_CODE_UNKNOWN_1PART_POTENTIAL)
      write(n_infor,'(a)') 'неизвестный одночастичный потенциал при расчете сил'
    case (ERR_CODE_UNKNOWN_2PART_POTENTIAL)
      write(n_infor,'(a)') 'неизвестный двухчастичный потенциал при расчете сил'
    case (ERR_CODE_UNKNOWN_3PART_POTENTIAL)
      write(n_infor,'(a)') 'неизвестный трехчастичный потенциал при расчете сил'
    case (ERR_CODE_READ_SPEED)
      write(n_infor,'(a)') 'ошибка при чтении смещений пар групп частиц'
    case (ERR_CODE_SPEED_WRONG_GROUP)
      write(n_infor,'(a)') 'некорректные номера в парах групп частиц'
    case (ERR_CODE_SPEED_ABROAD)
      write(n_infor,'(a)') 'недопустимое начальное расстояние между &
                            & центрами масс групп'
    case (ERR_CODE_SPEED_MIN_MAX)
      write(n_infor,'(a)') 'минимальное расстояние не может превышать &
                            & максимальное расстояние между центрами масс групп'
    case (ERR_CODE_UNKNOWN_SURFACE_PART_POTENTIAL)
      write(n_infor,'(a)') 'неизвестный потенциал поверхность-частица &
                            & при расчете сил'
    case (ERR_CODE_READ_GROUP)
      write(n_infor,'(a)') 'ошибка при чтении файла GROUP'   
    case (ERR_CODE_READ_GROUP_ELEM)
      write(n_infor,'(a)') 'ошибка при чтении элемента одной из групп &
                            & в файле GROUP' 
!     case(11)
!	write(n_infor,'(a)') 'не все'
    end select

    close(n_infor)
    stop
    return
    end subroutine ERROR
!*************************************************************

 subroutine TestCONTRParam(dt, step, kT, &
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
    
    ! /шаг моделирования/
    if (dt.le.0.0) then  
        call ERROR(ERR_CODE_CONTR_PARAM)
    endif
    ! /число шагов моделирования/
    if (step.lt.0) then 
        call ERROR(ERR_CODE_CONTR_PARAM)
    endif
    ! /энергия в дисперсии случайной силы/ 
    if (kT.lt.0) then 
        call ERROR(ERR_CODE_CONTR_PARAM)
    endif
    ! /параметр метода двойных списков/
    if (delta_r_sphere.lt.0) then 
        call ERROR(ERR_CODE_CONTR_PARAM)
    endif
    ! /периодичность вывода координат частиц/
    if (print_track.lt.0) then 
        call ERROR(ERR_CODE_CONTR_PARAM)
    endif
    ! /периодичность вывода картинок/
    if (print_picture.lt.0) then 
        call ERROR(ERR_CODE_CONTR_PARAM)
    endif 
    ! /периодичность вывода энергии системы/
    if (print_stats.lt.0) then 
        call ERROR(ERR_CODE_CONTR_PARAM)
    endif 
    if (print_tremor.lt.0) then 
        call ERROR(ERR_CODE_CONTR_PARAM)
    endif 
    !bondary_condition = 0.0   !  /размер коробки/
    !h = 0.0 ! /гидродинамический радиус/
    !hi_tensor_period =  1 ! /период перерасчета тензора инерции (только для метода cholesky)/
    !hi_method  = HI_METHOD_CHOLESKY ! /метод учета гидродинамических взаимодействий |cholesky | chebyshev | tea | krylov |/
    !n_anti_gyr = 0 ! /относительно какого атома будет происходить удаление вращения/
    
    return
end subroutine TestCONTRParam

      !>	@brief проверяем
      !! 	@n	1) не повторяются ли одинаковые связи 
      !! 	@n	2) есть ли петли
      !!
      !! 	@n	3) чтобы не было атомов с номером больше максимального, 
      !! 	@n	4) меньше 1 
      !!
      !!	@param	num_atom	[IN]	число атомов
      !!	@param	num_bond	[IN]	число связей
      !!	@param	bond1		[IN]	
      !!	@param	bond2		[IN]	массивы связей
      subroutine TestBonds(num_atom, num_bond, bond1, bond2)
	  ! входные
	  integer(4)	num_atom
	  integer(4)	num_bond
	  integer(4)	bond1(num_bond), bond2(num_bond)
	  ! для работы
	  integer(4)	i, j
	  
	  do i=1,num_bond
		! проверяем чтобы не было связей между атомами с одинаковыми номерами (типа 1 - 1)
		if (bond1(i).eq.bond2(i)) call ERROR(ERR_CODE_UNCORRECT_BONDS)
		
		do j=i+1,num_bond	
			if (bond1(i).eq.bond1(j).and.bond2(i).eq.bond2(j)) call ERROR(ERR_CODE_IDENTIC_BONDS)
			
			if (bond1(i).eq.bond2(j).and.bond2(i).eq.bond1(j)) call ERROR(ERR_CODE_IDENTIC_BONDS)
		enddo
	  enddo
	 
	  
	  do i=1,num_bond
 		! если есть меньше num_atom
 		if (bond1(i).gt.num_atom.or.bond2(i).gt.num_atom) call ERROR(ERR_CODE_BIG_NUM)
 		! если есть меньше 1
 		if (bond1(i).lt.1.or.bond2(i).lt.1) call ERROR(ERR_CODE_SMALL_NUM)
	  enddo

	 return
	  
      end subroutine TestBonds 
    
    
      !>	@brief проверяем
      !!	@n 1) не повторяются ли одинаковые углы
      !!	@n
      !!	@n 2) чтобы не было атомов с номером больше максимального, 
      !!	@n 3) меньше 1 
      !!
      !!	@param	num_atom	[IN]	число атомов
      !!	@param	num_bond	[IN]	число углов
      !!	@param	angles_l	[IN]	
      !!	@param	angles_c	[IN]	массивы углов
      !!	@param	angles_r	[IN]	
      subroutine TestAngles(num_atom, num_angles, angles_l, angles_c, angles_r)
	  integer(4)	num_atom
	  integer(4)	num_angles
	  integer(4)	angles_l(num_angles), angles_c(num_angles), angles_r(num_angles)
	  ! для работы
	  integer(4)	i, j
	  
	  do i=1,num_angles
		! проверяем чтобы не было одиннаковых вершин в угле:
		if (angles_l(i).eq.angles_r(i)) call ERROR(ERR_CODE_UNCORRECT_ANGLES)
		if (angles_c(i).eq.angles_r(i)) call ERROR(ERR_CODE_UNCORRECT_ANGLES)
		if (angles_c(i).eq.angles_l(i)) call ERROR(ERR_CODE_UNCORRECT_ANGLES)
		
		do j=i+1,num_angles	
			! если совпадает центр угла, то проверим краевые:
			if (angles_c(i).eq.angles_c(j)) then
			    if (angles_l(i).eq.angles_l(j).and.angles_r(i).eq.angles_r(j)) call ERROR(ERR_CODE_IDENTIC_ANGLES)
			
			    if (angles_l(i).eq.angles_r(j).and.angles_l(j).eq.angles_r(i)) call ERROR(ERR_CODE_IDENTIC_ANGLES)
			end if
		enddo
	  enddo
	 
      end subroutine TestAngles
      
      !
      !>	@brief	проверяем размер коробки
      !!
      !!	@param	num_atom	[IN]	число атомов
      !!	@param	num_bond	[IN]	число связей
      !!	@param	bond1		[IN]	
      !!	@param	bond2		[IN]	массивы связей
      !!	@param	box		[IN]	размер коробки
      !!	@param	rsphere		[IN]	радиус невалетного взаимодейтсвия
      !!	@param	lbond		[IN]	равновесная длина связи
      !
      subroutine TestBox(num_atom, num_bond, bond1, bond2, box, rsphere, lbond)
	  ! входные
	  integer(4)	num_atom, num_bond
	  integer(4)	bond1(num_atom), bond2(num_atom)
	  real(8)	box, rsphere, lbond
	  ! локальные
	  integer(4) list(num_atom)
	  integer(4) max_len
	  real(8) min_box
	  integer(4) i
	  
	  ! заполняем список степеней вершин атомов
	  list(:)=0
	  do i=1,num_bond
		list(bond1(i))=list(bond1(i))+1
		list(bond2(i))=list(bond2(i))+1
	  enddo

	  ! получаем максимальное расстояние между концами в абсолютных единицах:
	  call GetLongestWay(num_bond,bond1, bond2, num_atom, list, max_len)
	  ! находим минимальный размер коробки:
	  min_box = max_len*lbond + rsphere
	  
	  !write(*,*) min_box, box
	  
	  if (box.le.min_box) call ERROR(ERR_CODE_SMALL_BOX)
	  
	  return
      end subroutine

!!
!>	@brief	Ищем вес вершины
!!	@param	num_atom	[IN]	число атомов
!!	@paran	num_bond	[IN]	число связей
!!	@param	bond1		[IN]
!!	@param	bond2		[IN]	массивы связей
!!	@param	weight		[IN/OUT] массив весов 
!!	@param	curr_atom	{IN]	номер анализируемой вершины
!!
!!
recursive subroutine GetWeight(num_atom,num_bond, bond1, bond2, weight, curr_atom)
  ! у нас есть атом, массив пути
  ! длина пути (curr_len)
  integer(4) num_atom,num_bond
  integer(4) bond1(num_bond), bond2(num_bond), weight(num_atom)
  integer(4) curr_atom
  integer(4)	i, j
  
  integer(4)	min_w
  logical bfind
  
  ! определяем вес вершины:
  min_w = num_atom
  ! ищем те, в которых были:
  do i=1, num_bond
	! если нашли вес у соседней вершины
	if (curr_atom.eq.bond1(i)) then
	    if (weight(bond2(i)).lt.min_w) then
	       min_w = weight(bond2(i)) + 1
	    end if
	    
	end if
	! если нашли вес у соседней вершины
	if (curr_atom.eq.bond2(i)) then
	    if (weight(bond1(i)).lt.min_w) then
	       min_w = weight(bond1(i)) + 1
	    end if
	end if
  end do
  
  ! нашли минимальный:
  weight(curr_atom) = min_w
  
  ! смотрим куда пойти дальше:
  do i=1, num_bond
	
	if (curr_atom.eq.bond1(i)) then
	    ! если мы не обходили вершину, то посмотрим еЁ
	    if (weight(bond2(i)).ge.num_atom) then
		call GetWeight(num_atom,num_bond, bond1, bond2, weight, bond2(i))
	    end if
	end if
	
	if (curr_atom.eq.bond2(i)) then
	    ! если мы не обходили вершину, то посмотрим еЁ
	    if (weight(bond1(i)).ge.num_atom) then
		call GetWeight(num_atom,num_bond, bond1, bond2, weight, bond1(i))
	    end if
	end if
  end do
  
end subroutine GetWeight

!
!> 	@brief	находим самый длинный путь в дереве
!!	@param	num_bond	[IN]	число связей
!!	@param	bond1		[IN]	
!!	@param	bond2		[IN]	массивы связей
!!	@param	num_atom	[IN]	число атомов
!!	@param	list		[IN]	степени атомов
!!	@param	max_len		[OUT]	длина максимального пути в дереве
subroutine GetLongestWay(num_bond,bond1, bond2, num_atom, list, max_len)
    ! входные
    integer(4) max_len, num_bond, num_atom
    integer(4)	list(num_atom), bond1(num_bond), bond2(num_bond)
    ! локальные
    integer(4) i, j, k
    integer(4)	weight(num_atom)

     !BODY
    max_len = 0
				
    ! перебираем степени вершин
    do j=1, num_atom
	! если нашли краевой
	if (list(j).eq.1) then
	      ! пока у нас нет пути:
	      weight(:) = num_atom ! равносильно бесконечности 
	      weight(j) = 0 
	      
	      ! ищем атом с которым связан концевой:
	      do i=1, num_bond
		if (bond1(i).eq.j) then
		    call GetWeight(num_atom,num_bond, bond1, bond2, weight, bond2(i))
		end if
		
		if (bond2(i).eq.j) then
		    call GetWeight(num_atom,num_bond, bond1, bond2, weight, bond1(i))
		end if
	      end do

	      
	      ! тут уже все веса должны посчитаться
	      ! смотрим максимальный вес:
	      do k = 1, num_atom
		if (weight(k).gt.max_len.and.weight(i).lt.num_atom) then
		  max_len = weight(k)		  
		end if
	      end do
		
	endif ! обработка для одного краевого
    enddo  
	
    return
end subroutine GetLongestWay
      
end Module TestErrors