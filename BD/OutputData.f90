!> @brief	Работа с входными файлами и данными
!! @date		09.03.2016
!!
!! @author	Михайлов И.В., Шавыкин О.В.
Module OutputData
    Use CommonModule
    Use TestErrors 
    Use PeriodicCondition
Contains 

!> @brief проверяем существует ли файл с траеториями:
subroutine CheckTRACK(num_atom, x, y, z, nstep, &
                    & box, alpha_box)
    integer(4) num_atom
    real(8) x(num_atom), y(num_atom), z(num_atom)
    integer(8) nstep
    real(8) box, alpha_box
    
    ! локальные
    integer(4) ioer
    logical(1) file_exists
    integer(4) temp_int
    character(5) temp_word
    integer(4) i
                    
    ! ====================================================    
    ! Открываем выходные файлы, в которые будут записаны полученные результаты
    INQUIRE(FILE="TRACK", EXIST=file_exists) 
    if (file_exists.eqv..true.) then
      ! проверяем не пустой ли файл:
      open(n_track,file='TRACK')
      read(n_track,*, iostat=ioer) temp_word
      do i=1, num_atom
	read(n_track, *, iostat=ioer) temp_int
      end do
      close(n_track)
      
      ! если в файле есть хотя бы один снимок
      if (ioer.eq.0.0) then
	! читаем существующий файл
	open(n_track,file='TRACK', ACCESS='append')
	! отсчитываем координаты атомов, до записи номера отпечатка
	do i=1, num_atom + 1
	  backspace(n_track)
	end do
	
	! читаем число прочитанных отпечатков
	read(n_track,*) temp_word, nstep
	! прочитать последний отпечаток:
	do i=1, num_atom
	  read(n_track, *) temp_int, x(i), y(i), z(i)
	end do
	close(n_track)
      else
      !TODO: аналогичное со STATS и TREMOR
	open(n_track,file='TRACK')
#ifdef BOUNDARY_CONDITION_CUBIC
	write(n_track,'(a,1x,I0,1x,a,1x,F15.3)') "num_atom", num_atom, "box", box
#else
	write(n_track,'(a,1x,I0)') "num_atom", num_atom
#endif
	close(n_track)
	nstep = 0 ! действительное число принятых шагов
      end if 
    else
      open(n_track,file='TRACK')
#ifdef BOUNDARY_CONDITION_CUBIC
	write(n_track,'(a,1x,I0,1x,a,1x,F15.3)') "num_atom", num_atom, "box", box
#else
	write(n_track,'(a,1x,I0)') "num_atom", num_atom
#endif
      close(n_track)      
      nstep = 0 ! действительное число принятых шагов
    end if
    
    return
end subroutine CheckTRACK

subroutine CheckSTATS()    
    ! локальные
    logical(1) file_exists
                
    ! ====================================================    
    ! если нет файла, то создаём
    INQUIRE(FILE="STATS", EXIST=file_exists) 
    if (file_exists.eqv..false.) then
      ! файл с выводам энергии:
      open(n_stats,file='STATS')
      write(n_stats, '(a, a, a, a, a, a)')  '      step', &
				      '          E_v       ', &
				      '          E_lj      ', &
				      '          E_ang     ', &
				      '          sumE      '
      close(n_stats)
    end if
    
    return
end subroutine CheckSTATS

!
subroutine CheckTREMOR(num_group, natom_group, step, print_tremor)
    integer(4) num_group, natom_group
    integer(8) step
    integer(8) print_tremor
    ! локальные
    logical(1) file_exists
                
    ! ====================================================    
    ! если нет файла, то создаём
    INQUIRE(FILE="TREMOR", EXIST=file_exists) 
    if (file_exists.eqv..false.) then
      ! файл, куда выводятся силы, действующие на 
      ! атомы группы с фиксированным центром масс
      if (num_group.gt.0) then
	open(n_tremor,file='TREMOR')
	write(n_tremor, '(a,1x,I0)') "atom_group", natom_group
	close(n_tremor)
      end if
    end if
    
    if (num_group.eq.0) then
	! если нет групп частиц с фиксированном центром масс, 
	! то и нечего выводить
	print_tremor = step + 1
    end if
    
    return
end subroutine CheckTREMOR

! создаем файл с координатами конечной конфигурации
subroutine PrintCoordFinal(num_atom, X, Y, Z, &
              & num_type, type_name, atom_type)
    ! input
    integer(4) num_atom
    real(8) X(num_atom),Y(num_atom),Z(num_atom)
    integer(4) num_type
    character(2) type_name(num_type)
    integer(4) atom_type(num_atom)
    ! local
    integer(4) i

    open(n_coordf, file='COORDF')
    write(n_coordf,'(A, I0)') 'num_atom ', num_atom

    do i = 1, num_atom
        write(n_coordf,*) x(i), y(i), z(i), type_name(atom_type(i))
    enddo

    close(n_coordf)
    
end subroutine PrintCoordFinal


    !***********************************************************************
!>  @brief  Процедура создаёт "снимок" моделируемой системы в данный момент
!!       времени в виде ENT-файла и сохраняет его в каталоге с программой
!!  @param  num_atom [IN]   ЧИСЛО АТОМОВ В МОЛЕКУЛЕ
!!  @param  X    [IN]
!!  @param  Y    [IN}   МАССИВЫ КООРДИНАТ АТОМОВ
!!  @param  Z    [IN]
!!  @param  bond1    [IN]
!!  @param  bond2    [IN] массивы с номерами связанных атомов
!!  @param  bgroup_atom     [IN] для пометки групп атомов с фиксированным центром масс
!!  @param  num_pic  [IN] ЦЕЛОЕ ЧИСЛО, ХАРАКТЕРИЗУЮЩЕЕ "СНИМОК"
   subroutine PictureTypeAtom(num_atom, X, Y, Z, &
              & num_bond, bond1, bond2, &
              & num_type, type_name, atom_type, num_pic)
!***********************************************************************
    integer(8) num_pic
    integer(4) num_atom, num_bond
    integer(4) :: bond1(num_bond), bond2(num_bond)

    real(8) X(num_atom),Y(num_atom),Z(num_atom)
    integer(4) num_type
    character(2) type_name(num_type)
    integer(4) atom_type(num_atom)
    ! внутренние
    integer(4) i

    character(30) pach
    character(32) A0
    character(6) A1,A2
    character(3) A3
  
    parameter (A1='HETATM',A2='CONECT',A3='END')

    write(A0,*) num_pic
    A0=ADJUSTL(A0)
    write(pach,'(A,A,A)') 'PICTURE_',A0(1:LEN_TRIM(A0)),'.ENT'

    OPEN(2,FILE=pach)
    do i = 1, num_atom
        write(2,'(A,I5,2X,A2,6X,I5,4X,F8.3,F8.3,F8.3)')A1,I,type_name(atom_type(i)),I,X(I),Y(I),Z(I)
    end do

    do i=1, num_bond
       write(2,'(A,I5,I5)') A2,bond1(i),bond2(i)
    end do

    endfile(2)
    close(2)

    return
    end subroutine PictureTypeAtom

!*************************************************************
!>  @brief  Процедура создаёт "снимок" моделируемой системы в данный момент
!!       времени в виде ENT-файла и сохраняет его в каталоге с программой
!!  @param  num_atom [IN]   ЧИСЛО АТОМОВ В МОЛЕКУЛЕ
!!  @param  X    [IN]
!!  @param  Y    [IN}   МАССИВЫ КООРДИНАТ АТОМОВ
!!  @param  Z    [IN]
!!  @param  bond1    [IN]
!!  @param  bond2    [IN] массивы с номерами связанных атомов
!!  @param  bgroup_atom     [IN] для пометки групп атомов с фиксированным центром масс
!!  @param  num_pic  [IN] ЦЕЛОЕ ЧИСЛО, ХАРАКТЕРИЗУЮЩЕЕ "СНИМОК"
   subroutine Picture_Ex(num_atom, X, Y, Z, &
              & num_bond, bond1, bond2, &
              & bgroup_atom, box, alpha_box, num_pic)
    integer(8) num_pic
    integer(4) num_atom, num_bond
    integer(4) :: bond1(num_bond),bond2(num_bond)

    real(8) X(num_atom),Y(num_atom),Z(num_atom)
    logical(4) bgroup_atom(num_atom)
    real(8) box, alpha_box
    
#ifdef BOUNDARY_CONDITION_CUBIC
    call PictureBoxCubic(num_atom, X, Y, Z, num_bond, bond1, bond2, &
		      & bgroup_atom, box, num_pic)
    return
#else
    call Picture(num_atom, X, Y, Z, &
              & num_bond, bond1, bond2, &
              & bgroup_atom, num_pic)
#endif
    return
end subroutine Picture_Ex

    !***********************************************************************
!>  @brief  Процедура создаёт "снимок" моделируемой системы в данный момент
!!       времени в виде ENT-файла и сохраняет его в каталоге с программой
!!  @param  num_atom [IN]   ЧИСЛО АТОМОВ В МОЛЕКУЛЕ
!!  @param  X    [IN]
!!  @param  Y    [IN}   МАССИВЫ КООРДИНАТ АТОМОВ
!!  @param  Z    [IN]
!!  @param  bond1    [IN]
!!  @param  bond2    [IN] массивы с номерами связанных атомов
!!  @param  bgroup_atom     [IN] для пометки групп атомов с фиксированным центром масс
!!  @param  num_pic  [IN] ЦЕЛОЕ ЧИСЛО, ХАРАКТЕРИЗУЮЩЕЕ "СНИМОК"
   subroutine Picture(num_atom, X, Y, Z, &
              & num_bond, bond1, bond2, &
              & bgroup_atom, num_pic)
!***********************************************************************
    integer(8) num_pic
    integer(4) num_atom, num_bond
    integer(4) :: bond1(num_bond),bond2(num_bond)

    real(8) X(num_atom),Y(num_atom),Z(num_atom)
    logical(4) bgroup_atom(num_atom)
    ! внутренние
    integer(4) i

    character(30) pach
    character(32) A0
    character(6) A1,A2
    character(3) A3
    character(1) A4
    character(1) A_FIX
    parameter (A1='HETATM',A2='CONECT',A3='END',A4='C', A_FIX='O')

    write(A0,*) num_pic
    A0=ADJUSTL(A0)
    write(pach,'(A,A,A)') 'PICTURE_',A0(1:LEN_TRIM(A0)),'.ENT'

    OPEN(2,FILE=pach)
    do i=1, num_atom
      ! помечаем по-разному атомы
      if (bgroup_atom(i).eqv..true.) then
         write(2,'(A,I5,2X,A,7X,I5,4X,F8.3,F8.3,F8.3)')A1,I,A_FIX,I,X(I),Y(I),Z(I)
      else
         write(2,'(A,I5,2X,A,7X,I5,4X,F8.3,F8.3,F8.3)')A1,I,A4,I,X(I),Y(I),Z(I)
      end if
    end do

    do i=1, num_bond
       write(2,'(A,I5,I5)') A2,bond1(i),bond2(i)
    end do

    endfile(2)
    close(2)

    return
    end subroutine Picture
!***********************************************************************

#ifdef BOUNDARY_CONDITION_CUBIC
!***********************************************************************
!>  @brief  Процедура создаёт "снимок" моделируемой системы в данный момент
!!       времени в виде ENT-файла и сохраняет его в каталоге с программой
!!  @param  num_atom [IN]   ЧИСЛО АТОМОВ В МОЛЕКУЛЕ
!!  @param  X    [IN]
!!  @param  Y    [IN}   МАССИВЫ КООРДИНАТ АТОМОВ
!!  @param  Z    [IN]
!!  @param  bond1    [IN]
!!  @param  bond2    [IN] массивы с номерами связанных атомов
!!  @param  bgroup_atom     [IN] для пометки групп атомов с фиксированным центром масс
!!  @param  box [IN] кубические границы
!!  @param  num_pic  [IN] ЦЕЛОЕ ЧИСЛО, ХАРАКТЕРИЗУЮЩЕЕ "СНИМОК"
   subroutine PictureBoxCubic(num_atom, X, Y, Z, &
              & num_bond, bond1, bond2, &
              & bgroup_atom, box, num_pic)
!***********************************************************************

    integer(8) num_pic
    integer(4) num_atom, num_bond
    integer(4) :: bond1(num_bond),bond2(num_bond)

    real(8) X(num_atom),Y(num_atom),Z(num_atom)
    logical(4) bgroup_atom(num_atom)
    real(8) box
    ! внутренние
    real(8) Xt(num_atom), Yt(num_atom), Zt(num_atom)
    real(8) dx, dy, dz
    ! итератор
    integer(4) i

    character(30) pach
    character(32) A0
    character(6) A1,A2
    character(3) A3
    character(1) A4
    character(1) A_FIX
    parameter (A1='HETATM',A2='CONECT',A3='END',A4='C', A_FIX='O')

    write(A0,*) num_pic
    A0=ADJUSTL(A0)
    write(pach,'(A,A,A)') 'PICTURE_',A0(1:LEN_TRIM(A0)),'.ENT'

    Xt(:)=X(:)
    Yt(:)=Y(:)
    Zt(:)=Z(:)

    ! по связям:
    do i=1,num_bond
      dx = Xt(bond1(i))-Xt(bond2(i))
      dy = Yt(bond1(i))-Yt(bond2(i))
      dz = Zt(bond1(i))-Zt(bond2(i))

      ! если большое расстояние у связи:
      if (abs(dx).gt.box/2.0) then
        ! применяем граничные условия к частице с большим номером:
        if (bond1(i).gt.bond2(i)) then
           Xt(bond1(i)) = Xt(bond1(i)) + sign(box, dx)
        else
           Xt(bond2(i)) = Xt(bond2(i)) + sign(box, dx)
        end if
      end if

      ! если большое расстояние у связи:
      if (abs(dy).gt.box/2.0) then
        ! применяем граничные условия к частице с большим номером:
        if (bond1(i).gt.bond2(i)) then
           Yt(bond1(i)) = Yt(bond1(i)) + sign(box, dy)
        else
           Yt(bond2(i)) = Yt(bond2(i)) + sign(box, dy)
        end if
      end if

      ! если большое расстояние у связи:
      if (abs(dz).gt.box/2.0) then
        ! применяем граничные условия к частице с большим номером:
        if (bond1(i).gt.bond2(i)) then
           Zt(bond1(i)) = Zt(bond1(i)) + sign(box, dz)
        else
           Zt(bond2(i)) = Zt(bond2(i)) + sign(box, dz)
        end if
      end if
    end do

    OPEN(2,FILE=pach)

    do i=1, num_atom
      ! помечаем по-разному атомы
      if (bgroup_atom(i).eqv..true.) then
         write(2,'(A,I5,2X,A,7X,I5,4X,F8.3,F8.3,F8.3)')A1,I,A_FIX,I,Xt(I),Yt(I),Zt(I)
      else
         write(2,'(A,I5,2X,A,7X,I5,4X,F8.3,F8.3,F8.3)')A1,I,A4,I,Xt(I),Yt(I),Zt(I)
      end if
    end do

    do i=1, num_bond
       write(2,'(A,I5,I5)') A2,bond1(i),bond2(i)
    end do

    endfile(2)
    close(2)

    return
    end subroutine PictureBoxCubic
!***********************************************************************
#endif

End Module OutputData
