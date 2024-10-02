!
!>	@brief	Учёт гидродинамических взаимодействий
!!
!!	@n Подсчёт тензора диффузии
!!	@n Разложение Холецкого
!!	@n Метод Фиксмана
!!      @n Метод Гейера (TEA)
Module HydroDynamic
  Use CommonModule
  Use TestErrors
  
  implicit none

  
Contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!#define HYDRODYNAMIC_INTERACTION

#ifdef HYDRODYNAMIC_INTERACTION


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                   0. выделение памяти
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine HI_AllocateMemory(num_atom, D, hi_method, B)
    integer(4) num_atom
    ! тензор диффузии:
    real(8), dimension(:,:), allocatable :: D
    !
    integer(4) hi_method
    ! ТОЛЬКО в методе Холецкого:
    real(8), dimension(:,:), allocatable :: B

    !real(8), dimension(:,:), allocatable :: fixman_matrix
    !real(8), dimension(:), allocatable :: cheb_coeff

    !real(8), dimension(:), allocatable :: C_tea
    !>>>>>>
    ! нам понадобиться тензор диффузии
    allocate(D(1:3*num_atom, 1:3*num_atom))

    ! в зависимости от метода:
    select case(hi_method)
        case(HI_METHOD_CHOLESKY)
		allocate(B(1:3*num_atom, 1:3*num_atom))
        case(HI_METHOD_CHEBYSHEV)
        	! получаем оптимальный порядок разложения:
    		!call OptimOrderExpansion(h_hi, num_atom, L_fixman)
    		! выделяем память:
    		!allocate(cheb_coeff(1:L_fixman+1))
    		!allocate(fixman_matrix(1:3*num_atom, 1:3*num_atom))
        case(HI_METHOD_TEA)
        	!allocate(C_tea(1:3*num_atom))
	case(HI_METHOD_KRYLOV)
        	
    end select

    return
end subroutine HI_AllocateMemory
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                   1. расчет Тензора Ротне-Прагера
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! тензорное произведение векторов:
subroutine TensorProduct(n, a, b, M)
  ! IN
  integer(4)	n
  real(8) a(n), b(n)
  ! OUT
  real(8) M(n, n)
  ! work
  integer(4) i,j
  
  ! находим произведение:
  do i=1, n
    do j=1, n 
      M(i,j) = a(i)*a(j)
    end do
  end do
  
  return
end subroutine TensorProduct

!>	
!!	@brief	Получаем тензор Rotne-Prager-Yamakawa
!!
!!	@param	num_atom	[IN]
!!	@param	X		[IN]
!!	@param	Y		[IN]
!!	@param	Z		[IN]
!!	@param	kT		[IN]
!!	@param	friction	[IN]
!!	@param	a		[IN]
!!	@param	D		[OUT]
subroutine RotnePragerYamakawa(num_atom, X, Y, Z, kT, friction, a, D)
    ! IN
    integer(4) num_atom
    real(8) X(num_atom), Y(num_atom), Z(num_atom)
    real(8) kT, friction, a
    ! OUT
    real(8) D(3*num_atom, 3*num_atom)
    ! локальные
    integer(4)	i, j, k, m
    real(8)	E(3,3), tR(3,3) 
    real(8)	dr(3)
    real(8)	r

    D(:,:) = 0.0
    ! заполняем единичную матрицу:
    E(:,:) = 0
    do i=1, 3
      E(i,i) = 1.0
    end do
    
    ! считаем тензор:
    do i=1, num_atom
      D(3*(i-1)+1,3*(i-1)+1) = kt/friction 
      D(3*(i-1)+2,3*(i-1)+2) = kt/friction
      D(3*(i-1)+3,3*(i-1)+3) = kt/friction
      
      do j=i+1, num_atom
	dr(1) = x(j) - x(i) 
	dr(2) = y(j) - y(i)
	dr(3) = z(j) - z(i)
	
	r = sqrt(dr(1)**2 + dr(2)**2 + dr(3)**2)

	! тензорное произведение вектора между атомами:
	call TensorProduct(3, dr, dr, tR)
	
	if (r > 2*a) then
	  do k=1, 3
	    do m=1, 3
	      D(3*(i-1)+k,3*(j-1)+m) = kT*3.0*a/(friction*4.0*r)*(E(k,m) + tR(k,m)/r**2 + &
			      & (2.0*a**2)/(3.0*r**2)*(E(k,m) -3.0*tR(k,m)/r**2))
	      ! матрица симметричная:
	      D(3*(j-1)+k,3*(i-1)+m) = D(3*(i-1)+k,3*(j-1)+m)
	    end do
	  end do
	else  
	  do k=1, 3
	    do m=1, 3
	      D(3*(i-1)+k,3*(j-1)+m) = (kT/friction)*( (1.0 - 9.0*r/(32.0*a))*E(k,m) + &
			      & (3.0/(32.0*a))*tR(k,m)/r**2)
	      ! матрица симметричная:
	      D(3*(j-1)+m,3*(i-1)+k) = D(3*(i-1)+k,3*(j-1)+m)
	    end do
	  end do
	end if
	
      end do ! for j
    end do ! for i
    
    return
    
end subroutine RotnePragerYamakawa


!> @brief перерасчет тензора диффузии
subroutine HI_TensorRPY_Recalc(nstep, recalc_diffusion_tensor, num_atom, X, Y, Z, & 
			     & dt, kT, friction, hradius, D)
    ! input
    integer(8) nstep
    integer(4) recalc_diffusion_tensor
    
    integer(4) num_atom
    real(8) X(num_atom), Y(num_atom), Z(num_atom)
    ! шаг интегрирования:
    real(8) dt
    real(8) kT
    real(8) friction
    real(8) hradius

    ! output
    real(8) D(3*num_atom, 3*num_atom)
    
    ! проверяем наступило ли время расчета тензора диффузии:
    if (modulo(nstep-1,recalc_diffusion_tensor).eq.0) then
        ! считаем тензор диффузии:
        call RotnePragerYamakawa(num_atom, X, Y, Z, kT, friction, hradius, D)
    end if

    return
end subroutine HI_TensorRPY_Recalc


!>      @brief  Определитель матрицы
!!      @param  num     [IN]    размерность  
!!      @param  A       [IN]    исходная матрица
!!      @param  det     [OUT]   определитель матрицы
!!
subroutine Determinant(n, A, det)
    integer(4) n
    real(8) A(n, n)
    real(8) det
    ! work
    real(8) matrix(n,n)
    real(8) m, temp
    integer(4) i, j, k, l ! итераторы
    logical(1) DetExists
    ! начало
    DetExists = .true.
    l = 1
    matrix(:,:) = A(:,:)
    ! Переводим к верхнетреугольному виду
    do k = 1, n-1
        if (matrix(k,k).eq.0) then
            DetExists = .false.
            do i = k+1, n
                if (matrix(i,k).ne.0) then
                    do j = 1, n
                        temp = matrix(i,j)
                        matrix(i,j)= matrix(k,j)
                        matrix(k,j) = temp
                    end do
                    DetExists = .true.
                    l=-l
                    exit
                endif
            end do
            
            if (DetExists.eqv..false.) then
                det = 0
                return
            end if
        end if
        do j = k+1, n
            m = matrix(j,k)/matrix(k,k)
            do i = k+1, n
                matrix(j,i) = matrix(j,i) - m*matrix(k,i)
            end do
        end do
    end do
    
    ! подсчитываем определитель как произведение диагональных элементов:
    det = l
    do i = 1, n
        det = det * matrix(i,i)
    end do
    
    return
end subroutine Determinant

!>
subroutine TensorIsPositiveDefine_Test(n, A)
	integer(4) n
	real(8) A(n, n)
	! local
	integer(4) m
	real(8) det

	! проверяем все угловые миноры:
	if (A(1,1).le.0.0) then
		write(*, *) 'Error: Second law of Thermodynamic is unhappy:' 
		write(*, *) 'Tensor RPY is not a positive defined matrix'
		stop
	endif
	do m = 2, n
		call Determinant(m, A(1:m, 1:m), det)
		if (det.le.0.0) then
			write(*, *) 'Error: Second law of Thermodynamic is unhappy:' 
			write(*, *) 'Tensor RPY is not a positive defined matrix'
			stop
		endif
	enddo
	return
end subroutine TensorIsPositiveDefine_Test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                   2. Общий расчет для смещений положений частиц
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!!	@brief	Расчет смещений частиц с учетом гидродинамических взаимодействий
!!
!!	@param	num_atom [IN]	число частиц
!!	@param	Rx       [IN]	
!!	@param	Ry       [IN]	компоненты случайных сил
!!	@param	Rz       [IN]	
!!	@param	Fx       [IN]	
!!	@param	Fy       [IN]	силы действующие на частицы
!!	@param	Fz       [IN]	
!!      @param  dt       [IN]
!!	@param	D	 [IN]  тензор диффузии
!!      @param  delta_coord_x   [OUT]
!!      @param  delta_coord_y   [OUT]   компоненты итогового перемещения частиц
!!      @param  delta_coord_z   [OUT]
!!
subroutine HydrodynamicInteraction(num_atom, Rx, Ry, Rz, Fx, Fy, Fz, &
				  & dt, D, delta_coord_x, delta_coord_y, delta_coord_z)
      ! in
      integer(4) num_atom
      ! силы действующие на частицы:
      real(8) Rx(num_atom), Ry(num_atom), Rz(num_atom)
      real(8) Fx(num_atom), Fy(num_atom), Fz(num_atom)
      ! шаг интегрирования:
      real(8) dt
      real(8) D(3*num_atom, 3*num_atom)
      ! out
      real(8) delta_coord_x(num_atom), delta_coord_y(num_atom), delta_coord_z(num_atom)
      ! локальные
      integer(4) i,j
      real(8) sum_x, sum_y, sum_z
      
      ! считаем приращения:
      do i=1, num_atom
	sum_x = 0.0
	sum_y = 0.0
	sum_z = 0.0

	! основные силы
	do j=1, num_atom
	  sum_x = sum_x + D(3*(i-1)+1,3*(j-1) + 1)*Fx(j) 
	  sum_x = sum_x + D(3*(i-1)+1,3*(j-1) + 2)*Fy(j) 
	  sum_x = sum_x + D(3*(i-1)+1,3*(j-1) + 3)*Fz(j) 
	  
	  sum_y = sum_y + D(3*(i-1)+2,3*(j-1) + 1)*Fx(j) 
	  sum_y = sum_y + D(3*(i-1)+2,3*(j-1) + 2)*Fy(j) 
	  sum_y = sum_y + D(3*(i-1)+2,3*(j-1) + 3)*Fz(j) 
	  
	  sum_z = sum_z + D(3*(i-1)+3,3*(j-1) + 1)*Fx(j) 
	  sum_z = sum_z + D(3*(i-1)+3,3*(j-1) + 2)*Fy(j) 
	  sum_z = sum_z + D(3*(i-1)+3,3*(j-1) + 3)*Fz(j) 
	end do

        ! случайные:
	delta_coord_x(i) = sum_x*dt + Rx(i)
	delta_coord_y(i) = sum_y*dt + Ry(i)
	delta_coord_z(i) = sum_z*dt + Rz(i)
      end do


      return
end subroutine HydrodynamicInteraction


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                   3. функции для метода Холецкого
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>
!!	@brief	разложение Холецкого
!!	D = B B*, B* - транспонированная матрица B
!!
!!	@param	N	[IN]	размерность
!!	@param	D	[IN]	симметричная матрица
!!	@param	B	[OUT]	нижнетреугольная матрица
!!
!! https://algowiki-project.org/ru/%D0%A0%D0%B0%D0%B7%D0%BB%D0%BE%D0%B6%D0%B5%D0%BD%D0%B8%D0%B5_%D0%A5%D0%BE%D0%BB%D0%B5%D1%86%D0%BA%D0%BE%D0%B3%D0%BE_%28%D0%BC%D0%B5%D1%82%D0%BE%D0%B4_%D0%BA%D0%B2%D0%B0%D0%B4%D1%80%D0%B0%D1%82%D0%BD%D0%BE%D0%B3%D0%BE_%D0%BA%D0%BE%D1%80%D0%BD%D1%8F%29
!!
!!	@comment почему-то не работает в случае фантомных взаимодействий, возможно нарушается симметричность или положительная определенность
subroutine CholeskyDecomposition(N, D, B)
  ! IN
  integer(4)	N
  real(8)	D(N, N)
  ! OUT
  real(8)	B(N, N)
  
  ! вспомогательное
  integer(4)	i, j, k
  real(8) recur_sum
  
  B(:,:) = 0.0
  
  do i = 1, N
    recur_sum = 0.0
  
    do k = 1, i - 1
      recur_sum = recur_sum + B(i, k)**2
    end do
  
    ! диагональные элементы:
    B(i,i) = dsqrt(D(i,i) - recur_sum)
    if (D(i,i).lt.recur_sum) then
	open(n_infor,file='INFOR', access='append')
	write(n_infor, *) 'WARNING: HI: second law of thermodynamics'
	close(n_infor)

	write(*, *) 'WARNING: HI: second law of thermodynamics'
	stop
    endif
    
    do j = i + 1, N
      recur_sum = 0.0
      do k = 1, i -1
	recur_sum = recur_sum + B(i, k)*B(j,k)
      end do
      
      ! то что ниже диагонали:
      B(j,i) = (D(j,i) - recur_sum)/B(i,i)
    end do ! do j
 enddo ! do i

   do i = 1, N
    	do j = 1, N
      	   if (isnan(B(i,j)).eqv..true.) then
		open(n_infor,file='INFOR', access='append')
		write(n_infor, *) 'WARNING: HI: CholeskyDecomposition : nan'
		close(n_infor)

		write(*, *) 'WARNING: HI: CholeskyDecomposition : nan'
		stop
	   endif
    	enddo
    enddo

    return
end subroutine CholeskyDecomposition

!>
!!	@brief	разложение Холецкого
!!	D = B B*, B* - транспонированная матрица B
!!
!!	@param	N	[IN]	размерность
!!	@param	D	[IN]	симметричная матрица
!!	@param	B	[OUT]	нижнетреугольная матрица
!!
!!
!!	@comment также проблемы в случае фантомных моделей
subroutine CholeskyDecomposition2(N, D, B)
  ! IN
  integer(4)	N
  real(8)	D(N, N)
  ! OUT
  real(8)	B(N, N)
  
  ! вспомогательное
  integer(4)	i, j, k
  real(8) recur_sum
  
  B(:,:) = 0.0
  
  do i = 1, N
    recur_sum = 0.0
  
    do k = 1, i - 1
      recur_sum = recur_sum + B(i, k)**2
    end do
  
    ! диагональные элементы:
    B(i,i) = dsqrt(D(i,i) - recur_sum)
    if (D(i,i).lt.recur_sum) then
	open(n_infor,file='INFOR', access='append')
	write(n_infor, *) 'WARNING: HI: second law of thermodynamics'
	close(n_infor)

	write(*, *) 'WARNING: HI: second law of thermodynamics'
	stop
    endif
    
    do j = 1, i-1
      recur_sum = 0.0
      do k = 1, j -1
	recur_sum = recur_sum + B(i, k)*B(j,k)
      end do
      
      ! то что ниже диагонали:
      B(i,j) = (D(i,j) - recur_sum)/B(j,j)
    end do ! do j
    
 enddo ! do i

    do i = 1, N
    	do j = 1, N
      	   if (isnan(B(i,j)).eqv..true.) then
		open(n_infor,file='INFOR', access='append')
		write(n_infor, *) 'WARNING: HI: CholeskyDecomposition : nan'
		close(n_infor)

		write(*, *) 'WARNING: HI: CholeskyDecomposition : nan'
		stop
	   endif
    	enddo
    enddo

    return
end subroutine CholeskyDecomposition2
    
    
subroutine HI_RandomCholesky(num_atom, B, grand_x, grand_y, grand_z, dt, &
                            & Rx, Ry, Rz)
      ! input
      integer(4) num_atom
      real(8) B(3*num_atom, 3*num_atom)
      real(8) grand_x(num_atom), grand_y(num_atom), grand_z(num_atom)
      real(8) dt
      ! output
      real(8) Rx(num_atom), Ry(num_atom), Rz(num_atom)
      ! work
      integer(4) i, j
      
      Rx(:) = 0.0
      Ry(:) = 0.0
      Rz(:) = 0.0
      ! считаем приращения:
      do i=1, num_atom
	! основные силы
	do j=1, num_atom
 	  Rx(i) = Rx(i) + B(3*(i-1)+1, 3*(j-1) + 1)*grand_x(i)
	  Rx(i) = Rx(i) + B(3*(i-1)+1, 3*(j-1) + 2)*grand_y(i)
	  Rx(i) = Rx(i) + B(3*(i-1)+1, 3*(j-1) + 3)*grand_z(i)
	  
	  Ry(i) = Ry(i) + B(3*(i-1)+2, 3*(j-1) + 1)*grand_x(i)
	  Ry(i) = Ry(i) + B(3*(i-1)+2, 3*(j-1) + 2)*grand_y(i)
	  Ry(i) = Ry(i) + B(3*(i-1)+2, 3*(j-1) + 3)*grand_z(i)
	  
	  Rz(i) = Rz(i) + B(3*(i-1)+3, 3*(j-1) + 1)*grand_x(i)
	  Rz(i) = Rz(i) + B(3*(i-1)+3, 3*(j-1) + 2)*grand_y(i)
	  Rz(i) = Rz(i) + B(3*(i-1)+3, 3*(j-1) + 3)*grand_z(i)
	end do
     end do
    
    ! нормируем:
    Rx(:) =Rx(:)*sqrt(2.0*dt)
    Ry(:) =Ry(:)*sqrt(2.0*dt)
    Rz(:) =Rz(:)*sqrt(2.0*dt)
      
    return
end subroutine HI_RandomCholesky    

!> @brief перерасчет тензора для случайных сил в методе Холецкого
subroutine HI_CholeskyDecomposition_Recalc(nstep, recalc_diffusion_tensor, num_atom, D, B)
    ! input
    integer(8) nstep
    integer(4) recalc_diffusion_tensor
    integer(4) num_atom
    real(8) D(3*num_atom, 3*num_atom)
    ! output 
    real(8) B(3*num_atom, 3*num_atom)
    
    ! проверяем наступило ли время расчета тензора диффузии:
    if (modulo(nstep-1,recalc_diffusion_tensor).eq.0) then
        ! считаем тензор диффузии:
        call CholeskyDecomposition(3*num_atom, D, B)
    end if

    return
end subroutine HI_CholeskyDecomposition_Recalc



!>>
!  Расчет итогового перемещений
!>
subroutine HI_Displacement_Cholesky(nstep, recalc_diffusion_tensor, num_atom, X, Y, Z, &
    & grand_x, grand_y, grand_z, Fx, Fy, Fz, dt, kT, friction, hradius, D, &   
    & B, delta_coord_x, delta_coord_y, delta_coord_z)
    ! in
    integer(8) nstep
    integer(4) recalc_diffusion_tensor
    integer(4) num_atom
    real(8) X(num_atom), Y(num_atom), Z(num_atom)
    ! силы действующие на частицы:
    real(8) grand_x(num_atom), grand_y(num_atom), grand_z(num_atom)
    real(8) Fx(num_atom), Fy(num_atom), Fz(num_atom)
    ! шаг интегрирования:
    real(8) dt
    real(8) kT
    real(8) friction
    real(8) hradius
    real(8) D(3*num_atom, 3*num_atom)
       
    real(8) B(3*num_atom, 3*num_atom)

    ! out
    real(8) delta_coord_x(num_atom), delta_coord_y(num_atom), delta_coord_z(num_atom)
    ! work  
    real(8) Rx(num_atom), Ry(num_atom), Rz(num_atom)
    
    ! если нужно пересчитываем тензор диффузии:
    call HI_TensorRPY_Recalc(nstep, recalc_diffusion_tensor, num_atom, X, Y, Z, & 
			     & dt, kT, friction, hradius, D)

    ! если обновился тензор диффузии, то делаем разложение Холецкого:
    call HI_CholeskyDecomposition_Recalc(nstep, recalc_diffusion_tensor, num_atom, D, B)

    ! находим вектор случайных сил:
    call HI_RandomCholesky(num_atom, B, grand_x, grand_y, grand_z, dt, Rx, Ry, Rz)


    ! делаем итоговое перемещение:
    call HydrodynamicInteraction(num_atom, Rx, Ry, Rz, Fx, Fy, Fz, &
                            & dt, D, delta_coord_x, delta_coord_y, delta_coord_z)
    

    return
end subroutine HI_Displacement_Cholesky

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                   4. функции для метода Чебышева (Фиксмана)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
subroutine EigenValuesMaxMin(N, D, lambda_max, lambda_min)
    ! input
    integer(4) N
    real(8) D(N, N)
    ! output
    real(8) lambda_max, lambda_min
    
    call MaxEigenvalues(N, D, lambda_max)
    call MinEigenvalues(N, D, lambda_max, lambda_min)

    return
end subroutine EigenValuesMaxMin

!>
!!	@brief	
!!
!!	@param	N	[IN]	размерность
!!	@param	D	[IN]	симметричная матрица
!!	@param	v	[OUT]	собственные числа
!!
subroutine MaxEigenvalues(n, D, lambda_max)
   ! точность:
   real(8), parameter :: eps = 1e-7
   integer(4), parameter :: max_iter = 10000
   ! in
   integer(4) N
   real(8) D(N, N)
   ! out
   real(8) lambda_max
   ! work
   real(8) A(N,N)
   real(8) X(N)
   real(8) Y(N)
   real(8) lambda_s
   real(8) mult_yx, mult_xx
   real(8) normY

   
   integer(4) k
   
   ! чтобы не изменить входную матрицу:
   A(:,:) = D(:,:)
   
   X(:) = 1.0
   !
   k = 0
   do while (k.lt.max_iter)
      k = k + 1
      
      Y = matmul(A, X)
      
      if (k.ge.2) then
         lambda_s = lambda_max
      endif
      
      mult_yx = sum(Y(:)*X(:))
      mult_xx = sum(X(:)*X(:))
      
      lambda_max = mult_yx/mult_xx
      
      ! или сумма?
      !normY = maxval(Y(:))
      ! евклидова норма:
      normY = sqrt(sum(Y(:)**2))
      
      X = (1.0/normY)*Y
      if (k.ge.2) then
        if (abs(lambda_max - lambda_s)/abs(lambda_max).le.eps) then
            exit
        endif
      endif
   enddo
   
   if (k.ge.max_iter) then
        write(*, *) 'error: max_step-1!'
   endif
   
   return
end subroutine MaxEigenvalues 

!>
!!	@brief	
!!
!!	@param	N	[IN]	размерность
!!	@param	D	[IN]	симметричная матрица
!!	@param	v	[OUT]	собственные числа
!!
subroutine MinEigenvalues(n, D, lambda_max, lambda_min)
   ! точность:
   real(8), parameter :: eps = 1e-7
   integer(4), parameter :: max_iter = 10000
   real(8), parameter :: koff = 1.01 ! 1 процент
   ! in
   integer(4) N
   real(8) D(N, N)
   real(8) lambda_max
   ! out
   real(8) lambda_min
   ! work
   real(8) A(N,N)
   real(8) X(N)
   real(8) Y(N)
   real(8) lambda_s
   real(8) mult_yx, mult_xx
   integer(4) i
   real(8) normY

   integer(4) k
   
   ! чтобы не изменить входную матрицу:
   A(:,:) = -D(:,:)
   do i = 1, N
      A(i,i) = lambda_max*koff + A(i,i)
   enddo
   
   X(:) = 1.0
   !
   k = 0
   do while (k.lt.max_iter)
      k = k + 1
    
      Y = matmul(A, X)
      
      if (k.ge.2) then
         lambda_s = lambda_min
      endif
      
      mult_yx = sum(Y(:)*X(:))
      mult_xx = sum(X(:)*X(:))
      
      lambda_min = mult_yx/mult_xx
      
      ! или сумма?
      ! евклидова норма:
      normY = sqrt(sum(Y(:)**2))
      
      X = (1.0/normY)*Y
      if (k.ge.2) then
        if (abs(lambda_min - lambda_s)/abs(lambda_min).le.eps) then
            exit
        endif
      endif
   enddo
   
   lambda_min = koff*lambda_max - lambda_min
   
   if (k.ge.max_iter) then
        write(*, *) 'error: max_step-2!'
   endif
   
   return
end subroutine MinEigenvalues 

!>
!!	@brief	
!!	@param	N	[IN]	размерность
!!	@param	D	[IN]	симметричная матрица
!!      @param  L       [IN]
!!
subroutine ChebyshevCoef(N, L, lambda_min, lambda_max, a)
    ! IN
    integer(4)	N
    integer(4) L
    real(8) lambda_max, lambda_min ! наибольшее и наименьшее
    
    ! OUT
    real(8) a(L+1)
    
    real(8) a_plus, a_minus
    
    real(8) lambda(1:L+1)
    real(8) c(1:L+1, 1:L+1)
    real(8) normC(1:L+1)
    real(8) ee, ff
    !
    real(8), parameter :: PI = 3.1415926535897932
    ! итераторы:
    integer(4) i, k
    integer(4) j
    
   
    ! ===========================>
    a_plus =  (lambda_max + lambda_min)/2.0
    a_minus = (lambda_max - lambda_min)/2.0
      
    do k = 1, L + 1
        lambda(k) = cos(PI*(k-0.5)/(L+1))
    enddo
    
    !
    c(1, :) = 1
    c(2, :) = lambda(:)
    do j = 2, L
        do k = 1, L + 1
            c(j+1, k) = 2*lambda(k)*c(j, k) - c(j-1, k) 
        enddo
    enddo

    do i = 1, L + 1
        normC(i) = sum(c(i, :)**2) ! it's norm without square root!!!
    enddo
    a(:) = 0.0
    do i = 1, L + 1
        do k = 1, L + 1
            call ShiftedSquareRoot(a_minus, a_plus, lambda(k), ff)
            
            a(i) = a(i) + c(i, k)*ff/normC(i)
        enddo
    enddo
 
    return
end subroutine ChebyshevCoef

subroutine ShiftedSquareRoot(a_minus, a_plus, e, f)
    ! input
    real(8) a_plus, a_minus, e
    ! output
    real(8) f
    
    f = (a_plus + a_minus*e)**(1.0/2.0)

    return
end subroutine ShiftedSquareRoot



!>>>
subroutine HI_EigenValues_Recalc(nstep, recalc_diffusion_tensor, num_atom, D, &
			& L, lambda_max, lambda_min)
    ! IN
    integer(8) nstep
    integer(4) recalc_diffusion_tensor
    integer(4)	num_atom
    real(8)	D(3*num_atom, 3*num_atom)
    !  OUT
    integer(4) L
    real(8) lambda_max, lambda_min ! наибольшее и наименьшее

    !>>
    ! проверяем наступило ли время перерасчета:
    if (modulo(nstep-1,recalc_diffusion_tensor).eq.0) then
    	! рассчитываем собственный числа матрицы D:
    	call EigenValuesMaxMin(3*num_atom, D, lambda_max, lambda_min)
    	L = int(sqrt(lambda_max/lambda_min)) + 1
    endif

    return
end subroutine HI_EigenValues_Recalc

!>>
subroutine EigenValuesFixman(num_atom, D, lambda_max, lambda_min)
    integer(4)	num_atom
    real(8)	D(3*num_atom, 3*num_atom)
    !  OUT
    real(8) lambda_max, lambda_min ! наибольшее и наименьшее
   ! local
    integer(4) i, j

   !>>> 1. Максимальное собственное число
    ! - умножаем на единичный вектор с двух сторон: 
    lambda_max = sum(D(:, :)) / (3*num_atom)
 
    lambda_min = 0.0
    do i = 1, 3*num_atom
	do j = 1, 3*num_atom
	    lambda_min = lambda_min + (-1)**(i + j - 2)*D(i, j)
	enddo
    enddo
    lambda_min = lambda_min / (3*num_atom)

    return
end subroutine EigenValuesFixman
!>>>
subroutine HI_EigenValuesFixmanApprox_Recalc(nstep, recalc_diffusion_tensor, num_atom, D, &
			& L, lambda_max, lambda_min)
    ! IN
    integer(8) nstep
    integer(4) recalc_diffusion_tensor
    integer(4)	num_atom
    real(8)	D(3*num_atom, 3*num_atom)
    !  OUT
    integer(4) L
    real(8) lambda_max, lambda_min ! наибольшее и наименьшее
    !>>
    ! проверяем наступило ли время перерасчета:
    if (modulo(nstep-1,recalc_diffusion_tensor).eq.0) then
    	! рассчитываем собственный числа матрицы D:
    	call EigenValuesFixman(num_atom, D, lambda_max, lambda_min)
    	L = int(sqrt(lambda_max/lambda_min)) + 1
    endif

    return
end subroutine HI_EigenValuesFixmanApprox_Recalc

!>>>
subroutine HI_Chebyshev_Recalc(nstep, recalc_diffusion_tensor, num_atom, D, &
				& L, lambda_max, lambda_min, k1, k2, cheb_coeff)
    ! IN
    integer(8) nstep
    integer(4) recalc_diffusion_tensor
    integer(4)	num_atom
    real(8)	D(3*num_atom, 3*num_atom)
    integer(4) L
    real(8) lambda_max, lambda_min ! наибольшее и наименьшее
    !  OUT
    real(8) k1, k2
    real(8) cheb_coeff(L+1)
    !>>
    ! проверяем наступило ли время перерасчета:
    if (modulo(nstep-1,recalc_diffusion_tensor).eq.0) then
    	call ChebyshevCoef(3*num_atom, L, lambda_min, lambda_max, cheb_coeff)
 
    	k1 = 2.0/(lambda_max - lambda_min)
    	k2 = -(lambda_max + lambda_min)/(lambda_max - lambda_min)
    endif

    return
end subroutine HI_Chebyshev_Recalc

!>>>>>>>
subroutine ChebyshevApproximation(num_atom, D, k1, k2, cheb_coeff, L, dt, &
				& rand_x, rand_y, rand_z, Rx, Ry, Rz)
    ! IN
    integer(4)	num_atom
    real(8)	D(3*num_atom, 3*num_atom)
    integer(4) L
    real(8) lambda_max, lambda_min ! наибольшее и наименьшее
    real(8) dt
    real(8) rand_x(num_atom), rand_y(num_atom), rand_z(num_atom)
    real(8) k1, k2
    real(8) cheb_coeff(L+1)
    ! OUT
    real(8) Rx(num_atom), Ry(num_atom), Rz(num_atom)
    ! work
    
    real(8) z(0:L+1, 3*num_atom)
    integer(4) i, j
    
    !>>>>>>
    z(1, :) = 0.0
    do i = 1, num_atom
        z(1, 3*(i-1) + 1) = rand_x(i)
        z(1, 3*(i-1) + 2) = rand_y(i)
        z(1, 3*(i-1) + 3) = rand_z(i)
    enddo
    z(2, :) = k1*MATMUL(D, z(1, :)) + k2*z(1, :)
    
    do i = 2, L
        z(i+1, :) =  2.0*k1*MATMUL(D, z(i, :)) + 2.0*k2*z(i, :) - z(i-1, :)
    end do
    
    Rx(:) = 0.0
    Ry(:) = 0.0
    Rz(:) = 0.0
    do i = 1, L + 1
	do j = 1, num_atom
        	Rx(j) = Rx(j) + cheb_coeff(i)*z(i, 3*(j-1) + 1)
		Ry(j) = Ry(j) + cheb_coeff(i)*z(i, 3*(j-1) + 2)
		Rz(j) = Rz(j) + cheb_coeff(i)*z(i, 3*(j-1) + 3)
	enddo
    end do

   ! нормируем:
    Rx(:) =Rx(:)*sqrt(2.0*dt)
    Ry(:) =Ry(:)*sqrt(2.0*dt)
    Rz(:) =Rz(:)*sqrt(2.0*dt)

    
    return
end subroutine ChebyshevApproximation


!>>
!  Расчет итогового перемещений
!>
subroutine HI_Displacement_Chebyshev(nstep, recalc_diffusion_tensor, num_atom, X, Y, Z, &
    & grand_x, grand_y, grand_z, Fx, Fy, Fz, dt, kT, friction, hradius, D, &   
    & delta_coord_x, delta_coord_y, delta_coord_z)
    ! in
    integer(8) nstep
    integer(4) recalc_diffusion_tensor
    integer(4) num_atom
    real(8) X(num_atom), Y(num_atom), Z(num_atom)
    ! силы действующие на частицы:
    real(8) grand_x(num_atom), grand_y(num_atom), grand_z(num_atom)
    real(8) Fx(num_atom), Fy(num_atom), Fz(num_atom)
    ! шаг интегрирования:
    real(8) dt
    real(8) kT
    real(8) friction
    real(8) hradius
    real(8) D(3*num_atom, 3*num_atom)
    ! local
    real(8) lambda_max, lambda_min
    integer(4) L
    real(8), allocatable, dimension(:) :: cheb_coeff
    real(8) k1, k2

    ! out
    real(8) delta_coord_x(num_atom), delta_coord_y(num_atom), delta_coord_z(num_atom)
    ! work  
    real(8) Rx(num_atom), Ry(num_atom), Rz(num_atom)
    
    ! если нужно пересчитываем тензор диффузии:
    call HI_TensorRPY_Recalc(nstep, recalc_diffusion_tensor, num_atom, X, Y, Z, & 
			     & dt, kT, friction, hradius, D)

    call HI_EigenValues_Recalc(nstep, recalc_diffusion_tensor, num_atom, D, &
			& L, lambda_max, lambda_min)

    allocate(cheb_coeff(1:L+1))
    call HI_Chebyshev_Recalc(nstep, recalc_diffusion_tensor, num_atom, D, &
			& L, lambda_max, lambda_min, k1, k2, cheb_coeff)
    ! делаем приближение случайных сил методов Чебышева (Фиксмана)
    call ChebyshevApproximation(num_atom, D, k1, k2, cheb_coeff, L, dt, &
				& grand_x, grand_y, grand_z, Rx, Ry, Rz)
    deallocate(cheb_coeff)
    ! делаем итоговое перемещение:
    call HydrodynamicInteraction(num_atom, Rx, Ry, Rz, Fx, Fy, Fz, &
                            & dt, D, delta_coord_x, delta_coord_y, delta_coord_z)
    

    return
end subroutine HI_Displacement_Chebyshev


!>>
!  Расчет итогового перемещений
!>
subroutine HI_Displacement_Fixman(nstep, recalc_diffusion_tensor, num_atom, X, Y, Z, &
    & grand_x, grand_y, grand_z, Fx, Fy, Fz, dt, kT, friction, hradius, D, &   
    & delta_coord_x, delta_coord_y, delta_coord_z)
    ! in
    integer(8) nstep
    integer(4) recalc_diffusion_tensor
    integer(4) num_atom
    real(8) X(num_atom), Y(num_atom), Z(num_atom)
    ! силы действующие на частицы:
    real(8) grand_x(num_atom), grand_y(num_atom), grand_z(num_atom)
    real(8) Fx(num_atom), Fy(num_atom), Fz(num_atom)
    ! шаг интегрирования:
    real(8) dt
    real(8) kT
    real(8) friction
    real(8) hradius
    real(8) D(3*num_atom, 3*num_atom)
    ! local
    real(8) lambda_max, lambda_min
    integer(4) L
    real(8), allocatable, dimension(:) :: cheb_coeff
    real(8) k1, k2

    ! out
    real(8) delta_coord_x(num_atom), delta_coord_y(num_atom), delta_coord_z(num_atom)
    ! work  
    real(8) Rx(num_atom), Ry(num_atom), Rz(num_atom)
    
    ! если нужно пересчитываем тензор диффузии:
    call HI_TensorRPY_Recalc(nstep, recalc_diffusion_tensor, num_atom, X, Y, Z, & 
			     & dt, kT, friction, hradius, D)

    call HI_EigenValuesFixmanApprox_Recalc(nstep, recalc_diffusion_tensor, num_atom, D, &
			& L, lambda_max, lambda_min)

    allocate(cheb_coeff(1:L+1))
    call HI_Chebyshev_Recalc(nstep, recalc_diffusion_tensor, num_atom, D, &
			& L, lambda_max, lambda_min, k1, k2, cheb_coeff)
    ! делаем приближение случайных сил методов Чебышева (Фиксмана)
    call ChebyshevApproximation(num_atom, D, k1, k2, cheb_coeff, L, dt, &
				& grand_x, grand_y, grand_z, Rx, Ry, Rz)
    deallocate(cheb_coeff)

    ! делаем итоговое перемещение:
    call HydrodynamicInteraction(num_atom, Rx, Ry, Rz, Fx, Fy, Fz, &
                            & dt, D, delta_coord_x, delta_coord_y, delta_coord_z)
    

    return
end subroutine HI_Displacement_Fixman

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                   5. функции для метода TEA (Гейера)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!>>>
subroutine HI_TEA_Epsilon(n, D, eps)
    ! input
    integer(4) n
    real(8) D(n, n)
    ! output
    real(8) eps
    ! local
    integer(4) i, j
    !>>>
    eps = 0.0
    
    do i = 1, n
        do j = i + 1, n
            eps = eps + D(i, j)/D(i, i)
        enddo
    enddo
    eps = eps/(dble(n-1)*dble(n)/2.0)

    return
end subroutine HI_TEA_Epsilon

!>>
subroutine HI_TEA_BetaS(n, eps, beta_s)
    ! input
    integer(4) n
    real(8) eps
    ! output
    real(8) beta_s
    ! local
    real(8) temp

    temp = (n-1)*(eps**2) - (n-2)*eps
    
    beta_s = (1 - sqrt(1 - temp))/temp 

    return
end subroutine HI_TEA_BetaS

!>>
subroutine HI_TEA_NormalConst(n, D, beta_s, c)
    ! input
    integer(4) n
    real(8) D(n, n)
    real(8) beta_s
    ! output
    real(8) c(1:n)
    ! local
    integer(4) i, j
    real(8) temp

    c(:) = 0.0
    do i = 1, n
        temp = 0.0
        do j = 1, n
            if (j.ne.i) then
                temp = temp + (D(i, j)**2)/(D(i,i)*D(j,j))
            endif
        enddo
        c(i) = sqrt(1.0/(1 + (beta_s**2)*temp))
              
    enddo

   return
end subroutine HI_TEA_NormalConst

!>> получить выражения для случайного перемещения
subroutine TruncatedExpansionApproximation(num_atom, dt, D, c, beta_s, grand_x, grand_y, grand_z, Rx, Ry, Rz)
    ! IN
    integer(4)	num_atom
    real(8) dt
    real(8)	D(3*num_atom, 3*num_atom)
    real(8)	c(3*num_atom)
    real(8) beta_s
    real(8) grand_x(num_atom), grand_y(num_atom),  grand_z(num_atom)
    ! OUT
    real(8) Rx(num_atom), Ry(num_atom),  Rz(num_atom)
    ! work
    integer(4) i, j, k
    integer(4) ix, iy, iz
    integer(4) jx, jy, jz
   

    do i = 1, num_atom
	ix = 3*(i - 1) + 1
	iy = 3*(i - 1) + 2
	iz = 3*(i - 1) + 3

	Rx(i) = sqrt(D(ix, ix))*grand_x(i)
	Ry(i) = sqrt(D(iy, iy))*grand_y(i)
	Rz(i) = sqrt(D(iz, iz))*grand_z(i)
	do j = i + 1, num_atom
	    jx = 3*(j - 1) + 1
	    jy = 3*(j - 1) + 2
	    jz = 3*(j - 1) + 3
	    !>>
	    Rx(i) = Rx(i) + beta_s*D(ix, jx)/sqrt(D(ix, ix))*grand_x(j)
	    Rx(i) = Rx(i) + beta_s*D(ix, jy)/sqrt(D(ix, ix))*grand_y(j)
	    Rx(i) = Rx(i) + beta_s*D(ix, jz)/sqrt(D(ix, ix))*grand_z(j)
	    !>>
	    Ry(i) = Ry(i) + beta_s*D(iy, jx)/sqrt(D(iy, iy))*grand_x(j)
	    Ry(i) = Ry(i) + beta_s*D(iy, jy)/sqrt(D(iy, iy))*grand_y(j)
	    Ry(i) = Ry(i) + beta_s*D(iy, jz)/sqrt(D(iy, iy))*grand_z(j)
	    !>>
	    Rz(i) = Rz(i) + beta_s*D(iz, jx)/sqrt(D(iz, iz))*grand_x(j)
	    Rz(i) = Rz(i) + beta_s*D(iz, jy)/sqrt(D(iz, iz))*grand_y(j)
	    Rz(i) = Rz(i) + beta_s*D(iz, jz)/sqrt(D(iz, iz))*grand_z(j)
	enddo

	Rx(i) = c(ix)*Rx(i)
	Ry(i) = c(iy)*Ry(i)
	Rz(i) = c(iz)*Rz(i)
    enddo

    ! нормируем:
    Rx(:) =Rx(:)*sqrt(2.0*dt)
    Ry(:) =Ry(:)*sqrt(2.0*dt)
    Rz(:) =Rz(:)*sqrt(2.0*dt)

    
    return
end subroutine TruncatedExpansionApproximation

!> @brief перерасчет тензора диффузии
subroutine HI_factorsTEA_Recalc(nstep, recalc_diffusion_tensor, num_atom, D, &
				& beta_s, c)
    ! input
    integer(8) nstep
    integer(4) recalc_diffusion_tensor
    integer(4) num_atom
    real(8) D(3*num_atom, 3*num_atom)
    
    ! output
    real(8) beta_s
    real(8) c(3*num_atom)
    ! local
    real(8) eps
    
    ! проверяем наступило ли время расчета тензора диффузии:
    if (modulo(nstep-1,recalc_diffusion_tensor).eq.0) then
        ! считаем тензор диффузии:
        call HI_TEA_Epsilon(3*num_atom, D, eps)
        call HI_TEA_BetaS(3*num_atom, eps, beta_s)
	call HI_TEA_NormalConst(3*num_atom, D, beta_s, c)

    end if

    return
end subroutine HI_factorsTEA_Recalc

!>>
!  Расчет итогового перемещений
!>
subroutine HI_Displacement_TEA(nstep, recalc_diffusion_tensor, num_atom, X, Y, Z, &
    & grand_x, grand_y, grand_z, Fx, Fy, Fz, dt, kT, friction, hradius, D, &   
    & delta_coord_x, delta_coord_y, delta_coord_z)
    ! in
    integer(8) nstep
    integer(4) recalc_diffusion_tensor
    integer(4) num_atom
    real(8) X(num_atom), Y(num_atom), Z(num_atom)
    ! силы действующие на частицы:
    real(8) grand_x(num_atom), grand_y(num_atom), grand_z(num_atom)
    real(8) Fx(num_atom), Fy(num_atom), Fz(num_atom)
    ! шаг интегрирования:
    real(8) dt
    real(8) kT
    real(8) friction
    real(8) hradius
    real(8) D(3*num_atom, 3*num_atom)
    ! local
    real(8) beta_s
    real(8) c(3*num_atom)

    ! out
    real(8) delta_coord_x(num_atom), delta_coord_y(num_atom), delta_coord_z(num_atom)
    ! work  
    real(8) Rx(num_atom), Ry(num_atom), Rz(num_atom)
    
    ! если нужно пересчитываем тензор диффузии:
    call HI_TensorRPY_Recalc(nstep, recalc_diffusion_tensor, num_atom, X, Y, Z, & 
			     & dt, kT, friction, hradius, D)
    ! пересчитываем коэффициенты разложения TEA, если изменился D:
    call HI_factorsTEA_Recalc(nstep, recalc_diffusion_tensor, num_atom, D, beta_s, c)
    ! находим вектор случайных перемещений
    call TruncatedExpansionApproximation(num_atom, dt, D, c, beta_s, grand_x, grand_y, grand_z, Rx, Ry, Rz)

    ! делаем итоговое перемещение:
    call HydrodynamicInteraction(num_atom, Rx, Ry, Rz, Fx, Fy, Fz, &
                            & dt, D, delta_coord_x, delta_coord_y, delta_coord_z)
    

    return
end subroutine HI_Displacement_TEA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                   6. функции для метода Крылова
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> Фробениусова норма матрицы элементов ниже диагонали
subroutine MatrixNormBelowDiag(n, A, m_norm)
    integer(4) n
    real(8) A(n, n)
    ! 
    real(8) m_norm
    ! local
    integer(4) i, j
    
    m_norm = 0.0
    do i = 1, n-1
        do j = i+1, n
            !m_norm = m_norm + A(i, j)**2
            m_norm = m_norm + A(j, i)**2
        enddo
    enddo
    m_norm = sqrt(m_norm)
    
    return
end subroutine MatrixNormBelowDiag

!>
subroutine SchurDecompositionQR(n, A, Qr, B)
    integer(4) n
    real(8) A(n, n)
    ! out
    real(8) B(n, n) ! собственные значения:
    real(8) Qr(n, n)
    ! local
    real(8) Q(n, n)
    real(8) R(n, n)
    real(8) error
    integer(4) i
    !
    integer(4), parameter :: num_iter = 30
    real(8), parameter :: tolerance = 1e-5
    !>>>>>>>>>>
    Qr(:, :) = 0.0
    do i = 1, n
        Qr(i, i) = 1
    enddo
    
    B(:, :) = A(:, :)
    do i = 1, num_iter
        call QR_Decomposition_Householder(n, n, B, Q, R)
        B = matmul(R, Q)
        Qr = matmul(Qr, Q)
        call MatrixNormBelowDiag(n, B, error)
        if (error.lt.tolerance) exit
    enddo
    
    !write(*, '(A, I0, A)') 'QR solver Schur decompostion are used ', i, ' iterations'
    
    return
end subroutine SchurDecompositionQR


! пример:
! http://www.studfiles.ru/preview/4600713/
subroutine QR_Decomposition_Householder(n, m, A, Q, R)
    ! входные:
    integer(4) n, m
    real(8) A(n, n)
    ! выход:
    real(8) Q(n, n)
    real(8) R(n, n)
    ! local
    real(8) P(n, n)
    real(8) E(n, n)
    real(8) v(n)
    real(8) dX
    integer(4) i, j, k
    
    !!!!!!!!!!!!
    v(:) = 0.0
    E(:, :) = 0.0
    do j = 1, n
        E(j, j) = 1.0
    enddo
    
    Q(:, :) = E(:, :)
    R(:, :) = A(:, :)
    
    ! для каждого столбца Q находим матрицу отражений:
    do i = 1, n - 1
        dX = 0.0
        do j = i, n
            v(j) = R(j, i)
            dX = dX + v(j)**2
        enddo
        dX = sqrt(dX)
        
        ! поправка ?? + или -:
        v(i) = v(i) + dsign(dX, v(i))
        
        dX = 0.0
        do j = i, n
            dX = dX + v(j)**2
        enddo
    
        P(:, :) = E(:, :)
        ! находим отражение:
        do k = i, n
            do j = i, n
                P(k, j) = E(k, j) - 2*v(k)*v(j)/dX
            enddo
        enddo
        
        ! запоминаем результат в матрицу:
        Q = matmul(Q, P)
        R = matmul(P, R)       
    enddo
    

    return
end subroutine QR_Decomposition_Householder


!>>>>>>
subroutine KrylovSquareRoot(num_atom, dt, D, m_lanczos, grand_x, grand_y, grand_z, Rx, Ry, Rz)
    ! IN
    integer(4)	num_atom
    real(8) dt
    real(8)	D(3*num_atom, 3*num_atom)
    integer(4) m_lanczos
    real(8) grand_x(num_atom), grand_y(num_atom),  grand_z(num_atom)
    ! OUT
    real(8) Rx(num_atom), Ry(num_atom),  Rz(num_atom)
    ! local
    integer(4) i, j
    integer(4) ix, iy, iz
    !>>>
    real(8) z(3*num_atom)
    real(8) norm_z

    real(8) w(3*num_atom)

    real(8) v(m_lanczos+1, 3*num_atom)
    real(8) h(m_lanczos+1, m_lanczos)

    real(8) Bm(m_lanczos, m_lanczos)
    real(8) Qr(m_lanczos, m_lanczos)
    real(8) Lambda(m_lanczos, m_lanczos)
    
    
    real(8), parameter :: tolerance = 1e-10
    
    !x(:) = 0
    do i = 1, num_atom
	ix = 3*(i - 1) + 1
	iy = 3*(i - 1) + 2
	iz = 3*(i - 1) + 3
   	z(ix) = grand_x(i)
   	z(iy) = grand_y(i)
   	z(iz) = grand_z(i)
    enddo
    norm_z = sqrt(sum(z(:)**2))
    
    v(1, :) = z(:)/norm_z
    
    h(:, :) = 0.0
    do j = 1, m_lanczos
        w = matmul(D, v(j, :))
        
	if (j.gt.1) then
		w = w - h(j-1, j)*v(j-1, :)
	endif
        h(j, j) = sum(w(:)*v(j, :))
	
	
	if (j.lt.m_lanczos) then
		w(:) = w(:) - h(j, j)*v(j, :)
		h(j+1, j) = sqrt(sum(w(:)**2))
		h(j, j+1) = h(j+1, j) 

		if (h(j+1, j).eq.0) then
            		!write(*, *) j
            		write(*, *) 'WARNING : Krylov : square root!'
            		stop
       		endif
		v(j+1, :) = w(:)/h(j+1, j)
	endif
        
        
    enddo
   
    !>>> первый вариант нахождения корня из матрицы:
    call SchurDecompositionQR(m_lanczos, H(1:m_lanczos, :), Qr, Lambda)
     do j =1, m_lanczos
        !write(*, *) j, Lambda(j, j), sqrt(Lambda(j, j))
        Lambda(j, j) = sqrt(Lambda(j, j))
     enddo
    Bm = matmul(Qr, matmul(Lambda,  transpose(Qr) ))

    !>> второй способ:
    ! свяли задачу к меньшой размерности, в которой ищем корень из матрицы:
    !call CholeskyDecomposition(m_lanczos, H(1:m_lanczos, :),  Bm)
    
    ! нужен только первый столбец:
    ! matmul(transpose(v(1:m_lanczos, 1:n)), Qr)
    Rx(:) = 0.0
    Ry(:) = 0.0
    Rz(:) = 0.0
    do i = 1, num_atom
	ix = 3*(i - 1) + 1
	iy = 3*(i - 1) + 2
	iz = 3*(i - 1) + 3

	do j = 1, m_lanczos
	    	Rx(i) = Rx(i) + v(j, ix)*Bm(1, j)
		Ry(i) = Ry(i) + v(j, iy)*Bm(1, j)
		Rz(i) = Rz(i) + v(j, iz)*Bm(1, j)
	enddo

	Rx(i) = Rx(i)*norm_z
	Ry(i) = Ry(i)*norm_z
	Rz(i) = Rz(i)*norm_z
    enddo
    ! нормируем:
    Rx(:) =Rx(:)*sqrt(2.0*dt)
    Ry(:) =Ry(:)*sqrt(2.0*dt)
    Rz(:) =Rz(:)*sqrt(2.0*dt)

    
    return
end subroutine KrylovSquareRoot
    
!>>
!  Расчет итогового перемещений
!>
subroutine HI_Displacement_Krylov(nstep, recalc_diffusion_tensor, num_atom, X, Y, Z, &
    & grand_x, grand_y, grand_z, Fx, Fy, Fz, dt, kT, friction, hradius, D, m_lanczos, &   
    & delta_coord_x, delta_coord_y, delta_coord_z)
    ! in
    integer(8) nstep
    integer(4) recalc_diffusion_tensor
    integer(4) num_atom
    real(8) X(num_atom), Y(num_atom), Z(num_atom)
    ! силы действующие на частицы:
    real(8) grand_x(num_atom), grand_y(num_atom), grand_z(num_atom)
    real(8) Fx(num_atom), Fy(num_atom), Fz(num_atom)
    ! шаг интегрирования:
    real(8) dt
    real(8) kT
    real(8) friction
    real(8) hradius
    real(8) D(3*num_atom, 3*num_atom)
    integer(4) m_lanczos
    ! out
    real(8) delta_coord_x(num_atom), delta_coord_y(num_atom), delta_coord_z(num_atom)
    ! work  
    real(8) Rx(num_atom), Ry(num_atom), Rz(num_atom)
    ! local
     
    ! если нужно пересчитываем тензор диффузии:
    call HI_TensorRPY_Recalc(nstep, recalc_diffusion_tensor, num_atom, X, Y, Z, & 
			     & dt, kT, friction, hradius, D)
    
    ! находим вектор случайных перемещений
    call KrylovSquareRoot(num_atom, dt, D, m_lanczos, grand_x, grand_y, grand_z, Rx, Ry, Rz)

    ! делаем итоговое перемещение:
    call HydrodynamicInteraction(num_atom, Rx, Ry, Rz, Fx, Fy, Fz, &
                            & dt, D, delta_coord_x, delta_coord_y, delta_coord_z)
    

    return
end subroutine HI_Displacement_Krylov

!>>
!  Расчет итогового перемещений
!>
subroutine HI_Displacement_Test(nstep, recalc_diffusion_tensor, num_atom, X, Y, Z, &
    & grand_x, grand_y, grand_z, Fx, Fy, Fz, dt, kT, friction, hradius, lambda, &   
    & delta_coord_x, delta_coord_y, delta_coord_z)
    ! in
    integer(8) nstep
    integer(4) recalc_diffusion_tensor
    integer(4) num_atom
    real(8) X(num_atom), Y(num_atom), Z(num_atom)
    ! силы действующие на частицы:
    real(8) grand_x(num_atom), grand_y(num_atom), grand_z(num_atom)
    real(8) Fx(num_atom), Fy(num_atom), Fz(num_atom)
    ! шаг интегрирования:
    real(8) dt
    real(8) kT
    real(8) friction
    real(8) hradius
    real(8) lambda
    ! out
    real(8) delta_coord_x(num_atom), delta_coord_y(num_atom), delta_coord_z(num_atom)
    ! local
    real(8) D(3*num_atom, 3*num_atom)
    integer(4) i

	! если нужно пересчитываем тензор диффузии:
    	call HI_TensorRPY_Recalc(nstep, recalc_diffusion_tensor, num_atom, X, Y, Z, & 
			     & dt, kT, friction, hradius, D)
	! но нужен он только для теста:
	call TensorIsPositiveDefine_Test(3*num_atom, D)
	do i=1, num_atom
    	    ! находим перемещение
 	    delta_coord_x(i) = Fx(i)*dt/friction + lambda*grand_x(i)
  	    delta_coord_y(i) = Fy(i)*dt/friction + lambda*grand_y(i)
   	    delta_coord_z(i) = Fz(i)*dt/friction + lambda*grand_z(i)
 	  enddo 

	return
end subroutine HI_Displacement_Test
#endif


end Module HydroDynamic