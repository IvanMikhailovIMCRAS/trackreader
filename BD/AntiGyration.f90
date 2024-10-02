Module AntiGyration
  Use CommonModule
  Use TestErrors
  
  implicit none

  
Contains

!>
!!	@brief	Убираем из системы вращение
!!
!!	@param	num_atom [IN]	        число атомов
!!	@param	X,Y,Z 	 [IN]	        координаты атомов
!!	@param	dX,dY,dZ [IN/OUT]	перемещения атомов
!!
subroutine AntiRotation(num_atom, X, Y, Z, dX, dY, dZ, &
                        & center_x, center_y, center_z)
	! in / out
	integer(4) num_atom
	real(8) X(num_atom), Y(num_atom), Z(num_atom)
	real(8) dX(num_atom), dY(num_atom), dZ(num_atom)
	real(8) center_x, center_y, center_z
	! work
	integer(4) i, j
	real(8) Iij(3,3),Li(3),W(3)
	real(8) WX,WY,WZ,sum_x,sum_y,sum_z,XI,YI,ZI

	DATA Iij,Li,W /9*0, 3*0, 3*0/

	sum_x = center_x
	sum_y = center_y
	sum_z = center_z

	do i=1, num_atom
	   sum_x = sum_x + X(i)
	   sum_y = sum_y + Y(i)
	   sum_z = sum_z + Z(i)
	end do

	
	do i=1, num_atom
	   XI = X(i)- sum_x/num_atom
	   YI = Y(i)- sum_y/num_atom
	   ZI = Z(i)- sum_z/num_atom

	   Iij(1,1) = Iij(1,1) + YI**2 + ZI**2
	   Iij(1,2) = Iij(1,2) - XI*YI
	   Iij(1,3) = Iij(1,3) - XI*ZI
	   Iij(2,2) = Iij(2,2) + XI**2 + ZI**2
	   Iij(2,3) = Iij(2,3) - YI*ZI
	   Iij(3,3) = Iij(3,3) + XI**2 + YI**2

	   Li(1) = Li(1) + YI*dZ(I) - ZI*dY(I)
	   Li(2) = Li(2) + ZI*dX(I) - XI*dZ(I)
	   Li(3) = Li(3) + XI*dY(I) - YI*dX(I)
        end do

	Iij(2,1) = Iij(1,2)
	Iij(3,1) = Iij(1,3)
	Iij(3,2) = Iij(2,3)

	call GAUSS(3, Iij, Li, W)

	do i=1, num_atom
	   XI = X(I) - sum_x/num_atom
	   YI = Y(I) - sum_y/num_atom
	   ZI = Z(I) - sum_z/num_atom

	   dX(I) = dX(I) + (YI*W(3) - ZI*W(2))
	   dY(I) = dY(I) + (ZI*W(1) - XI*W(3))
	   dZ(I) = dZ(I) + (XI*W(2) - YI*W(1))
	enddo
	
	return
end subroutine AntiRotation

!>
!!	@brief	Решение системы уравнений методом Гаусса
!!
!!	@param	N [IN] размерность
!!	@param	A [IN] симметричная матрица
!!	@param	Y [IN] свободные члены
!!	@param	X [OUT]	вектор неизвестных
!!
!***********************************************************************
subroutine GAUSS(N, A, Y, X)
    ! in
    integer(4) N
    real(8) A(N,N), Y(N)
    ! out
    real(8) X(N)
    ! work
    integer(4) i, j, k
    real(8) B(N,N+1), SUMY, MULT

    X(:) = 0.0
    SUMY = 0.0

    ! заполняем систему:
    do i = 1, N
        do j= 1, N
	   B(I,J)=A(I,J)
	      
           if (i.eq.j.and.B(i,j).eq.0.0) then
	      write(*,'(A,I3,A,A)') 'ERROR! - Aij=0, (i=j=',J,') - PRESS', &
                            &	  '<ENTER>'
     	      read(*,*)
	      
	      stop
           end if
        end do
	
	B(i, N+1) = Y(i)
	SUMY = SUMY + Y(i) ! dabs(Y(i))
    end do

    ! если все свободные члены равны нулю:
    if (SUMY.eq.0.0) then
        write(*,'(A,I3,A)') 'ERROR! - Yi=0,(i=1,',N,') - PRESS <ENTER>'
        read(*,*)
        stop
    end if

    ! зануляем элементы ниже главной диагонали:
    do i = 1, N
        do k = i+1, N
            MULT = B(K,I)/B(I,I)
            do j=1, N+1
	      B(K,J)=B(K,J)-MULT*B(I,J)
            end do
        end do
    end do
			
    ! зануляем строки выше главной диагонали и находим икс-ы
    do i=N,1,-1
        do j=1, N
	   B(i, N+1) = B(i, N+1) - B(i,j)*X(j)
        end do
	X(i) = B(i,N+1)/B(i,i)
    end do
	
    return
end subroutine GAUSS

end Module AntiGyration