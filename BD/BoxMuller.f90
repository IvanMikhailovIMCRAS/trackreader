!> @brief	Преобразование Бокса-Мюллера
!! @date		18.12.2014
!!
!! @author	Михайлов И.В., Шавыкин О.В.
Module BoxMuller  

Contains 
    !> @brief получаем вектор случайных сил для одного шага интегрирования
    !!
    !! @param 
    subroutine GetRandomForce(num_atom, grand_x, grand_y, grand_z) 
      Implicit None
      !INPUT
      integer(4) num_atom
      !OUTPUT
      real(8)   grand_x(num_atom), grand_y(num_atom), grand_z(num_atom)    
      ! WORK
      integer(4) i ! итерататор
      real(8) sum_delta_x, sum_delta_y, sum_delta_z
      
      ! случайная сила
	  sum_delta_x = 0.0
	  sum_delta_y = 0.0
	  sum_delta_z = 0.0
	  ! заготавливаем массив перемещений
	  do i=1, num_atom
	      call GaussRand(grand_x(i))
	      call GaussRand(grand_y(i))
	      call GaussRand(grand_z(i))
	    
	      ! находим сумму перемещений:
	      sum_delta_x = sum_delta_x + grand_x(i) 
	      sum_delta_y = sum_delta_y + grand_y(i)
	      sum_delta_z = sum_delta_z + grand_z(i)
	  enddo
	  
	  ! среднее значение:
	  sum_delta_x = sum_delta_x / num_atom
	  sum_delta_y = sum_delta_y / num_atom
	  sum_delta_z = sum_delta_z / num_atom
 	  
 	  ! убираем импульс:
	  do i=1,num_atom
	      grand_x(i) = grand_x(i) - sum_delta_x
	      grand_y(i) = grand_y(i) - sum_delta_y
	      grand_z(i) = grand_z(i) - sum_delta_z  
	  enddo
	  
    end subroutine GetRandomForce

    !>	@brief возвращает случайное число из Гауссового распределения. 
    !!		Было использовано преобразование Бокса-Мюллера
    !!	@param	gauss_rand	[OUT]	случайное нормальнораспределённое число
    subroutine GaussRand(gauss_rand) 
      Implicit None
      !OUTPUT
      real(8)   gauss_rand
      !WORK
      real(8)   PI
      parameter( PI=3.1415926535897934)

      save     compute,store ! делаем переменные статическими

      real(8)   phi, r  ! параметры преобразования

      real(8)   store
      logical  compute
      data     compute  / .true. /
      !BODY
      if (compute) then
        r = dble(randm())
        ! чтобы не получился NaN - 
        ! случайное число не должно быть очень маленьким
        do while (r < 0.01)
            r = dble(randm())
        enddo

        phi = dble(randm())

        ! преобразование Бокса-Мюллера (1 вариант)
        gauss_rand = cos(2.0*PI*phi)*sqrt(-2.0*log(r))

        store = sin(2.0*PI*phi)*sqrt(-2.0*log(r))

        compute = .false.
      else
         gauss_rand  = store
         compute  = .true.
      end if

      return
      end subroutine GaussRand
!*****************************************


!=======================================================================
!>	@brief Генератор случайных чисел равномеро распределённых на (0,1]
!!	@return	Случайное число.
!!	@author П.Г. Халатур
    real*4 function randm()
!======================================================================
! randmon number generator
!======================================================================
      save
      real*4 temp(2048)
      data   key /2049/, len /2048/

      key = key + 1
      if ( key .gt. len ) then
         if ( key .eq. 2050 ) then
            call ctime00 ( tx, ty )
            ij=1 + int( 31327. * ( ty-aint( ty ) ) )
            do i = 1, 10000
               call ctime00 ( tx, ty )
            end do
            kl = 1 + int(30080.*(10000.*(tx+ty)-aint(10000.*(tx+ty))))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     ij = 1
!     kl = 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           if ( ij .lt. 0  .or.  ij .gt. 31328 ) ij = 1
           if ( kl .lt. 0  .or.  kl .gt. 30081 ) kl = 1
           call rmarin ( ij, kl )
!          write(*,*)' initialization of the ranmar()'
!          write(*,*)' ij & kl: ',ij, kl
         end if
         call ranmar ( temp, len )
         key = 1
      end if
      randm = temp ( key )
      return
      end function randm
!=============================================================================

!=============================================================================
!>	@brief Создания инициализирующего массива для генератора случайных чисел
!!	@param	ij	[IN]	от 0 до 31328 
!!	@param	kl	[IN}	от 0 до 30081
!!	@author П.Г. Халатур
!!
    subroutine rmarin (ij, kl)
!=============================================================================
! this is the initialization routine for the randmom number generator ranmar()
! note: the seed variables can have values between:    0 <= ij <= 31328
!                                                      0 <= kl <= 30081
! the randmom number sequences created by these two seeds are of sufficient
! length to complete an entire calculation with. for example, if sveral
! different groups are working on different parts of the same calculation,
! each group could be assigned its own ij seed. this would leave each group
! with 30000 choices for the second seed. that is to say, this randmom
! number generator can create 900 million different subsequences -- with
! each subsequence having a length of approximately 10^30.
!
! use ij = 1802 & kl = 9373 to test the randmom number generator. the
! subroutine ranmar should be used to generate 20000 randmom numbers.
! then display the next six randmom numbers generated multiplied by 4096*4096
! if the randmom number generator is working properly, the randmom numbers
! should be:
!           6533892.0  14220222.0  7275067.0
!           6172232.0  8354498.0   10633180.0
!==========================================================================
      real u(97), c, cd, cm
      integer i97, j97
      logical test
      common /raset1/ u, c, cd, cm, i97, j97, test

      test=.false.


      i = mod(ij/177, 177) + 2
      j = mod(ij    , 177) + 2
      k = mod(kl/169, 178) + 1
      l = mod(kl,     169)

      do ii = 1, 97
         s = 0.0
         t = 0.5
         do jj = 1, 24
            m = mod(mod(i*j, 179)*k, 179)
            i = j
            j = k
            k = m
            l = mod(53*l+1, 169)
            if (mod(l*m, 64) .ge. 32) then
               s = s + t
            endif
            t = 0.5 * t
ENDDO
         u(ii) = s
ENDDO

      c = 362436.0 / 16777216.0
      cd = 7654321.0 / 16777216.0
      cm = 16777213.0 /16777216.0

      i97 = 97
      j97 = 33

      test = .true.
      return
      end
!==========================================================================

!==========================================================================
!>	@brief Конгруэнтный генератор случайных чисел предложенный G. Marsaglia
!!	@param	rvec	[OUT]	массив случайных чисел
!!	@param	len 	[IN]	длина выходного массива
!!	@author П.Г. Халатур
!!
      subroutine ranmar ( rvec, len )
!==========================================================================
! this is the randmom number generator proposed by george marsaglia in
! florida state university report: fsu-scri-87-50
! it was slightly modified by f. james to produce an array of pseudorandmom
! numbers.
!==========================================================================
      real rvec(*)
      real u(97), c, cd, cm
      integer i97, j97
      logical test
      common /raset1/ u, c, cd, cm, i97, j97, test
      integer ivec

      do ivec = 1, len
         uni = u(i97) - u(j97)
         if( uni .lt. 0.0 ) uni = uni + 1.0
         u(i97) = uni
         i97 = i97 - 1
         if(i97 .eq. 0) i97 = 97
         j97 = j97 - 1
         if(j97 .eq. 0) j97 = 97
         c = c - cd
         if( c .lt. 0.0 ) c = c + cm
         uni = uni - c
         if( uni .lt. 0.0 ) uni = uni + 1.0
         rvec(ivec) = uni
ENDDO
      return
      end
!============================================================================
!>	@brief возвращает число секунд и миллисекунд прошедших с полуночи
!!	@param	tx	[OUT]	число секунд
!!	@param	ty	[OUT]	число миллисекунд
!!	@author П.Г. Халатур
!!
subroutine ctime00 ( tx, ty )
!******************************************************************
!  returns the number of seconds and hundredths of seconds elapsed
!  since midnight
!******************************************************************
      save
      integer ih, im, is, ihu
      integer hms(3)
      data    key /0/

!     use the Unix-style "itime" and "idate" intrinsic functions,
!     this code works for all compilers except those noted below

!     call gettim (ih, im, is, ihu)
!     tx = ih*3600.0 + im*60 + is + ihu/100.0

!    call System_Clock (ih, im, is)
!    tx = real (ih) / real (im)

      call itime (hms)
      ih = hms(1)
      im = hms(2)
      is = hms(3)

      tx = (ih * 3600.0 + im * 60.0 + is) / 60.

      if (key .gt. 0) go to 10
      key = 1
      tr  = tx
      ihu = 0
   10 tx  = tx - tr
      ty  = tr
      return
      end
      
End Module BoxMuller
