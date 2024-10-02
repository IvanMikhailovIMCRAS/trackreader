!>	@brief	Функции для расчета энергии и сил частиц и системы в целом
!!
! Расскоментировать одну из следующих строк, если нужны граничные условия:
!#define BOUNDARY_CONDITION_CUBIC 
!#define BOUNDARY_CONDITION_CYLINDRICAL
Module ForceEnergy
  Use CommonModule
  Use TestErrors
  use Potential
  use PeriodicCondition ! граничные условия
  implicit none

  
Contains

!>>> составляем список сил, действующих на одну частицу
subroutine OneParticleInteractionList(num_atom, atom_type, num_potential, forces, num_onepart, onepart, onepart_p)
    ! input
    integer(4) num_atom
    integer(4) atom_type(num_atom)
    integer(4) num_potential
    type(force_info) forces(num_potential)
    ! output
    integer(4) num_onepart
    integer(4), dimension(:), allocatable :: onepart
    integer(4), dimension(:), allocatable :: onepart_p
    ! local
    integer(4) i, j
    integer(4) num_onepart_forces
    integer(4), dimension(:), allocatable :: mem_onepart_force
    integer(4), dimension(:), allocatable :: mem_onepart_type
    
    ! сначала сканируем потенциалы, узнаем какие из них одночастичные
    num_onepart_forces = 0
    do i = 1, num_potential
        if (forces(i)%num_type.eq.1) then
            num_onepart_forces = num_onepart_forces + 1
        endif
    enddo   
    
    ! если нет таких потенциалов, то извиняйте - завершаем работу процедуры!
    if (num_onepart_forces.eq.0) then
        num_onepart = 0
        return
    endif
    
    ! памяти под это:
    allocate(mem_onepart_force(1:num_onepart_forces))
    allocate(mem_onepart_type(1:num_onepart_forces))
    
    ! снова ищем и уже запоминаем необходимые индетифакаторы потенциалов и типов частиц:
    num_onepart_forces = 0
    do i = 1, num_potential
        if (forces(i)%num_type.eq.1) then
            num_onepart_forces = num_onepart_forces + 1
            ! запоминаем необхомые номера потенциалов и типов частиц:
            mem_onepart_force(num_onepart_forces) = i
            mem_onepart_type(num_onepart_forces) = forces(i)%types(1)
        endif
    enddo   
    
    ! сканируем доступные типы частиц и сравниваем со списком частиц на которые действуют одночастичные потенциалы
    num_onepart = 0
    do i = 1, num_atom
        do j = 1, num_onepart_forces
            ! наши совпадение, увеличиваем счетчик:
            if (mem_onepart_type(j).eq.atom_type(i)) then
                num_onepart = num_onepart + 1
            endif
        enddo
    enddo
    
    ! память:
    allocate(onepart(1:num_onepart))
    allocate(onepart_p(1:num_onepart))
    
    ! заполняем
    num_onepart = 0
    do i = 1, num_atom
        do j = 1, num_onepart_forces
            ! наши совпадение, увеличиваем счетчик:
            if (mem_onepart_type(j).eq.atom_type(i)) then
                num_onepart = num_onepart + 1
                    
                onepart(num_onepart) = i
                onepart_p(num_onepart) = mem_onepart_force(j)
            endif
        enddo
    enddo
   
    open(n_infor,file='INFOR', access='append')
    write(n_infor,'(a,1x,I0)') 'число постоянных силовых полей, &
                            & действующих на одну частицу:', num_onepart
    write(n_infor, *)
    close(n_infor)
   
    return
end subroutine OneParticleInteractionList

!>>> составляем список сил, действующих на одну частицу
subroutine ParticleParticleBondInteractionList(num_bond, bond1, bond2, num_atom, atom_type, &
                                & num_potential, forces, bond_p)
    ! input
    integer(4) num_bond
    integer(4) bond1(num_bond)
    integer(4) bond2(num_bond)
    integer(4) num_atom
    integer(4) atom_type(num_atom)
    integer(4) num_potential
    type(force_info) forces(num_potential)
    ! output
    integer(4), dimension(:), allocatable :: bond_p
    ! local
    integer(4) i, j
    logical(1) bfind
    
    allocate(bond_p(1:num_bond))
    
    ! для каждой связи ищем потенциал:
    do i = 1, num_bond
        bfind = .false.
        do j = 1, num_potential
            ! если это парный потенциал:
            if (forces(j)%num_type.eq.2) then
                ! и если потенциал без радиуса обрезки:
                if (forces(j)%bunbond.eqv..false.) then
                    ! то проверяем подходит ли этот потенциал для связи:
                    if (forces(j)%types(1).eq.atom_type(bond1(i)).and.forces(j)%types(2).eq.atom_type(bond2(i))) then
                        bond_p(i) = j
                        
                        bfind = .true.
                    endif
                    if (forces(j)%types(2).eq.atom_type(bond1(i)).and.forces(j)%types(1).eq.atom_type(bond2(i))) then
                        bond_p(i) = j
                        
                        bfind = .true.
                    endif
                endif
            endif
        enddo ! j
        
        if (bfind.eqv..false.) then
            call ERROR(ERR_CODE_NO_BOND_POTENTIAL)
        endif
    enddo
    
    
    return
end subroutine ParticleParticleBondInteractionList


!>>> составляем список сил, действующих на одну частицу
subroutine ParticleParticleAnglInteractionList(num_angles, angles_r, angles_c, angles_l, num_atom, atom_type, &
                                & num_potential, forces, angles_p)
    ! input
    integer(4) num_angles
    integer(4) angles_r(num_angles)
    integer(4) angles_c(num_angles)
    integer(4) angles_l(num_angles)
    integer(4) num_atom
    integer(4) atom_type(num_atom)
    integer(4) num_potential
    type(force_info) forces(num_potential)
    ! output
    integer(4), dimension(:), allocatable :: angles_p
    ! local
    integer(4) i, j
    logical(1) bfind
    
    allocate(angles_p(1:num_angles))
    
    ! для каждой связи ищем потенциал:
    do i = 1, num_angles
        bfind = .false.
        do j = 1, num_potential
            ! если это тройной потенциал:
            if (forces(j)%num_type.eq.3) then
                ! и если потенциал без радиуса обрезки (проверка на всякий случай):
                if (forces(j)%bunbond.eqv..false.) then
                    ! то проверяем подходит ли этот потенциал для связи:
                    if (forces(j)%types(1).eq.atom_type(angles_l(i)).and.forces(j)%types(2).eq.atom_type(angles_c(i)) &
                                       & .and.forces(j)%types(3).eq.atom_type(angles_r(i))) then
                        angles_p(i) = j
                        
                        bfind = .true.
                    endif
                    if (forces(j)%types(1).eq.atom_type(angles_r(i)).and.forces(j)%types(2).eq.atom_type(angles_c(i)) &
                                       & .and.forces(j)%types(3).eq.atom_type(angles_l(i))) then
                        angles_p(i) = j
                        
                        bfind = .true.
                    endif
                endif
            endif
        enddo ! j
        
        if (bfind.eqv..false.) then
            call ERROR(ERR_CODE_NO_ANGL_POTENTIAL)
        endif
    enddo
    
    
    return
end subroutine ParticleParticleAnglInteractionList

!>>
subroutine ParticleParticleUnBondInteractionMatrix(num_atom, atom_type, &
                                & num_potential, forces, num_type, P_aa_matrix)
    integer(4) num_atom
    integer(4) atom_type(num_atom)
    integer(4) num_potential
    type(force_info) forces(num_potential)
    integer(4) num_type
    ! output
    integer(4), dimension(:,:), allocatable :: P_aa_matrix
    ! local
    integer(4) i, j, p
    logical(1) bfind
    
    allocate(P_aa_matrix(1:num_type,1:num_type))
    
    ! для каждой связи ищем потенциал:
    do i = 1, num_type
        do j = i, num_type
            bfind = .false.
            do p = 1, num_potential
                ! если это двойной потенциал:
                if (forces(p)%num_type.eq.2) then
                    ! (проверка на всякий случай):
                    if (forces(p)%bunbond.eqv..true.) then
                    
                        ! то проверяем подходит ли этот потенциал для связи:
                        if (forces(p)%types(1).eq.i.and.forces(p)%types(2).eq.j) then
                            P_aa_matrix(i, j) = p
                            P_aa_matrix(j, i) = p
                        
                            bfind = .true.
                        endif
                        if (forces(p)%types(2).eq.i.and.forces(p)%types(1).eq.j) then
                            P_aa_matrix(i, j) = p
                            P_aa_matrix(j, i) = p
                        
                            bfind = .true.
                        endif
                    endif
                endif
            enddo ! p
            
            if (bfind.eqv..false.) then
                P_aa_matrix(i,j) = 0
                P_aa_matrix(j,i) = 0
            endif
        enddo ! j
    enddo ! i
    
    return
end subroutine ParticleParticleUnBondInteractionMatrix



!>>
subroutine SurfaceParticleUnBondInteractionMatrix(num_atom, atom_type, &
                                & num_potential, forces, num_surface, num_type, P_sa_matrix)
    integer(4) num_atom
    integer(4) atom_type(num_atom)
    integer(4) num_potential
    type(force_info) forces(num_potential)
    integer(4) num_surface
    integer(4) num_type
    ! output
    integer(4), dimension(:,:), allocatable :: P_sa_matrix
    ! local
    integer(4) i, j, p
    logical(1) bfind
    
    allocate(P_sa_matrix(1:num_surface,1:num_type))
    
    ! для каждой связи ищем потенциал:
    do i = 1, num_surface
        do j = 1, num_type
            bfind = .false.
            do p = 1, num_potential
                ! если это двойной потенциал:
                if (forces(p)%num_type.eq.2) then
                    ! (проверка на всякий случай):
                    if (forces(p)%bunbond.eqv..true.) then
                    
                        ! то проверяем подходит ли этот потенциал для связи:
                        if (forces(p)%types(1).eq.i.and.forces(p)%types(2).eq.j) then
                            P_sa_matrix(i, j) = p
                            bfind = .true.
                        endif
                    endif
                endif
            enddo ! p
            
            if (bfind.eqv..false.) then
                P_sa_matrix(i,j) = 0
            endif
        enddo ! j
    enddo ! i

    return
end subroutine SurfaceParticleUnBondInteractionMatrix


!>	@brief	Подсчет энергии объёмных взаимодействий
!!	@param	num_atom	[IN]	число атомов
!!	@param	X		[IN]	
!!	@param	Y		[IN]	координаты атомов
!!	@param	Z		[IN]
!!	@param	num_unbond	[IN]	число невалентных взаимодействий
!!	@param	unbond1		[IN]
!!	@param	unbond2		[IN]	невалентные взаимодействия
!!	@param	eps		[IN]	энергия невалентного взаимодействия
!!	@param	sig 		[IN]	диаметр атомов
!!	@param	rcut		[IN]	радиус обрезки
!!	@param	solvent		[IN]	качество расворителя
!!	@param	box		[IN]	размер коробки
!!	@param	U		[OUT]	энергия невалентного взаимодействия
subroutine Calculation_Lennard_Jones_Energy(num_atom, X, Y, Z, &
                                    & num_unbond, unbond1, unbond2, &
                                    & eps, sig, rcut,solvent, box, alpha_box, U)
      
	integer (4) num_atom,  num_unbond
	integer(4) unbond1(num_unbond), unbond2(num_unbond)
        real(8) X(num_atom), Y(num_atom), Z(num_atom)
        real(8)  eps, sig, rcut, U
        real(8) solvent, box, alpha_box
        ! WORK
        real(8) r, dx, dy, dz
        integer(4)  i
        real(8) energy
        
        U = 0.0
	do i=1,num_unbond
            ! получаем разницу между атомами:

            dx = X(unbond2(i)) - X(unbond1(i))
            dy = Y(unbond2(i)) - Y(unbond1(i))
            dz = Z(unbond2(i)) - Z(unbond1(i))
#ifdef BOUNDARY_CONDITION_CUBIC           
            ! граничные условия
            call BoundaryConditionCubic(dx, dy, dz, box)
#endif     
            r = sqrt(dx**2 + dy**2 + dz**2)
            
            call Lennard_Jones_Energy(r,energy,sig,rcut,eps,solvent)
            
            U = U + 2*energy
        end do
      
        return
end subroutine Calculation_Lennard_Jones_Energy


!> 	@brief Подсчёт энергии валентных взаимодействий
!!	@param	num_atom	[IN]	число атомов
!!	@param	X		[IN]	
!!	@param	Y		[IN]	координаты атомов
!!	@param	Z		[IN]
!!	@param	num_bond	[IN]	число связей
!!	@param	bond1		[IN]
!!	@param	bond2		[IN]	связи
!!	@param	kbond		[IN]	жесткость связи
!!	@param	lbond		[IN]	равновесная длина связи
!!	@param	box		[IN]	размер коробки
!!	@param	U		[OUT]	энергия валентного взаимодействия
subroutine Calculation_Harmonic_Energy(num_atom, X, Y, Z, &
                                & num_bond, bond1, bond2, &
                                & kbond, lbond, box, alpha_box, U)
      integer(4) num_bond, num_atom
      real(8) X(num_atom), Y(num_atom), Z(num_atom)
      integer(4) bond1(num_bond), bond2(num_bond)
      real(8) box, alpha_box, U, kbond, lbond
      
      ! локальные
      integer(4) i
      real(8) dx, dy, dz, r
      real(8) energy
      
      U = 0.0
      do i=1,num_bond
! получаем разницу между атомами:

            dx = X(bond2(i)) - X(bond1(i))
            dy = Y(bond2(i)) - Y(bond1(i))
            dz = Z(bond2(i)) - Z(bond1(i))
#ifdef BOUNDARY_CONDITION_CUBIC           
            ! граничные условия
            call BoundaryConditionCubic(dx, dy, dz, box)
#endif            
            r = sqrt(dx**2 + dy**2 + dz**2)
            
            call Harmonic_Energy(r, energy, kbond, lbond)
            U = U + 2*energy
     enddo

end subroutine Calculation_Harmonic_Energy



!> @brief	Подсчёт энергии угловых взаимодействий
!!	@param	num_atom	[IN]	число атомов
!!	@param	X		[IN]	
!!	@param	Y		[IN]	координаты атомов
!!	@param	Z		[IN]
!!	@param	num_bond	[IN]	число углов
!!	@param	angles_l	[IN]
!!	@param	angles_c	[IN]	углы
!!	@param	angles_r	[IN]
!!	@param	k_ang		[IN]	жесткость угла
!!	@param	av_ang		[IN]	равновесное значение угла
!!	@param	box		[IN]	размер коробки
!!	@param	U		[OUT]	энергия валентного взаимодействия
subroutine Calculation_Angle_Energy(num_atom, X, Y, Z, &
                                & num_angles, angles_r, angles_c, angles_l,&
                                & k_ang, av_ang, box, alpha_box, U)

      integer(4) num_angles, num_atom
      real(8) x(num_atom), y(num_atom), z(num_atom)
      integer(4) angles_l(num_atom), angles_c(num_atom), angles_r(num_atom)
      real(8) k_ang, av_ang, box, alpha_box, U
      
      ! локальные
      integer(4) i
      real(8) dx, dy, dz, r
      real(8) ang_val
      real(8) energy
                                
      U = 0.0
      do i=1, num_angles
! 	! расстояние между двумя точками:
! 	dx = x(angles_l(i)) - x(angles_r(i))
! 	dy = y(angles_l(i)) - y(angles_r(i))
! 	dz = z(angles_l(i)) - z(angles_r(i))
! 
! #ifdef BOUNDARY_CONDITION_CUBIC            
!         ! граничные условия
!         call BoundaryConditionCubic(dx, dy, dz, box)
! #endif            
!             
! 	r = sqrt(dx**2 + dy**2 + dz**2)
	
	! получаем угол
	call GetAngle(x(angles_l(i)), y(angles_l(i)), z(angles_l(i)), &
		  & x(angles_c(i)), y(angles_c(i)), z(angles_c(i)), &
		  & x(angles_l(i)), y(angles_l(i)), z(angles_l(i)), ang_val)
	
                                
                                
      call Angle_Energy(ang_val, energy, k_ang, av_ang)

      U = 2*energy
      
      end do
end subroutine Calculation_Angle_Energy      

!************************************************************

!>	@brief	Подсчет энергии объёмного взаимодействия
!!	@param	r		[IN]	расстояние между атомами
!!	@param	U		[OUT]	энергия невалентного взаимодействия
!!	@param	sig 		[IN]	диаметр атомов
!!	@param	rcut		[IN]	радиус обрезки
!!	@param	eps		[IN]	энергия невалентного взаимодействия
!!	@param	solvent		[IN]	качество расворителя
    subroutine Lennard_Jones_Energy(r,U,sig,rcut,eps, solvent)
!************************************************************
    implicit none
    real(8) r,U
    real(8) sig,rcut,eps
    real(8) rinf
    real(8) solvent

    rinf=0.5*sig
    !write(*,*) sig, r

    U=0.0

    ! если слишком маленькое расстояние
    if (r.le.rinf) call ERROR(ERR_CODE_DISTANCE)			  
    
    ! 0 < r <= sigma*2^(1/6)
    if (r.le.sig*2**(1/6)) U=4.0*eps*((sig/r)**12-(sig/r)**6) + eps - &
                            ! с этим слагаемым функция становится непрерывной в точке r=sigma*2^(1/6)
                            & ( solvent/2**(1/6) ) * (1 - (sig*2**(1/6)/rcut)**2)**2

    ! sigma*2^(1/6) < r <= rcut
    if (r.gt.sig*2**(1/6).and.r.le.rcut) U= -solvent*sig/r*(1.0 - (r/rcut)**2)**2


    return
    end subroutine Lennard_Jones_Energy
!**************************************************************

!**************************************************************
!>	@brief	Энергия связи
!!	@param	r	[IN]	расстояние между атомами
!!	@param	U	[OUT]	энергия связи
!!	@param	kbond	[IN]	жесткость связи
!!	@param	lbond	[IN]	равновесная длина связи
!!
    subroutine Harmonic_Energy(r, U, kbond, lbond)
!**************************************************************
        implicit none
        real(8) r, U, kbond, lbond
        

        U=0.0

        if (r.le.0.3*lbond.or.r.ge.1.8*lbond) then
            !write(n_infor,*) r, lbond, rinf
            call ERROR(ERR_CODE_LENGTH)
        end if
        U = 0.5 * kbond * (r - lbond)**2

        return
    end subroutine Harmonic_Energy
!**************************************************************

!
!>	@brief	Энергия углового взаимодействия
!!	@param	ang_val	[IN]	угол между векторами
!!	@param	U	[OUT]	энергия углового взаимодействия
!!	@param	k_ang	[IN]	жесткость угла
!!	@param	av_ang	[IN]	равновесный угол
subroutine Angle_Energy(ang_val, U, k_ang, av_ang)
      ! входные
      real(8) k_ang, av_ang
      
      real(8) ang_val, U
      
      U = 0.5 * k_ang*(ang_val - av_ang)**2
	    
      return
end subroutine Angle_Energy

!
!>	@brief	считаем кинетическую энергию системы
!!
!!	@param	nfree		[IN]	число степеней свободы
!!	@param	num_atom	[IN]	число атомов
!!	@param	dx		[IN]	
!!	@param	dy		[IN]	перемещения атомов
!!	@param	dz		[IN]
!!	@param	dt		[IN]	шаг интегрирования
!!	@param	ekin		[OUT]	кинетическая энергия
!!	
!!	@todo	разобраться со степенями свободы
!!	@todo	нужно ли делить на шаг интегрирования?
!!	@todo	не рассчитывается температура
subroutine KineticEnergy(nfree, num_atom, dx, dy, dz, dt, ekin)
	integer(4) nfree, num_bond, num_atom
	real(8) dx(num_atom), dy(num_atom), dz(num_atom)
	real(8) dt
	real(8) temp
	! вспомогательные:
	real(8) vv_sum, ekin
	real(8) vx, vy, vz
	integer(4) i
	
	! считаем средний квадрат скоростей
	vv_sum = 0
	do i = 1, num_atom
		vx = dx(i)  !/dt
		vy = dy(i) !/dt
		vz = dz(i) !/dt
		vv_sum = vv_sum + vx**2 + vy**2 + vz**2
	enddo

	! усреднённая кинетическая энергия системы после преобразования
	ekin   = 0.5 * vv_sum / num_atom
	! температура системы после преобразования
	!temp   = 2.0 * ekin / nfree

	return
end subroutine KineticEnergy




!>>>
subroutine GetRsphere(delta_r_sphere, num_potential, forces, num_type, P_aa_matrix, r_sphere2)
    !input
    real(8) delta_r_sphere
    integer(4) num_potential
    type(force_info) forces(num_potential)
    integer(4) num_type
    integer(4) P_aa_matrix(num_type, num_type)
    ! output
    real(8), dimension(:, :), allocatable :: r_sphere2
    !local
    integer(4) i, j
    real(8) r_cut
    
    
    allocate(r_sphere2(num_type,num_type))

    do i = 1, num_type
        do j = i, num_type
            if (P_aa_matrix(i, j).gt.0) then
                call PotentialRcut(forces(P_aa_matrix(i, j)), r_cut)
            else
                r_cut = 0.0
            endif
            
            r_sphere2(i, j) = (r_cut + delta_r_sphere)**2
            r_sphere2(j, i) = r_sphere2(i, j)
        enddo
    enddo

    return
end subroutine GetRsphere

!************************************************************
!> 	@brief Подсчет сил 
!!	@param	num_atom	[IN]	число атомов
!!	@param	X		[IN]	
!!	@param	Y		[IN]	координаты атомов
!!	@param	Z		[IN]
!!	@param	num_bond	[IN]	число связей
!!	@param	bond1		[IN]
!!	@param	bond2		[IN]	связи
!!	@param	bond_p		[IN]	потенциалы
!!	@param	num_unbond	[IN]	число невалентных взаимодействий
!!	@param	unbond1		[IN]
!!	@param	unbond2		[IN]	невалентные взаимодействия
!!	@param	P_aa_matrix	[IN]	потенциалы
!!	@param	num_angles	[IN]	число угловых взаимодействий
!!	@param	angles_l	[IN]
!!	@param	angles_c	[IN]	углы
!!	@param	angles_r	[IN]
!!	@param	angles_p	[IN]    потенциалы
!!	@param	box		[IN]	размер коробки
!!	@param	Fx		[OUT]
!!	@param	Fy		[OUT]	проекции сил действующих на частицы
!!	@param	Fz		[OUT]	
  subroutine Calculation_Force(num_atom, X, Y, Z, atom_type, num_type, &
                                & num_onepart, onepart, onepart_p, &
                                & num_bond, bond1, bond2, bond_p, &
                                & num_unbond, unbond1, unbond2, P_aa_matrix, &
                                & num_angles, angles_r, angles_c, angles_l, angles_p, &
                                & num_unbond_sa, unbond_sa1, unbond_sa2, P_sa_matrix, &
                                & num_potential, forces, &
                                & num_surface, surf, &
                                & box, &
                                & Fx, Fy, Fz)
!************************************************************
        ! входные
        integer(4) num_atom, num_bond, num_unbond, num_angles
        integer(4) atom_type(num_atom)
        integer(4) num_type
        real(8) X(num_atom), Y(num_atom), Z(num_atom)
        !>
        integer(4) num_onepart
        integer(4) onepart(num_onepart)
        integer(4) onepart_p(num_onepart)
        !>
        integer(4) bond1(num_bond), bond2(num_bond)
        integer(4) bond_p(num_bond)
        !>
        integer(4) unbond1(num_unbond), unbond2(num_unbond)
        integer(4) P_aa_matrix(num_type,num_type)
        !>
        integer(4) angles_r(num_angles), angles_c(num_angles), angles_l(num_angles)
        integer(4) angles_p(num_angles)
        !>
        integer(4) num_unbond_sa
        integer(4) unbond_sa1(num_unbond_sa), unbond_sa2(num_unbond_sa)
        integer(4) P_sa_matrix(num_surface, num_type)
        !>
        integer(4) num_potential
        type(force_info) forces(num_potential)
        !>
        integer(4) num_surface
        type(surface_info) surf(num_surface)
        !>
        real(8) box
        ! выходные
        real(8) Fx(num_atom), Fy(num_atom), Fz(num_atom)
        ! локальные:
        real(8) r, F
        real(8) dx, dy, dz
        !
        real(8) F_x, F_y, F_z
        integer(4)  i
        
        Fx(:) = 0.0
        Fy(:) = 0.0
        Fz(:) = 0.0
        
        ! 1. поля и потоки действующие на одну частицу:
        do i = 1, num_onepart
            call PotentialForceOneParticle(X(onepart(i)), Y(onepart(i)), Z(onepart(i)), & 
                                    & forces(onepart_p(i)), F_x, F_y, F_z)
                                    
            Fx(onepart(i)) = Fx(onepart(i)) + F_x
            Fy(onepart(i)) = Fy(onepart(i)) + F_y
            Fz(onepart(i)) = Fz(onepart(i)) + F_z
        enddo
        ! 2. парные валентные взаимодействия
        do i=1, num_bond
            call PotentialForceTwoParticle(X(bond1(i)), Y(bond1(i)), Z(bond1(i)), &
                                         & X(bond2(i)), Y(bond2(i)), Z(bond2(i)), forces(bond_p(i)), &
                                         & F_x, F_y, F_z)
                                         
            Fx(bond1(i)) = Fx(bond1(i)) + F_x
            Fy(bond1(i)) = Fy(bond1(i)) + F_y
            Fz(bond1(i)) = Fz(bond1(i)) + F_z
	  
	    Fx(bond2(i)) = Fx(bond2(i)) - F_x
            Fy(bond2(i)) = Fy(bond2(i)) - F_y
            Fz(bond2(i)) = Fz(bond2(i)) - F_z
        enddo
        
        ! 3. тройные валентные взаимодействия (углы)
        do i=1, num_angles
            call PotentialForceThreeParticle( x(angles_l(i)), y(angles_l(i)), z(angles_l(i)), &
                                         & x(angles_c(i)), y(angles_c(i)), z(angles_c(i)), &
                                         & x(angles_r(i)), y(angles_r(i)), z(angles_r(i)), &
                                         & forces(angles_p(i)), &
                                         & F_x, F_y, F_z)
        
            Fx(angles_l(i)) = Fx(angles_l(i)) + F_x
            Fy(angles_l(i)) = Fy(angles_l(i)) + F_y
            Fz(angles_l(i)) = Fz(angles_l(i)) + F_z
	  
	    Fx(angles_r(i)) = Fx(angles_r(i)) - F_x
            Fy(angles_r(i)) = Fy(angles_r(i)) - F_y
            Fz(angles_r(i)) = Fz(angles_r(i)) - F_z
        enddo
    
        ! 4. парные невалентные взаимодействия 
        do i=1, num_unbond
            call PotentialForceTwoParticle(X(unbond1(i)), Y(unbond1(i)), Z(unbond1(i)), &
                                         & X(unbond2(i)), Y(unbond2(i)), Z(unbond2(i)), &
                                         & forces(P_aa_matrix( atom_type(unbond1(i)),atom_type(unbond2(i)) )), &
                                         & F_x, F_y, F_z)
                                         
            Fx(unbond1(i)) = Fx(unbond1(i)) + F_x
            Fy(unbond1(i)) = Fy(unbond1(i)) + F_y
            Fz(unbond1(i)) = Fz(unbond1(i)) + F_z
	  
	    Fx(unbond2(i)) = Fx(unbond2(i)) - F_x
            Fy(unbond2(i)) = Fy(unbond2(i)) - F_y
            Fz(unbond2(i)) = Fz(unbond2(i)) - F_z
        enddo
          
        ! TODO: поверхности
        do i=1, num_unbond_sa
            call PotentialForceTwoParticle(X(unbond_sa1(i)), Y(unbond_sa1(i)), Z(unbond_sa1(i)), &
                                         & X(unbond_sa2(i)), Y(unbond_sa2(i)), Z(unbond_sa2(i)), &
                                         & forces(P_sa_matrix( atom_type(unbond_sa1(i)),atom_type(unbond_sa2(i)) )), &
                                         & F_x, F_y, F_z)
                  
	    Fx(unbond_sa2(i)) = Fx(unbond_sa2(i)) - F_x
            Fy(unbond_sa2(i)) = Fy(unbond_sa2(i)) - F_y
            Fz(unbond_sa2(i)) = Fz(unbond_sa2(i)) - F_z
        enddo
 	
        return
    end subroutine Calculation_Force
!************************************************************




!***********************************************************
!>	@brief	обновляет список несвязанных атомов
!!	@param	num_atom	[IN]	число атомов
!!	@param	x		[IN]	
!!	@param	y		[IN]	массивы координат атомов
!!	@param	z		[IN]
!!	@param	max_unbond_num	[IN]	максимальное число шариков в сфере
!!	@param	num_unbond	[OUT]	число невалетных взаимодействий
!!	@param	unbond1		[OUT]	
!!	@param	unbond2		[OUT]	массивы невалетных взаимодействий
!!	@param	box		[IN]	размер коробки
!!	@param	r_sphere	[IN]	радиус сферы
!!	@param	label_taboo	[IN]	метка "нет объёмных взаимодействий"
!!
    subroutine UnBond_List(num_atom,x,y,z,max_unbond_num, &
			  & num_unbond, unbond1, unbond2, &
			  & box, &
			  & atom_type, num_type, r_sphere2,label_taboo)
!************************************************************
    ! входные
    integer(4) num_atom
    integer(4) max_unbond_num, num_unbond
    real(8) x(num_atom),y(num_atom),z(num_atom)
    integer(4) unbond1(max_unbond_num), unbond2(max_unbond_num)
    logical label_taboo(num_atom,num_atom)
    real(8)	box
    integer(4) atom_type(num_atom)
    integer(4) num_type
    real(8) r_sphere2(num_type,num_type)
    
    ! рабочие
    real(8) distance
    integer(4) i,j, k
    real(8) delta_x, delta_y, delta_z
    
!     n_atom_sphere(:)=0
!     contact_unbond(:,:)=0

    num_unbond = 0
    
    ! просматриваем попарно атомы:
    do i=1,num_atom-1
        do j=i+1,num_atom
            if (label_taboo(i,j).eqv..false.) then
		delta_x = x(i)-x(j)
		delta_y = y(i)-y(j)
		delta_z = z(i)-z(j)
		! проверка граничных условий
#ifdef BOUNDARY_CONDITION_CUBIC           
                call BoundaryConditionCubic(delta_x, delta_y, delta_z, box)
#endif                
                ! находим расстояние между атомами:
                distance = delta_x**2 + delta_y**2 +delta_z**2
                
                if (distance.lt.r_sphere2(atom_type(i),atom_type(j))) then		    
                   ! если попали сюда, то значит добавляем невалентное взаимодействие:
		   num_unbond = num_unbond + 1
		      
		   unbond1(num_unbond) = j
		   unbond2(num_unbond) = i
               
                endif
            endif
        enddo
    enddo


    return
    end subroutine UnBond_List
!************************************************************

end Module ForceEnergy