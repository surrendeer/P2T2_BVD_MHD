!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModSetInitialCondition

  use BATL_lib, ONLY: &
       test_start, test_stop, iBlockTest, iTest, jTest, kTest

  use ModVarIndexes, ONLY: nVar
  use ModMain, ONLY: NamePrimitive_V
  use ModPhysics, ONLY: &
       UseShocktube, ShockLeftDim_V, ShockRightDim_V, ShockLeft_V, &
       ShockRight_V, ShockPosition, ShockSlope, nVectorVar, iVectorVar_I
  use ModBatsrusUtility, ONLY: stop_mpi, get_ivar
  use ModNumConst, ONLY: cTwoPi,cPi, cDegToRad

  implicit none

  private ! except

  public:: read_initial_cond_param ! read parameters for initial condition
  public:: set_initial_condition   ! set initial condition for one block
  public:: add_rotational_velocity ! transform between rotating/inertial frames

  ! Local variables

  ! Entropy constant for isentropic initial condition. Only used if positive.
  real :: EntropyConstant = -1.0

  ! Use 2D initial state algorithm
  logical:: UseInitialStateDefinition = .false.

  ! Wave/Tophat/Gaussian
  logical:: UseWave = .false.
  logical:: DoRotateWave = .true.

  ! Wave parameters
  real:: Width, Amplitude, Phase, LambdaX, LambdaY, LambdaZ
  real, dimension(nVar):: Width_V=0.0, Ampl_V=0.0, Phase_V=0.0, &
       x_V=0.0, y_V=0.0, z_V=0.0, KxWave_V=0.0, KyWave_V=0.0, KzWave_V=0.0
  integer:: iPower_V(nVar)=1


   !dzy modifed 2024/2/29    :)
  logical::testAdvectVortex =.false.!.true.!
 
  logical::testDecayAlfven  =.false.!.true.!

  logical::test2DAcousticWave  =.false.!.true.!  for Eular Equation

  logical::test30DegreesAlfven=.false.!.true.!  Old pi/6
  logical::testsmooothAlfven=.false.! .true.!   new pi/4
  
  logical::testRotor=.false.!.true.!

  logical::testCloud=.false.!.true.!.true.!

  logical::testBlast=.false.!.true.!

  logical::testTorsionalAlfvenPulse=.false.  !negative pressure :(
  
  logical::testRayleighTaylorHD=.false. !.true.!

  logical::test2DRiemann =.false.!.true.! ! for Eular Equation

  logical::testSphericalExplosion=.false.!.true.!
  logical::testOrszarg=.false.!.true.!
contains
  !============================================================================
  subroutine read_initial_cond_param(NameCommand)

    use ModReadParam, ONLY: read_var
    use ModUtilities, ONLY: join_string
    use ModInitialState, ONLY: init_initial_state, read_initial_state_param

    character(len=*), intent(in):: NameCommand

    character(len=500):: StringVar
    character(len=20):: NameVar

    integer:: iVar

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'read_initial_cond_param'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    select case(NameCommand)
    case("#STATEDEFINITION")
       UseInitialStateDefinition = .true.
       call join_string(nVar, NamePrimitive_V(1:nVar), StringVar)
       call init_initial_state(StringVar)
       call read_initial_state_param(NameCommand)

    case("#STATEINTERFACE")
       call read_initial_state_param(NameCommand)

    case("#UNIFORMSTATE")
       UseShockTube = .true.
       do iVar = 1, nVar
          call read_var(NamePrimitive_V(iVar), ShockLeftDim_V(iVar))
       end do
       ShockRightDim_V = ShockLeftDim_V

    case("#SHOCKTUBE")
       UseShockTube = .true.
       do iVar = 1, nVar
          call read_var(NamePrimitive_V(iVar)//' left', ShockLeftDim_V(iVar))
       end do
       do iVar = 1, nVar
          call read_var(NamePrimitive_V(iVar)//' right', &
               ShockRightDim_V(iVar))
       end do

    case("#SHOCKPOSITION")
       call read_var('ShockPosition', ShockPosition)
       call read_var('ShockSlope', ShockSlope)

    case("#ENTROPYCONSTANT")
       call read_var('EntropyConstant', EntropyConstant)

    case("#WAVE", "#WAVE2", "#WAVE4", "#WAVE6")
       UseWave = .true.
       call read_var('NameVar', NameVar)
       call get_ivar(NameVar, iVar)
       call read_var('Width', Width)
       call read_var('Amplitude', Amplitude)
       call read_var('LambdaX', LambdaX)
       call read_var('LambdaY', LambdaY)
       call read_var('LambdaZ', LambdaZ)
       call read_var('Phase', Phase)
       Width_V(iVar) = Width
       Ampl_V(iVar)  = Amplitude
       Phase_V(iVar) = Phase*cDegToRad

       if(NameCommand == '#WAVE6')then
          iPower_V(iVar) = 6
       elseif(NameCommand == '#WAVE4')then
          iPower_V(iVar) = 4
       elseif(NameCommand == '#WAVE2')then
          iPower_V(iVar) = 2
       else
          iPower_V(iVar) = 1
       end if

       ! if wavelength is smaller than 0, then the wave number is set to 0
       KxWave_V(iVar) = max(0.0, cTwoPi/LambdaX)
       KyWave_V(iVar) = max(0.0, cTwoPi/LambdaY)
       KzWave_V(iVar) = max(0.0, cTwoPi/LambdaZ)

    case("#BUMP")
       UseWave = .true.
       call read_var('NameVar',   NameVar)
       call get_ivar(NameVar, iVar)
       call read_var('Amplitude', Ampl_V(iVar))
       call read_var('WidthX',    LambdaX)
       call read_var('WidthY',    LambdaY)
       call read_var('WidthZ',    LambdaZ)
       call read_var('CenterX',   x_V(iVar))
       call read_var('CenterY',   y_V(iVar))
       call read_var('CenterZ',   z_V(iVar))
       call read_var('nPower',    iPower_V(iVar))
       iPower_V(iVar) = -iPower_V(iVar)

       ! Negative width sets 0 for 1/width (unlimited in that direction)
       KxWave_V(iVar) = max(0.0, 1/LambdaX)
       KyWave_V(iVar) = max(0.0, 1/LambdaY)
       KzWave_V(iVar) = max(0.0, 1/LambdaZ)

    case default
       call stop_mpi(NameSub//': unknown command='//NameCommand)
    end select

    call test_stop(NameSub, DoTest)
  end subroutine read_initial_cond_param
  !============================================================================
  subroutine set_initial_condition(iBlock)

    use ModMain
    use ModAdvance
    use ModB0, ONLY: B0_DGB, set_b0_cell, subtract_b0
    use ModGeometry, ONLY: Used_GB
    use ModIO, ONLY: IsRestart
    use ModPhysics, ONLY: FaceState_VI, CellState_VI, ShockSlope, &
         UseShockTube, ShockPosition, iUnitPrim_V, Io2No_V, Gamma_I
    use ModUserInterface ! user_set_ics
    use ModSaMhd,          ONLY: UseSaMhd, init_samhd
    use ModConstrainDivB, ONLY: constrain_ics
    use ModMultiFluid
    use ModRestartFile, ONLY: UseRestartWithFullB
    use ModBoundaryGeometry, ONLY: iBoundary_GB
    use ModInitialState, ONLY: get_initial_state
    use ModIonElectron,   ONLY: &
         correct_electronfluid_efield , DoCorrectElectronFluid, DoCorrectEfield
    use BATL_lib, ONLY: Xyz_DGB, IsPeriodic_D

    integer, intent(in) :: iBlock

    real   :: SinSlope, CosSlope, Rot_II(2,2), x, y
    integer:: i, j, k, iVar, iBoundary, iFluid, iGang
    real:: cPhi,eta_x,eta_y
    real::cos_alpha,sin_alpha,cbeta,v_parallel,v_perpendicular,B_parallel,B_perpendicular
    real::radius ,func_radius

    logical:: DoTestCell
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'set_initial_condition'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)

    iGang = 1
#ifdef _OPENACC
    iGang = iBlock
#endif

    DtMax_CB(:,:,:,iBlock) = 0.0

    Flux_VXI(:,:,:,:,iGang) = 0.0
    Flux_VYI(:,:,:,:,iGang) = 0.0
    Flux_VZI(:,:,:,:,iGang) = 0.0

    if(Unused_B(iBlock))then
       do iVar = 1, nVar
          State_VGB(iVar,:,:,:,iBlock) = DefaultState_V(iVar)
       end do
    else
       ! If used, initialize solution variables and parameters.
       if(UseB0) call set_b0_cell(iBlock)

       ! Subtract B0 from Full B0+B1 from restart to obtain B1
       if(UseB0 .and. IsRestart .and. UseRestartWithFullB) &
            call subtract_b0(iBlock)

       if(.not.IsRestart)then

          if(UseShockTube)then
             ! Calculate sin and cos from the tangent = ShockSlope
             CosSlope = 1/sqrt(1 + ShockSlope**2)
             SinSlope = ShockSlope*CosSlope

             ! Set rotational matrix
             Rot_II = reshape([CosSlope, SinSlope, -SinSlope, CosSlope],[2,2])

             if(ShockSlope /= 0.0 .and. UseWave .and. DoRotateWave) &
                  call rotate_wave

          end if  ! UseShockTube

          ! Loop through all the cells
          do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI

             DoTestCell = DoTest .and. i==iTest .and. j==jTest .and. k==kTest
             ! Default state
             State_VGB(:,i,j,k,iBlock) = CellState_VI(:,Coord1MaxBc_)

             if(DoTestCell)write(*,*) NameSub,': default state=', &
                  State_VGB(:,i,j,k,iBlock)

             if(.not.Used_GB(i,j,k,iBlock))then
                ! Cells outside the domain
                iBoundary = iBoundary_GB(i,j,k,iBlock)

                State_VGB(1:nVar,i,j,k,iBlock) = FaceState_VI(1:nVar,iBoundary)

                if(DoTestCell)write(*,*) NameSub,': face state=', &
                     State_VGB(:,i,j,k,iBlock)

                ! Convert velocity to momentum
                do iFluid = 1, nFluid
                   if(nFluid > 1) call select_fluid(iFluid)
                   State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock) = &
                        FaceState_VI(iUx:iUz,iBoundary) &
                        *FaceState_VI(iRho,iBoundary)
                end do

                if(DoTestCell)write(*,*) NameSub,': face cons.=', &
                     State_VGB(:,i,j,k,iBlock)

             elseif(UseInitialStateDefinition)then
                ! Cells defined by the #STATEDEFINITION command
                x = Xyz_DGB(x_,i,j,k,iBlock)
                y = Xyz_DGB(y_,i,j,k,iBlock)
                call get_initial_state( [x, y], State_VGB(:,i,j,k,iBlock) )

            !==================================================================================================
            !------------------the initialCondition of MHD test examples----------
                !!!!!           State_VGB(rho,rho*Ux,rho*Uy,rho*Uz,Bx,By,Bz,P)

            !==================================================================================================

                !testAdvectVortex=.true.
             elseif(testAdvectVortex)then
               x = Xyz_DGB(x_,i,j,k,iBlock)
               y = Xyz_DGB(y_,i,j,k,iBlock)                         
                        !!case("#AdvectVortex")   
                     State_VGB(1,i,j,k,iBlock)=1.0
                        !set Veloity
                     State_VGB(2,i,j,k,iBlock)=1.0-y*((1.0/cTwoPi)*(exp(0.5*(1.0-x**2-y**2))))
                                                              
                     State_VGB(3,i,j,k,iBlock)=1.0+x*((1.0/cTwoPi)*(exp(0.5*(1.0-x**2-y**2))))

                     State_VGB(4,i,j,k,iBlock)=0.0
                        !set Magnetic field

                     State_VGB(5,i,j,k,iBlock)=0.0-y*((1.0/cTwoPi)*(exp(0.5*(1.0-x**2-y**2))))
                                                
                     State_VGB(6,i,j,k,iBlock)=0.0+x*((1.0/cTwoPi)*(exp(0.5*(1.0-x**2-y**2))))

                     State_VGB(7,i,j,k,iBlock)= 0.0      
                        !set pressure
                     State_VGB(8,i,j,k,iBlock)=1.0-(((x**2+y**2)/(8.0*cPi**2))*exp(1.0-x**2-y**2))
                                                     
                     !write(*,*)"AdvectionVortex 'bump' has initialized"
             !--------------------------------------------------------------------------------------
             elseif(testDecayAlfven)then     !(BALSARA,2003)
               x = Xyz_DGB(x_,i,j,k,iBlock)
               y = Xyz_DGB(y_,i,j,k,iBlock) 

               eta_x=1.0/(1.0+x**2+y**2)**0.5         !eq(7.1)
               eta_y=eta_x*(x**2+y**2)**0.5

               cPhi=(x*cTwoPi/(x**2+y**2)**0.5) +cTwoPi*y      !eq(7.2)                
                        
                     State_VGB(1,i,j,k,iBlock)=1.0

                     State_VGB(2,i,j,k,iBlock)= -0.2*eta_y*cos(cPhi)      !eq(7.3)
                     State_VGB(3,i,j,k,iBlock)= +0.2*eta_x*cos(cPhi) 
                     State_VGB(4,i,j,k,iBlock)= +0.2*sin(cPhi)

                     State_VGB(5,i,j,k,iBlock)= 1.0*eta_x+0.2*eta_y*2*sqrt(cPI)*cos(cPhi)   !eq(7.4)
                     State_VGB(6,i,j,k,iBlock)= 1.0*eta_y-0.2*eta_x*2*sqrt(cPI)*cos(cPhi)
                     State_VGB(7,i,j,k,iBlock)= -0.2*2*sqrt(cPI)*sin(cPhi)

                     State_VGB(8,i,j,k,iBlock)= 1.0

                 
            elseif(test2DAcousticWave)then !Eular equation used!
               x = Xyz_DGB(x_,i,j,k,iBlock)
               y = Xyz_DGB(y_,i,j,k,iBlock) 
                  State_VGB(1,i,j,k,iBlock)=1.0
                  State_VGB(2,i,j,k,iBlock)=0.0
                  State_VGB(3,i,j,k,iBlock)=0.0
                  State_VGB(4,i,j,k,iBlock)=0.0
                     if(sqrt((x-0.5)**2+(y-0.5)**2) >0.3)then
                        State_VGB(5,i,j,k,iBlock)=0.6
                     else
                        State_VGB(5,i,j,k,iBlock)=&
                                 0.6+0.1*exp(-((sqrt((x-0.5)**2+(y-0.5)**2)/0.15)**2 ))*&
                                 cos(0.25*cPi*(sqrt((x-0.5)**2+(y-0.5)**2)/0.15))**6
                     endif

            elseif(test30DegreesAlfven)then !pi/6 is used
               x = Xyz_DGB(x_,i,j,k,iBlock)
               y = Xyz_DGB(y_,i,j,k,iBlock)
               cos_alpha=cos(cPi/6.0)
               sin_alpha=sin(cPi/6.0)

                State_VGB(1,i,j,k,iBlock)=1.0

                State_VGB(2,i,j,k,iBlock)=&
                           -0.5*(0.1*sin(cTwoPi*(x*cos_alpha+y*sin_alpha))  )
                State_VGB(3,i,j,k,iBlock)=&
                           -1.0*sqrt(3.0)*State_VGB(2,i,j,k,iBlock)
                State_VGB(4,i,j,k,iBlock)=&
                           0.1*cos(cTwoPi*(x*cos_alpha+y*sin_alpha))

                State_VGB(5,i,j,k,iBlock)=&
                           sqrt(3.0)/2.0-0.5*(0.1*sin(cTwoPi*(x*cos_alpha+y*sin_alpha)))
                State_VGB(6,i,j,k,iBlock)=&
                           2.0-sqrt(3.0)*State_VGB(5,i,j,k,iBlock)
                State_VGB(7,i,j,k,iBlock)=&
                           0.1*cos(cTwoPi*(x*cos_alpha+y*sin_alpha))

                State_VGB(8,i,j,k,iBlock)=0.1

                !Liu, M., Zhang, M., Li, C., Shen, F., 2021. A new locally divergence-free wls-eno scheme based on 
                !the positivity-preserving finite volume method for ideal mhd equations
            elseif(testsmooothAlfven)then !pi/4 is used
                  x = Xyz_DGB(x_,i,j,k,iBlock)
                  y = Xyz_DGB(y_,i,j,k,iBlock)
                  cos_alpha=cos(cPi/4.0)
                  sin_alpha=sin(cPi/4.0)

                  cbeta=(x*cos_alpha+y*sin_alpha)
                  v_parallel=0.0
                  v_perpendicular=0.1*sin(cTwoPi*cbeta)
                  B_parallel=1.0
                  B_perpendicular=v_perpendicular

   
                   State_VGB(1,i,j,k,iBlock)=1.0
   
                   State_VGB(2,i,j,k,iBlock)=&
                                             cos_alpha*v_parallel-sin_alpha*v_perpendicular
                   State_VGB(3,i,j,k,iBlock)=&
                                             cos_alpha*v_perpendicular+sin_alpha*v_parallel
                   State_VGB(4,i,j,k,iBlock)=&
                                             0.1*cos(cTwoPi*cbeta)
   
                   State_VGB(5,i,j,k,iBlock)=&
                                             cos_alpha*B_parallel-sin_alpha*B_perpendicular
                              
                   State_VGB(6,i,j,k,iBlock)=&
                                             cos_alpha*B_perpendicular+sin_alpha*B_parallel
                             
                   State_VGB(7,i,j,k,iBlock)=0.1*cos(cTwoPi*cbeta)
                        
                   State_VGB(8,i,j,k,iBlock)=0.1

            elseif(testRotor)then 
               x = Xyz_DGB(x_,i,j,k,iBlock)
               y = Xyz_DGB(y_,i,j,k,iBlock) 

                  radius = sqrt((x-0.5)**2+(y-0.5)**2)
                  func_radius=(0.115-radius)/0.015
                  if (radius<0.1)then

                     State_VGB(1,i,j,k,iBlock)= 10.0
                     State_VGB(8,i,j,k,iBlock)= 0.5
                     State_VGB(5,i,j,k,iBlock)= 2.50/(sqrt(4*cPi))
                     State_VGB(6,i,j,k,iBlock)= 0.0
                     State_VGB(7,i,j,k,iBlock)= 0.0
     
                     State_VGB(2,i,j,k,iBlock)= (-(y-0.5)/0.1)*State_VGB(1,i,j,k,iBlock)
                     State_VGB(3,i,j,k,iBlock)= ( (x-0.5)/0.1)*State_VGB(1,i,j,k,iBlock)
                     State_VGB(4,i,j,k,iBlock)=    0.0

                  elseif(radius.ge.0.1 .and. radius < 0.115)then

                     State_VGB(1,i,j,k,iBlock)= 1.0+9.0*func_radius
                     State_VGB(8,i,j,k,iBlock)= 0.5
                     State_VGB(5,i,j,k,iBlock)= 2.50/(sqrt(4*cPi))
                     State_VGB(6,i,j,k,iBlock)= 0.0
                     State_VGB(7,i,j,k,iBlock)= 0.0
     
                     State_VGB(2,i,j,k,iBlock)= -1.0*(y-0.5)*func_radius/radius*State_VGB(1,i,j,k,iBlock)
                     State_VGB(3,i,j,k,iBlock)=  1.0*(x-0.5)*func_radius/radius*State_VGB(1,i,j,k,iBlock)
                     State_VGB(4,i,j,k,iBlock)=  0.0
                  elseif(radius .ge. 0.115)then
                     State_VGB(1,i,j,k,iBlock)= 1.0
                     State_VGB(8,i,j,k,iBlock)= 0.5
                     State_VGB(5,i,j,k,iBlock)= 2.50/(sqrt(4*cPi))
                     State_VGB(6,i,j,k,iBlock)= 0.0
                     State_VGB(7,i,j,k,iBlock)= 0.0
     
                     State_VGB(2,i,j,k,iBlock)=  0.0
                     State_VGB(3,i,j,k,iBlock)=  0.0
                     State_VGB(4,i,j,k,iBlock)=  0.0
                  endif

            elseif(testCloud)then 
               x = Xyz_DGB(x_,i,j,k,iBlock)
               y = Xyz_DGB(y_,i,j,k,iBlock)  
               radius = sqrt((x-0.25)**2+(y-0.5)**2)  
               
                  if(x<0.05)then!leftstate
                     State_VGB(1,i,j,k,iBlock)= 3.86859
                     State_VGB(2,i,j,k,iBlock)= 11.2536*3.86859!rho*ux
                     State_VGB(3,i,j,k,iBlock)= 0.0
                     State_VGB(4,i,j,k,iBlock)= 0.0      
                     State_VGB(5,i,j,k,iBlock)= 0.0
                     State_VGB(6,i,j,k,iBlock)= 2.1826182
                     State_VGB(7,i,j,k,iBlock)= -2.1826182
                     State_VGB(8,i,j,k,iBlock)= 167.345

                  elseif(x>0.05 .and. radius>0.15)then   !rightstate
                     State_VGB(1,i,j,k,iBlock)= 1.0

                     State_VGB(2,i,j,k,iBlock)= 0.0
                     State_VGB(3,i,j,k,iBlock)=  0.0
                     State_VGB(4,i,j,k,iBlock)=  0.0
         
                     State_VGB(5,i,j,k,iBlock)= 0.0
                     State_VGB(6,i,j,k,iBlock)= 0.56418958
                     State_VGB(7,i,j,k,iBlock)= 0.56418958
                     State_VGB(8,i,j,k,iBlock)= 1.0

                  else
                     State_VGB(1,i,j,k,iBlock)= 10.0

                     State_VGB(2,i,j,k,iBlock)= 0.0
                     State_VGB(3,i,j,k,iBlock)=  0.0
                     State_VGB(4,i,j,k,iBlock)=  0.0
         
                     State_VGB(5,i,j,k,iBlock)= 0.0
                     State_VGB(6,i,j,k,iBlock)= 0.56418958
                     State_VGB(7,i,j,k,iBlock)= 0.56418958
                     State_VGB(8,i,j,k,iBlock)= 1.0

                  
                  endif
            elseif(testBlast)then
               x = Xyz_DGB(x_,i,j,k,iBlock)
               y = Xyz_DGB(y_,i,j,k,iBlock)  
               radius = sqrt((x-0.5)**2+(y-0.5)**2) 
               if(radius>0.1-1.d-15)then
                  State_VGB(1,i,j,k,iBlock)= 1.0
                  State_VGB(8,i,j,k,iBlock)= 0.1
                  State_VGB(5,i,j,k,iBlock)= 100.0/(sqrt(4*cPi))!sqrt(2.0)/2.0   !
                  State_VGB(6,i,j,k,iBlock)= 0.0!sqrt(2.0)/2.0   !100.0/(sqrt(4*cPi))!0.0     ! 
                  State_VGB(7,i,j,k,iBlock)= 0.0

                  State_VGB(2,i,j,k,iBlock)=  0.0
                  State_VGB(3,i,j,k,iBlock)=  0.0
                  State_VGB(4,i,j,k,iBlock)=  0.0
               else 
                  State_VGB(1,i,j,k,iBlock)= 1.0
                  State_VGB(8,i,j,k,iBlock)= 1000.0!10.0
                  State_VGB(5,i,j,k,iBlock)= 100.0/(sqrt(4*cPi))!sqrt(2.0)/2.0        !100.0/(sqrt(8*cPi))
                  State_VGB(6,i,j,k,iBlock)= 0.0!sqrt(2.0)/2.0        !100.0/(sqrt(8*cPi))!0.0         !100.0/(sqrt(8*cPi))
                  State_VGB(7,i,j,k,iBlock)= 0.0
   
                  State_VGB(2,i,j,k,iBlock)=  0.0
                  State_VGB(3,i,j,k,iBlock)=  0.0
                  State_VGB(4,i,j,k,iBlock)=  0.0
               endif


            
            elseif(testTorsionalAlfvenPulse)then
               x = Xyz_DGB(x_,i,j,k,iBlock)
               y = Xyz_DGB(y_,i,j,k,iBlock)

               cPhi=(cPi/8.0)*(tanh((x+0.25)/0.005)+1.0)*(tanh((0.25-x)/0.005)+1.0)

               State_VGB(1,i,j,k,iBlock)= 1.0

               State_VGB(2,i,j,k,iBlock)= 10.0
               State_VGB(3,i,j,k,iBlock)=  10.0*cos(cPhi)
               State_VGB(4,i,j,k,iBlock)=  10.0*sin(cPhi)
   
               State_VGB(5,i,j,k,iBlock)= 10.0
               State_VGB(6,i,j,k,iBlock)= -20.0*sqrt(cPi)*cos(cPhi)
               State_VGB(7,i,j,k,iBlock)= -20.0*sqrt(cPi)*sin(cPhi)
               State_VGB(8,i,j,k,iBlock)= 0.01

            elseif(testRayleighTaylorHD)then
                  x = Xyz_DGB(x_,i,j,k,iBlock)
                  y = Xyz_DGB(y_,i,j,k,iBlock)


                  if(y<0.5)then
                     State_VGB(1,i,j,k,iBlock)= 2.0
   
                     State_VGB(2,i,j,k,iBlock)= 0.0


                   
                     State_VGB(4,i,j,k,iBlock)=  0.0
         
                     State_VGB(5,i,j,k,iBlock)= 0.0
                     State_VGB(6,i,j,k,iBlock)= 0.0
                     State_VGB(7,i,j,k,iBlock)= 0.0
   
                     State_VGB(8,i,j,k,iBlock)= 2.0*y+1
                     if(x<0.125)then
                        State_VGB(3,i,j,k,iBlock)=-0.025*cos(8*cPi*x)*2.0*sqrt((5.0/3.0)*(2.0*y+1)/2.0)
                     else
                        State_VGB(3,i,j,k,iBlock)=-0.025*cos(0.25-8*cPi*x)*2.0*sqrt((5.0/3.0)*(2.0*y+1)/2.0)     
                     endif
                  else
                     State_VGB(1,i,j,k,iBlock)= 1.0
   
                     State_VGB(2,i,j,k,iBlock)= 0.0

                     
                     State_VGB(4,i,j,k,iBlock)=  0.0
         
                     State_VGB(5,i,j,k,iBlock)= 0.0
                     State_VGB(6,i,j,k,iBlock)= 0.0
                     State_VGB(7,i,j,k,iBlock)= 0.0
   
                     State_VGB(8,i,j,k,iBlock)= y+1.5
                     if(x<0.125)then
                        State_VGB(3,i,j,k,iBlock)=-0.025*cos(8*cPi*x)*1.0*sqrt((5.0/3.0)*(y+1.5)/1.0)
                     else
                        State_VGB(3,i,j,k,iBlock)=-0.025*cos(0.25-8*cPi*x)*1.0*sqrt((5.0/3.0)*(y+1.5)/1.0)
      
                     endif
                  endif

            elseif(test2DRiemann)then
                     x = Xyz_DGB(x_,i,j,k,iBlock)
                     y = Xyz_DGB(y_,i,j,k,iBlock)
   
   
                     if(x>0.3d0-1.d-15 .and. y>0.3d0-1.d-15)then
                        State_VGB(1,i,j,k,iBlock)= 1.5d0
      
                        State_VGB(2,i,j,k,iBlock)= 0.d0
                        State_VGB(3,i,j,k,iBlock)= 0.d0
                        State_VGB(4,i,j,k,iBlock)= 0.d0     

                        State_VGB(5,i,j,k,iBlock)= 0.d0
                        State_VGB(6,i,j,k,iBlock)= 0.d0
                        State_VGB(7,i,j,k,iBlock)= 0.d0
      
                        State_VGB(8,i,j,k,iBlock)= 1.5d0
                     elseif(x<0.3d0-1.d-15 .and. y>0.3d0+1.d-15)then

                        State_VGB(1,i,j,k,iBlock)= 0.5323d0
      
                        State_VGB(2,i,j,k,iBlock)= 1.206d0*0.5323d0
                        State_VGB(3,i,j,k,iBlock)= 0.d0                       
                        State_VGB(4,i,j,k,iBlock)=  0.d0
            
                        State_VGB(5,i,j,k,iBlock)= 0.d0
                        State_VGB(6,i,j,k,iBlock)= 0.d0
                        State_VGB(7,i,j,k,iBlock)= 0.d0
      
                        State_VGB(8,i,j,k,iBlock)= 0.3d0
                     elseif(x<0.3d0+1.d-15 .and. y<0.3d0+1.d-15)then

                        State_VGB(1,i,j,k,iBlock)= 0.138d0
      
                        State_VGB(2,i,j,k,iBlock)= 1.206d0*0.138d0
                        State_VGB(3,i,j,k,iBlock)= 1.206d0*0.138d0                     
                        State_VGB(4,i,j,k,iBlock)=  0.0d0
            
                        State_VGB(5,i,j,k,iBlock)= 0.0d0
                        State_VGB(6,i,j,k,iBlock)= 0.0d0
                        State_VGB(7,i,j,k,iBlock)= 0.0d0
      
                        State_VGB(8,i,j,k,iBlock)= 0.029d0
                     elseif(x>0.3d0+1.d-15 .and. y<0.3d0-1.d-15)then

                        State_VGB(1,i,j,k,iBlock)= 0.5323d0
      
                        State_VGB(2,i,j,k,iBlock)= 0.0d0
                        State_VGB(3,i,j,k,iBlock)= 1.206d0*0.5323d0                     
                        State_VGB(4,i,j,k,iBlock)=  0.0d0
            
                        State_VGB(5,i,j,k,iBlock)= 0.0d0
                        State_VGB(6,i,j,k,iBlock)= 0.0d0
                        State_VGB(7,i,j,k,iBlock)= 0.0d0
      
                        State_VGB(8,i,j,k,iBlock)= 0.3d0
                     endif
   

            elseif(testSphericalExplosion)then
                        x = Xyz_DGB(x_,i,j,k,iBlock)
                        y = Xyz_DGB(y_,i,j,k,iBlock)
                        radius = sqrt((x-0.0)**2+(y-0.0)**2) 
      
      
                        if(radius<10.0)then
                           State_VGB(1,i,j,k,iBlock)= 1.0
         
                           State_VGB(2,i,j,k,iBlock)= 0.0     
                           State_VGB(3,i,j,k,iBlock)= 0.0
                           State_VGB(4,i,j,k,iBlock)= 0.0
               
                           State_VGB(5,i,j,k,iBlock)= 0.0
                           State_VGB(6,i,j,k,iBlock)= 50.0/(sqrt(cPi))
                           State_VGB(7,i,j,k,iBlock)= 0.0
         
                           State_VGB(8,i,j,k,iBlock)= 100.0
                        elseif(radius<12.0 .and. radius>10.0)then
                           State_VGB(1,i,j,k,iBlock)= 1.0

                           State_VGB(2,i,j,k,iBlock)= 0.0     
                           State_VGB(3,i,j,k,iBlock)= 0.0
                           State_VGB(4,i,j,k,iBlock)= 0.0
               
                           State_VGB(5,i,j,k,iBlock)= 0.0
                           State_VGB(6,i,j,k,iBlock)= 50.0/(sqrt(cPi))
                           State_VGB(7,i,j,k,iBlock)= 0.0
         
                           State_VGB(8,i,j,k,iBlock)= 100.0-(99.0/2)*(radius-10.0)
                        else
                           State_VGB(1,i,j,k,iBlock)= 1.0

                           State_VGB(2,i,j,k,iBlock)= 0.0     
                           State_VGB(3,i,j,k,iBlock)= 0.0
                           State_VGB(4,i,j,k,iBlock)= 0.0
               
                           State_VGB(5,i,j,k,iBlock)= 0.0
                           State_VGB(6,i,j,k,iBlock)= 50.0/(sqrt(cPi))
                           State_VGB(7,i,j,k,iBlock)= 0.0
         
                           State_VGB(8,i,j,k,iBlock)= 1.0

                        endif
            elseif(testOrszarg)then
               x = Xyz_DGB(x_,i,j,k,iBlock)
               y = Xyz_DGB(y_,i,j,k,iBlock)

               State_VGB(1,i,j,k,iBlock)= 25.0/9.0
      
               State_VGB(2,i,j,k,iBlock)= -sin(y)*25.0/9.0
               State_VGB(3,i,j,k,iBlock)=  sin(x)*25.0/9.0                     
               State_VGB(4,i,j,k,iBlock)=  0.0
    
               State_VGB(5,i,j,k,iBlock)= -sin(y)
               State_VGB(6,i,j,k,iBlock)=  sin(2.0*x)
               State_VGB(7,i,j,k,iBlock)=  0.0

               State_VGB(8,i,j,k,iBlock)= 5.0/3.0








                  


 






               






















                if(DoTestCell)write(*,*) NameSub,': InitalState=', &
                     State_VGB(:,i,j,k,iBlock)
             
                     elseif(UseShocktube)then
                if( (Xyz_DGB(x_,i,j,k,iBlock) - ShockPosition) &
                     < -ShockSlope*Xyz_DGB(y_,i,j,k,iBlock))then
                   ! Set all variables first
                   State_VGB(:,i,j,k,iBlock) = ShockLeft_V

                   if(DoTestCell)write(*,*) NameSub,': left State=', &
                        State_VGB(:,i,j,k,iBlock)
#ifndef SCALAR
                   ! Rotate vector variables
                   do iFluid = 1, nFluid
                      if(nFluid > 1) call select_fluid(iFluid)
                      State_VGB(iUx:iUy,i,j,k,iBlock) = &
                           matmul(Rot_II, ShockLeft_V(iUx:iUy))
                   end do
                   if(UseB) State_VGB(Bx_:By_,i,j,k,iBlock) = &
                        matmul(Rot_II, ShockLeft_V(Bx_:By_))
#endif
                   if(DoTestCell)write(*,*) NameSub,': final left State=', &
                        State_VGB(:,i,j,k,iBlock)
                else
                   ! Set all variables first
                   State_VGB(:,i,j,k,iBlock) = ShockRight_V
                   if(DoTestCell)write(*,*) NameSub,': right State=', &
                        State_VGB(:,i,j,k,iBlock)
#ifndef SCALAR
                   ! Set vector variables
                   do iFluid = 1, nFluid
                      if(nFluid > 1) call select_fluid(iFluid)
                      State_VGB(iUx:iUy,i,j,k,iBlock) = &
                           matmul(Rot_II, ShockRight_V(iUx:iUy))
                   end do
                   if(UseB) State_VGB(Bx_:By_,i,j,k,iBlock) = &
                        matmul(Rot_II, ShockRight_V(Bx_:By_))
#endif
                   if(DoTestCell)write(*,*) NameSub,': final right State=', &
                        State_VGB(:,i,j,k,iBlock)
                end if

                ! Apply "wave" perturbations
                if(UseWave)then
                   call apply_wave(i, j, k, iBlock)
                   if(DoTestCell)write(*,*) NameSub,': wave State=', &
                        State_VGB(:,i,j,k,iBlock)
                end if

#ifndef SCALAR
                ! Convert velocity to momentum
                do iFluid = 1, nFluid
                   if(nFluid > 1) call select_fluid(iFluid)
                   State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock) = &
                        State_VGB(iRho,i,j,k,iBlock) * &
                        State_VGB(iUx:iUz,i,j,k,iBlock)
                end do
#endif
                if(.not.UseB0)CYCLE
                ! Remove B0 from B (if any)
                State_VGB(Bx_:Bz_,i,j,k,iBlock) = &
                     State_VGB(Bx_:Bz_,i,j,k,iBlock) - B0_DGB(:,i,j,k,iBlock)
                if(DoTestCell)write(*,*) NameSub,': cons shock State=', &
                     State_VGB(:,i,j,k,iBlock)

                ! end if UseShockTube

             elseif(UseWave)then
                ! Apply waves/bump without shocktube
                call apply_wave(i, j, k, iBlock)
                if(DoTestCell)write(*,*) NameSub,': wave State=', &
                     State_VGB(:,i,j,k,iBlock)
#ifndef SCALAR
                ! Convert velocity to momentum
                do iFluid = 1, nFluid
                   if(nFluid > 1) call select_fluid(iFluid)
                   State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock) = &
                        State_VGB(iRho,i,j,k,iBlock) * &
                        State_VGB(iUx:iUz,i,j,k,iBlock)
                end do
#endif
                if(.not.UseB0)CYCLE
                ! Remove B0 from B (if any)
                State_VGB(Bx_:Bz_,i,j,k,iBlock) = &
                     State_VGB(Bx_:Bz_,i,j,k,iBlock) - B0_DGB(:,i,j,k,iBlock)
                if(DoTestCell)write(*,*) NameSub,': cons wave State=', &
                     State_VGB(:,i,j,k,iBlock)

             end if

          end do; end do; end do

          if(EntropyConstant > 0.0)then
             do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
                State_VGB(iP_I,i,j,k,iBlock) = &
                     EntropyConstant*State_VGB(iRho_I,i,j,k,iBlock)**Gamma_I
             end do; end do; end do
             if(DoTest)write(*,*) NameSub,': entropy State=', &
                  State_VGB(:,iTest,jTest,kTest,iBlockTest)
          end if

          ! Correct electron fluid for wave perturbations
          if(UseEfield .and. DoCorrectElectronFluid .and. UseWave) then
             call correct_electronfluid_efield(State_VGB(:,:,:,:,iBlock), &
                  1, nI, 1, nJ, 1, nK, iBlock, DoHallCurrentIn=.true.,    &
                  DoGradPeIn=.false., DoCorrectEfieldIn=DoCorrectEfield)
             if(DoTest)write(*,*) NameSub,': electric State=', &
                  State_VGB(:,iTest,jTest,kTest,iBlockTest)
          end if
        !  if(UseConstrainB) call constrain_ics(iBlock)
         ! dzy this line will set B(t=0) to 0. 

          ! Apply user defined initial conditions
          if(UseUserICs)then
             call user_set_ics(iBlock)
             if(DoTest)write(*,*) NameSub,': user State=', &
                  State_VGB(:,iTest,jTest,kTest,iBlockTest)
          end if

          if(iSignRotationIC /= 0)then
             call add_rotational_velocity(iSignRotationIC, iBlock)
             if(DoTest)write(*,*) NameSub,': rot vel State=', &
                  State_VGB(:,iTest,jTest,kTest,iBlockTest)
          end if
          if(UseSaMhd)then
             call init_samhd(iBlock)
             if(DoTest)write(*,*) NameSub,': samhd State=', &
                  State_VGB(:,iTest,jTest,kTest,iBlockTest)
          end if
       end if ! not IsRestart

    end if ! Unused_B

    if(DoTest)write(*,*) &
         NameSub,': final State=',State_VGB(:,iTest,jTest,kTest,iBlockTest)

    call test_stop(NameSub, DoTest, iBlock)

  contains
    !==========================================================================
    subroutine rotate_wave

      ! Rotate wave parameters with the angle of the shock slope

      integer:: iVar, iVectorVar, iVarX, iVarY
      real:: x, y
      !------------------------------------------------------------------------
      DoRotateWave = .false.

      ! Rotate vector variables
      do iVectorVar = 1, nVectorVar
         ! X and Y indexes of the vector variables
         iVarX = iVectorVar_I(iVectorVar)
         iVarY = iVarX + 1

         ! Make sure that the X and Y components of vector variables
         ! have consistent wave parameters
         if(Width_V(iVarX) > 0.0 .and. Width_V(iVarY) == 0.0) &
              call copy_wave(iVarX, iVarY)
         if(Width_V(iVarY) > 0.0 .and. Width_V(iVarX) == 0.0) &
              call copy_wave(iVarY, iVarX)

         ! Rotate amplitudes with the angle of the shock slope
         x = Ampl_V(iVarX); y = Ampl_V(iVarY)
         Ampl_V(iVarX) = CosSlope*x - SinSlope*y
         Ampl_V(iVarY) = SinSlope*x + CosSlope*y
      end do

      ! Rotate wave vectors
      do iVar = 1, nVar
         x = KxWave_V(iVar); y = KyWave_V(iVar)
         KxWave_V(iVar) = CosSlope*x - SinSlope*y
         KyWave_V(iVar) = SinSlope*x + CosSlope*y
      end do

    end subroutine rotate_wave
    !==========================================================================
    subroutine apply_wave(i, j, k, iBlock)

      ! Apply wave/Gaussian/tophat perturbations at a given grid cell

      use ModGeometry, ONLY: &
           xMinBox, xMaxBox, yMinBox, yMaxBox, zMinBox, zMaxBox

      integer, intent(in):: i, j, k, iBlock

      integer:: iVar
      real:: x, y, z, r, Ampl
      character(len=*), parameter:: NameSub = 'apply_wave'
      !------------------------------------------------------------------------
      do iVar = 1, nVar

         ! Normalize amplitude
         Ampl = Ampl_V(iVar)*Io2No_V(iUnitPrim_V(iVar))

         if(Ampl == 0.0) CYCLE

         if(DoTestCell) &
              write(*,*) NameSub, 'iVar, Ampl, kX, kY, kZ, iPower=', &
              iVar, Ampl, KxWave_V(iVar), KyWave_V(iVar), KzWave_V(iVar), &
              iPower_V(iVar)

         if(iPower_V(iVar) <= 0)then
            ! Bump
            x = Xyz_DGB(x_,i,j,k,iBlock) - x_V(iVar)
            y = Xyz_DGB(y_,i,j,k,iBlock) - y_V(iVar)
            z = Xyz_DGB(z_,i,j,k,iBlock) - z_V(iVar)
            if(IsPeriodic_D(1))then
               if(x > +(xMaxBox-xMinBox)/2) x = x - (xMaxBox-xMinBox)
               if(x < -(xMaxBox-xMinBox)/2) x = x + (xMaxBox-xMinBox)
            end if
            if(IsPeriodic_D(2))then
               if(y > +(yMaxBox-yMinBox)/2) y = y - (yMaxBox-yMinBox)
               if(y < -(yMaxBox-yMinBox)/2) y = y + (yMaxBox-yMinBox)
            end if
            if(IsPeriodic_D(3))then
               if(z > +(zMaxBox-zMinBox)/2) z = z - (zMaxBox-zMinBox)
               if(z < -(zMaxBox-zMinBox)/2) z = z + (zMaxBox-zMinBox)
            end if
            ! normalize
            x = x*KxWave_V(iVar)
            y = y*KyWave_V(iVar)
            z = z*KzWave_V(iVar)
            r = sqrt(x**2 + y**2 + z**2)
            if(r < 0.5) State_VGB(iVar,i,j,k,iBlock) = &
                 State_VGB(iVar,i,j,k,iBlock) &
                 + Ampl * cos(cPi*r)**abs(iPower_V(iVar))
         else
            ! cos^n profile
            if(Width_V(iVar) > 0)then
               if(DoTestCell)write(*,*) NameSub,': Width=', Width_V(iVar)
               if(KxWave_V(iVar) > 0.0)then
                  if(abs(Xyz_DGB(x_,i,j,k,iBlock) &
                       + ShockSlope*Xyz_DGB(y_,i,j,k,iBlock)) &
                       > Width_V(iVar) ) CYCLE
               elseif(KyWave_V(iVar) > 0.0)then
                  if(abs(Xyz_DGB(y_,i,j,k,iBlock)) > Width_V(iVar) ) CYCLE
               elseif(KzWave_V(iVar) > 0.0)then
                  if(abs(Xyz_DGB(z_,i,j,k,iBlock)) > Width_V(iVar) ) CYCLE
               end if
            end if

            State_VGB(iVar,i,j,k,iBlock) =        &
                 State_VGB(iVar,i,j,k,iBlock)          &
                 + Ampl*cos(Phase_V(iVar)      &
                 + KxWave_V(iVar)*Xyz_DGB(x_,i,j,k,iBlock)  &
                 + KyWave_V(iVar)*Xyz_DGB(y_,i,j,k,iBlock)  &
                 + KzWave_V(iVar)*Xyz_DGB(z_,i,j,k,iBlock))**iPower_V(iVar)
         end if
      end do

    end subroutine apply_wave
    !==========================================================================
    subroutine copy_wave(iVar, jVar)

      ! Copy wave parameters from iVar to jVar for rotated problems

      integer, intent(in):: iVar, jVar
      logical:: DoTest
      character(len=*), parameter:: NameSub = 'copy_wave'
      !------------------------------------------------------------------------
      call test_start(NameSub, DoTest)
      Width_V(jVar)  = Width_V(iVar)
      KxWave_V(jVar) = KxWave_V(iVar)
      KyWave_V(jVar) = KyWave_V(iVar)
      KzWave_V(jVar) = KzWave_V(iVar)
      Phase_V(jVar)  = Phase_V(iVar)
      iPower_V(jVar) = iPower_V(iVar)

      call test_stop(NameSub, DoTest)
    end subroutine copy_wave
    !==========================================================================
  end subroutine set_initial_condition
  !============================================================================
  subroutine add_rotational_velocity(iSign, iBlockIn)

    use ModSize,     ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK, nBlock, x_, y_
    use ModMain,     ONLY: Unused_B, NameThisComp
    use ModAdvance,  ONLY: State_VGB
    use ModGeometry, ONLY: Used_GB
    use ModPhysics,  ONLY: OmegaBody
    use ModMultiFluid, ONLY: iRho_I, iRhoUx_I, iRhoUy_I
    use BATL_lib,    ONLY: Xyz_DGB, iProc

    integer, intent(in):: iSign
    integer, optional, intent(in):: iBlockIn

    ! Transform velocities between inertial and rotating frames
    ! where Omega is the angular velocity of the rotating frame
    ! Since Omega = (0,0,OmegaBody)
    ! ux = ux - iSign*OmegaBody*y
    ! uy = uy + iSign*OmegaBody*x
    ! iSign=+1: from rotating to inertial frame
    ! iSign=-1: from inertial to rotating frame
    !
    ! If iBlockIn is present, do that block, otherwise do all blocks.

    integer :: i, j, k, iBlock, iBlockFirst, iBlockLast
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'add_rotational_velocity'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    if(present(iBlockIn))then
       iBlockFirst = iBlockIn; iBlockLast = iBlockIn
    else
       iBlockFirst = 1; iBlockLast = nBlock
       if(iProc==0)write(*,'(a)')&
            NameThisComp//': add rotational velocity to convert coords'
    end if

    do iBlock = iBlockFirst, iBlockLast
       if(Unused_B(iBlock))CYCLE
       do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI
          if(.not.Used_GB(i,j,k,iBlock)) CYCLE
          State_VGB(iRhoUx_I,i,j,k,iBlock) = State_VGB(iRhoUx_I,i,j,k,iBlock) &
               - iSign*State_VGB(iRho_I,i,j,k,iBlock) &
	       *OmegaBody*Xyz_DGB(y_,i,j,k,iBlock)

          State_VGB(iRhoUy_I,i,j,k,iBlock) = State_VGB(iRhoUy_I,i,j,k,iBlock) &
	       + iSign*State_VGB(iRho_I,i,j,k,iBlock) &
	       *OmegaBody*Xyz_DGB(x_,i,j,k,iBlock)

       end do; end do; end do
    end do

    call test_stop(NameSub, DoTest)
  end subroutine add_rotational_velocity
  !============================================================================

end module ModSetInitialCondition
!==============================================================================
