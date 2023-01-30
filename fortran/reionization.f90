
    module Reionization
    use Precision
    use MiscUtils
    use MathUtils
    use classes
    use constants
    use results
    implicit none
    private
    !This module puts smooth asymmetric reionization of specified mid-point (z_{re})
    !and end-point (z_{end})

    !Rionization_zexp=1.5 has the property that for the same z_{re}
    !the optical depth agrees with infinitely sharp model for matter domination
    !So tau and zre can be mapped into each other easily (for any symmetric window)
    !However for generality the module maps tau into z_{re} using a binary search
    !so could be easily modified for other monatonic parameterizations.

    !The ionization history must be twice differentiable.

    !AL March 2008
    !AL July 2008 - added trap for setting optical depth without use_optical_depth
    !AG modified 2022 for asymmetric reionization

    !See CAMB notes for further discussion: http://cosmologist.info/notes/CAMB.pdf

    real(dl), parameter :: TTanhReionization_DefFraction = -1._dl
    !if -1 set from YHe assuming Hydrogen and first ionization of Helium follow each other

    real(dl) :: Tanh_zexp = 1.5_dl


    type, extends(TReionizationModel) :: TTanhReionization
        logical    :: use_optical_depth = .false.
        real(dl)   :: redshift = 0._dl
        real(dl)   :: optical_depth = 0._dl
        real(dl)   :: fraction = TTanhReionization_DefFraction
        !Parameters for the second reionization of Helium
        logical    :: include_helium_fullreion  = .true.
        real(dl)   :: helium_redshift  = 3.5_dl
        real(dl)   :: helium_delta_redshift  = 0.4_dl
        real(dl)   :: helium_redshiftstart  = 5.5_dl
        real(dl)   :: tau_solve_accuracy_boost = 1._dl
        real(dl)   :: timestep_boost =  1._dl
        real(dl)   :: z_end = 0._dl
        real(dl)   :: z_early = 20._dl
        !The rest are internal to this module.
        real(dl), private :: reion_tol = 1d-5
        real(dl), private ::  fHe
        real(dl), private :: alpha = 0_dl
        class(CAMBdata), pointer :: State
    contains
    procedure :: ReadParams => TTanhReionization_ReadParams
    procedure :: Init => TTanhReionization_Init
    procedure :: x_e => TTanhReionization_xe
    procedure :: get_timesteps => TTanhReionization_get_timesteps
    procedure, nopass :: SelfPointer => TTanhReionization_SelfPointer
    procedure, nopass ::  GetZreFromTau => TTanhReionization_GetZreFromTau
    procedure, nopass ::  GetTauFromXe => TTanhReionization_GetTauFromXe
    procedure, private :: SetParamsForAlpha => TTanhReionization_SetParamsForAlpha
    procedure, private :: SetAlphaForZre => TTanhReionization_SetAlphaForZre
    procedure, private :: zreFromOptDepth => TTanhReionization_zreFromOptDepth
    end type TTanhReionization

    public TTanhReionization
    contains


    function TTanhReionization_xe(this, z, tau, xe_recomb)
    !a and time tau and redundant, both provided for convenience
    !xe_recomb is xe(tau_start) from recombination (typically very small, ~2e-4)
    !xe should map smoothly onto xe_recomb
    class(TTanhReionization) :: this
    real(dl), intent(in) :: z
    real(dl), intent(in), optional :: tau, xe_recomb
    real(dl) TTanhReionization_xe
    real(dl) tgh, xod
    real(dl) xstart
    write(*, *) 'test2', xe_recomb, tau

    xstart = PresentDefault( 0._dl, xe_recomb)

    if (z < this%z_end) then
        TTanhReionization_xe = this%fraction-xstart
    else if (z > this%z_early) then
        TTanhReionization_xe = xstart
    else 
        TTanhReionization_xe = ((this%z_early-z)/(this%z_early - this%z_end ))**(this%alpha)
        TTanhReionization_xe = (this%fraction-xstart)*TTanhReionization_xe
    end if

    if (this%include_helium_fullreion .and. z < this%helium_redshiftstart) then

        !Effect of Helium becoming fully ionized is small so details not important
        xod = (this%helium_redshift - z)/this%helium_delta_redshift
        if (xod > 100) then
            tgh=1.d0
        else
            tgh=tanh(xod)
        end if

        TTanhReionization_xe =  TTanhReionization_xe + this%fHe*(tgh+1._dl)/2._dl

    end if

    end function TTanhReionization_xe

    subroutine TTanhReionization_get_timesteps(this, n_steps, z_start, z_complete)
    !minimum number of time steps to use between tau_start and tau_complete
    !Scaled by AccuracyBoost later
    !steps may be set smaller than this anyway
    class(TTanhReionization) :: this
    integer, intent(out) :: n_steps
    real(dl), intent(out):: z_start, z_Complete

    n_steps = nint(50 * this%timestep_boost)
    z_start = this%z_early
    z_complete = 0.d0

    end subroutine TTanhReionization_get_timesteps

    subroutine TTanhReionization_ReadParams(this, Ini)
    use IniObjects
    class(TTanhReionization) :: this
    class(TIniFile), intent(in) :: Ini

    this%Reionization = Ini%Read_Logical('reionization')
    if (this%Reionization) then

        this%use_optical_depth = Ini%Read_Logical('re_use_optical_depth')
        if (this%use_optical_depth) then
            this%optical_depth = Ini%Read_Double('re_optical_depth')
        else
            this%redshift = Ini%Read_Double('re_redshift')
            this%z_end = Ini%Read_Double('re_zend')
        end if

        call Ini%Read('re_ionization_frac',this%fraction)
        call Ini%Read('re_helium_redshift',this%helium_redshift)
        call Ini%Read('re_helium_delta_redshift',this%helium_delta_redshift)

        this%helium_redshiftstart  = Ini%Read_Double('re_helium_redshiftstart', &
            this%helium_redshift + 5*this%helium_delta_redshift)

    end if

    end subroutine TTanhReionization_ReadParams

    subroutine TTanhReionization_SetParamsForAlpha(this)
    class(TTanhReionization) :: this

    if (this%alpha < 1.) then
        this%redshift = this%z_early
    else
        this%redshift = this%z_early - (this%z_early - this%z_end)/2.**(1./this%alpha)
    end if

    end subroutine TTanhReionization_SetParamsForAlpha

    subroutine TTanhReionization_SetAlphaForZre(this)
    class(TTanhReionization) :: this

    this%alpha = log(1./2.) / log((this%z_early-this%redshift)/(this%z_early- this%z_end))
    if (this%alpha <= 1) then
        write (*,*) 'WARNING: there is an issue with your zre value:', this%redshift
        this%alpha = 1.
    end if

    end subroutine TTanhReionization_SetAlphaForZre


    subroutine TTanhReionization_Init(this, State)
    use constants
    use MathUtils
    class(TTanhReionization) :: this
    class(TCAMBdata), target :: State
    procedure(obj_function) :: dtauda

    select type (State)
    class is (CAMBdata)
        this%State => State

        this%fHe =  State%CP%YHe/(mass_ratio_He_H*(1.d0-State%CP%YHe))

        if (this%Reionization) then

            if (this%optical_depth /= 0._dl .and. .not. this%use_optical_depth) &
                write (*,*) 'WARNING: You seem to have set the optical depth, but use_optical_depth = F'

            if (this%use_optical_depth.and.this%optical_depth<0.001 &
                .or. .not.this%use_optical_depth .and. this%Redshift<0.001) then
                this%Reionization = .false.
            end if

            if (this%fraction==TTanhReionization_DefFraction) then
                this%fraction = 1._dl + this%fHe  !H + singly ionized He
            end if


            if (this%use_optical_depth) then 
                this%redshift = 0._dl
                this%alpha = 0._dl
                this%z_end = 5.5_dl
                call this%zreFromOptDepth()
                call this%SetAlphaForZre()
                if (global_error_flag/=0) return
                if (FeedbackLevel > 1) then
                    write(*,'("Reion redshift       =  ",f6.3)') this%redshift
                    write(*,'("Reion endpoint       =  ",f6.3)') this%z_end
                end if
            else
                if (this%redshift == 0._dl) write (*,*) 'WARNING: You seem to have set use_optical_depth = F but redshift = 0'
                if (this%z_end == 0._dl) then
                    write (*,*) 'WARNING: z_end not set. Automatically setting it to 5.5'
                    this%z_end = 5.5_dl
                end if
                call this%SetAlphaForZre()
                this%optical_depth = this%State%GetReionizationOptDepth()
                if (FeedbackLevel > 1) write(*,'("Optical depth       =  ",f6.3)') this%optical_depth
            end if

            !this is a check, agrees very well in default parameterization
            ! if (FeedbackLevel > 1) write(*,'("Integrated opt depth = ",f7.4)') this%State%GetReionizationOptDepth()

        end if
    end select
    end subroutine TTanhReionization_Init

    subroutine TTanhReionization_Validate(this, OK)
    class(TTanhReionization),intent(in) :: this
    logical, intent(inout) :: OK

    if (this%Reionization) then
        if (this%use_optical_depth) then
            if (this%optical_depth<0.04 .or. this%optical_depth > 0.9  .or. &
                this%include_helium_fullreion .and. this%optical_depth<0.01) then
                OK = .false.
                write(*,*) 'Optical depth is strange. You have:', this%optical_depth
            end if
        else
            if (this%redshift < this%z_end .or. this%redshift > this%z_early .or. &
                this%include_helium_fullreion .and. this%redshift < this%helium_redshift) then
                OK = .false.
                write(*,*) 'Reionization redshift strange. You have: ',this%redshift
            end if
        end if
        if (this%fraction/= TTanhReionization_DefFraction .and. (this%fraction < 0 .or. this%fraction > 1.5)) then
            OK = .false.
            write(*,*) 'Reionization fraction strange. You have: ',this%fraction
        end if
        if (this%alpha > 30 .or. this%alpha < 1.1) then
            !Extreme values of alpha lead to similar reionization histories
            OK = .false.
            write(*,*) 'Reionization alpha is strange. You have: ',this%alpha
        end if
    end if

    end subroutine TTanhReionization_Validate

    subroutine TTanhReionization_zreFromOptDepth(this)
    !General routine to find zre parameter given optical depth
    class(TTanhReionization) :: this
    real(dl) try_b, try_t
    real(dl) tau, last_top, last_bot
    real(dl) a, b, step, criterium
    integer i

    criterium=0.002_dl
    step=1._dl

    !coefficients for linear relation between ln(alpha) and ln(tau) 
    a=-18.
    b=-52.8
    !rough estimate of alpha for given tau to start the binary search
    this%alpha= (this%optical_depth**a) * exp(b) !min value of tau is 0.044 so max value of alpha is about 33

    i=0
    do
        call this%SetParamsForAlpha() ! computes z_re for given alpha
        tau = this%State%GetReionizationOptDepth()
        if (abs(tau - this%optical_depth) <= criterium) exit !success

        if (this%alpha >= 30._dl) then
            step=10._dl !increase the step in alpha because for large alphas a small jump in tau correponds to a large jump in tau (see alpha vs. tau plot)
        else if (this%alpha <  10._dl) then !and reciprocally
            step=0.5_dl
            if (i>50) step=0.1_dl
            if (i>100) step=0.05_dl
            if (i>150) step=0.01_dl
            if (i>150 .and. this%alpha <= 1.1_dl) exit !for these values a small jump in alpha is a huge jump in tau so wont converge
        else 
            step=1._dl
        end if

        if (tau>this%optical_depth) then
            this%alpha= this%alpha + step
        else 
            this%alpha= this%alpha - step
        end if

        if (i>95) write(*,*) i, criterium, step, tau, this%optical_depth, this%alpha, this%redshift

        i=i+1

    end do

    end subroutine TTanhReionization_zreFromOptDepth

    real(dl) function TTanhReionization_GetZreFromTau(P, tau)
    type(CAMBparams) :: P, P2
    real(dl) tau
    integer error
    type(CAMBdata) :: State

    P2 = P

    select type(Reion=>P2%Reion)
    class is (TTanhReionization)
        Reion%Reionization = .true.
        Reion%use_optical_depth = .true.
        Reion%optical_depth = tau
    end select
    call State%SetParams(P2,error)
    if (error/=0)  then
        TTanhReionization_GetZreFromTau = -1
    else
        select type(Reion=>State%CP%Reion)
        class is (TTanhReionization)
            TTanhReionization_GetZreFromTau = Reion%redshift
        end select
    end if

    end function TTanhReionization_GetZreFromTau

    real(dl) function TTanhReionization_GetTauFromXe(P, zre, zend)
    type(CAMBparams) :: P, P2
    real(dl) zre, zend
    integer error
    type(CAMBdata) :: State


    P2 = P

    select type(Reion=>P2%Reion)
    class is (TTanhReionization)
        Reion%Reionization = .true.
        Reion%use_optical_depth = .false.
        Reion%redshift = zre
        Reion%z_end = zend
    end select


    call State%SetParams(P2,error)
    if (error/=0)  then
        TTanhReionization_GetTauFromXe = -1
    else
        select type(Reion=>State%CP%Reion)
        class is (TTanhReionization)
            call Reion%SetAlphaForZre()
            TTanhReionization_GetTauFromXe = Reion%State%GetReionizationOptDepth()

        end select
    end if

    end function TTanhReionization_GetTauFromXe

    subroutine TTanhReionization_SelfPointer(cptr,P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (TTanhReionization), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine TTanhReionization_SelfPointer

    end module Reionization
