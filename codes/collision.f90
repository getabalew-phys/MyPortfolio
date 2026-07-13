PROGRAM collision
    IMPLICIT NONE

    !=======================================================================
    ! parameters
    !=======================================================================

    ! simulation control
    INTEGER,      PARAMETER :: Np     = 150           ! number of particles
    INTEGER,      PARAMETER :: D      = 3             ! dimension of the grid
    REAL(KIND=4), PARAMETER :: t_span = 1e7       ! total time span
    REAL(KIND=4), PARAMETER :: dt     = 1e2          ! time step
    INTEGER,      PARAMETER :: Nloop  = t_span/dt     ! number of time steps

    ! physical constants
    REAL(KIND=4), PARAMETER :: kb = 1.38e-6           ! boltzmann constant
    REAL(KIND=4), PARAMETER :: PI = 4.0*ATAN(1.0)    ! pi

    ! grid and particle size
    REAL(KIND=4), PARAMETER :: l  = 100.0             ! grid size
    REAL(KIND=4), PARAMETER :: dl = 4.0              ! particle size

    !=======================================================================
    ! variables
    !=======================================================================

    ! position, velocity, momentum of each particle
    REAL(KIND=4), DIMENSION(Np, D) :: X, V, P
    REAL(KIND=4), DIMENSION(Np)    :: M              ! mass of each particle

    ! relative velocity, CM velocity, reduced mass
    REAL(KIND=4) :: v_r(D), v_cm(D), m_reduced
    REAL(KIND=4) :: r_hat(D)                         ! unit separation vector
    REAL(KIND=4) :: rand(D)                          ! normalized random vector workspace

    ! initial temperature and maxwellian sampling variables
    REAL(KIND=4) :: T_init = 300.0                   ! initial temperature
    REAL(KIND=4) :: sigma                            ! sqrt(kb*T/m), maxwellian width
    REAL(KIND=4) :: u1, u2                           ! box-muller uniform samples

    ! time and collision counters
    REAL(KIND=4) :: t = 0.0
    INTEGER      :: count        = 0                 ! timestep count
    INTEGER      :: boundary_hit = 0                 ! boundary crossing count
    INTEGER      :: Ncol         = 0                 ! collision count

    ! tracer particle arrays (particle 1, x-component)
    REAL, DIMENSION(Nloop) :: Vx_tracer, Vy_tracer, Vz_tracer, X_tracer, Y_tracer, Z_tracer

    ! output parameters
    REAL :: Vin, Vf                                  ! initial and final mean speed
    REAL :: Temp_sim                                 ! measured temperature
    REAL :: vol_box, Nd                              ! box volume, number density
    REAL :: cross_section_area                       ! pi*dl^2
    REAL :: lambda_theory, lambda_sim, lambda_analysis
    REAL :: coll_freq_sim, coll_freq_theory

    ! loop indices
    INTEGER :: i, j

    ! cpu timing
    REAL :: begin, endt

    !=======================================================================
    ! initialisation
    !=======================================================================

    call cpu_time(begin)

    X = 0.0;  V = 0.0;  P = 0.0;  M = 10.0

    sigma = SQRT(kb * T_init / M(1))

    ! randomly initializing particle positions and velocities
    loop1: DO i = 1, Np
        loop2: DO j = 1, D
            CALL RANDOM_NUMBER(X(i,j))
            CALL RANDOM_NUMBER(u1)
            CALL RANDOM_NUMBER(u2)
            X(i,j) = -l + 2*l*X(i,j)                         ! uniform in [-l, l]
            V(i,j) = sigma * SQRT(-2*LOG(u1)) * COS(2*PI*u2)  ! box-muller maxwellian
        END DO loop2
    END DO loop1

    !V(1,:) = 0.9; M(1)=0.1
    Vin = SQRT(SUM(V**2)) / Np                               ! initial average speed
    

    ! output of initial velocity distribution
    OPEN(UNIT=9, FILE="initial_v.csv", STATUS="REPLACE")
    WRITE(9, '(3A20)') "vx", "vy", "vz"
    DO i = 1, Np
        WRITE(9, '(3F20.10)') V(i,1), V(i,2), V(i,3)
    END DO
    CLOSE(9)

    !=======================================================================
    ! collision kernel
    !=======================================================================

    loop5: DO WHILE (t < t_span)

        ! updating time and positions
        t     = t + dt
        X     = X + V*dt
        count = count + 1

        ! recording tracer particle
        Vx_tracer(count) = V(1,1)
        X_tracer(count) = X(1,1)
        Y_tracer(count) = X(1,2)


        loop6: DO i = 1, Np
            loop7: DO j = i+1, Np

                ! looking for approaching particle pairs
                IF (SUM((X(i,:) - X(j,:))**2) < dl**2) THEN

                    v_r       = V(i,:) - V(j,:)                              ! relative velocity
                    v_cm      = (M(i)*V(i,:) + M(j)*V(j,:)) / (M(i) + M(j)) ! center of mass velocity
                    m_reduced = M(i)*M(j) / (M(i) + M(j))                   ! reduced mass

                    ! generating normalized random 3D vector for post-collision direction
                    CALL RANDOM_NUMBER(rand)
                    rand = 2*rand - 1
                    rand = rand / SQRT(SUM(rand**2))
                    v_r = SQRT(SUM(v_r**2)) * rand

                    r_hat = (X(i,:) - X(j,:)) / SQRT(SUM((X(i,:) - X(j,:))**2)) ! unit separation vector

                    ! updating velocities and resolving overlap
                    V(i,:) = v_cm + (m_reduced/M(i)) * v_r
                    V(j,:) = v_cm - (m_reduced/M(j)) * v_r
                    X(i,:) = X(i,:) + 0.6*dl * r_hat
                    X(j,:) = X(j,:) - 0.6*dl * r_hat

                    Ncol = Ncol + 1
                END IF

            END DO loop7

            ! imposing periodic boundary condition
            loop8: DO j = 1, D
                IF (ABS(X(i,j)) > l) THEN
                    X(i,j)       = X(i,j) - 2*l * NINT(X(i,j)/(2*l))
                    boundary_hit = boundary_hit + 1
                END IF
            END DO loop8

        END DO loop6

    END DO loop5

    call cpu_time(endt)

    !=======================================================================
    ! computing output parameters
    !=======================================================================

    Vf                 = SQRT(SUM(V**2)) / Np
    Temp_sim           = M(1) * SUM(V**2) / (3.0 * Np * kb)  ! equipartition temperature
    vol_box            = (2*l)**3
    Nd                 = Np / vol_box                          ! number density
    cross_section_area = PI * dl**2
    lambda_theory      = 0.707 / (cross_section_area * Nd)
    coll_freq_theory   = Vf / lambda_theory
    coll_freq_sim      = REAL(Ncol) / (Np * t_span)
    lambda_sim         = Vin / coll_freq_sim
    lambda_analysis    = (vol_box / (4.0*PI*(dl/2)**3/3.0)) &
                         * dl / (Np-1) * 0.707 / Vf

    !=======================================================================
    ! output
    !=======================================================================

    ! writing tracer particle trajectory
    OPEN(UNIT=9, FILE="tracer.csv", STATUS="REPLACE")
    WRITE(9, '(3A20)') "vx_tracer", "x_tracer", "y_tracer"
    DO i = 1, Nloop
        WRITE(9, '(3F20.10)') Vx_tracer(i), X_tracer(i), Y_tracer(i)
    END DO
    CLOSE(9)

    ! output of final velocity distribution
    OPEN(UNIT=9, FILE="final_v.csv", STATUS="REPLACE")
    WRITE(9, '(3A20)') "vfx", "vfy", "vfz"
    DO i = 1, Np
        WRITE(9, '(3F20.10)') V(i,1), V(i,2), V(i,3)
    END DO
    CLOSE(9)

    ! printing simulation summary
    WRITE(*, '(A)') "===== simulation summary ====="
    WRITE(*, '(A, I10)')    "  timesteps          : ", count
    WRITE(*, '(A, F25.1)')  "  final time[s]         : ", t
    WRITE(*, '(A, I10)')    "  boundary crossings : ", boundary_hit
    WRITE(*, '(A, I10)')    "  total collisions   : ", Ncol
    WRITE(*, '(A, F12.6)')  "  Vin, Vf/Vin             : ",Vin, "  ", Vf/Vin
    WRITE(*, '(A, F12.4)')  "  T_sim [K]          : ", Temp_sim
    WRITE(*, '(A, ES12.4)') "  lambda (theory)    : ", lambda_theory
    WRITE(*, '(A, ES12.4)') "  lambda (sim)       : ", lambda_sim
    WRITE(*, '(A, ES12.4)') "  coll freq (sim)    : ", coll_freq_sim
    WRITE(*, '(A, ES12.4)') "  coll freq (theory) : ", coll_freq_theory
    WRITE(*, '(A, F12.4)')  "  wall hits/collision: ", REAL(boundary_hit)/Ncol
    WRITE(*, '(A, F12.4)')  "  collisions/particle: ", REAL(Ncol)/Np
    WRITE(*, '(A, F10.3)')  "  cpu time [s]       : ", endt - begin
    WRITE(*, '(A)') "=============================="

END PROGRAM collision