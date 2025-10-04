!===============================================================================
! PROGRAM: CosmologyMCMC.f90
!
! DESCRIPTION:
!     This program performs parallel Markov Chain Monte Carlo (MCMC) sampling
!     using the Metropolis-Hastings algorithm to estimate cosmological parameters.
!     It incorporates observational data from Type Ia Supernovae (SN), Baryon
!     Acoustic Oscillations (BAO), and the Cosmic Microwave Background (CMB).
!
!     The program uses MPI (Message Passing Interface) to run multiple MCMC chains
!     concurrently, each on a separate processor. Chains are synchronized 
!     periodically to check for convergence using the Gelman-Rubin Rhat diagnostic.
!
! FEATURES:
!     - Supports flat, open, and closed LCDM models
!     - Uses DL(z), DA(z), and H(z) observables from SN, BAO, and CMB
!     - Modular structure: separate routines for DL, DA, chi-squareds, and priors
!     - Adaptive stopping criterion based on Gelman-Rubin Rhat statistic
!     - MPI parallelization using MPI_Bcast, MPI_Gatherv, etc.
!
! INPUT:
!     - Initial parameter values (initpoint)
!     - Prior bounds (priormin, priormax)
!     - Jump sizes (jumpsize)
!     - Observational data files (SN, BAO, CMB), assumed to be preloaded in modules
!
! OUTPUT:
!     - Accepted parameter samples (MCMC chains) compatible with
!       the GetDist analysis tool
!     - Convergence statistics
! 
! Dependencies:
!   - LAPACK and BLAS libraries for linear algebra operations
!
! USAGE:
!     - Compile with MPI Fortran compiler, e.g.:
!           mpiifort -O2 CosmologyMCMC.f90 -o CosmologyMCMC -llapack -lblas
!     - Run with desired number of chains/processes:
!           mpirun -np <nchains> ./CosmologyMCMC
!
! NOTES:
!     - Chains run independently and share convergence data via MPI
!
! Repository:
!   https://github.com/krezazadeh/CosmologyMCMC
!
! Reference:
!   https://arxiv.org/abs/????.?????
!
! Author:
!   Kazem Rezazadeh
!   School of Astronomy,
!   Institute for Research in Fundamental Sciences (IPM)
!   kazem.rezazadeh@ipm.ir
!   https://github.com/krezazadeh
!
! Date:
!   21 May 2025
!=======================================================================

program CosmologyMCMC

    use mpi

    implicit none
    
    ! Constants
    real(8), parameter :: pi = acos(-1.0d0)
    real(8), parameter :: TCMB = 2.7255d0
    real(8), parameter :: c = 299792.458d0
    real(8), parameter :: KtoeV = 8.621738d-5
    real(8), parameter :: cmtoeV = 1.0d0/8065.544005d0
    real(8), parameter :: rhocroh2 = 8.098d-11 ! eV**4
    real(8), parameter :: ggamma = 2.0d0
    real(8), parameter :: Omegagammah2 = (pi**2/30.0d0)*ggamma*(TCMB*KtoeV)**4/rhocroh2
    real(8), parameter :: Neff = 3.046d0

    ! Model parameters

    ! Constant parameters of the model
    ! real(8), parameter :: ok = 0.0d0
    
    ! Number of varying parameters of the model
    integer, parameter :: nparams = 3
    
    ! Varying parameters of the model
    real(8) :: h, ob, oc
    real(8) :: ok
    
    ! DLsol
    integer, parameter :: n_array_z = 100000
    real(8) :: array_z(n_array_z)
    real(8) :: array_dL(n_array_z)
    real(8) :: array_dL_buffer(n_array_z)

    ! CMB_Planck2018_Zhai2018
    integer, parameter :: ndataCMB = 4
    real(8), dimension(ndataCMB) :: dataCMB
    real(8), dimension(ndataCMB, ndataCMB) :: covCMB, invcovCMB
    ! Curvature of the universe
    logical, parameter :: is_flat = .true.

!     ! CMB_Planck2018_Chen2018
!     integer, parameter :: ndataCMB = 3
!     real(8), dimension(ndataCMB) :: dataCMB
!     real(8), dimension(ndataCMB, ndataCMB) :: covCMB, invcovCMB
!     ! Curvature of the universe
!     logical, parameter :: is_flat = .false.

    ! BAO_DESI_DR2
    integer, parameter :: ndataBAO = 13
    real(8), dimension(ndataBAO) :: zBAO, valueBAO
    character(len=20), dimension(ndataBAO) :: quantityBAO
    real(8), dimension(ndataBAO, ndataBAO) :: covBAO, invcovBAO

    ! SN_DESY5
    integer, parameter :: ndataSN = 1829
    real(8) :: dataSN(ndataSN, 7)
    real(8) :: Dij(ndataSN, ndataSN)
    real(8) :: CijSN(ndataSN, ndataSN), invCijSN(ndataSN, ndataSN)
    real(8) :: Id(ndataSN)

!     ! SN_Pantheon
!     integer, parameter :: ndataSN = 1048
!     real(8) :: dataSN(ndataSN, 18)
!     real(8) :: Dij(ndataSN, ndataSN)
!     real(8) :: CijSN(ndataSN, ndataSN), invCijSN(ndataSN, ndataSN)
!     real(8) :: Id(ndataSN)

    ! MCMC
    integer, parameter :: max_try = 10000
    integer, parameter :: max_iter = 1000
    real(8), parameter :: Rhat_minus_1_tol = 0.01d0
    integer :: ierr, rank, nprocess
    integer :: try
    real(8), allocatable :: points(:, :, :) ! points(i_process, i_point, i_param)
    integer :: npoints
    real(8), allocatable :: points_local(:, :, :) ! points_local(rank + 1, i_iter, i_param) 
    real(8), allocatable :: points_temp(:, :, :) ! points_local(rank + 1, i_iter, i_param) 
    integer :: i, j, k
    character(len=30) :: filename
    integer :: unit
    real(8) :: params(nparams), params_new(nparams)
    real(8), dimension(nparams) :: initpoint
    real(8), dimension(nparams) :: jumpsize
    real(8), dimension(nparams) :: priormin, priormax
    real(8) :: chi2params, chi2params_new
    real(8) :: alpha
    real(8) :: rand
    integer :: NN
    integer :: nchains
    real(8), allocatable :: mean_chain(:, :) ! (i_process, i_param)
    real(8), allocatable :: mean_param(:), B(:), W(:), Rhat(:) ! (i_param)
    real(8) :: Rhat_minus_1
    logical :: converged
    integer, parameter :: sendcount = max_iter*nparams
    integer, allocatable :: recvcounts(:), displs(:)
    real(8), allocatable :: sendbuf(:), recvbuf(:)

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nprocess, ierr)

    call CMB_Planck2018()

    call BAO_DESI_DR2()
    
    call SN_DESY5()
!     call SN_Pantheon()

    allocate(points_local(nprocess, max_iter, nparams))
    allocate(sendbuf(sendcount))

    if (rank == 0) then
        allocate(points(nprocess, max_iter*max_try, nparams))
        allocate(points_temp(nprocess, max_iter, nparams))
        allocate(mean_chain(nprocess, nparams))
        allocate(mean_param(nparams))
        allocate(B(nparams))
        allocate(W(nparams))
        allocate(Rhat(nparams))
        allocate(recvbuf(nprocess*sendcount))
        allocate(recvcounts(nprocess))
        allocate(displs(nprocess))
    end if

    unit = 10 + rank
    write(filename, '(A,I0,A)') 'chains_', rank, '.txt'
    open(unit=unit, file=filename, status='replace')
    
    call random_seed()

    initpoint = (/ 0.688993d0, 0.2511d0, 0.003d0 /)
    priormin = (/ 0.6d0, 0.1d0, -0.1d0 /)
    priormax = (/ 0.8d0, 0.5d0, 0.1d0 /)
    jumpsize = (/ 0.001d0, 0.001d0, 0.001d0 /)

    ! start from the given initial point.
    do i = 1, nparams
        params(i) = initpoint(i)
    end do

    ! start from an arbitrary point in the prior intervals.
    ! do i = 1, nparams
    !    call random_number(rand)
    !    params(i) = priormin(i) + (priormax(i) - priormin(i))*rand
    ! end do

    h = params(1)
    ob = 0.02218d0/h**2
    oc = params(2)
    ok = params(3)
    call DLsol()
    chi2params = chi2total()

    converged = .false.

    if (rank == 0) then
        try = 0
    end if

    do while (.not. converged)
    
    if (rank == 0) then
        try = try + 1
        print *, "Attempt: ", try
        if (try > max_try) then
            print *, "WARNING: MCMC completed all iterations but failed to converge."
            print *, "Consider increasing the number of trys or adjusting initial chain parameters."
            stop
        end if
    end if

    do i = 1, max_iter

        call propose_jump(params, jumpsize, params_new)

        do j = 1, nparams
            if (params_new(j) < priormin(j)) params_new(j) = priormin(j)
            if (params_new(j) > priormax(j)) params_new(j) = priormax(j)
        end do

        h = params_new(1)
        ob = 0.02218d0/h**2
        oc = params_new(2)
        ok = params_new(3)
        call DLsol()
        chi2params_new = chi2total()

        points_local(rank + 1, i, :) = params_new

        write(10 + rank, "(10e25.16)") 1.0d0, chi2params_new/2.0d0, params_new, &
        ob, H0(), om(), ol(), age()

        alpha = min(1.0d0, exp(-0.5d0 * (chi2params_new - chi2params)))
        call random_number(rand)
        if (rand <= alpha) then
            params = params_new
            chi2params = chi2params_new
        end if

    end do

    ! Computation of Gelman-Rubin parameter (Rhat)
    ! The convergence step (Gelman-Rubin) is performed at rank zero.

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    call pack_2D(points_local(rank+1, :, :), sendbuf)

    if (rank == 0) then
        do i = 0, nprocess - 1
            recvcounts(i+1) = sendcount
            displs(i+1) = i*sendcount
        end do
    end if
    
    call MPI_Gatherv(sendbuf, sendcount, MPI_DOUBLE_PRECISION, &
    recvbuf, recvcounts, displs, MPI_DOUBLE_PRECISION, &
    0, MPI_COMM_WORLD, ierr)

    if (rank == 0) then

        call unpack_3D(recvbuf, points_temp)

        npoints = nprocess*max_iter*try
        nchains = nprocess
        NN = max_iter*try
        
        do i = 1, nprocess
            do j = 1, max_iter
                do k = 1, nparams
                    points(i, (try - 1)*max_iter + j, k) = points_temp(i, j, k)
                end do
            end do
        end do

        mean_chain = 0.0
        do k = 1, nparams
            do i = 1, nprocess
                do j = 1, NN
                    mean_chain(i, k) = mean_chain(i, k) + points(i, j, k)
                end do
                mean_chain(i, k) = mean_chain(i, k) / real(NN, 8)
            end do
        end do

        mean_param = 0.0
        do k = 1, nparams
            do i = 1, nprocess
                do j = 1, NN
                    mean_param(k) = mean_param(k) + points(i, j, k)
                end do
            end do
            mean_param(k) = mean_param(k) / real(nprocess*NN, 8)
        end do

        B = 0.0
        W = 0.0
        Rhat = 0.0d0
        do k = 1, nparams
            do i = 1, nprocess
                B(k) = B(k) + (mean_chain(i, k) - mean_param(k))**2
                do j = 1, NN
                    W(k) = W(k) + (points(i, j, k) - mean_chain(i, k))**2
                end do
            end do
            B(k) = B(k)/real(nprocess - 1, 8)
            W(k) = W(k) / real(nprocess*(NN - 1), 8)
            Rhat(k) = (((real(NN - 1, 8))/real(NN, 8))*W(k) + &
            (1.0d0 + 1.0d0/real(nprocess, 8))*B(k))/W(k)
        end do

        Rhat_minus_1 = maxval(abs(Rhat - 1.0d0))
        converged = (Rhat_minus_1 < Rhat_minus_1_tol)

        if (converged) then
            print *, "Converged: Rhat - 1 = ", Rhat_minus_1
        else
            print *, "Not converged: Rhat - 1 = ", Rhat_minus_1
        end if
    
    end if

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    call MPI_BCAST(converged, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

    end do ! while

    deallocate(points_local)
    deallocate(sendbuf)

    if (rank == 0) then
        deallocate(points)
        deallocate(points_temp)
        deallocate(mean_chain)
        deallocate(mean_param)
        deallocate(B)
        deallocate(W)
        deallocate(Rhat)
        deallocate(recvbuf)
        deallocate(recvcounts)
        deallocate(displs)
    end if

    close(unit)

    if (rank == 0) print *, "Done!"

    call MPI_Finalize(ierr)

contains

function Omegagamma()

    implicit none

    real(8) :: Omegagamma

    Omegagamma = Omegagammah2/h**2

end function Omegagamma

function Omeganu()

    implicit none

    real(8) :: Omeganu

    Omeganu = Neff*(7.0d0/8.0d0)*(4.0d0/11.0d0)**(4.0d0/3.0d0)*Omegagamma()

end function Omeganu

function Omegar()

    implicit none

    real(8) :: Omegar

    Omegar = Omegagamma() + Omeganu()

end function Omegar

function Ht(a)

    implicit none

    real(8) :: Ht
    real(8), intent(in) :: a

    Ht = sqrt(1.0d0 - ob - oc + (ob + oc)/a**3 + &
    ok/a**2 - Omegar() + Omegar()/a**4)

end function Ht

function HH(a)

    implicit none

    real(8) :: HH
    real(8), intent(in) :: a

    HH = 100.0*h*Ht(a)

end function HH

function om()

    implicit none

    real(8) :: om

    om = ob + oc

end function om

function obh2()

    implicit none
    
    real(8) :: obh2

    obh2 = ob*h**2

end function obh2

function cs(a)

    implicit none

    real(8) :: cs
    real(8), intent(in) :: a

    cs = c/(sqrt(3.0d0)*sqrt(1.0d0 + (3.0d0*a*obh2())/ &
    (4.0d0*Omegagammah2)))

end function cs

function zCMBHE()

    implicit none

    real(8) :: zCMBHE

    zCMBHE = 1048.0d0*(1.0d0 + 0.00124d0/obh2()**0.738d0)* &
    (1.0d0 + (0.0783d0*(h**2*om())** &
    (0.56d0/(1.0d0 + 21.1d0*obh2()**1.81d0)))/ &
    ((1.0d0 + 39.5d0*obh2()**0.763d0)* &
    obh2()**0.238d0))

end function zCMBHE

function zdragHE()

    implicit none

    real(8) :: zdragHE

    zdragHE = (1291.0d0*(h**2*om())**0.251d0* &
    (1.0d0 + (0.313d0*obh2()** &
    (0.238d0*(h**2*om())**0.223d0)* &
    (1.0d0 + 0.607d0*(h**2*om())**0.674d0))/ &
    (h**2*om())**0.419d0))/ &
    (1.0d0 + 0.659d0*(h**2*om())**0.828d0)

end function zdragHE

function zCMBGA()

    implicit none

    real(8) :: zCMBGA

    zCMBGA = (h**2*om())**(-0.731631d0) + &
    obh2()**0.93681d0*(h**2*om())**0.0192951d0* &
    (937.422d0/obh2()**0.97966d0 + &
    391.672d0/(h**2*om())**0.372296d0)

end function zCMBGA

function zdragGA()

    implicit none

    real(8) :: zdragGA

    zdragGA = (1.0d0 + 428.169d0*obh2()**0.256459d0* &
    (h**2*om())**0.616388d0 + &
    925.56d0*(h**2*om())**0.751615d0)/ &
    (h**2*om())**0.714129d0

end function zdragGA

function zCMB()

    implicit none

    real(8) :: zCMB

    zCMB = zCMBGA()

end function zCMB

function zdrag()

    implicit none

    real(8) :: zdrag

    zdrag = zdragGA()

end function zdrag

subroutine DLsol()

    implicit none
    
    real(8) :: zitial, z_final, dz
    real(8) :: dLitial
    integer :: i
    real(8) :: z_i, dL_i
    real(8) :: z_im1, dL_im1
    real(8) :: k_dL_1, k_dL_2, k_dL_3, k_dL_4
    real(8) :: d1cubicspline, d2cubicspline

    zitial = 0.0d0
    z_final = 1500.0d0
    dz = (z_final - zitial)/real(n_array_z - 1, 8)

    dLitial = 0.0d0

    ! i = 1
    z_i = zitial
    dL_i = dLitial
    
    array_z(1) = z_i
    array_dL(1) = dL_i

    do i = 2, n_array_z
    
        z_im1 = z_i
        dL_im1 = dL_i
    
        k_dL_1 = dz*ddL(z_im1, dL_im1)
        k_dL_2 = dz*ddL(z_im1 + dz/2.0d0, dL_im1 + k_dL_1/2.0d0)
        k_dL_3 = dz*ddL(z_im1 + dz/2.0d0, dL_im1 + k_dL_2/2.0d0)
        k_dL_4 = dz*ddL(z_im1 + dz, dL_im1 + k_dL_3)
    
        z_i = z_im1 + dz
        dL_i = dL_im1 + (1.0d0/6.0d0)*(k_dL_1 + 2.0d0*k_dL_2 + 2.0d0*k_dL_3 + &
        k_dL_4)
    
        array_z(i) = z_i
        array_dL(i) = dL_i

    end do

    d1cubicspline = (array_dL(2) - array_dL(1))/(array_z(2) - array_z(1))
    d2cubicspline = (array_dL(n_array_z) - array_dL(n_array_z - 1))/ &
    (array_z(n_array_z) - array_z(n_array_z - 1))

    call spline(array_z(1:n_array_z), array_dL(1:n_array_z), n_array_z, &
    d1cubicspline, d2cubicspline, array_dL_buffer(1:n_array_z))

end subroutine DLsol

function ddL(z, dL)

    implicit none

    real(8) :: ddL
    real(8), intent(in) :: z, dL

    if (ok == 0.0d0) then
        ddL = dL/(1.0d0 + z) + (1.0d0 + z)/Ht(1.0d0/(1.0d0 + z))
    else if (ok > 0.0d0) then
        ddL = ((1.0d0 + z)**2*Sqrt(1.0d0 + &
        (ok*dL**2)/(1.0d0 + z)**2) + &
        dL*Ht(1.0d0/(1.0d0 + z)))/((1.0d0 + z)*Ht(1.0d0/(1.0d0 + z)))
    else
        ddL = ((1.0d0 + z)**2*Sqrt(Abs(1.0d0 - &
        (Abs(ok)*dL**2)/(1.0d0 + z)**2)) + &
        dL*Ht(1.0d0/(1.0d0 + z)))/((1.0d0 + z)*Ht(1.0d0/(1.0d0 + z)))
    end if

end function ddL

function DL(z)

    implicit none

    real(8) :: DL
    real(8), intent(in) :: z
    real(8) :: dL_cubicspline

    call spline_out(array_z(1:n_array_z), array_dL(1:n_array_z), &
    array_dL_buffer(1:n_array_z), &
    n_array_z, z, dL_cubicspline)

    DL = (c/(100.0d0*h))*dL_cubicspline

end function DL

subroutine spline(x, y, n, d11, d1n, d2)
    integer, intent(in) :: n
    real(8), intent(in) :: x(n), y(n), d11, d1n
    real(8), intent(out) :: d2(n)
    integer i
    real(8) xp, qn, sig, un, xxdiv, u(n-1), d1l, d1r
    
    d1r= (y(2)-y(1))/(x(2)-x(1))
    if (d11>0.99d30) then
        d2(1)=0.0d0
        u(1)=0.0d0
    else
        d2(1)=-0.5d0
        u(1)=(3.0d0/(x(2)-x(1)))*(d1r-d11)
    endif
    
    do i=2,n-1
        d1l=d1r
        d1r=(y(i+1)-y(i))/(x(i+1)-x(i))
        xxdiv=1.0d0/(x(i+1)-x(i-1))
        sig=(x(i)-x(i-1))*xxdiv
        xp=1.0d0/(sig*d2(i-1)+2.0d0)
    
        d2(i)=(sig-1.q0)*xp
    
        u(i)=(6.0d0*(d1r-d1l)*xxdiv-sig*u(i-1))*xp
    end do
    d1l=d1r
    
    if (d1n>0.99d30) then
        qn=0.0d0
        un=0.0d0
    else
        qn=0.5d0
        un=(3.0d0/(x(n)-x(n-1)))*(d1n-d1l)
    endif
    
    d2(n)=(un-qn*u(n-1))/(qn*d2(n-1)+1.d0)
    do i=n-1,1,-1
        d2(i)=d2(i)*d2(i+1)+u(i)
    end do
end subroutine spline
    
subroutine spline_out(xarr,yarr,yarr_buff,n,x,y)
    integer n,llo_out,lhi_out,midp
    real(8) xarr(n),yarr(n),yarr_buff(n)
    real(8) x,y,a0_out,b0_out,ho_out

    llo_out=1
    lhi_out=n
    do while((lhi_out-llo_out).gt.1.0d0)
        midp=(llo_out+lhi_out)/2.0d0
        if (xarr(midp).gt.x) then
            lhi_out=midp
        else
            llo_out=midp
        endif
    enddo
    
    ho_out=xarr(lhi_out)-xarr(llo_out)
    a0_out=(xarr(lhi_out)-x)/ho_out
    b0_out=(x-xarr(llo_out))/ho_out
    y=a0_out*yarr(llo_out)+b0_out*yarr(lhi_out)+((a0_out**3.0d0-a0_out)*&
            yarr_buff(llo_out)+(b0_out**3.0d0-b0_out)*&
            yarr_buff(lhi_out))*ho_out*ho_out /6.0d0
end subroutine spline_out

function DA(z)
    
    implicit none

    real(8) :: DA
    real(8) :: z

    DA = (1.0d0/(1.0d0 + z)**2)*DL(z)

end function DA

function rombint(f, a, b, tol)

    !  Rombint returns the integral from a to b of using Romberg integration.
    !  The method converges provided that f(x) is continuous in (a,b).
    !  f must be real(dl) and must be declared external in the calling
    !  routine. tol indicates the desired relative accuracy in the integral.

    implicit none
    integer, parameter :: MAXITER=20
    integer, parameter :: MAXJ=5
    dimension g(MAXJ+1)
    real(8) f
    ! external f
    real(8) :: rombint
    real(8), intent(in) :: a, b, tol
    integer :: nint, i, k, jmax, j
    real(8) :: h_rombint, gmax, error, g, g0, g1, fourj

    h_rombint=0.5d0*(b-a)
    gmax=h_rombint*(f(a)+f(b))
    g(1)=gmax
    nint=1
    error=1.0d20
    i=0
10  i=i+1
    if (i.gt.MAXITER.or.(i.gt.5.and.abs(error).lt.tol)) &
        go to 40
    !  Calculate next trapezoidal rule approximation to integral.
    g0=0.0d0
    do 20 k=1,nint
        g0=g0+f(a+(k+k-1)*h_rombint)
20  continue
    g0=0.5d0*g(1)+h_rombint*g0
    h_rombint=0.5d0*h_rombint
    nint=nint+nint
    jmax=min(i,MAXJ)
    fourj=1.0d0
    do 30 j=1,jmax
        !  Use Richardson extrapolation.
        fourj=4.0d0*fourj
        g1=g0+(g0-g(j))/(fourj-1.0d0)
        g(j)=g0
        g0=g1
30  continue
    if (abs(g0).gt.tol) then
        error=1.0d0-gmax/g0
    else
        error=gmax
    end if
    gmax=g0
    g(jmax+1)=g0
    go to 10
40  rombint=g0
    if (i.gt.MAXITER.and.abs(error).gt.tol)  then
        write(*,*) 'Warning: Rombint failed to converge; '
        write (*,*)'integral, error, tol:', rombint, error, tol
    end if

end function rombint

function rs_zCMB()

    implicit none

    real(8) :: rs_zCMB

    rs_zCMB = rombint(drsdx, 1.0d-10, 1.0d0/(1.0d0 + zCMB()), 1.0d-6)

end function rs_zCMB

function rs_zdrag()

    implicit none

    real(8) :: rs_zdrag
    
    rs_zdrag = rombint(drsdx, 1.0d-10, 1.0d0/(1.0d0 + zdrag()), 1.0d-6)

end function rs_zdrag

function drsdx(x)

    implicit none

    real(8) :: drsdx
    real(8), intent(in) :: x

    drsdx = cs(x)/(x**2*HH(x))

end function drsdx

function DV(zbao)

    implicit none

    real(8) :: DV
    real(8), intent(in) :: zbao

    DV = ((c*zbao*DL(zbao)**2)/ &
    ((1.0d0 + zbao)**2*HH(1.0d0/(1.0d0 + zbao))))**0.3333333333333333d0

end function DV

function dz(zbao)

    implicit none

    real(8) :: dz
    real(8), intent(in) :: zbao

    dz = rs_zdrag()/DV(zbao)

end function dz

function DH(z)

    implicit none

    real(8) :: DH
    real(8), intent(in) :: z

    DH = c/HH(1.0d0/(1.0d0 + z))

end function DH

function DM(z)

    implicit none

    real(8) :: DM
    real(8) :: z

    DM = (1.0d0 + z)*DA(z)

end function DM

function R()

    implicit none

    real(8) :: R

    R = (100.0d0*DA(zCMB())* &
    Sqrt(h**2*om())*(1.0d0 + zCMB()))/c

end function R

function la()

    implicit none

    real(8) :: la

    la = (pi*DA(zCMB())*(1.0d0 + zCMB()))/rs_zCMB()

end function la

! Returns the inverse of a matrix calculated by finding the LU
! decomposition. Depends on LAPACK.

function inv(A) result(Ainv)

    real(8), dimension(:,:), intent(in) :: A
    real(8), dimension(size(A,1),size(A,2)) :: Ainv
    
    real(8), dimension(size(A,1)) :: work  ! work array for LAPACK
    integer, dimension(size(A,1)) :: ipiv   ! pivot indices
    integer :: n, info
    
    ! External procedures defined in LAPACK
    external DGETRF
    external DGETRI
    
    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    n = size(A,1)
    
    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call DGETRF(n, n, Ainv, n, ipiv, info)
    
    if (info /= 0) then
        stop 'Matrix is numerically singular!'
    end if
    
    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
    call DGETRI(n, Ainv, n, ipiv, work, n, info)
    
    if (info /= 0) then
        stop 'Matrix inversion failed!'
    end if
    
end function inv

! CMB_Planck18
! Zhai2018
! arXiv:1811.07425
! Central values for {R, la, obh2, ns} from Planck 2018 likelihood
! base plikHM TTTEEE lowl lowE lensing

subroutine CMB_Planck2018()

    implicit none

    integer :: i, j

    if (is_flat == .true.) then
        open (unit = 11, file = './data/CMB_Planck2018_Zhai2018/CMB_Planck2018-flat.txt', status = 'old')
    else
        open (unit = 11, file = './data/CMB_Planck2018_Zhai2018/CMB_Planck2018-non-flat.txt', status = 'old')
    end if

    do i = 1, ndataCMB
        read(11, *) dataCMB(i)
    end do

    close(11)

    if (is_flat == .true.) then
        open (unit = 11, file = './data/CMB_Planck2018_Zhai2018/CMB_Planck2018-flat-cov.txt', status = 'old')
    else
        open (unit = 11, file = './data/CMB_Planck2018_Zhai2018/CMB_Planck2018-non-flat-cov.txt', status = 'old')
    end if

    do i = 1, ndataCMB
        do j = 1, ndataCMB
            read(11, *) covCMB(i, j)
        end do
    end do

    close(11)

    invcovCMB = inv(covCMB)

end subroutine CMB_Planck2018

function vecCMB()

    implicit none

    real(8), dimension(ndataCMB) :: vecCMB

    vecCMB(1) = R() - dataCMB(1)
    vecCMB(2) = la() - dataCMB(2)
    vecCMB(3) = obh2() - dataCMB(3)
    vecCMB(4) = 0.0d0

end function vecCMB

function chi2_CMB_Planck2018()

    implicit none

    real(8) :: chi2_CMB_Planck2018

    real(8) :: temp(4)

    temp = matmul(invcovCMB, vecCMB())
    chi2_CMB_Planck2018 = dot_product(vecCMB(), temp)

end function chi2_CMB_Planck2018

! ! CMB_Planck18
! ! Chen2018
! ! arXiv:1808.05724
! ! Central values for {R, la, obh2} from Planck 2018 likelihood
! ! base Planck 2018 TT,TE,EE + lowE
!
! subroutine CMB_Planck2018()
!
!     implicit none
!
!     integer :: i, j
!
!     if (is_flat == .true.) then
!         open (unit = 11, file = './data/CMB_Planck2018_Chen2018/CMB_Planck2018-flat.txt', status = 'old')
!     else
!         open (unit = 11, file = './data/CMB_Planck2018_Chen2018/CMB_Planck2018-non-flat.txt', status = 'old')
!     end if
!
!     do i = 1, ndataCMB
!         read(11, *) dataCMB(i)
!     end do
!
!     close(11)
!
!     if (is_flat == .true.) then
!         open (unit = 11, file = './data/CMB_Planck2018_Chen2018/CMB_Planck2018-flat-invcov.txt', status = 'old')
!     else
!         open (unit = 11, file = './data/CMB_Planck2018_Chen2018/CMB_Planck2018-non-flat-invcov.txt', status = 'old')
!     end if
!
!     do i = 1, ndataCMB
!         do j = 1, ndataCMB
!             read(11, *) invcovCMB(i, j)
!         end do
!     end do
!
!     close(11)
!
! end subroutine CMB_Planck2018
!
! function vecCMB()
!
!     implicit none
!
!     real(8), dimension(ndataCMB) :: vecCMB
!
!     vecCMB(1) = R() - dataCMB(1)
!     vecCMB(2) = la() - dataCMB(2)
!     vecCMB(3) = obh2() - dataCMB(3)
!
! end function vecCMB
!
! function chi2_CMB_Planck2018()
!
!     implicit none
!
!     real(8) :: chi2_CMB_Planck2018
!
!     real(8) :: temp(ndataCMB)
!
!     temp = matmul(invcovCMB, vecCMB())
!     chi2_CMB_Planck2018 = dot_product(vecCMB(), temp)
!
! end function chi2_CMB_Planck2018

! BAO_DESI_DR2
! arXiv:2503.14738

subroutine BAO_DESI_DR2()

    implicit none

    integer :: i, j

    open (unit = 11, file = './data/DESI_BAO_DR2/desi_gaussian_bao_ALL_GCcomb_mean.txt', status = 'old')

    do i = 1, ndataBAO
        read(11, *) zBAO(i), valueBAO(i), quantityBAO(i)
    end do

    close(11)

    open (unit = 11, file = './data/DESI_BAO_DR2/desi_gaussian_bao_ALL_GCcomb_cov.txt', status = 'old')

    do i = 1, ndataBAO
        read(11, *) covBAO(i, 1:ndataBAO)
    end do

    close(11)

    invcovBAO = inv(covBAO)

end subroutine BAO_DESI_DR2

function vecBAO()

    implicit none

    integer, parameter :: ndataBAO = 13
    real(8), dimension(ndataBAO) :: vecBAO

    integer :: i

    do i = 1, ndataBAO
        if (quantityBAO(i) == "DV_over_rs") then
            vecBAO(i) = valueBAO(i) - DV(zBAO(i))/ &
            rs_zdrag()
        else if (quantityBAO(i) == "DM_over_rs") then
            vecBAO(i) = valueBAO(i) - DM(zBAO(i))/ &
            rs_zdrag()
        else if (quantityBAO(i) == "DH_over_rs") then
            vecBAO(i) = valueBAO(i) - DH(zBAO(i))/ &
            rs_zdrag()
        end if
    end do

end function vecBAO

function chi2_BAO_DESI_DR2()

    implicit none

    real(8) :: chi2_BAO_DESI_DR2
    
    real(8) :: temp(ndataBAO)

    temp = matmul(invcovBAO, vecBAO())
    chi2_BAO_DESI_DR2 = dot_product(vecBAO(), temp)

end function chi2_BAO_DESI_DR2

! SN_DESY5
! arXiv:2401.02929

subroutine SN_DESY5()

    implicit none

    integer :: i, j
    integer :: ios
    character(len=200) :: line

    character(len=20) :: CID(ndataSN)
    integer :: IDSURVEY(ndataSN)
    real(8) :: zCMB(ndataSN), zHD(ndataSN), zHEL(ndataSN), MU(ndataSN), MUERR(ndataSN)
    real(8) :: Cij(ndataSN, ndataSN), invCij(ndataSN, ndataSN)

    open (unit = 11, file = './data/DESY5_SN/DES-SN5YR_HD.txt', status = 'old')

    read(11,'(A)', iostat=ios) line

    do i = 1, ndataSN
        read(11, *) CID(i), IDSURVEY(i), zCMB(i), zHD(i), zHEL(i), MU(i), MUERR(i)
    end do

    close(11)

    ! We use this to prevent compiler error.
    do i = 1, ndataSN
        dataSN(i, 1) = 0.0d0
        dataSN(i, 2) = 0.0d0
        dataSN(i, 3) = zCMB(i)
        dataSN(i, 4) = zHD(i)
        dataSN(i, 5) = zHEL(i)
        dataSN(i, 6) = MU(i)
        dataSN(i, 7) = MUERR(i)
    end do

    open(11, file = "./data/DESY5_SN/covsys_000.txt", status='old')

        read(11, *) ! skip dimension line

        do i = 1, ndataSN
            read(11, *) (Dij(i,j), j = 1, ndataSN)
        end do

    close(11)

    do i = 1, ndataSN
        do j = 1, ndataSN
            Cij(i,j) = Dij(i,j)
        end do
        Cij(i,i) = Cij(i,i) + dataSN(i, 7)**2
        Id(i) = 1.0d0
    end do

    CijSN = Cij

    invCij = inv(Cij)

    invCijSN = invCij

end subroutine SN_DESY5

function chi2_SN_DESY5()

    real(8) :: chi2_SN_DESY5

    real(8) :: z, mu_obs, mu_model
    real(8) :: delta_m(ndataSN)
    integer :: i, j

    real(8) :: A

    do i = 1, ndataSN
        z = dataSN(i, 3)
        mu_obs = dataSN(i, 6)
        mu_model = 5.d0*log10(DL_DESY5(z)) + 25.0d0
        delta_m(i) = mu_obs - mu_model
    end do

    A = 0.0d0
    do i = 1, ndataSN
        do j = 1, ndataSN
            A = A + delta_m(i)*invCijSN(i,j)*delta_m(j)
        end do
    end do

    chi2_SN_DESY5 = A

end function chi2_SN_DESY5

function M0_SN()

    real(8) :: M0_SN

    real(8) :: B, E
    real(8) :: z, mu_model, flux
    real(8) :: delta_m(ndataSN: ndataSN)
    integer :: i, j

    do i = 1, ndataSN
        z = dataSN(i, 3)
        flux = dataSN(i, 5)
        mu_model = 5.d0 * log10((1 + dataSN(i, 4))/(1 + z) * DL(z)) + 25.0d0
        delta_m(i) = flux - mu_model
    end do

    B = 0.d0
    E = 0.d0
    do i = 1, ndataSN
        do j = 1, ndataSN
            B = B + delta_m(i) * invCijSN(i,j) * Id(j)
            E = E + Id(i) * invCijSN(i,j) * Id(j)
        end do
    end do

    M0_SN = B/E

end function M0_SN

function DL_DESY5(z)

    implicit none

    real(8) :: DL_DESY5
    real(8) :: z

    if (ok == 0.0d0) then
        DL_DESY5 = (c/(100.0*h))*(z + 1.0d0)* &
        rombint(DL_DESY5_integrand, 0.0d0, z, 1.0d-6)
    else if (ok > 0.0d0) then
        DL_DESY5 = (c/(100.0*h))*(z + 1.0d0)* &
        (1.0d0/sqrt(ok))*sinh(sqrt(ok)* &
        rombint(DL_DESY5_integrand, 0.0d0, z, 1.0d-6))
    else
        DL_DESY5 = (c/(100.0*h))*(z + 1.0d0)* &
        (1.0d0/sqrt(abs(ok)))*sin(sqrt(abs(ok))* &
        rombint(DL_DESY5_integrand, 0.0d0, z, 1.0d-6))
    end if

end function DL_DESY5

function DL_DESY5_integrand(z)

    implicit none

    real(8) :: DL_DESY5_integrand
    real(8) :: z

    DL_DESY5_integrand = 1.0d0/Ht(1.0d0/(z + 1.0d0))

end function DL_DESY5_integrand

! SN_Pantheon
! D. M. Scolnic et al., Astrophys. J. 859, 101 (2018) [arXiv:1710.00845].

subroutine SN_Pantheon()

    implicit none
    
    integer :: i, j
    integer :: ios
    character(len=200) :: line

    character(len=20) :: name_array(ndataSN)
    real(8) :: zCMB_array(ndataSN), zhel_array(ndataSN), dz_array(ndataSN)
    real(8) :: mb_array(ndataSN), dmb_array(ndataSN)
    real(8) :: x1_array(ndataSN)
    real(8) :: dx1_array(ndataSN), color_array(ndataSN), dcolor_array(ndataSN)
    real(8) :: thirdvar_array(ndataSN), d3rdvar_array(ndataSN)
    real(8) :: cov_m_s_array(ndataSN), cov_m_c_array(ndataSN)
    real(8) :: cov_s_c_array(ndataSN), set_array(ndataSN), ra_array(ndataSN)
    real(8) :: dec_array(ndataSN), biascor_array(ndataSN)

    real(8) :: Cij(ndataSN, ndataSN), invCij(ndataSN, ndataSN)

    character(len=20) :: name_temp

    open (unit = 11, file = './data/data_Pantheon/lcparam_full_long_zhel.txt', status = 'old')

    read(11,'(A)', iostat=ios) line

    dataSN = 0.0d0

    do i = 1, ndataSN
        read(11, *) name_array(i), zCMB_array(i), zhel_array(i), dz_array(i), &
        mb_array(i), dmb_array(i), x1_array(i), &
        dx1_array(i), color_array(i), dcolor_array(i), &
        thirdvar_array(i), d3rdvar_array(i), &
        cov_m_s_array(i), cov_m_c_array(i), &
        cov_s_c_array(i), set_array(i), ra_array(i), &
        dec_array(i)
    end do

    ! do i = 1, ndataSN
    !     read(11, *) name_temp
    !     read(11, *) dataSN(i, 2)
    !     read(11, *) dataSN(i, 3)
    !     read(11, *) dataSN(i, 4)
    !     read(11, *) dataSN(i, 5)
    !     read(11, *) dataSN(i, 6)
    !     read(11, *) dataSN(i, 7)
    !     read(11, *) dataSN(i, 8)
    !     read(11, *) dataSN(i, 9)
    !     read(11, *) dataSN(i, 10)
    !     read(11, *) dataSN(i, 11)
    !     read(11, *) dataSN(i, 12)
    !     read(11, *) dataSN(i, 13)
    !     read(11, *) dataSN(i, 14)
    !     read(11, *) dataSN(i, 15)
    !     read(11, *) dataSN(i, 16)
    !     read(11, *) dataSN(i, 17)
    !     read(11, *) dataSN(i, 18)
    ! end do

    close(11)

    ! We use this to prevent compiler error.    
    do i = 1, ndataSN
        dataSN(i, 2) = zCMB_array(i)
        dataSN(i, 3) = zhel_array(i)
        dataSN(i, 4) = dz_array(i)
        dataSN(i, 5) = mb_array(i)
        dataSN(i, 6) = dmb_array(i)
    end do

    open(11, file = "./data/data_Pantheon/sys_full_long.txt", status='old')

        read(11, *) ! skip dimension line
    
        do i = 1, ndataSN
            read(11, *) (Dij(i,j), j = 1, ndataSN)
        end do
    
    close(11)

    do i = 1, ndataSN
        do j = 1, ndataSN
            Cij(i,j) = Dij(i,j)
        end do
        Cij(i,i) = Cij(i,i) + dataSN(i, 6)**2
        Id(i) = 1.0d0
    end do

    CijSN = Cij

    invCij = inv(Cij)

    invCijSN = invCij

end subroutine SN_Pantheon

function chi2_SN_Pantheon()
    
    real(8) :: chi2_SN_Pantheon
    
    real(8) :: z, mu_obs, mu_model
    real(8) :: delta_m(ndataSN)
    integer :: i, j

    real(8) :: A, B, E

    do i = 1, ndataSN
        z = dataSN(i, 2)
        mu_obs = dataSN(i, 5)
        mu_model = 5.d0 * log10(((1 + dataSN(i, 3))/(1 + z)) * DL(z)) + 25.0d0
        delta_m(i) = mu_obs - mu_model
    end do

    A = 0.0d0
    B = 0.0d0
    E = 0.0d0
    do i = 1, ndataSN
        do j = 1, ndataSN
            A = A + delta_m(i) * invCijSN(i,j) * delta_m(j)
            B = B + delta_m(i) * invCijSN(i,j) * Id(j)
            E = E + Id(i) * invCijSN(i,j) * Id(j)
        end do
    end do

    chi2_SN_Pantheon = A + log(E/(2.0d0 * 3.141592653589793d0)) - (B**2)/E

end function chi2_SN_Pantheon

! H0 : SH0ES
! arXiv : 2112.04510
! dataH0 = {73.04, 1.04};

function chi2H0()

    implicit none

    real(8) :: chi2H0

    chi2H0 = (100.0d0*h - 73.04d0)**2/1.04d0**2

end function chi2H0

! BBN
! arXiv : 2401.15054
! dataBBN = {0.02218, 0.00055};

function chi2BBN()

    implicit none

    real(8) :: chi2BBN

    chi2BBN = (obh2() - 0.02218d0)**2/0.00055d0**2

end function chi2BBN

function chi2total()

    implicit none

    real(8) :: chi2total

    chi2total = chi2_SN_DESY5() + chi2_CMB_Planck2018() + chi2_BAO_DESI_DR2()

end function

function H0()

    real(8) :: H0

    H0 = 100.0*h

end function H0

function ol()

    real(8) :: ol

    ol = 1.0d0 - (ob + oc + Omegar())

end function ol

function age()

    real(8) :: age

    ! age = 977.8d0*rombint(age_integrand, 0.0d0, 1.0d6, 1.0d-6)
    age = 977.8d0*simpson(age_integrand, 0.0d0, 1.0d5, 1000000)

end function age

function age_integrand(z)

    real(8) :: age_integrand
    real(8) :: z

    age_integrand = 1.0d0/((z + 1.0d0)*HH(1.0d0/(z + 1.0d0)))

end function age_integrand

function simpson(f, a, b, n)

    real(8) :: simpson
    real(8) :: f, a, b
    integer :: n
    ! n must be even.
    integer :: i
    real(8) :: h_simpson, integral, x

    if (mod(n, 2) /= 0) then
    ! print *, 'Error: n must be even.'
    ! stop
    n = n + 1
    end if

    h_simpson = (b - a) / n
    integral = f(a) + f(b)

    do i = 1, n - 1
    x = a + i * h_simpson
    if (mod(i, 2) == 0) then
    integral = integral + 2.0d0 * f(x)
    else
    integral = integral + 4.0d0 * f(x)
    end if
    end do

    simpson = integral * h_simpson / 3.0d0

end function simpson

subroutine propose_jump(center, jumpsize, proposal)
    real(8), intent(in) :: center(nparams)
    real(8), intent(in) :: jumpsize(nparams)
    real(8), intent(out) :: proposal(nparams)
    real(8) :: randn(nparams), z
    integer :: i
    do i = 1, nparams
        call random_number(z)
        randn(i) = sqrt(-2.d0 * log(z)) * cos(2.d0 * 3.141592653d0 * z)
        proposal(i) = center(i) + randn(i) * jumpsize(i)
    end do
end subroutine propose_jump

subroutine pack_2D(array2D, array1D)
    real(8), intent(in) :: array2D(:, :)
    real(8), intent(out) :: array1D(:)
    integer :: i, j, k
    k = 1
    do i = 1, size(array2D, 1)
        do j = 1, size(array2D, 2)
            array1D(k) = array2D(i, j)
            k = k + 1
        end do
    end do
end subroutine

subroutine unpack_3D(array1D, array3D)
    real(8), intent(in) :: array1D(:)
    real(8), intent(out) :: array3D(:, :, :)
    integer :: i, j, k, idx
    idx = 1
    do i = 1, size(array3D, 1)
        do j = 1, size(array3D, 2)
            do k = 1, size(array3D, 3)
                array3D(i, j, k) = array1D(idx)
                idx = idx + 1
            end do
        end do
    end do
end subroutine

end program CosmologyMCMC
