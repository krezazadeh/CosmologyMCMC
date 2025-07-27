program test_background

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

    integer, parameter :: n_array_background = 1000000
    real(8), dimension(n_array_background) :: arrayNe, arraya, arrayHt
    real(8), dimension(n_array_background) :: arrayphit, arrayphitp
    integer :: imax_array_background
    real(8) :: Omegam0, Omegar0
    real(8) :: ainitial
    real(8) :: Vt0

    integer, parameter :: n_data = 100000
    integer :: i

    real(8) :: h, ob, oc, ok, log10phitinitial

    real(8) :: a

    open(unit = 11, file = "test_background.txt", status = 'replace')

    h = 0.68d0
    ob = 0.04d0
    oc = 0.26d0
    ok = 0.0d0

!     log10phitinitial = -4.0d0
!     log10phitinitial = -0.46d0
    log10phitinitial = -3.0d0

    call compute_background()

    do i = 1, n_data

        a = 1.0d-6 + (1.0d0 - 1.0d-8)*real(i - 1, 8)/real(n_data - 1, 8)

        write(11, "(2e25.16)") a, Ht(a)

    end do

    close(11)

contains

subroutine compute_background()

real(8) :: Vt0old
real(8) :: phitinitial

real(8) :: Neinitial, Htinitial, phitpinitial
real(8) :: Ne0

real(8) :: phit0, phitp0

real(8) :: dNe

real(8) :: Nei, Hti, phiti, phitpi
real(8) :: Neim1, Htim1, phitim1, phitpim1
real(8) :: kphit1, kphit2, kphit3, kphit4
real(8) :: kphitp1, kphitp2, kphitp3, kphitp4

integer :: i

Omegam0 = om()
Omegar0 = Omegar()

ainitial = 1.0d-6

Neinitial = 0.0d0
Ne0 = -log(ainitial)
dNe = 1.0d-4

phitinitial = 10.0d0**log10phitinitial

Vt0old = 0.0d0

Vt0 = (6.0d0 - 6.0d0*ok - 6.0d0*Omegam0 - 6.0d0*Omegar0)/ &
(2.0d0*fphit(phitinitial))

do while (abs((Vt0 - Vt0old)/Vt0) >= 1.0d-3)

phitpinitial = -((ainitial**4*Vt0*dfphit(phitinitial))/ &
(3.0d0*ainitial**2*ok + 3.0d0*ainitial*Omegam0 + &
3.0d0*Omegar0 + ainitial**4*Vt0*fphit(phitinitial)))

i = 1

Nei = Neinitial
phiti = phitinitial
phitpi = phitpinitial
Hti = Ht_Ne(Nei, phiti, phitpi)

arrayNe(i) = Nei
arrayHt(i) = Hti

arrayphit(i) = phiti
arrayphitp(i) = phitpi

do while (Hti >= 1.0d0)

    i = i + 1

    Neim1 = Nei
    phitim1 = phiti
    phitpim1 = phitpi

    kphit1 = dNe*dphit(Neim1, phitim1, phitpim1)
    kphitp1 = dNe*dphitp(Neim1, phitim1, phitpim1)

    kphit2 = dNe*dphit(Neim1 + dNe/2.0d0, phitim1 + kphit1/2.0d0, phitpim1 + kphitp1/2.0d0)
    kphitp2 = dNe*dphitp(Neim1 + dNe/2.0d0, phitim1 + kphit1/2.0d0, phitpim1 + kphitp1/2.0d0)

    kphit3 = dNe*dphit(Neim1 + dNe/2.0d0, phitim1 + kphit2/2.0d0, phitpim1 + kphitp2/2.0d0)
    kphitp3 = dNe*dphitp(Neim1 + dNe/2.0d0, phitim1 + kphit2/2.0d0, phitpim1 + kphitp2/2.0d0)

    kphit4 = dNe*dphit(Neim1 + dNe, phitim1 + kphit3, phitpim1 + kphitp3)
    kphitp4 = dNe*dphitp(Neim1 + dNe, phitim1 + kphit3, phitpim1 + kphitp3)

    Nei = Neim1 + dNe
    phiti = phitim1 + (1.0d0/6.0d0)*(kphit1 + 2.0d0*kphit2 + 2.0d0*kphit3 + kphit4)
    phitpi = phitpim1 + (1.0d0/6.0d0)*(kphitp1 + 2.0d0*kphitp2 + 2.0d0*kphitp3 + kphitp4)

    Hti = Ht_Ne(Nei, phiti, phitpi)

    arrayNe(i) = Nei
    arrayHt(i) = Hti

    arrayphit(i) = phiti
    arrayphitp(i) = phitpi

end do ! while (Hti >= 1.0d0)

phit0 = phiti
phitp0 = phitpi

Vt0old = Vt0

Vt0 = (6.0d0 - 6.0d0*ok - 6.0d0*Omegam0 - 6.0d0*Omegar0 - phitp0**2)/ &
(2.0d0*fphit(phit0))

! print *, "Vt0", Vt0old, Vt0, abs((Vt0 - Vt0old)/Vt0)
! print *, "Ne", Nei, Ne0, Nei - Ne0

end do ! while (abs((Vt0 - Vt0old)/Vt0) >= 1.0d-3)

imax_array_background = i

! store background
open(unit = 12, file = "compute_background.txt", status = 'replace')

do i = 1, imax_array_background

arraya(i) = exp(arrayNe(i) - Nei)
arrayHt(i) = arrayHt(i)/Hti

! store background
write(12, "(5e25.16)") arraya(i), arrayHt(i), arrayphit(i), arrayphitp(i), arrayphitp(i)**2/6.0d0

end do

! store background
close(12)

end subroutine compute_background

function dphit(Ne, phit, phitp)

    real(8) :: dphit
    real(8) :: Ne, phit, phitp

    dphit = phitp

end function dphit

function dphitp(Ne, phit, phitp)

    real(8) :: dphitp
    real(8) :: Ne, phit, phitp

    dphitp = -3.0q0*phitp + (phitp* &
    (2.0q0*ainitial**2*exp(2.0q0*Ne)*ok + &
    3.0q0*ainitial*exp(Ne)*Omegam0 + 4.0q0*Omegar0 + &
    ainitial**4*exp(4.0q0*Ne)*Ht_Ne(Ne, phit, phitp)**2*phitp**2))/ &
    (2.0q0*ainitial**4*exp(4.0q0*Ne)*Ht_Ne(Ne, phit, phitp)**2) - &
    (Vt0*dfphit(phit))/Ht_Ne(Ne, phit, phitp)**2

end function dphitp

function Ht_Ne(Ne, phit, phitp)

    real(8) :: Ht_Ne
    real(8) :: Ne, phit, phitp

    Ht_Ne = sqrt((6.0d0*(ainitial**2*exp(2.0d0*Ne)*ok + &
    ainitial*exp(Ne)*Omegam0 + Omegar0))/ &
    exp(4.0d0*Ne) + 2.0d0*ainitial**4*Vt0*fphit(phit))/ &
    (ainitial**2*sqrt(6.0d0 - phitp**2))

end function Ht_Ne

function Ht(ainput)

real(8) :: Ht
real(8) :: ainput

integer :: leftpoint, rightpoint, midpoint
integer :: i

if (ainput <= arraya(1)) then

! Radiation dominated: Ht*a**2 = const
Ht = arrayHt(1)*arraya(1)**2/ainput**2

else if ((ainput > arraya(1)) .and. (ainput <= arraya(imax_array_background))) then

! bisection method
leftpoint = 1
rightpoint = imax_array_background
do while ((rightpoint - leftpoint) > 1)
   midpoint = (rightpoint + leftpoint)/2
   if (arraya(midpoint) > ainput) then
      rightpoint = midpoint
   else
      leftpoint = midpoint
   endif
enddo
i = leftpoint

Ht = (ainput*arrayHt(i) - arraya(i + 1)*arrayHt(i) - ainput*arrayHt(i + 1) + &
arraya(i)*arrayHt(i + 1))/(arraya(i) - arraya(i + 1))

else if (ainput > arraya(imax_array_background)) then

Ht = ((ainput - arraya(imax_array_background))* &
arrayHt(imax_array_background - 1) + &
(-ainput + arraya(imax_array_background - 1))* &
arrayHt(imax_array_background))/ &
(arraya(imax_array_background - 1) - arraya(imax_array_background))

end if

end function Ht

function fphit(phit)

    real(8) :: fphit
    real(8) :: phit

    fphit = 1.0d0 - phit**2

end function fphit

function dfphit(phit)

    real(8) :: dfphit
    real(8) :: phit

    dfphit = -2.0d0*phit

end function dfphit

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

function HH(a)

    real(8) :: HH
    real(8), intent(in) :: a

    HH = 100.0*h*Ht(a)

end function HH

function om()

    implicit none

    real(8) :: om

    om = ob + oc

end function om

end program test_background
