
        ! initpackage burstatmos lmodel.dat /home/UTU/jmjkuu/heasoft-6.21/local
        ! lmod burstatmos /home/UTU/jmjkuu/heasoft-6.21/local

        ! TODO:
        !  1. The paths for the datafiles for atmosphere models and for fluxes, ells and teffs
        !     are now given explicitly. The could should read them automatically from the same 
        !     directory where this script is located
        !  2. The file handler integers should be determined automatically


        subroutine burstatmos_READ(ear, ne2, param, ifl, photar)

            implicit none

            integer, intent(in) :: ne2, ifl
            real, intent(in) :: ear(0:ne2)
            real, intent(in), dimension(5) :: param
            real, intent(out), dimension(ne2) :: photar

            integer, parameter :: Ne = 300
            integer, parameter :: Ng = 50
            integer, parameter :: Nlg = 3
            integer, parameter :: Na = 3

            real, dimension(Ne, Ng, Nlg, Na) :: flx
            real, dimension(Ng, Nlg, Na) :: ells
            real, dimension(Ng, Nlg, Na) :: tefs

            real, dimension(3) :: As, gs
            real, dimension(Ng) :: gradg_grid
            real, dimension(Ne) :: eteff_grid, fluxes, energy
            real :: grg, logg, comp, ell, teff, M, R, D, z
            real :: flux_low, flux_mid, flux_high, dist_const, flux_int
            real :: ener0, ener1, ener2, ergkev, constants
            integer, dimension(4) :: dimm = (/Ne, Ng, Nlg, Na/)    
            integer, dimension(3) :: dimm2 = (/Ng, Nlg, Na/)    
            integer :: i

            !write(*,*) ifl

            ! Fit parameters
            ! param = g_rad/g, hydrogen mass fraction, mass, radius, distance
            ! mass in solar units, radius in km, distance in units of 10 kpc
            grg = param(1)
            comp = param(2)
            M = param(3)
            R = param(4)
            D = param(5)

            write(*,*) param

            call calculate_logg(M, R, logg)
            call calculate_redshift(M, R, z)


            ! Construct the flux etc arrays from the models for the given parameters
            call make_atms(eteff_grid, gradg_grid, gs, As, &
             & flx, ells, tefs)

            ! Interpolate the flux in the given energy range
            do i = 1, Ne
                call interp4d(dimm, eteff_grid, gradg_grid, gs, As, &
                 & flx, eteff_grid(i), grg, logg, comp, fluxes(i))
            enddo

            ! interpolate the ell and teff
            call interp3d(dimm2, gradg_grid, gs, As, ells, grg, &
             & logg, comp, ell)
            call interp3d(dimm2, gradg_grid, gs, As, tefs, grg, &
             & logg, comp, teff)
            energy = eteff_grid * teff


            ! kpc to km
            ! D is in units of kpc, R in units of km
            dist_const = 3.086e16
            ! erg to keV
            ergkev = 6.2415e8
            ! Combine constants
            constants = ergkev*(R/(D*dist_const))**2
            
            i = 1
            do i=1, ne2

                ener0 = ear(i-1)*z
                ener2 = ear(i)*z
                ener1 = (ener0 + ener2)/2.0

                ! interpolate the flux for the given energies
                call interp1d(Ne, energy, fluxes, ener0, flux_low)
                call interp1d(Ne, energy, fluxes, ener1, flux_mid)
                call interp1d(Ne, energy, fluxes, ener2, flux_high)

                ! Convert the energy flux at NS surface to the comoving photon flux at distance D
                ! N_e = F'(E*redshift) * erg_to_keV * (R/D)^2 / (E * redshift)
                ! here redshift = z as in code = 1/(1+z) 


                flux_low = 10**(flux_low) * constants / ener0
                flux_mid = 10**(flux_mid) * constants / ener1
                flux_high = 10**(flux_high) * constants / ener2


                ! integrate over the bin size using simpsons rule
                flux_int = flux_low + flux_mid*4.0 + flux_high

                photar(i) = ((ear(i)-ear(i-1))/6.0)*flux_int

            enddo

            return
        end subroutine burstatmos





        subroutine make_atms(eteff_grid, gradg_grid, gs, As, &
            & flx, ells, tefs)

        implicit none
        integer, parameter :: Ne = 300
        integer, parameter :: Ng = 50
        integer, parameter :: Nlg = 3
        integer, parameter :: Na = 3

        real, dimension(Ne, Ng, Nlg, Na) :: flx
        real, dimension(Ng, Nlg, Na) :: ells
        real, dimension(Ng, Nlg, Na) :: tefs

        real, dimension(3) :: As, gs
        real, dimension(Ng) :: gradg_grid
        real, dimension(Ne) :: eteff_grid
        integer :: j, i, k, l
        character(128) :: ell_file, teff_file, flux_file


        As = (/ 0.0, 0.738, 1.0 /)    ! Array of model hydrogen fractions
        gs = (/ 14.0, 14.3, 14.6 /)   ! Array of model log g values
 
        ! Create a uniform g_rad/g grid of 50 values between 0.1 and 1.0
        do j = 1, Ng
            gradg_grid(j) = 0.1 + (j-1)*0.9/49.0
        enddo

        ! Create a uniform (in log scale) eteff grid of 300 values between -2 and 2
        do j = 1, Ne
            eteff_grid(j) = 10**(-2.0 + (j-1)*4/299.0)
        enddo


        ell_file = "/home/UTU/jmjkuu/heasoft-6.21/local/ells.dat"
        teff_file = "/home/UTU/jmjkuu/heasoft-6.21/local/teffs.dat"
        flux_file = "/home/UTU/jmjkuu/heasoft-6.21/local/fluxes.dat"

        open(unit=10, file=ell_file, status='old', form='formatted')
        open(unit=11, file=teff_file, status='old', form='formatted')
        open(unit=12, file=flux_file, status='old', form='formatted')
        do i=1, 3, 1
            do j=1, 3, 1
                do k=1, 50, 1
                    read(10,"(F10.8)") ells(k, j, i)
                    read(11,"(F10.8)") tefs(k, j, i)
                    do l=1, 300, 1
                        read(12,"(F10.7)") flx(l, k, j, i)
                    enddo
                enddo
            enddo
        enddo
        close(10)
        close(11)
        close(12)

        ! Arrays returned:
        ! eteff_grid, grad_grid, gs, As
        ! flx, ells, tefs

        return

        end subroutine make_atms





        subroutine calculate_redshift(M, R, z)
            ! Calculate redshift based on given mass and radius
            ! Mass in units of solar mass and radius in km
            ! here z means 1/(1+z) 
            implicit none
            real, intent(in) :: M, R
            real, intent(out) :: z
            real :: M_sol = 1.9891e30
            real :: G = 6.67408e-11
            real :: c = 299792458
            real :: rs
            rs = 2*G*M*M_sol/c**2
            z = 1.0/sqrt(1.0 - rs/(R*1e3))
            return
        end subroutine calculate_redshift


        subroutine calculate_logg(M, R, logg)
            ! Calculate log of surface gravity based on given mass and radius
            ! Mass in units of solar mass and radius in km            
            implicit none
            real, intent(in) :: M, R
            real, intent(out) :: logg
            real :: M_sol = 1.9891e30
            real :: G = 6.67408e-11
            real :: z
            call calculate_redshift(M, R, z)
            logg = G * M * M_sol * z / (R*1e3)**2
            logg = log10(logg*100.0)
            return
        end subroutine calculate_logg




        subroutine read_atms_files(ich, ig, models, energ, lums, fluxes)
            implicit none
            integer, intent(in) :: ich, ig
            integer :: lf, var1, j, k
            integer, intent(in) :: models
            character(128) :: datafile
            integer, parameter :: lines = 5400    ! Number of lines in datafile
            real, dimension(29) :: vars
            real, intent(out), dimension(300) :: energ
            real, intent(out), dimension(models, 300) :: fluxes
            real, intent(out), dimension(models) :: lums
            lums(1:26) = (/0.001,0.003,0.01,0.03,0.05,0.07, &
                            0.1,0.15,0.2,0.3,0.4,0.5,0.55,0.6,0.65,0.7,& 
                            0.75,0.8,0.85,0.9,0.95,0.98,1.00,1.02,1.04,1.06/)
            if (ig .ge. 1) lums(27) = 1.08
            if (ig .eq. 2) lums(28) = 1.10
            lf = ich*3 + ig + 1
            j = 1

            datafile = "/home/UTU/jmjkuu/heasoft-6.21/local/J_A+A_545_A120/tabled2.csv"

            open(unit=11, file=datafile, status='old')
            do k=1, lines
                read(11,*) var1, vars(1:29)
                if(var1 .eq. lf) then
                    energ(j) = vars(1)
                    fluxes(:,j) = vars(2:models+1)
                    j = j + 1
                endif
            enddo
            close(11)
            return
        end subroutine


        subroutine readmo(ich, ig, models, lledd, gradg, fcol1, norm1, teff)
            implicit none
            integer, intent(in) :: ich, ig
            integer :: lf, j, k, var1
            integer, intent(in) :: models
            character(128) :: datafile
            integer, parameter :: lines = 484    ! Number of lines in datafile
            real, intent(out), dimension(models) :: lledd, gradg, teff, norm1, fcol1
            real :: var2, var3, var4, var5, var6, var7, var8, var9, var10, var11, var12
            real, dimension(11) :: vars
            lf = ich*3 + ig + 1
            j = 1

            datafile = "/home/UTU/jmjkuu/heasoft-6.21/local/J_A+A_545_A120/tabled1.dat"

            open(unit=10, file=datafile, status='old')
            do k=1, lines
                read(10,*) var1, var2, var3, var4, var5, var6, var7, var8, var9, var10, var11, var12
                read(10,*) var1, vars(1:11)
                if (var1 .eq. lf) then
                    lledd(j) = vars(1)
                    gradg(j) = vars(2)
                    teff(j) = vars(3)
                    fcol1(j) = vars(4)
                    norm1(j) = vars(9) / vars(4)**4.0
                    j = j + 1
                endif
            enddo
            close(10)
            return
        end subroutine
        


        subroutine interp1d(dimm, Xgrid, Ygrid, X0, Y)

                ! 1D Piecewise multilinear interpolator
                ! See Weiser & Zarantonello (1988) and subroutine interp4d

                implicit none
                integer, intent(in) :: dimm
                real, intent(in), dimension(dimm) :: Xgrid, Ygrid
                real, intent(in) :: X0
                real, intent(out) :: Y
                real :: x
                integer :: a       

                call interk(X0, Xgrid, dimm, a, x)

                if (abs(x) .lt. 0.001) then
                    Y = Ygrid(a)
                else
                    Y = Ygrid(a) + x*(Ygrid(a+1) - Ygrid(a))
                endif

                return
        end subroutine interp1d





        subroutine interp2d(dimm, Xgrid, Ygrid, Zgrid, X0, Y0, Z)

                ! 2D Piecewise multilinear interpolator
                ! See Weiser & Zarantonello (1988) and subroutine interp4d

                implicit none
                integer, intent(in), dimension(2) :: dimm
                real, intent(in), dimension(dimm(1)) :: Xgrid
                real, intent(in), dimension(dimm(2)) :: Ygrid
                real, intent(in), dimension(dimm(1), dimm(2)) :: Zgrid 
                real, intent(in) :: X0, Y0
                real, intent(out) :: Z
                real :: x, y, Gx0, Gx1
                integer :: a, b, nx, ny

                nx = dimm(1)
                ny = dimm(2)

                call interk(X0, Xgrid, nx, a, x)
                call interk(Y0, Ygrid, ny, b, y)

                if (abs(x) .lt. 0.001) then
                    Gx0 = Zgrid(a, b)
                    Gx1 = Zgrid(a, b+1)
                else
                    Gx0 = Zgrid(a, b) + x*(Zgrid(a+1, b) - Zgrid(a, b))
                    Gx1 = Zgrid(a, b+1) + x*(Zgrid(a+1, b+1) - Zgrid(a, b+1))
                endif

                Z = Gx0 + y*(Gx1 - Gx0)

                return
        end subroutine interp2d



        subroutine interp3d(dimm, Xgrid, Ygrid, Zgrid, Rgrid, X0, Y0, Z0, R)

                ! 3D Piecewise multilinear interpolator
                ! See Weiser & Zarantonello (1988) and subroutine interp4d

                implicit none

                integer, intent(in), dimension(3) :: dimm
                real, intent(in), dimension(dimm(1)) :: Xgrid
                real, intent(in), dimension(dimm(2)) :: Ygrid
                real, intent(in), dimension(dimm(3)) :: Zgrid
                real, intent(in), dimension(dimm(1), dimm(2), dimm(3)) :: Rgrid
                real, intent(in) :: X0, Y0, Z0
                real, intent(out) :: R
                real :: x, y, z
                real :: Gx11, Gx10, Gx01, Gx00, Gxy1, Gxy0
                integer :: a, b, c, nx, ny, nz

                nx = dimm(1)
                ny = dimm(2)
                nz = dimm(3)

                call interk(X0, Xgrid, nx, a, x)
                call interk(Y0, Ygrid, ny, b, y)
                call interk(Z0, Zgrid, nz, c, z)


                if (abs(x) .lt. 0.001) then
                    Gx11 = Rgrid(a, b+1, c+1)
                    Gx01 = Rgrid(a, b, c+1)
                    Gx10 = Rgrid(a, b+1, c)
                    Gx00 = Rgrid(a, b, c)
                else
                    Gx11 = Rgrid(a, b+1, c+1) + x*(Rgrid(a+1, b+1, c+1) - Rgrid(a, b+1, c+1))
                    Gx01 = Rgrid(a, b, c+1) + x*(Rgrid(a+1, b, c+1) - Rgrid(a, b, c+1))
                    Gx10 = Rgrid(a, b+1, c) + x*(Rgrid(a+1, b+1, c) - Rgrid(a, b+1, c))
                    Gx00 = Rgrid(a, b, c) + x*(Rgrid(a+1, b, c) - Rgrid(a, b, c))
                endif

                Gxy1 = Gx01 + y*(Gx11 - Gx01)
                Gxy0 = Gx00 + y*(Gx10 - Gx00)
                R = Gxy0 + z*(Gxy1 - Gxy0)

                return
        end subroutine interp3d






        subroutine interp4d(dimm, Xgrid, Ygrid, Zgrid, Wgrid, Rgrid, X0, Y0, Z0, W0, R)
                
                ! 4D Piecewise multilinear interpolator
                ! the interpolant F_m can be calculated recursively
                ! F_m = G(y_1, ..., y_N) = G(y_1, ..., y_[N-1], 0) + y_N*{G(y_1, ..., y_[N-1], 1) - G(y_1, ..., y_[N-1], 0)}
                ! See Weiser & Zarantonello (1988)

                implicit none

                integer, intent(in), dimension(4) :: dimm
                real, intent(in), dimension(dimm(1)) :: Xgrid
                real, intent(in), dimension(dimm(2)) :: Ygrid
                real, intent(in), dimension(dimm(3)) :: Zgrid
                real, intent(in), dimension(dimm(4)) :: Wgrid
                real, intent(in), dimension(dimm(1), dimm(2), dimm(3), dimm(4)) :: Rgrid
                real, intent(in) :: X0, Y0, Z0, W0
                real, intent(out) :: R
                real :: x, y, z, w
                real :: Gx111, Gx011, Gx101, Gx001, Gx110, Gx010, Gx100, Gx000
                real :: Gxy11, Gxy01, Gxy10, Gxy00, Gxyz1, Gxyz0
                integer :: a, b, c, d, nx, ny, nz, nw

                nx = dimm(1)
                ny = dimm(2)
                nz = dimm(3)
                nw = dimm(4)

                ! Find the grid points Xgrid(i) and Xgrid(i+1) so that Xgrid(i) < X0 < Xgrid(i+1)
                ! a is the index i of Xgrid(i), x is the distance of X0 from Xgrid(i) when Xgrid(i)=0 and Xgrid(i+1)=1
                call interk(X0, Xgrid, nx, a, x)
                call interk(Y0, Ygrid, ny, b, y)
                call interk(Z0, Zgrid, nz, c, z)
                call interk(W0, Wgrid, nw, d, w)
 

                ! Interpolate in the x-direction
                if (abs(x) .lt. 0.001) then 
                    Gx111 = Rgrid(a, b+1, c+1, d+1)
                    Gx011 = Rgrid(a, b, c+1, d+1)
                    Gx101 = Rgrid(a, b+1, c, d+1)
                    Gx001 = Rgrid(a, b, c, d+1)
                    Gx110 = Rgrid(a, b+1, c+1, d)
                    Gx010 = Rgrid(a, b, c+1, d)
                    Gx100 = Rgrid(a, b+1, c, d)
                    Gx000 = Rgrid(a, b, c, d)
                else
                    Gx111 = Rgrid(a, b+1, c+1, d+1) + x*(Rgrid(a+1, b+1, c+1, d+1) - Rgrid(a, b+1, c+1, d+1))
                    Gx011 = Rgrid(a, b, c+1, d+1) + x*(Rgrid(a+1, b, c+1, d+1) - Rgrid(a, b, c+1, d+1))
                    Gx101 = Rgrid(a, b+1, c, d+1) + x*(Rgrid(a+1, b+1, c, d+1) - Rgrid(a, b+1, c, d+1))
                    Gx001 = Rgrid(a, b, c, d+1) + x*(Rgrid(a+1, b, c, d+1) - Rgrid(a, b, c, d+1))
                    Gx110 = Rgrid(a, b+1, c+1, d) + x*(Rgrid(a+1, b+1, c+1, d) - Rgrid(a, b+1, c+1, d))
                    Gx010 = Rgrid(a, b, c+1, d) + x*(Rgrid(a+1, b, c+1, d) - Rgrid(a, b, c+1, d))
                    Gx100 = Rgrid(a, b+1, c, d) + x*(Rgrid(a+1, b+1, c, d) - Rgrid(a, b+1, c, d))
                    Gx000 = Rgrid(a, b, c, d) + x*(Rgrid(a+1, b, c, d) - Rgrid(a, b, c, d))
                endif

                ! Interpolate in the y-direction
                Gxy11 = Gx011 + y*(Gx111 - Gx011)
                Gxy01 = Gx001 + y*(Gx101 - Gx001)
                Gxy10 = Gx010 + y*(Gx110 - Gx010)
                Gxy00 = Gx000 + y*(Gx100 - Gx000)

                ! Interpolate in the z-direction
                Gxyz1 = Gxy01 + z*(Gxy11 - Gxy01)
                Gxyz0 = Gxy00 + z*(Gxy10 - Gxy00)

                ! Interpolate in the w-direction
                R = Gxyz0 + w*(Gxyz1 - Gxyz0)

                return
        end subroutine interp4d



        subroutine interk(X0, Xgrid, nx, a, x)

            implicit none
            integer, intent(in) :: nx
            real, intent(in) :: X0
            real, intent(in), dimension(nx) :: Xgrid
            integer, intent(out) :: a
            real, intent(out) :: x
            real :: dx = 0.0
            integer :: mid, start, finish, rangee

            a = 0

            ! If X0 is outside the grid, then do linear extrapolation
            if (X0 .lt. Xgrid(1)) then
                a = 1
                dx = abs(Xgrid(1) - Xgrid(2))
                x = (X0 - Xgrid(1))/dx
            else if (X0 .gt. Xgrid(nx)) then
                a = nx-1
                dx = abs(Xgrid(nx) - Xgrid(nx-1))
                x = (X0 - Xgrid(nx-1))/dx

            ! Otherwise interpolate
            else
                start = 1
                finish = nx
                rangee = finish - start
                mid = (start + finish)/2
                do while ( rangee > 1 .and. Xgrid(mid) /= X0 )
                    if (X0 > Xgrid(mid)) then
                        start = mid
                    else
                        finish = mid
                    end if
                    rangee = finish - start
                    mid = (start + finish)/2
                enddo

                a = mid
                dx = Xgrid(mid+1) - Xgrid(mid)
                x = (X0 - Xgrid(mid))/dx
                
            endif

            return
        end subroutine interk


