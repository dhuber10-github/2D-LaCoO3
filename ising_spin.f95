!Versionskommentar
!3 Spinzustände für tau: -1 | 0 | 1


module global_variables
	integer, parameter :: NMax=48
	double precision :: acc1, acc0, accm1, mu
	double precision, parameter :: JJ=-1.74d0, V1=64.d0, V2=2.06d0
	!double precision, parameter :: JJ=0d0, V1=64d0, V2=0d0
end module global_variables



Program allinone
	use global_variables
	implicit none

	double precision CalcE
	external CalcE
	double precision NNEnergy
	external NNEnergy

	integer :: mag(Nmax,Nmax)=0			!/* 2D Ising Lattice */
	integer :: i, j, k, m, n			!/* Loop counters */
	integer :: s, d					!/* Lattice spin variables */
	integer :: a, b, c				!L/3 und 2*L/3 und L/2
	integer :: counter=0
	double precision :: Energy=0d0			!/* Total lattice energy */
	double precision :: spinsum=0d0			!Mittelwertbildung pro sweep der Summe über mag(:,:)
	double precision :: T=0d0			!/* temperature (in units of J/k_B) */
	double precision :: beta=0d0
	double precision :: start, finish
	double precision :: tstep=0.001d0			!Temperaturschrittweite zum dynamischen Sampling
	double precision, parameter :: tmax=200d0		!Temperatur bei der die Simulation beginnt
	double precision, parameter :: mumax=3d0	!mumax-0.1 ist der maximale wert von mu in der Simulation, da auch 0
	integer, parameter :: sweeps=1000		!/* number of measurement sweeps */
	integer, parameter :: warm=10d4			!/* number of warm-up sweeps */
	integer, parameter :: L=48			!/* lattice dimension */

	!Fouriertransformation
	double precision :: fsum=0d0, peaksum=0d0, centerpeak=0d0	!Summe über f(:,:), Summe der Peaks, für Schachbrettmuster
	double precision :: wsave(4*L+15)
	double precision :: f(L,L)=0d0, back(L,L)=0d0, fsym(L,L)=0d0, fabs2sum(L,L)=0d0
	double complex :: fcomplex(L,L)=(0d0,0d0), f1d(L,L)=(0d0,0d0), fcomplexsum(L,L)=(0d0,0d0)
	double complex :: column(L)=(0d0,0d0), line(L)=(0d0,0d0)



 ! /***************************
 !  * Initialization          * 
 !  ***************************/
	call cpu_time(start)
	call srand(time())
	call zffti(L, wsave)

	a=int(L/3)					!L/3...für die Peakwerte
	b=int(2*L/3)					!2*L/3...für die Peakwerte
	c=int(L/2)					!L/2...für die Peakwerte
	T=tmax
	mu=0d0

	!Standardwerte: L=48, fft: abs2, Mittelwertbildung MW
!	open(10, file='data_w10e3s10e3_T10-01_mu0.dat')
!	open(20, file='fft_w10e3s10e3_T10-01_mu0.dat')
!	open(30, file='fftabs2b_w10e3s10e3_T10-01_mu0.dat')
!	open(40, file='fftsym_w10e3s10e3_T10-01_mu-10-3.dat')
	open(50, file='snapshots_ising_spin.dat')



!	do while (mu > -10.09d0)				!Bereich für chemisches Potential

	call random_init(mag,L)				!Start jedes unterschiedlichen mu mit zufälligem Anfangsgitter
	call outputmag(mag,L)

	do while (T > 0.791d0)				!Bereich für Temperatur

!		if(T >= (1.0) .and. T <= (1.5)) then

!			if(T >= (1.15) .and. T <= (1.25)) then
!				tstep=0.001
!			else
!				tstep=0.01
!			end if

!		else
!			tstep=0.1
!		end if
		beta=1.d0/T


 
! /***************************
! /* warum up sweeps */
! /***************************
	acc1=0d0
	acc0=0d0
	accm1=0d0
	do i=1, warm
		call sweep(mag,L,beta)
	end do
	call outputmag(mag,L)
	write(*,*) "1 =", acc1/(L**2 * warm)
	write(*,*) "0 =", acc0/(L**2 * warm)
	write(*,*) "-1 =", accm1/(L**2 * warm)
	stop
 ! /***************************
! /* end warum up sweeps */ 
! /***************************

 

! /***************************
! /*  sweeps */
! /***************************
	do i=1, sweeps

!		if(i==sweeps .and. T<=0.1) then
!			call outputmag(mag,L)
!			counter=counter+1
!		end if

		call sweep(mag,L,beta)
		spinsum=spinsum+sum(mag)

	!++++++++++++++++++++++++++++++++++++++
	!FFT
	!++++++++++++++++++++++++++++++++++++++
		do n=1, L				!1D-FFT of columns
			column=mag(n,:)
			call zfftf(L, column, wsave)
			f1d(n,:)=column
		end do

		do m=1, L				!Fouriertransformierte Spalten werden zeilenweise fouriertransformiert
			line=f1d(:,m)
			call zfftf(L, line, wsave)
			fcomplex(:,m)=line
		end do
	!++++++++++++++++++++++++++++++++++++++
	!ende FFT
	!++++++++++++++++++++++++++++++++++++++

		fcomplex=fcomplex/L
		f(:,:) = (abs(fcomplex(:,:)))**2
		fabs2sum = fabs2sum + f
!		fcomplexsum=fcomplexsum+fcomplex	!Aufsummieren über anzahl der snapshots

	end do						!anzahl der sweeps
! /***************************
! /*  end sweeps */ 
! /***************************



	spinsum = spinsum/real(sweeps)/real(L)**2
	Energy = CalcE(mag,L)/real(L)**2
	fabs2sum = fabs2sum/sweeps
!	fcomplexsum = fcomplexsum/sweeps
!	f(:,:) = (abs(fcomplexsum(:,:)))**2		!Absolutbetrag der komplexen Fourierkoeffizienten
!	fsum = sum(f)/real(L)**2
!	peaksum = (f(a+1,a+1)+f(a+1,b+1)+f(b+1,a+1)+f(b+1,b+1))/real(L)**2
!	centerpeak = f(c+1,c+1)/real(L)**2

	write(10,*) mu, T, spinsum, Energy!, fsum, peaksum, centerpeak



!--------------------------------------
!Ausgabe der Fouriertransformierten
!--------------------------------------
	do m=2, L
		write(20,*) fabs2sum(:,m)
	end do

!	do n=1, L
!		do m=1, L
!			fsym(n,m)=0.5*(f(n,m)+f(L-m+1,n))
!		end do
!	end do
!	do m=2, L-1
!		write(40,*) fsym(:,m)
!	end do
!--------------------------------------
!ende Ausgabe der Fouriertransformierten
!--------------------------------------



!######################################
!FFT rücktransformation
!######################################
	do n=1, L
		column=fabs2sum(n,:)
		call zfftb(L, column, wsave)
		f1d(n,:)=column
	end do

	do m=1, L
		line=f1d(:,m)
		call zfftb(L, line, wsave)
		back(:,m)=line
	end do
!######################################
!ende FFT rücktransformation
!######################################

	back=back/L

!--------------------------------------
!Ausgabe der Rücktransformierten
!--------------------------------------
	do m=1, L
		write(30,*) back(:,m)
	end do
!--------------------------------------
!ende Ausgabe der Rücktransformierten
!--------------------------------------



		T=T-tstep
		spinsum=0d0

	end do						!Bereich für Temperatur

!		T=tmax
!		mu=mu-0.1d0

!	end do						!Bereich für chemisches Potential

	close(10)
	close(20)
	close(30)
!	close(40)
	close(50)

	write(*,*) "counter:", counter

	call cpu_time(finish)
	write(*,*) (finish-start)/60, 'minutes'

end Program allinone







subroutine random_init(mag,L)
use global_variables
implicit none
  integer mag(Nmax,Nmax), L
  double precision :: r
  integer i,j


  do i=1,L
     do j=1,L

	r = rand()

       	if (r .lt. (1d0/3d0)) then
	   mag(i,j) = -1
        else if (r .gt. (2d0/3d0)) then
           mag(i,j) = 1
	else 
	   mag(i,j) = 0
        endif
     enddo
  enddo

  return
end subroutine random_init



subroutine sweep(mag,L,beta)
use global_variables
 implicit none 
  integer mag(Nmax,Nmax), L, magtemp
  integer i,j
  double precision r, beta, option
  double precision :: Etemp=0,deltaE=0
  double precision NNEnergy

do i=1, L
	do j=1, L
		magtemp = mag(i,j)
		Etemp=NNEnergy(mag,L,i,j)
		option = rand()

		if (mag(i,j) == -1) then

			if(option .gt. 0.5) then
				mag(i,j) = 1
			else
				mag(i,j) = 0
			end if

		else if (mag(i,j) == 0) then

			if(option .gt. 0.5) then
				mag(i,j) = 1
			else
				mag(i,j) = -1
			end if

		else

			if(option .gt. 0.5) then
				mag(i,j) = 0
			else
				mag(i,j) = -1
			end if

		end if


		deltaE=NNEnergy(mag,L,i,j)-Etemp		!der energieunterschied vor und nach der spinänderung wird berechnet
!deltaE niedriger als bei ising
		r=dexp((-deltaE)*beta)				!Boltzmann Gewicht

		if (rand() > (r/(1d0+r))) then			!Wärmebad Verfahren | neue konfiguration ist energetisch günstiger
			mag(i,j)=magtemp
		end if

		if (mag(i,j) == 1) acc1 = acc1 + 1
		if (mag(i,j) == 0) acc0 = acc0 + 1
		if (mag(i,j) == -1) accm1 = accm1 + 1

	end do
end do
return
end subroutine sweep



subroutine outputmag(mag,L)
use global_variables
implicit none
  integer mag(Nmax,Nmax), L
 
  integer i,j
  write (50,*) 'Spin snapshot'

  do i=1,L
     do j=1,L
        if(mag(i,j) .eq. 1) then
	   write(50,'(a)',ADVANCE='NO') '↑ '
        else if (mag(i,j) .eq. 0) then
           write(50,'(a)',ADVANCE='NO') '0 '
	else
	   write(50,'(a)',ADVANCE='NO') '↓ '
        endif
      
     enddo
     !     Newline to complete
     write (50,*) 
!     write (50,*)
  enddo
  return 
end subroutine outputmag



double precision function CalcE(mag,L)
use global_variables
implicit none
  integer mag(Nmax,Nmax),L
  
  integer i,j
  integer m,n,o,p
  integer i1, i2, j1, j2
  double precision Energy, V1part, V2part, Jpart

  !/* Determine the initial energy */
  Energy = 0d0
  
  do i=1,L
     do j=1,L


	    !/* Check periodic boundary conditions */

	p=i-1
	m=i-2
	n=j-1
	o=j-2

	if (p .eq. 0) then
		p=L
		m=L-1
	endif
	if (m .eq. 0) then
		m=L
	endif

	if (n .eq. 0) then
		n=L
		o=L-1
	endif
	if (o .eq. 0) then
		o=L
	endif



	i1=i+1
	i2=i+2
	j1=j+1
	j2=j+2

	if (i .eq. L) then
		i1=1
		i2=2
	endif
	if (i1 .eq. L) then
		i2=1
	endif
	if (j .eq. L) then
		j1=1
		j2=2
	endif
	if (j1 .eq. L) then
		j2=1
	endif



	V1part = V1*mag(i,j)*(mag(p,j) + mag(i1,j) + mag(i,n) + mag(i,j1))
	V2part = V2*mag(i,j)*(mag(m,j) + mag(i2,j) + mag(i,o) + mag(i,j2))
	Jpart = JJ/2d0*(4+mag(i,j)*(mag(p,n) + mag(i1,n) + mag(p,j1) + mag(i1,j1)))

	    !/* Update the energy */
	Energy = Energy + V1part + V2part + Jpart


	!/*Calculate the contribution from the field H */	... HH-->mu
	Energy = Energy + 2.d0*mu*mag(i,j)
     enddo
   enddo


   !/* Account for double counting */
   Energy = Energy/2.d0

   CalcE=Energy
   return
 end function CalcE



double precision function NNEnergy(mag,L,i,j)
use global_variables
implicit none
  integer i,j
  integer m,n,o,p
  integer i1, i2, j1, j2
  integer mag(Nmax,Nmax),L
  double precision E, V1part, V2part, Jpart
E=0d0
	    !/* Check periodic boundary conditions */

	p=i-1
	m=i-2
	n=j-1
	o=j-2

	if (p .eq. 0) then
		p=L
		m=L-1
	endif
	if (m .eq. 0) then
		m=L
	endif

	if (n .eq. 0) then
		n=L
		o=L-1
	endif
	if (o .eq. 0) then
		o=L
	endif



	i1=i+1
	i2=i+2
	j1=j+1
	j2=j+2

	if (i .eq. L) then
		i1=1
		i2=2
	endif
	if (i1 .eq. L) then
		i2=1
	endif
	if (j .eq. L) then
		j1=1
		j2=2
	endif
	if (j1 .eq. L) then
		j2=1
	endif



	V1part = V1*mag(i,j)*(mag(p,j) + mag(i1,j) + mag(i,n) + mag(i,j1))
	V2part = V2*mag(i,j)*(mag(m,j) + mag(i2,j) + mag(i,o) + mag(i,j2))
	Jpart = JJ/2d0*(4+mag(i,j)*(mag(p,n) + mag(i1,n) + mag(p,j1) + mag(i1,j1)))


	    !/* Update the energy */
	E = V1part + V2part + Jpart

	E = E + mu*mag(i,j)


NNEnergy=E
	
end function NNEnergy
