! Author: Sofiane Achour
! Description: multiphase flash Calculation with capillary pressure
  ! To run this code on a mac, follow the follwing steps:
    ! 1. open the terminal
    ! 2. paste the following: cd <filepath> where filepath is the project directory
    ! 3. paste the following: ifort filename.f90 && ./a.out

  ! To compile into an executable, follow the following steps:
    ! 1. open the terminal
    ! 2. paste the following: cd <filepath> where filepath is the project directory
    ! 3. paste the following: ifort -debug -O0 -c ./<program>.f90 -o ./<program>.o
    ! 4. paste the following: ifort -debug -O0 ./<program>.o -o ./<program>
    ! 5. Look at the results
! Miscellaneous:
  ! 1. The convention for the K matrix is the following: Kij, where i is the
  ! compt1nt and column number, and j is the phase number and the row number
  ! 2. The initial guess for beta may be non-physical, but it will still converge.
  ! Physical beta may allow it to converge faster
  ! 3. The Binary interaction coefficient must be at least lower triangular
  ! 4. V2 contains a better formatted input file with more options for the units of inputs
! This code attempts to include the effects of capillary pressure
! 03/26/2018
! 03/27/2018 improved the convergence criterion for stability analysis to match logic in SS flash
! + also corrected a mistake in the cardano code
program SS_flash

  implicit none
  character(len=100) :: fil1,fil2,text
  character(len=4) :: Tunit,Punit,zunit,Pcapstr
  integer ::  Nc , Npmax  ! Number of compt1nts and the max number of phases

  Real*8, allocatable, dimension(:):: Tc,Pc,w,z,Beta,&
                                        & xi,Z_z,Vm,logPcap

  Real*8, allocatable, dimension(:,:):: BIP,K,x,lnFug,Fugacity
  Real *8 :: residuals(1000),dummy,temp,T,P,Radius,&
  !       Variables below are conversion factors
  &         convT1,convz,Pcap(2),R,MaxTolRR,MaxTolSA,MaxTolPTSS

  integer :: i,j,ios,iter,Ns,Np,err,linsrch,MaxitSA,MaxitRR,MaxitPTSS
  logical :: unstable
  Pcap = (/0.d0,0.d0/)
!________________________________________________________________OPEN FILES_____________________________________________
  fil1 = 'input_files/input.txt'
  open(unit=1, iostat = ios,file=fil1, action= 'read', status = 'old')
  if ( ios /= 0 ) stop "Error opening the input data file"

  fil2 = 'output_files/output_data.txt'
  open(unit=2, iostat = ios,file=fil2, action = 'write', status = 'replace')
  if ( ios /= 0 ) stop "Error opening the output data file"


!______________________________________________IMPORT SIMULATION FLAGS AND ITERATIVE PARAMETERS_____________________________________________
!This flag is one when line_search is activated, 0 otherwise
  read(1,'(/,I1,/)',iostat = ios) linsrch

  read(1,'(/,ES20.10,/)',iostat = ios) MaxTolRR

  read(1,'(/,I4,/)',iostat = ios) MaxitRR

  read(1,'(/,ES20.10,/)',iostat = ios) MaxTolSA

  read(1,'(/,I4,/)',iostat = ios) MaxitSA

  read(1,'(/,ES20.10,/)',iostat = ios) MaxTolPTSS

  read(1,'(/,I4,/)',iostat = ios) MaxitPTSS

!________________________________________________________________IMPORT PHYSICAL PARAMETER_____________________________________________
  read(1,'(/,T9,ES20.8,/)',iostat = ios) Radius

  read(1,'(/,T9,A4,/,T9,ES20.8,/)',iostat = ios) Tunit,T
!
  read(1,'(/,T9,A4,/,T9,ES20.8,/)',iostat = ios) Punit,P
!
  read(1,'(T24,I2,/,T24,I1,/)',iostat = ios) Nc,Npmax

! Parameters to be used for converting T,P,Tc,Pc
  print *, 'radius = ',radius
  convT1 = 0d0

!
! ! Format of conversion is T = (T+convT1) to be used again for converting Tc,Pc
  if ( Tunit == "degF" ) then
    T = T+459.67d0
    R = 10.7316d0
    convT1 = 459.67d0
  else if ( Tunit == "degC" ) THEN
    T = T+273.15d0
    R = 8.3144598d0
    convT1 = 273.15d0

  end if

  if (Tunit == "degK") R = 8.3144598d0
  if (Tunit == "degR") R = 10.7316d0


! ________________________________________________________________ALLOCATE____________________________________________


  allocate(Vm(Npmax), stat=err)
  if ( err/= 0) print *, "Vm: Allocation request denied"

  allocate(Tc(Nc), stat=err)
  if ( err/= 0) print *, "Tc: Allocation request denied"

  allocate(Pc(Nc), stat=err)
  if ( err /= 0) print *, "Pc(Nc): Allocation request denied"

  allocate(w(Nc), stat=err)
  if ( err /= 0) print *, "w(Nc): Allocation request denied"
  allocate(BIP(Nc,Nc), stat=err)
  if ( err /= 0) print *, "BIP(Nc,Nc): Allocation request denied"

  allocate(z(Nc), stat=err)
  if ( err /= 0) print *, "z(Nc): Allocation request denied"

  allocate(K(Nc,Npmax), stat=err)
  if ( err /= 0) print *, "K(Nc,Npmax-1): Allocation request denied"

  allocate(x(Nc,Npmax), stat=err)
  if ( err /= 0) print *, "x(Nc,Npmax): Allocation request denied"

  allocate(Beta(Npmax), stat=err)
  if ( err /= 0) print *, "Beta(Npmax): Allocation request denied"

  allocate(lnFug(Nc,Npmax), stat=err)
  if ( err /= 0) print *, "lnFug(Nc,Npmax): Allocation request denied"

  allocate(fugacity(Nc,Npmax), stat=err)
  if ( err /= 0) print *, "fugacity(Nc,Npmax): Allocation request denied"

  allocate(logPcap(Npmax), stat=err)
  if ( err /= 0) print *, "logPcap(Npmax): Allocation request denied"

! ________________________________________________________________READ____________________________________________
  do j = 1, Npmax, 1
    logPcap(j) = DLOG(1d0+Pcap(j)/P)
  end do

  read(1,"(A100)",iostat = ios) text
  read(1,'(T9,A4)',iostat = ios) zunit

  if (zunit == "perc") THEN
    convz = 1d-2
  else
    convz = 1d0

  end if

  read(1,'(<Nc>(ES20.8,/))',iostat = ios) (z(i),i = 1,Nc)

  read(1,"(A100)",iostat = ios) text
  read(1,'(<Nc>(ES20.8,/))',iostat = ios) (Tc(i),i = 1,Nc)

  read(1,"(A100)",iostat = ios) text
  read(1,'(<Nc>(ES20.8,/))',iostat = ios) (Pc(i),i = 1,Nc)

  read(1,"(A100)",iostat = ios) text
  read(1,'(<Nc>(ES20.8,/))',iostat = ios) (w(i),i = 1,Nc)

  read(1,"(A100)",iostat = ios) text
  read(1,*,iostat = ios) ((BIP(i,j),j = 1,Nc),i = 1,Nc)

  if ( ios /= 0 ) stop "Error reading file input_data Tc,Pc,w,BIP"


  ! Data conversion
print *, Nc, Npmax
  do i = 1, Nc, 1
    z(i) = z(i) * convz
    Tc(i) = Tc(i)+convT1
  end do



print *, 'z'
print *, z
print *, 'Tc'
print *, Tc
print *, 'Pc'
print *, Pc
print *, 'omega'
print *, w
print *, "BIP"
call print_matrix(BIP,Nc,Nc)
! !________________________________________________________________
! Stab Calc

! First, calculate the a,b,A,B,Z of the main phase' composition

call PT_flash(Nc, Npmax,linsrch,MaxitSA,MaxitRR,MaxitPTSS,&
                & Tc,Pc,w,z,BIP,&
                & Radius,Pcap,R,MaxTolRR,MaxTolSA,MaxTolPTSS,T,P,&
                ! outputs are below
                & x,lnFug,Fugacity,Beta,Vm,residuals,Np,iter)


  print *, 'converged after',iter-1,'iterations'
  print *, 'residuals is',residuals(iter)
  print *, 'molar volume:',Vm
!________________________________________________________________OUTPUT FILE CREATION_______________________________________________________

    print *, Np
    Write(2,'(1X,A4,I2,A22)') "Your",Np," phase model is stable"

    write(unit=2, fmt = '(/,1X,A22/,1X,I2)', iostat=ios) "Number of components: ",Nc


    write(unit=2, fmt = '(/,1X,A55)', iostat=ios) "The original composition the mixture is the following: "
    write(unit=2, fmt = '(/,1X,(A20))', iostat=ios) "z"
    if ( ios /= 0 ) stop "Write error in output file for original composition data titles"

    write(unit=2, fmt = '(<Nc>(1X,(F20.16)/))', iostat=ios) (z(i),i=1,Nc)
    if ( ios /= 0 ) stop "Write error in output file for original composition data"

! print *, beta(1)
    write(unit = 2,fmt = '(/,1X,A24)', iostat = ios) "Phase molar fractions"
    if ( ios /= 0 ) stop "Write error in file unit 2"
    write(2,'(1X,<Np>(1X,F22.16))') (Beta(j),j=1,Np)

    write(unit = 2,fmt = '(/,1X,A23)', iostat = ios) "Molar volume"
    if ( ios /= 0 ) stop "Write error in file unit 2"
    write(2,'(1X,<Np>(1X,F20.16))') (Vm(j),j=1,Np)


    write(unit=2, fmt = '(/,1X,A48)', iostat=ios) "The composition of each phase is the following: "
    write(unit=2, fmt = '(/,1X,<Np>(A18,I2))', iostat=ios) ("xi",j,j=1,Np)
    if ( ios /= 0 ) stop "Write error in output file for composition data titles"

    write(unit=2, fmt = '(<Nc>(1X,<Np>(F20.16)/))', iostat=ios) ((x(i,j),j=1,Np),i=1,Nc)
    if ( ios /= 0 ) stop "Write error in output file for composition data"

    write(unit=2, fmt = '(/,1X,A45)', iostat=ios) "The Fugacity of each phase is the following: "
    write(unit=2, fmt = '(/,1X,<Np>(A18,I2,A5))', iostat=ios) ("f",j," (psia)",j=1,Np)
    if ( ios /= 0 ) stop "Write error in output file for composition data titles"

    write(unit=2, fmt = '(<Nc>(1X,<Np>(ES25.16)/))', iostat=ios) ((Fugacity(i,j),j=1,Np),i=1,Nc)
    if ( ios /= 0 ) stop "Write error in output file for composition data"

    write(unit=2, fmt = '(/,1X,A54)', iostat=ios) "The value of the objective function is the following: "
    if ( ios /= 0 ) stop "Write error in output file for composition data titles"

    write(unit=2, fmt = '(<iter>(1X,(ES25.16)/))', iostat=ios) (residuals(i),i=1,iter)
    ! write(unit=2, fmt = '(<iter>(1X,(I5,ES25.16)/))', iostat=ios) (i,residuals(i),i=1,iter)
    if ( ios /= 0 ) stop "Write error in output file for composition data"
  if (.not.(unstable)) then
    write(unit=2, fmt = '(/,1X,A37,I4,A12)', iostat=ios) "The flash calculation converged after",iter," iterations."
  else

    write(unit=2, fmt = '(/,1X,A46,I4,A12)', iostat=ios) "The flash calculation failed to converge after",iter," iterations."

  end if
! !

  do j = 1, Np, 1
    if ( beta(j).lt.0d0.or.beta(j).gt.1d0 ) then
      print *, 'Capillary pressure suppressed the bubbles from forming',Beta
      exit
    end if
  end do


  close (2)
  close (1)


! ____________________________________________________________________DEALLOCATION___________________________________________________

  if (allocated(Vm)) deallocate(Vm, stat = err)
  if ( err /= 0) print *, "Vm: Deallocation request denied"

  if (allocated(fugacity)) deallocate(fugacity, stat = err)
  if ( err /= 0) print *, "fugacity: Deallocation request denied"

  if (allocated(lnFug)) deallocate(lnFug, stat = err)
  if ( err /= 0) print *, "lnFug: Deallocation request denied"

  if (allocated(Beta)) deallocate(Beta, stat = err)
  if ( err /= 0) print *, "Beta: Deallocation request denied"

  if (allocated(x)) deallocate(x, stat = err)
  if ( err /= 0) print *, "x: Deallocation request denied"

  if (allocated(K)) deallocate(K, stat = err)
  if ( err /= 0) print *, "K: Deallocation request denied"

  if (allocated(z)) deallocate(z, stat = err)
  if ( err /= 0) print *, "z: Deallocation request denied"

  if (allocated(BIP)) deallocate(BIP, stat = err)
  if ( err /= 0) print *, "BIP: Deallocation request denied"

  if (allocated(w)) deallocate(w, stat = err)
  if ( err /= 0) print *, "w: Deallocation request denied"

  if (allocated(Pc)) deallocate(Pc, stat = err)
  if ( err /= 0) print *, "Pc: Deallocation request denied"

  if (allocated(Tc)) deallocate(Tc, stat=err)
  if ( err/= 0) print *, "Tc: Deallocation request denied"

  if (allocated(logPcap)) deallocate(logPcap, stat=err)
  if ( err/= 0) print *, "logPcap: Deallocation request denied"

end program SS_flash

subroutine PT_flash(Nc, Npmax,linsrch,MaxitSA,MaxitRR,MaxitPTSS,&
                & Tc,Pc,w,z,BIP,&
                & Radius,Pcap,R,MaxTolRR,MaxTolSA,MaxTolPTSS,T,P,&
                ! outputs are below
                & x,lnFug,Fugacity,Beta,Vm,residuals,Np,iter)
  implicit none
  integer,intent(in) :: Nc, Npmax,linsrch,MaxitSA,MaxitRR,MaxitPTSS

  Real*8, intent(in):: Tc(Nc),Pc(Nc),w(Nc),z(Nc),BIP(Nc,Nc)

  Real*8, intent(in) :: Pcap(Npmax),Radius,R,MaxTolRR,MaxTolSA,MaxTolPTSS,T,P

  Real*8, intent(out) :: x(Nc,Npmax),lnFug(Nc,Npmax),Fugacity(Nc,Npmax),&
              &   Beta(Npmax),Vm(Npmax),residuals(1000)

  integer,intent(out) :: Np,iter

  Real*8 :: rtAi(Nc),bi(Nc),Par(Nc),&
            & xi(Nc),Z_z(Npmax),logPcap(Npmax),K(Nc,Npmax-1)

  Real *8 :: dummy,temp

  integer :: i,j,ios,err,Ns
  logical :: unstable
  ! Print *, 'Radius = ',Radius
  do j = 1, Npmax, 1
    logPcap(j) = DLOG(1d0+Pcap(j)/P)
  end do
  iter = 0
  unstable = .true.
  Ns = 1

  j = 1
  call Calc_Mixture_Ln_Fug_coef(Tc,Pc,w,T,P,BIP,z,lnFug(:,j),Z_z(j),Nc,0)
  ! call Calc_Mixture_Ln_Fug_coef(rtAi,bi,T,P,BIP,x(:,j),lnFug(:,j),Z_z(j),Nc)
  print *, 'ln(phi) = '
  print *, lnFug(:,j)

  call Parachor(Tc,Pc,w,Nc,Par)

  call SinglePhase_Stability_check(z,T,P,Tc,Pc,w,BIP,Par,Radius,lnFug(:,j),Z_z(j),R,&
                                  &, unstable,x,Nc,MaxitSA,MaxTolSA)

  ! print *,  '__________________________________'
  if ( unstable ) then

    Write(*,'(1X,A5,I2,A8,I2,A24)') "Phase",j," of your",Ns," phase model is unstable"
    print *, 'Initiating 2 phase PT flash calculation'
    Ns = 2
  else
    x(:,1) = z
    Np = 1
    print *, "Only one phase forms for the given Temperature and Pressure conditions"
  end if

  Ns = 2
residuals(1) = 1d0
! 2 phase flash K-values:
  if ( unstable ) then

    call Wilson(T,P,Tc,Pc,w,Nc,K(:,1))

  end if
  j = 1
  print *, Ns
  !________________________________________________________________
  !________________REMOVE LINE BELOW_______________________________
  !________________________________________________________________
  ! unstable = .False.
! initiate multiphase flash
  stability_loop: do while ( unstable .and. Ns <= Npmax)

    residuals(iter) = 1d0

    SS_residuals_loop : Do while ((residuals(iter).ge.MaxTolPTSS).and.(iter-1 .le. maxitPTSS))
      iter = iter+1
  ! Using the K values previously estimated, find the x values

      call RR_flash(Nc,Ns,z,K(:,1:Ns-1),x(:,1:Ns),beta(1:Ns-1),linsrch,MaxitRR,MaxTolRR)

  !
      calc_ln_fug_coeff_for_NcCompenents_and_NsPhases : do j = 1, Ns, 1

        call Calc_Mixture_Ln_Fug_coef(Tc,Pc,w,T,P+Pcap(j),BIP,x(:,j),lnFug(:,j),Z_z(j),Nc,0)

  !
      end do calc_ln_fug_coeff_for_NcCompenents_and_NsPhases

!
      residuals(iter) = 0d0
!     Calaculate the new Kvalues for all Nc comp and Np - 1 phases
      Loop_over_phases : do j = 1, Ns-1, 1
        Loop_over_components : do i = 1, Nc, 1

          dummy = lnFug(i,Ns)-lnFug(i,j)-logPcap(j)

          temp =  abs(LOG(K(i,j)) - dummy)
          if ( iter.lt.2 ) then
            ! print *, temp
          end if
          if ( residuals(iter).lt.temp ) then
            residuals(iter) = temp
          end if
          ! if ( iter.lt.2 ) then
            ! print *, 'residuals',residuals(iter)
          ! end if

          K(i,j) = exp(dummy)
        end do Loop_over_components

      end do Loop_over_phases

    end do SS_residuals_loop
!
!
    print *, Ns,"phase stability check"
    call MultiPhase_Stability_check(x(:,1),T,P+Pcap(1),Tc,Pc,w,BIP,lnFug(:,1), &
    & unstable,xi,Nc,MaxitSA,MaxTolSA)

    if ( unstable ) then
      Calc_new_K: do i = 1, Nc, 1
        K(i,Ns) = xi(i)/x(i,Ns)
        ! print *, K(i,Ns)
      end do Calc_new_K
      exit
    end if

    if ( unstable ) then
      Write(*,'(1X,A5,I2,A8,I2,A24)') "Phase",j," of your",Ns," phase model is unstable"
      print *, 'Initiating phase-split calculations'
      Ns = Ns+1
    else
      Np = Ns
    end if

  end do stability_loop
  ! print *, Ns


  if (unstable) then
    Np = Ns-1
  end if

  Calc_Fugacity: do j = 1, Np, 1
    Vm(j) = R*T*Z_z(j)/(P+Pcap(j))
    do i = 1, Nc, 1
      Fugacity(i,j) = x(i,j) * (P+Pcap(j)) * exp(lnFug(i,j))
    end do

  end do Calc_Fugacity

  dummy = 0d0
  Calc_beta_Np: do j = 1, Np-1, 1
    dummy = dummy + beta(j)

  end do Calc_beta_Np

  beta(Np) = 1d0 - dummy

  print *, 'Final phase mole fractions',beta
end subroutine PT_flash

subroutine PReos_AB(T,P,w,Tc,Pc,bi,rtAi,Nc)
  Integer, intent(in) :: Nc
  real*8, intent(in) :: T,P,w(Nc),Tc(Nc),Pc(Nc)
  real*8, intent(out):: rtAi(Nc),bi(Nc)
  real*8 alpha(Nc),kappa(Nc),Tr(Nc),invTr(Nc),Pr(Nc),PT,PT2,Ai(Nc)
  real*8::smalla,smallb
  integer :: i
!

  ! write(unit=*, fmt='(1X,A27)') 'PR_EOS_AB'
  ! write(unit=*, fmt='(1X,2ES27.18,I3)') T,P,Nc
PT = P/T
PT2 = P/(T*T)
! print *, 'PT',PT
! print *, 'PT2',PT2
! print *, 'Tc',Tc
! print *, 'Pc',Pc

  where ( w .le. 0.49d0 )
    kappa=0.37464d0+1.54226d0*w-0.26992d0*w*w
  elsewhere
    kappa=0.37964d0+w*(1.48503d0+w*(-0.164423d0+0.016666d0*w))
  end where

  do i = 1, Nc, 1

    ! write(unit=*, fmt='(1X,ES27.18)')w(i)
    invTr(i) = Tc(i)/T
    Tr(i) = T/Tc(i)
    Pr(i) = P/Pc(i)
    ! print *, "Tr =",Tr
    Ai(i) = 0.457236d0*Tc(i)*Tc(i)*PT2/(Pc(i)) &
      &      *(1.d0+kappa(i)*(1.d0-dsqrt(Tr(i))))&
      &      *(1.d0+kappa(i)*(1.d0-dsqrt(Tr(i))))


    smalla = 0.457236d0*Tc(i)*Tc(i)*8.3144598d-5**2.d0/(Pc(i)*1.d-5) &
      &      *(1.d0+kappa(i)*(1.d0-dsqrt(Tr(i))))&
      &      *(1.d0+kappa(i)*(1.d0-dsqrt(Tr(i))))
    rtAi(i) = DSQRT(Ai(i))

    Bi(i) = 0.0778D0 * Tc(i)*PT/(Pc(i))
    smallb = 0.0778D0 * Tc(i)*8.3144598d0/(Pc(i))
    ! write(*, fmt="(1X,A27,3(/,1X,ES27.18))") ,'Ai, Bi, kappa',smalla,smallb,kappa(i)
    ! write(*, fmt="(1X,A27,4(/,1X,ES27.18))") ,'Ai, rtAi,Bi, kappa',Ai(i),rtAi(i),Bi(i),kappa(i)
  end do

  ! write (*,'(1X,A32,<Nc>(/,ES20.8))') "ai's have the following values: ", ai

end


subroutine Wilson(T,P,Tc,Pc,w,Nc,K)
  implicit none
  integer ,intent(in):: Nc
  Real*8, intent(in) ::  T,P,Tc(Nc),Pc(Nc),w(Nc)
  Real*8, intent(out) :: K(Nc)
  Real *8 :: Pin,Tin
  integer :: i

  Tin = 1.d0/T
  Pin = 1.d0/P
  ! print *, 'T = ',T
  ! print *, 'P = ',P
  Wilson_corel: do i = 1, Nc, 1

    K(i) = Pc(i)*Pin*dexp(5.37d0*(1.d0+w(i))*(1.d0-Tc(i)*Tin))

  end do Wilson_corel
end subroutine Wilson

subroutine Parachor(Tc,Pc,w,Nc,Par)
  implicit none
  integer ,intent(in):: Nc
  Real*8, intent(in) ::  Tc(Nc),Pc(Nc),w(Nc)
  Real*8, intent(out) :: Par(Nc)
  integer :: i

  ! print *, 'T = ',T
  ! print *, 'P = ',P
  Parachor_corel: do i = 1, Nc, 1
    Par(i) = (8.21307d0+1.97473d0*w(i))*(Tc(i)**1.03406d0)*(Pc(i)**(-0.82636d0))*1.d-6
  end do Parachor_corel

end subroutine Parachor


subroutine Calc_Pc(x,rhox,Nc,Par,rhoz,zi,radius,Pcap)
  implicit none
  integer,intent(in) :: Nc
  Real*8, intent(in) :: x(Nc),rhox,Par(Nc),rhoz,zi(Nc),radius
  Real*8,intent(out) :: Pcap
  Real*8:: IFT

  integer :: n,i,iter

  ! print *, 'x',x
  ! print *, 'parachors',par
  ! print *, 'bulk density',rhoz
  ! print *, 'bulk composition',zi
  ! print *, 'radius',radius
  ! print *, 'This is the value of the ideal gas constant',R
  ! Pjk = (/P,P/)

  ! print *, 'Phase input pressures',P,ID

    calc_IFT: do i = 1, Nc, 1
      IFT = IFT+Par(i)*(zi(i)*rhoz-x(i)*rhox)
    end do calc_IFT
    ! IFT = IFT**(3.6d0)
    ! print *, 'IFT',IFT**3.6d0,IFT
    if ( IFT.gt.0d0 ) then
      IFT = IFT**(3.88d0)
      ! print *, 'IFT was gt 0'
    else
      IFT = -(-IFT)**(3.88d0)
      ! print *, 'IFTIFT was lt 0'
    end if
    Pcap = 2.d0*IFT/radius


end subroutine Calc_Pc
!

subroutine VDW_mixing_rules(Nc,rtAi,Bi,T,P,BIP,x,Amix,Bmix,Aik)
  implicit none
  integer, intent(in) :: Nc
  Real*8, intent(in) :: rtAi(Nc),Bi(Nc),BIP(Nc,Nc),T,P,x(Nc)
  Real*8, intent(out) :: Amix,Aik(Nc,Nc),Bmix
  Real*8 :: dummy


  integer :: i,k

  Amix = 0.d0
  Bmix = 0.d0

  do i = 1, Nc, 1
    Aik(i,i) = rtAi(i)*rtAi(i)

    Bmix = Bmix + Bi(i)*x(i)
    do k = 1, i-1, 1
      Aik(i,k) = rtAi(i)*rtAi(k)*(1.d0-BIP(i,k))

      Aik(k,i) = Aik(i,k)

      Amix = Amix + x(i)*x(k)*Aik(i,k)+x(i)*x(k)*Aik(i,k)
    end do

    Amix = Amix + Aik(i,i)*x(i)*x(i)
  end do


end subroutine VDW_mixing_rules


subroutine CALC_Z(A,B,Zl,Zu,n)
  implicit none
  Real*8, intent(in) :: A,B
  Real*8, intent(out) :: Zl , Zu
  integer, intent(out) :: n ! n is th enumber of real roots
  Real*8 :: D,E,Del,al,be,ga,theta,Z3,dummy,dummy2,Ecubed,rtnE
  Real*8, parameter :: Pi = 3.141592653589793238d0


  al = B-1.d0
  be = A-3.d0*B*B-2.d0*B
  ga = B*B*B + B*B - A*B
  ! print *, 'cubic par',al,be,ga
  ! This subroutine solves the cubic equation:
  ! Z^3 + al* Z^2 + be * Z + ga = 0

  D = (2.d0*al*al*al-al*be*9.d0+27.d0*ga)/54.d0
  E = (be*3.d0-al*al)/9.d0
  Ecubed = E*E*E
  Del = D*D + Ecubed

  if (del.gt.0d0) then ! One Real solution only
    n = 1
    dummy = (abs(D)+DSQRT(Del))**(1.d0/3.d0)
    if ( D.gt.0d0 ) dummy = -dummy
    Zl =  dummy - E/dummy- al/3.d0
    ! print *, "screamo"
    Zu= 0d0
  else if ( del.lt.0d0 ) then ! Three real solutions
    n = 2
    theta = DACOS(-1d0*D/DSQRT(-1d0*Ecubed))
    rtnE = DSQRT(-E)

    Zl = 2.d0*rtnE * DCOS((theta+2.d0*Pi)/3.d0)-al/3.d0
    Zu = 2.d0*rtnE * DCOS((theta+4.d0*Pi)/3.d0)-al/3.d0
    Z3 = 2.d0*rtnE * DCOS(theta/3.d0)-al/3.d0
    ! print *, "Z:"
    ! print *, Zl, Zu, Z3
! Re-order Zl and Zu such that Zl = min(Zi), Zu = max(Zi), i =1,2,3
    if ( Zl<Zu ) then
      dummy = Zl
      dummy2 = Zu
    else
      dummy = Zu
      dummy2 = Zl
    end if
    if ( dummy>Z3 ) then
      dummy = Z3
    end if
    if ( dummy2<Z3 ) then
      dummy2 = Z3
    end if

    Zl = dummy
    Zu = dummy2
  else if ( abs(del).lt.1d-16.and.abs(D).lt.1d-16 ) then ! Three real repeated solutions
    n = 1
    Zl = -al/3d0
    Zu= 0d0
  else  ! Two repeated real solutions plus t1 unique solution
    n = 2

    Zl = 2d0*(-1d0*D)**(1.d0/3.d0) - al*1d0/3d0
    Zu = -1d0*(-1d0*D)**(1.d0/3.d0) - al*1d0/3d0
    ! Zl is not necessarily less than Zu
  end if

end subroutine

pure function cubroot(x)
  implicit none
  Real*8, intent(in) :: x
  Real*8 :: cubroot

  if (x.lt.0) then
    cubroot = -DEXP(DLOG(-1d0*x)*1d0/3d0)
  else if (x.gt.0) then
    cubroot = DEXP(DLOG(1d0*x)*1d0/3d0)
  else
    cubroot = 0d0
  end if

end function cubroot


subroutine CALC_lnFugacityi(Aik,B,Bmix,Amix,x,Z,lnFug,Nc)
  implicit none
  integer,intent(in) :: Nc
  Real*8,intent(in) :: Amix,Bmix,Aik(Nc,Nc),B(Nc),Z,x(Nc)
  Real*8,intent(out) :: lnFug(Nc)
  Real*8 :: sumAikxk,logterm,term1,term2,term3,&
          & t2rt2,t1prt2,t1mrt2,rt2
  integer :: i,k
  Data t2rt2,rt2,t1prt2,t1mrt2 &
   &/ 2.82842712474619029095d0,1.41421356237309514547d0,&
   &  2.41421356237309514547d0,-0.41421356237309514547d0/

  logterm = DLOG((Z+t1prt2*Bmix)/(Z+t1mrt2*Bmix))
  term1 = -DLOG(Z-Bmix)
  term2 = -logterm/(rt2*Bmix)
  term3 = (Z-1.d0)/Bmix+Amix*logterm/(t2rt2*Bmix*Bmix)

  lnFugacity_of_ij: do i = 1, Nc, 1
    sumAikxk = 0.d0
    do k = 1,Nc, 1
      sumAikxk = sumAikxk + x(k)*Aik(i,k)
    end do
    lnFug(i) = term1 + term2*sumAikxk + B(i)*term3
  end do lnFugacity_of_ij

end subroutine CALC_lnFugacityi



subroutine Calc_Mixture_Ln_Fug_coef(Tc,Pc,w,T,P,BIP,x,lnFug,Z,Nc,ID)
  implicit none
  integer,intent(in) :: Nc,ID
  ! ID is 0 for conventional root selection, 2 for metastable liquid,
  ! and 1 for metastable gas
  Real*8, intent(in) :: Tc(Nc),Pc(Nc),w(Nc),BIP(Nc,Nc),x(Nc),T,P
  Real*8,intent(out) :: lnFug(Nc),Z
  Real*8 :: Amix,Bmix,Aik(Nc,Nc),A(Nc),B(Nc),Zl,Zu,temp(Nc),dummy,&
            & rtAi(Nc),bi(Nc)
  integer :: n,i

  call PReos_AB(T,P,w,Tc,Pc,bi,rtAi,Nc)


  call VDW_mixing_rules(Nc,rtAi,bi,T,P,BIP,x,Amix,Bmix,Aik)
  ! write (*,'(1X,A31,<Nc>(/,ES20.8),/)') "For the following composition: ", x

  ! print *, 'Amix',Amix,'Bmix',Bmix
  call CALC_Z(Amix,Bmix,Zl,Zu,n)
  ! print *, 'n,Zl,Zu',n,Zl,Zu
  ! print *, ID, (ID==0)
  if ( ID==0.or.ID.gt.2 ) then
    if ( n==2 ) then
      call CALC_lnFugacityi(Aik,Bi,Bmix,Amix,x ,Zl,lnFug,Nc)
      call CALC_lnFugacityi(Aik,Bi,Bmix,Amix,x ,Zu,temp,Nc)

      dummy = 0d0
      Calc_Delta_G : do i = 1, Nc, 1
        dummy = dummy + x(i)*(lnFug(i)-temp(i))
      end do Calc_Delta_G
      if ( dummy.gt.0d0 ) then

        lnFug  = temp
        Zl = Zu
      end if

    else
      call CALC_lnFugacityi(Aik,Bi,Bmix,Amix,x,Zl,lnFug ,Nc)
    end if
    Z = Zl
  else
    if ( n.eq.2.and.ID.eq.1 ) then ! not in VDW loop and Gas
      Z = Zu
    else ! it is a liquid or not in VDW loop
      Z = Zl
    end if
    call CALC_lnFugacityi(Aik,Bi,Bmix,Amix,x,Z,lnFug ,Nc)
  end if

end subroutine Calc_Mixture_Ln_Fug_coef

subroutine SinglePhase_Stability_check(z,T,P,Tc,Pc,w,BIP,Par,Radius,lnPhi,Z_CF_b,R, &
                                        & unstable,xi,Nc,MaxitSA,MaxTolSA)
  implicit none

  integer, parameter :: Ng = 2! number of guesses
  integer , intent(in) :: Nc,MaxitSA
  Real*8, intent(in) :: T,P,Tc(Nc),Pc(Nc),w(Nc),BIP(Nc,Nc),z(Nc),Par(Nc),&
                      & Radius,lnPhi(Nc),Z_CF_b,R,MaxTolSA
  Real*8,intent(out) :: xi(Nc)
  logical,intent(out) :: unstable

  Real*8 :: lnFug(Nc,Ng),x(Nc,Ng),yi(Nc,Ng),residuals,sumyi(Ng),K(Nc),Z_CF,&
            & rtAi(Nc),bi(Nc),temp1,temp2,di(Nc),Pcap(Ng),rhox,rhoz
  integer :: i,j,iter
  iter = 0
  sumyi = (/(0d0,j = 1,Ng)/)
  rhoz = P/(Z_CF_b*R*T)
  Wilson_corel: do i = 1, Nc, 1 ! Use two of Michelsen's initial guesses

    K(i) = (Pc(i)/P)*exp(5.37d0*(1.d0+w(i))*(1.d0-Tc(i)/T))

    yi(i,1) = z(i)/K(i)
    yi(i,2) = z(i)*K(i)
    ! print *, K(i),yi(i,1),yi(i,2)
    do j = 1, Ng, 1
      sumyi(j) = sumyi(j) + yi(i,j)
    end do
    ! print *, 'phase',j

    ! di(i) = lnPhi(i) + DLOG(z(i))
    di(i) = lnPhi(i) + DLOG(z(i)*P)
  end do Wilson_corel


  normalizex: do j = 1, Ng, 1

    do i = 1, Nc, 1
      x(i,j) = yi(i,j) / sumyi(j)
    end do

  end do normalizex

  calc_Pcapj: do i = 1, Ng, 1
    call Calc_Mixture_Ln_Fug_coef(Tc,Pc,w,T,P,BIP,x(:,j),lnFug(:,j),Z_CF,Nc,0)
    rhox = P/(Z_CF*R*T)
    ! print *, rhox
    call Calc_Pc(x(:,j),rhox,Nc,Par,rhoz,z,radius,Pcap(j))

  end do calc_Pcapj
  ! print *, 'Capillary pressure', Pcap
  residuals = 1d0

  do while ((residuals.ge.MaxTolSA).and.(iter.lt.MaxitSA))

    do j = 1, Ng, 1

      call Calc_Mixture_Ln_Fug_coef(Tc,Pc,w,T,P,BIP,x(:,j),lnFug(:,j),Z_CF,Nc,0)
    end do

    residuals = 1d-16
    update_tol: do j = 1, Ng, 1
      temp2 = 1d0
      do i = 1, Nc, 1

        temp1 = abs(DLOG(yi(i,j)) + lnFug(i,j)- di(i))
        if ( temp1.lt.temp2 ) then
          temp2 = temp1
        end if
        end do
      ! print *, temp2
      if ( temp2.gt.residuals ) then
        residuals = temp2
      end if
    end do update_tol

    iter = iter+1
    sumyi = (/(0d0,j = 1,Ng)/)

    y_tp1: do j = 1, Ng, 1
      do i = 1, Nc, 1
        yi(i,j) = z(i) * Dexp(lnPhi(i) - lnFug(i,j))
        sumyi(j) = sumyi(j) + yi(i,j)
      end do
    end do y_tp1


    calc_yi: do j = 1, Ng, 1
      do i = 1, Nc, 1
        x(i,j) = yi(i,j)/sumyi(j)
      end do
    end do calc_yi

      ! print *, residuals
      ! print *, sumyi
  end do
  print *, 'SinglePhase Stability check'
  do j = 1, 2, 1
    Print *, 'converged guess',j
    print *, x(:,j)
    print *,'sum is ',sumyi(j)
  end do
  unstable = .false.
  do j = 1, Ng, 1
    if ( sumyi(j)-1.d0 > 1d-05 ) then
      unstable = .true.
      xi = x(:,j)
    end if
  end do

end subroutine SinglePhase_Stability_check


subroutine MultiPhase_Stability_check(z,T,P,Tc,Pc,w,BIP,lnPhi,unstable,xi,Nc,MaxitSA,MaxTolSA)
  implicit none

  integer, parameter :: Ng = 4! number of guesses
  integer , intent(in) :: Nc,MaxitSA
  Real*8, intent(in) :: T,P,Tc(Nc),Pc(Nc),w(Nc),BIP(Nc,Nc),z(Nc),&
                      & lnPhi(Nc),MaxTolSA
  Real*8,intent(out) :: xi(Nc)
  logical,intent(out) :: unstable

  Real*8 :: lnFug(Nc,Ng),x(Nc,Ng),yi(Nc,Ng),residuals,sumyi(Ng),K(Nc),Z_CF,&
            & rtAi(Nc),bi(Nc),di(Nc),temp1,temp2
  integer :: i,j,iter
print *, "New multi-phase stability check"
print *, z
  iter = 0
  sumyi = (/(0d0,j = 1,Ng)/)
  ! yi stands for X in Michelsen's paper

  initial_yi: do i = 1, Nc, 1 ! For multiphase Flash, use guesses described in Okuno Thesis

    yi(i,3) = z(i)*exp(lnPhi(i))
    yi(i,4) = 1d0/Nc

    do j = 3, Ng, 1
      sumyi(j) = sumyi(j) + yi(i,j)
    end do
    di(i) = lnPhi(i) + DLOG(z(i))
  end do initial_yi

  normalizex: do j = 3, Ng, 1

    do i = 1, Nc, 1
      x(i,j) = yi(i,j) / sumyi(j)
    end do
    print *, x(:,j)

  end do normalizex
  ! print *, MaxTolSA,MaxitSA
  ! print *, residuals,MaxTolSA,iter,MaxitSA
  x(2:Nc,1) = (/(0d0,i = 1,Nc-1)/)
  x(1:Nc-1,2) = (/(0d0,i = 1,Nc-1)/)
  x(1,1) = 1.d0
  x(Nc,2) = 1d0

  residuals = 1d0
  iter = 1

  do while ( (residuals.ge.MaxTolSA).and.(iter.lt.MaxitSA) )
  iter = iter+1
    do j = 1, Ng, 1
      call Calc_Mixture_Ln_Fug_coef(Tc,Pc,w,T,P,BIP,x(:,j),lnFug(:,j),Z_CF,Nc,0)
    end do

    residuals = 1d-16
    update_tol: do j = 1, Ng, 1
      temp2 = 1d0
      do i = 1, Nc, 1

        temp1 = abs(DLOG(yi(i,j)) + lnFug(i,j)- di(i))
        if ( temp1.lt.temp2 ) then
          temp2 = temp1
        end if
        end do
      ! print *, 'temp2',temp2
      if ( temp2.gt.residuals ) then
        residuals = temp2
      end if
    end do update_tol

    sumyi = (/(0d0,j = 1,Ng)/)

    y_tp1: do j = 1, Ng, 1
      do i = 1, Nc, 1
        yi(i,j) = z(i) * Dexp(lnPhi(i) - lnFug(i,j))
        sumyi(j) = sumyi(j) + yi(i,j)
      end do
    end do y_tp1


    calc_yi: do j = 1, Ng, 1
      do i = 1, Nc, 1
        x(i,j) = yi(i,j)/sumyi(j)
      end do
    end do calc_yi

  end do

! print *, 'sumyi is' ,sumyi
  ! write(unit=*, fmt='(A25,7(/,ES21.13))') 'The original phas comp is',z
  unstable = .false.
  do j = 1, Ng, 1

    if ( sumyi(j).gt.1.d0+1d-5 ) then
      unstable = .true.
      xi = x(:,j)
      write(unit=*, fmt='(A25,7(/,ES21.13))') 'The distabilizing comp is',xi
    end if
  end do

end subroutine MultiPhase_Stability_check






!________________________________________________________________
!________________________________________________________________
!___________________K-FLASH CALCULATION__________________________
!________________________________________________________________

subroutine RR_flash(Nc,Np,z,K,x,Beta,linsrch,MaxitRR,MaxTolRR)
  implicit none
!
  integer, intent(in) ::  linsrch,Nc,Np,MaxitRR  ! Linesearch flag, number of components, the number of phases
  Real*8, intent(in) :: z(Nc),K(Nc,Np-1),MaxTolRR
  Real*8, intent(out) :: x(Nc,Np),Beta(Np-1)

  Real*8 :: b(Nc),Beta_est(Np-1),dummy,f_vector(Np-1),&
          & invt(Nc),Hessian(Np-1,Np-1),NewtonDir(Np-1),&
          & lmax,s,residuals,betaNp,sumx(Np),sumxin(Np)
  Real*8, allocatable, dimension(:,:) :: Beta_inter

  integer :: i,j,n,l,err,fact,length


! First, allocate the allocatables

  dummy = fact(Np - 1)*fact(Nc-Np + 1)
  length = fact(Nc)/dummy

  allocate(Beta_inter(Np-1,length), stat=err)
  if ( err/= 0) print *, "Beta_inter: Allocation request denied"



! Calculate b
  calcb: do i = 1, Nc, 1
    b(i) = 1d0-z(i)
    do j = 1, Np-1, 1
      if (1.D0-k(i,j)*z(i).lt.b(i)) then
        b(i) = 1.D0-K(i,j)*z(i)
      end if
    end do
  end do calcb
  close (1)

! Find the intersection of the equations in set S, and copmute initial estimate for beta

  call find_intersec_vertices(K,b,Beta_inter,Nc,Np,l,length)



  Beta = (/(0.D0,i = 1,Np-1)/)
  dummy = 1.D0/DBLE(l)

  calc_beta_init: do i = 1,l, 1
    beta = beta+beta_inter(:,i)*dummy
  end do calc_beta_init
  Beta_est = Beta



  n = 0
  residuals = 1.D0

  ! do while (residuals >= 1d-18)

  do while ((n.lt.MaxitRR).and.(residuals.ge.MaxTolRR))

    n = n+1
    ! Calculate the f, vector consisting of Ratchford-Rice equations
    call RR_gradient_calc(z,K,Beta,f_vector,invt,Nc,Np)


  ! Calculate the Hessian of f

    call Hessian_calc(K,z,invt,Hessian,Nc,Np)


    ! Calculate Newton's Direction
    call solve_matmul_by_LU(Hessian,-1.D0*f_vector,Np-1,NewtonDir)



    call lambda_max(K,b,Beta,NewtonDir,lmax,Nc,Np)

!   if lnsrch is 1, then use armijo's rule to determine s otherwise, s = 1
    if ( linsrch.eq.1 ) then
      call line_search(z,K,beta,NewtonDir,lmax,s,Nc,Np)
    else
      s = 1.d0
    end if



    do i = 1, Np-1, 1
      beta(i) = beta(i) + s*lmax*NewtonDir(i)
    end do

    ! print *, residuals
    residuals = 0.D0
    ! residuals is calculated as the maximum norm of the gradient containing
    ! the RR equations
    do i = 1, Np-1, 1
      ! print *, f_vector(i),s*lmax*NewtonDir(i)
      if (abs(f_vector(i)).gt.residuals) THEN
        residuals = abs(f_vector(i))
      end if
    end do

  end do

  calc_invti: do i = 1, Nc, 1
    dummy= 1.D0
    calc_ti : do j = 1, Np-1, 1
      dummy = dummy+(K(i,j)-1.d0)*beta(j)
    end do calc_ti
    invt(i) = 1.D0/dummy
  end do calc_invti
  do i = 1, Nc, 1
  end do

  sumx(Np) = 0.d0
  do i = 1, Nc, 1
    x(i,Np) = z(i)*invt(i)

    sumx(Np) = sumx(Np) + x(i,Np)
  end do

  betaNp = 1.D0
  do j = 1, Np-1, 1
    sumx(j) = 0.d0
    betaNp = betaNp - beta(j)
  end do

  do i = 1, Nc, 1
    do j = 1, Np-1, 1
      x(i,j) = x(i,Np)*K(i,j)
      sumx(j) = sumx(j)+x(i,j)
    end do
  end do

  do j = 1, Np, 1
    sumxin(j) = 1.d0/sumx(j)
  end do
  do i = 1, Nc, 1
    do j = 1, Np, 1
      x(i,j) = x(i,j)*sumxin(j)
    end do
  end do


  if (allocated(beta_inter)) deallocate(beta_inter, stat=err)
  if ( err/= 0) print *, "indices: Deallocation request denied"


end subroutine RR_flash








!____________________________________________________________________________________________________
! These are the subroutines used for step 1
function fact(n)
  implicit none
  integer, intent(in) :: n
  integer :: i, fact

  fact = 1

  do i = 1, n, 1
    fact = fact*i
  end do
  return
end

! This recursive subrouting provides the indices of the operations that must eb applied to find intersection of sets
recursive subroutine recursive_test(Nc,Np,indices,length,i,m,pos)
  implicit none
  integer, intent(in) :: Nc,Np,length,m,pos
  integer, intent(inout) :: indices(length,Np),i
  integer :: j,n,elements(Nc-m+1),counter

  elements = (/(n, n = m,Nc,1)/)

  counter = 1

  do j = 1,Nc-(Np-pos)-m+1,1

    if (pos == Np.and.counter == 1) then
      indices(i,Np) = elements(j)
      counter = counter+1

      i = i+1
    else if  (pos == Np.and.counter /= 1) then
      indices(i,1:Np-1) = indices(i-1,1:Np-1)
      indices(i,Np) = elements(j)

      i = i+1
    else if (pos /=Np .and. counter == 1) then

      indices(i,pos) = elements(j)

      call recursive_test(Nc,Np,indices,length,i,j+m,pos+1)

      counter = counter+1

    else if (pos /=Np .and. counter /= 1) then

      indices(i,1:pos-1) = indices(i-1,1:pos-1)
      indices(i,pos) = elements(j)

      call recursive_test(Nc,Np,indices,length,i,j+m,pos+1)


    end if

  end do

end subroutine recursive_test


subroutine find_intersec_vertices(K,b,beta_inter,Nc,Np,l,length)
  implicit none
  integer, intent(in) :: Np,Nc,length
  Real*8, intent(in) :: K(Nc,Np-1),b(Nc)
  real*8, intent(out) :: beta_inter(Np-1,length)
  integer, intent(out) :: l
  logical :: test
  real*8 :: A(Np-1,Np-1),bi(Np-1),temp(Np-1)
  integer :: i,j,indices(length,Np)


  l=0
  i = 1
  call recursive_test(Nc,Np-1,indices,length,i,1,1)


  do i = 1, length, 1
    do j = 1, Np-1, 1
      A(j,:) = 1.D0-K(indices(i,j),:)
      bi(j) =b(indices(i,j))

    end do


    call solve_matmul_by_LU(A,bi,Np-1,temp)

    call check_intersec_point_range(1.D0-K,temp,b,test,Nc,Np)

    if (test) then
      ! write(*,*) "The intersection point number",i, "does not belong to sets S and P"
    else
      ! write(*,'(1X,A39,I3,A28)') "Success : The intersection point number",i, "does belong to sets S and P."
      l = l+1
      beta_inter(:,l) = temp

    end if


  end do

end subroutine

subroutine check_intersec_point_range(A,Beta,b,result,Nc,Np)
  implicit none
  integer, intent(in) :: Nc,Np
  real*8,intent(in) :: A(Nc,Np-1),Beta(Np-1),b(Nc)
  logical,intent(out) :: result
  real*8 :: dummy
  integer :: i,j

  result = .false.

  do i = 1, Nc, 1
    dummy = 0.D0
    do j = 1, Np-1, 1
      dummy = dummy + A(i,j)*beta(j)
    end do


    if (dummy .gt. b(i)+1D-12) then
      result = .true.

      return
    end if

  end do

end subroutine check_intersec_point_range
!
!
subroutine solve_matmul_by_LU(A,b,Nc,sol)
  ! A is the coeff matrix, Nc by Nc
  ! b is the variable vector
  ! sol is the solution vector
  ! Nc is the size of the square matrix and vectors
  ! This solves the eqation Asol=b, where sol is the unknown
  implicit none
  integer, intent(in):: Nc ! This is the size of the matrix
  real*8, intent(in) :: A(Nc,Nc),b(Nc)
  real*8, intent(out) :: sol(Nc)
  real*8 :: d(Nc),L(Nc,Nc),U(Nc,Nc)
  integer :: i,j
!
  U = A
  L = 0D0
  forall(i = 1:Nc,j=1:Nc,i==j)
        L(i,j)=1
  end forall
  d = 0D0
!
! First, create the inverse L and U matrices such that A = L*U
  do j = 1, Nc
    do i = j+1, Nc, 1
      L(i,j) = U(i,j)/U(j,j)
      U(i,:) = U(i,:)-U(j,:)*L(i,j)
    end do
  end do
!
! Then, create vector d such that L*d = b
  do i = 1,Nc
    d(i) = b(i)/L(i,i)
    do j = 1,i-1,1
      d(i) = d(i)-d(j)*L(i,j)/L(i,i)
    end do
  end do
!
! Finally, just solve the equation U*x = d
  do i = Nc,1,-1
    sol(i) = d(i)/U(i,i)
    do j = i+1,Nc
      sol(i) = sol(i)-sol(j)*U(i,j)/U(i,i)
    end do
  end do
!
  return
!
end subroutine


!____________________________________________________________________________________________________
! These are the subroutines used for step 2
!
!
subroutine RR_gradient_calc(z,K,Beta,f_vector,invt,Nc,Np)
!            This subroutine calculates the RR rice factor, the gradient of F
  implicit none

  integer, intent(in) :: Np,Nc
  real*8, intent(in) :: z(Nc),K(Nc,Np-1),Beta(Np-1)
  real*8, intent(out) :: f_vector(Np-1),invt(Nc)
  real*8 :: dummy
  integer*8 :: i,j

  calc_invti: do i = 1, Nc, 1
    dummy= 1.D0

    calc_ti : do j = 1, Np-1, 1
      dummy = dummy-(1.D0-K(i,j))*beta(j)
    end do calc_ti

    invt(i) = 1.D0/dummy
  end do calc_invti

  calc_fvect: do j = 1, Np-1, 1
  f_vector(j) = 0.D0

    calc_fj: do i = 1, Nc, 1
      f_vector(j) = f_vector(j) + (1.D0 - K(i,j))*z(i)*invt(i)
    end do calc_fj

  end do calc_fvect

  return
end subroutine RR_gradient_calc

!____________________________________________________________________________________________________
! These are the subroutines used for step 3
!

subroutine Hessian_calc(Km,z,invt,Hessian,Nc,Np)
  implicit none
  integer, intent(in) :: Np,Nc
  real*8, intent(in) :: Km(Nc,Np-1),z(Nc),invt(Nc)
  real*8, intent(out) :: Hessian(Np-1,Np-1)
  integer :: i,j,k

  j_loop: do j = 1, Np-1, 1

    k_loop: do k = 1, Np-1, 1

      Hessian(j,k) = 0.D0
      iloop: do i = 1, Nc, 1
        Hessian(j,k) = Hessian(j,k)+(1.D0-Km(i,j))*(1.D0-Km(i,k))*z(i)*invt(i)*invt(i)
      end do iloop
    end do k_loop
  end do j_loop


end subroutine Hessian_calc
!
!
!
!
!____________________________________________________________________________________________________
! These are the subroutines used for step 4
!
!
subroutine lambda_max(K,b,Beta,d,lmax,Nc,Np)
  implicit none
  Integer, intent(in) :: Nc,Np
  Real*8, intent(in):: K(Nc,Np-1),b(Nc),Beta(Np-1),d(Np-1)
  real*8, intent(out) :: lmax
  integer :: i,j,n
  real*8 :: newl,A(Nc,Np-1),dot1,dot2,invdot2

  A = 1.D0-K
  n = 0
  dot1 = 0.D0
  dot2 = 0.D0

  ! start by initiating the value of lambda
  do j = 1, Np-1, 1
    dot1 = dot1 + A(1,j)*beta(j)
    dot2 = dot2 + A(1,j)*d(j)
  end do

  invdot2 = 1.D0 / dot2
  newl = (b(1) - dot1)*invdot2

  if (n == 0 .and. (invdot2>0)) then
    n = n+1
    lmax = newl

  else if (n == 0 .and. (invdot2<0)) then

    lmax = 1.D30
    n = n+1
  end if

  !Now, initiate the loop
  Main_loop: do i = 2, Nc, 1
    dot1 = 0.D0
    dot2 = 0.D0
    do j = 1, Np-1, 1
      dot1 = dot1 + A(i,j)*beta(j)
      dot2 = dot2 + A(i,j)*d(j)
    end do

    invdot2 = 1.D0 / dot2
    newl = (b(i) - dot1)*invdot2


    if (newl < lmax .and. (invdot2>0)) then
      lmax = newl

    end if
  end do Main_loop


  if (lmax >=1.D0) then
    lmax = 1.D0
  end if

end subroutine lambda_max
!
!
!
!
!____________________________________________________________________________________________________
! These are the subroutines used for step 5
!
!
!
subroutine line_search(z,K,beta,d,lmax,s,Nc,Np)
  implicit none
  integer, Intent(in) :: Np,Nc
  real*8, intent(in) :: z(Nc),K(Nc,Np-1),Beta(Np-1),d(Np-1),lmax
  real*8, intent(out) :: s
  integer :: i,j,iter
  real*8 :: f(Np-1),Hessian(Np-1,Np-1),invt(Nc),newbeta(Np-1),d1,d2,invd2,temp(Np-1)

  s = 1.D0

  do iter = 1, 10, 1

    do i = 1, Np-1, 1
      newbeta(i) = beta(i)+s*lmax*d(i)
    end do

    call RR_gradient_calc(z,K,newbeta,f,invt,Nc,Np)
    call Hessian_calc(K,z,invt,Hessian,Nc,Np)

    d1 = 0.D0
    do i = 1, Np-1, 1
      d1 = (d1 + f(i)*d(i))*lmax
    end do

    if(iter ==1 .and. d1 <=0) return

    d2 = 0.D0
    do i = 1, Np-1, 1

      temp(i) = 0.D0

      do j = 1, Np-1, 1

        temp (i) = temp(i) + Hessian(i,j)*d(j)
      end do
    end do

    do i = 1, Np-1, 1
      d2 = d2 + lmax*lmax*d(i)*temp(i)
    end do

    invd2 = 1.D0 / d2

    s = s - d1*invd2

  end do

end subroutine line_search
!
!
!
!
subroutine print_matrix(matrix,m,n)
  implicit none
  integer, intent(in) :: m,n
  Real*8, intent(in) :: matrix(m,n)
  integer :: i

  do i = 1, m, 1
    print *, "( " ,matrix(i,:), " )"
  end do

end subroutine print_matrix

subroutine print_matrixi(matrix,m,n)
  implicit none
  integer, intent(in) :: m,n
  integer, intent(in) :: matrix(m,n)
  integer :: i

  do i = 1, m, 1
    print *, "( " ,matrix(i,:), " )"
  end do

end subroutine print_matrixi
