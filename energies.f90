!------------------------------------------------------------------------------
! MODULE: Modulename
!------------------------------------------------------------------------------
! DESCRIPTION: 
!> @brief
!!This module computes the total energy and the various contributions to
!!it in two ways. The first method evaluates the density functional, 
!!by direct integration to compute the <em> integrated
!!energy</em> \c ehfint. The second method uses the sum of
!!single-particle energies plus the rearrangement energies, \f$ E_{3,\rm
!!corr} \f$ for the density dependent part and \f$ E_{C,\rm corr} \f$ for
!!Coulomb exchange.
!>
!>@details 
!!The two ways of calculating the energy are assigned to the subroutines
!!\c integ_energy, which also calculates the rearrangement energies,
!!and \c sum_energy. Note that since \c integ_energy is always
!!called briefly before \c sum_energy, the rearrangement energies
!!are correctly available.
!!
!!In addition subroutine \c sum_energies also calculates the summed
!!spin, orbital, and total angular momenta.
!!
!!TODO: SHOULD WE ADD THE EQUATIONS FOR THE ENERGIES HERE?
!------------------------------------------------------------------------------
MODULE Energies
   USE Params, ONLY: db,tcoul,b_d17,C_0,C_1,a0,a1,a2,a3,alpha,a4,h2ml,maxiter,serr,sumflu,iter,prho
  USE Forces
  USE Densities
  USE Levels
  USE Pairs, ONLY: epair
  IMPLICIT NONE
  SAVE
  ! total energies calculated in subroutine "energy" and "hfenergy"
  REAL(db) :: ehft            !< kinetic energy
  REAL(db) :: ehf0            !< t0 contribution
  REAL(db) :: ehf1,ehf1_even,ehf1_odd            !< b1 contribution (current part)
  REAL(db) :: ehf2            !< b2 contribution (Laplacian part)
  REAL(db) :: ehf3            !< t3 contribution which models density dependence
  REAL(db) :: ehfls           !< spin-orbit contribution (time even)
  REAL(db) :: ehflsodd        !< spin-orbit contribution (odd-odd)
  REAL(db) :: ehfc            !< Coulomb contribution
  REAL(db) :: ecorc           !< Slater & Koopman exchange
  REAL(db) :: ehfint          !< integrated total energy
  REAL(db) :: efluct1         !< energyfluctuation \f$ \hat{h}^2 \f$
  REAL(db) :: efluct1prev     !< energyfluctuation \f$ \hat{h}^2 \f$ of previous iteration
  REAL(db) :: efluct2         !< fluctuation \f$ \hat{h}\cdot {\tt efluct}\f$
  REAL(db) :: efluct2prev     !< fluctuation \f$ \hat{h}\cdot {\tt efluct}\f$ of previous iteration
  REAL(db) :: tke             !< kinetic energy summed
  REAL(db) :: ehf             !< Hartree-Fock energy from s.p. levels
  REAL(db) :: ehfprev         !< Hartree-Fock energy from s.p. levels of previous iteration
  REAL(db) :: e3corr          !< rearrangement energy
  REAL(db) :: orbital(3)      !< the three components of the total orbital
                              !!angular momentum in units of \f$ \hbar \f$.
  REAL(db) :: spin(3)         !< the three components of the total spin in
                              !!units of \f$ \hbar \f$.
  REAL(db) :: total_angmom(3) !< the three components of the total
                              !!angular momentum in units of \f$ \hbar \f$.
  REAL(db) :: Zeeman_spin, Zeeman_orb,Zeeman_A2 !Zeeman energy (a) spin part, (b) orbital part, (c) diamagnetic part
 ! REAL(db) :: zeeman__spin_2!the second method to calculate the spin part of Zeeman energy
   REAL(db):: e0odd,e3odd,efinodd,Eeffodd,egradsodd,etodd !time-odd contributions to the energy
   REAL(db):: ehf_lambda,Zeeman_lambda ! Lambda functional contribution
CONTAINS
!---------------------------------------------------------------------------  
! DESCRIPTION: integ_energy
!> @brief
!!The purpose of this subroutine is to calculate the integrated energy.
!>
!> @details
!!This is implemented pretty straightforwardly. The only
!!programming technique worth noting is that intermediate variables such
!!as \c rhot for the total density are used to avoid repeating the
!!lengthy index lists. Compilers will eliminate these by optimization.
!!
!!In principle the integration loops in the subroutine could be
!!combined, but some space is saved by using the array \c worka for
!!different purposes in different loops.
!!
!!The calculation proceeds in the following steps:
!!  - <b> Step 1:</b> the Laplacian of the densities is calculated in
!!    \c worka, then the integrals for 
!!    \c ehf0,\c ehf2, and \c ehf3 are performed.
!!    After the loop the result for \c ehf3 is also used to calculate
!!    \c e3corr.
!!  - <b> Step 2:</b> the integral for \c ehf1 is evaluated
!!    using \c worka for the \f$ \vec\jmath_q{}^2 \f$ term.
!!  - <b> Step 3:</b> the spin-orbit contribution of 
!!    \c ehfls is calculated using \c worka as storage for
!!    \f$ \nabla\cdot\vec J_q \f$.
!!  - <b> Step 4:</b> the Coulomb energy \c ehfc is evaluated 
!!    with the Slater correction taken into account if the
!!    force's \c ex is nonzero. At the same time the Coulomb correction
!!    for the summed energy is calculated
!!    and stored in \c ecorc. It will be used in
!!    the subroutine \c sum_energy.
!!  - <b> Step 5: </b> the kinetic energy is integrated for 
!!    \c ehft. Note that only at this point the correct prefactor
!!    \f$ \hbar^2/2m \f$ is added; the use of \c tau in other expressions
!!    assumes its absence.
!!  - <b> Step 6: </b> Finally all terms are added to produce the total
!!  energy, \c efundet, from which the pairing energies are subtracted.
!--------------------------------------------------------------------------- 
  SUBROUTINE integ_energy
    USE Trivial, ONLY: rmulx,rmuly,rmulz
    USE Grids, ONLY: wxyz,der1x,der2x,der1y,der2y,der1z,der2z,x,y,z
    USE Coulomb, ONLY: wcoul
    USE Moment, ONLY:cmtot
   INTEGER :: ix,iy,iz,iq,i,ic
    REAL(db) :: rhot,rhon,rhop,d2rho,d2rhon,d2rhop,sc
    REAL(db) :: worka(nx,ny,nz,3),current_t(nx,ny,nz,2)
   REAL(db) :: div2rho(nx,ny,nz,3)
    REAL(db) :: workb(nx,ny,nz,3,2)
    REAL(db) ::factor,muN,gq(3)
    REAL(db) :: lzcoef,a2coef,xx(nx),yy(ny)
    REAL(db) :: d2s(nx,ny,nz,3,2)
    REAL(db) :: s2t, s2n, s2p
   REAL(db) :: A_fin, B_fin, A_st, B_st
    REAL(db) :: tmp_sn_d2sn, tmp_sp_d2sp, tmp_sn_d2sp, tmp_sp_d2sn
    REAL(db) :: jtot2, jn2, jp2

    real(db) :: dr, r_max,V_r(100), r_count(100),r_bin,r_center(100),V_lambda(nx,ny,nz)
    integer :: n_r=100, bin
    CHARACTER(len=200) :: fname
    ! Step 1: compute laplacian of densities, then ehf0, ehf2, and ehf3
    DO iq=1,3
       CALL rmulx(der2x,rho(:,:,:,iq),div2rho(:,:,:,iq),0)
       CALL rmuly(der2y,rho(:,:,:,iq),div2rho(:,:,:,iq),1)
       CALL rmulz(der2z,rho(:,:,:,iq),div2rho(:,:,:,iq),1)
    ENDDO
    ehf0=0.0D0
    ehf3=0.0D0
    ehf2=0.0D0
    DO iz=1,nz
       DO iy=1,ny
          DO ix=1,nx
             rhot=rho(ix,iy,iz,1)+rho(ix,iy,iz,2)
             rhon=rho(ix,iy,iz,1)
             rhop=rho(ix,iy,iz,2)
             ehf0=ehf0+wxyz*(b0*rhot**2-b0p*(rhop**2+rhon**2))/2.D0
             ehf3=ehf3+wxyz*rhot**f%power*(b3*rhot**2 &
                  -b3p*(rhop**2+rhon**2))/3.D0
             d2rho=div2rho(ix,iy,iz,1)+div2rho(ix,iy,iz,2)
             d2rhon=div2rho(ix,iy,iz,1)
             d2rhop=div2rho(ix,iy,iz,2)
             ehf2=ehf2+wxyz*(-b2*rhot*d2rho+b2p*(rhop* &
                  d2rhop+rhon*d2rhon))/2.D0
          ENDDO
       ENDDO
    ENDDO
    e3corr=-f%power*ehf3/2.D0
    ! Step 2: b1-contribution. worka=square of the current vector
    current_t=current(:,:,:,1,:)**2+current(:,:,:,2,:)**2+current(:,:,:,3,:)**2
    ehf1_even=wxyz*SUM(b1*(rho(:,:,:,1)+rho(:,:,:,2))*(tau(:,:,:,1)+tau(:,:,:,2))&
         -b1p*(rho(:,:,:,1)*tau(:,:,:,1)+rho(:,:,:,2)*tau(:,:,:,2)))

    ehf1_odd=wxyz*SUM(-b1*((current(:,:,:,1,1)+current(:,:,:,1,2))**2 &
         + (current(:,:,:,2,1)+current(:,:,:,2,2))**2 &
         + (current(:,:,:,3,1)+current(:,:,:,3,2))**2)&
         -b1p*(-current_t(:,:,:,1)-current_t(:,:,:,2)))
    ehf1=ehf1_even+ehf1_odd

    ! Step 3: the spin-orbit contribution
    !         (a) Time even part: worka=div J
    DO iq=1,3
       CALL rmulx(der1x,sodens(:,:,:,1,iq),worka(:,:,:,iq),0)
       CALL rmuly(der1y,sodens(:,:,:,2,iq),worka(:,:,:,iq),1)
       CALL rmulz(der1z,sodens(:,:,:,3,iq),worka(:,:,:,iq),1)
    ENDDO
    ehfls=wxyz*SUM(-b4*(rho(:,:,:,1)+rho(:,:,:,2)) &
         *(worka(:,:,:,1)+worka(:,:,:,2)) &
         -b4p *(rho(:,:,:,2)*worka(:,:,:,2)+rho(:,:,:,1)*worka(:,:,:,1)))
    ! Step 3: the spin-orbit contribution
    !         (b) odd-odd part: workb = curl J
    DO iq = 1,2
       CALL rmuly(der1y,current(:,:,:,3,iq),workb(:,:,:,1,iq),0)
       CALL rmulz(der1z,current(:,:,:,2,iq),workb(:,:,:,1,iq),-1)
       CALL rmulz(der1z,current(:,:,:,1,iq),workb(:,:,:,2,iq),0)
       CALL rmulx(der1x,current(:,:,:,3,iq),workb(:,:,:,2,iq),-1)
       CALL rmulx(der1x,current(:,:,:,2,iq),workb(:,:,:,3,iq),0)
       CALL rmuly(der1y,current(:,:,:,1,iq),workb(:,:,:,3,iq),-1)
    ENDDO
    ehflsodd=wxyz*SUM(-b4*(sdens(:,:,:,:,1)+sdens(:,:,:,:,2))* &
         (workb(:,:,:,:,1)+workb(:,:,:,:,2))-b4p*(sdens(:,:,:,:,1)* &
         workb(:,:,:,:,1)+sdens(:,:,:,:,2)*workb(:,:,:,:,2)))
    ehfls = ehfls + ehflsodd
    !
    ! Step 4: Coulomb energy with Slater term, also correction
    ! term to Koopman formula
    ehfc=0.0D0
    ecorc=0.0D0
    IF(tcoul) THEN
       IF(f%ex/=0) THEN
          sc=-3.0D0/4.0D0*slate
       ELSE
          sc=0.D0
       END IF
       DO iz=1,nz
          DO iy=1,ny
             DO ix=1,nx
                rhop=rho(ix,iy,iz,2)
                ehfc=ehfc+wxyz *(0.5D0*rhop*wcoul(ix,iy,iz) &
                     +sc*rhop**(4.0D0/3.0D0))
                ecorc=ecorc+wxyz*sc/3.0D0*rhop**(4.0D0/3.0D0)
             ENDDO
          ENDDO
       ENDDO
    ENDIF

       ! Step 4: Zeeman energy
       !         (a) spin part
       !         (b) orbital part
       !         (c) diamagnetic part
       !         (d) lambda spin part

   !(a) spin part
   zeeman_spin = 0.0_db
   gq(1)=-3.8263
   gq(2)=5.5856
   gq(3)=-1.226
   muN=0.31524512
   DO iq = 1, 2
         factor = -gq(iq) * muN * b_d17
         DO iz = 1, nz
            DO iy = 1, ny
               DO ix = 1, nx
                  zeeman_spin = zeeman_spin + factor * sdens(ix,iy,iz,3,iq) * wxyz/2.0
               END DO
            END DO
         END DO
   END DO

   !zeeman_spin_2=0d0
   ! DO i = 1, npsi(1)              ! pnrtot 是总轨道数或类似的自旋数组长度
   !   zeeman_spin_2 = zeeman_spin_2 + ( -gq(1) * muN * b_d17 )  * wocc(i) * sp_spin(3,i)
   ! END DO
    !DO i = npsi(1)+1, nstmax              ! pnrtot 是总轨道数或类似的自旋数组长度
   !   zeeman_spin_2 = zeeman_spin_2  + ( -gq(2) * muN * b_d17 )  * wocc(i) * sp_spin(3,i)
   ! END DO
   !(b) orbital part
   Zeeman_orb = 0.0_db
   DO i = npsi(1)+1, nstmax
      Zeeman_orb = Zeeman_orb - 0.31523022d0 * b_d17 * wocc(i) * sp_orbital(3,i)
   END DO

   !(c) diamagnetic part
   xx = x - cmtot(1)
   yy = y - cmtot(2)
   lzcoef = muN * b_d17
   a2coef = (lzcoef*lzcoef) / (4.D0*f%h2m(2))

   zeeman_A2 = 0.D0
    DO iz=1,nz
      DO iy=1,ny
        DO ix=1,nx
          zeeman_A2 = zeeman_A2 + wxyz * a2coef * &
               (xx(ix)*xx(ix)+yy(iy)*yy(iy)) * rho(ix,iy,iz,2)
        ENDDO
      ENDDO
    ENDDO
!(d) lambda spin part
zeeman_lambda = 0.0
    DO i = npsi(2)+1, nstmax            
      zeeman_lambda = zeeman_lambda  + (-gq(3) * muN * b_d17 )  * wocc(i) * sp_spin(3,i)
    END DO



    ! Step 5: kinetic energy contribution
    ehft=wxyz*SUM(f%h2m(1)*tau(:,:,:,1)+f%h2m(2)*tau(:,:,:,2))



 !=========================================================
    !  Step 0: initialize
    !=========================================================
    e0odd    = 0.0D0
    e3odd    = 0.0D0
    efinodd  = 0.0D0
    Eeffodd  = 0.0D0
   egradsodd= 0.0D0

    !=========================================================
    !  Step 1: E0^odd and E3^odd (keep your original forms)
    !    E0odd = ∫ d^3r  (1/4)t0 [ x0*s^2 - s_n^2 - s_p^2 ]
    !    E3odd = ∫ d^3r (1/24)t3 rho^gamma [ x3*s^2 - s_n^2 - s_p^2 ]
    !=========================================================
    DO iz=1,nz
       DO iy=1,ny
          DO ix=1,nx
             rhot = rho(ix,iy,iz,1) + rho(ix,iy,iz,2)
             s2n = sdens(ix,iy,iz,1,1)**2 + sdens(ix,iy,iz,2,1)**2 + sdens(ix,iy,iz,3,1)**2
             s2p = sdens(ix,iy,iz,1,2)**2 + sdens(ix,iy,iz,2,2)**2 + sdens(ix,iy,iz,3,2)**2
             s2t = (sdens(ix,iy,iz,1,1)+sdens(ix,iy,iz,1,2))**2 &
                 + (sdens(ix,iy,iz,2,1)+sdens(ix,iy,iz,2,2))**2 &
                 + (sdens(ix,iy,iz,3,1)+sdens(ix,iy,iz,3,2))**2

             e0odd = e0odd + wxyz * 0.25D0*f%t0 * ( f%x0*s2t - s2n - s2p )
             e3odd = e3odd + wxyz * (1.D0/24.D0)*f%t3 * rhot**f%power * ( f%x3*s2t - s2n - s2p )
          ENDDO
       ENDDO
    ENDDO



    !=========================================================
    !  Step 3: E_fin^odd (keep your given t,x form)
    !    H_fin^odd =
    !      A_fin*( s_n·∇^2 s_n + s_p·∇^2 s_p )
    !    + B_fin*( s_n·∇^2 s_p + s_p·∇^2 s_n )
    !=========================================================
   ! DO iq=1,2
   !    DO ic=1,3
   !       CALL rmulx(der2x,sdens(:,:,:,ic,iq),d2s(:,:,:,ic,iq),0)
   !       CALL rmuly(der2y,sdens(:,:,:,ic,iq),d2s(:,:,:,ic,iq),1)
   !       CALL rmulz(der2z,sdens(:,:,:,ic,iq),d2s(:,:,:,ic,iq),1)
   !    ENDDO
   ! ENDDO

    !A_fin = (1.D0/32.D0) * ( 3.D0*f%t1*(1.D0-f%x1) + f%t2*(1.D0+f%x2) )
    !B_fin = (1.D0/32.D0) * ( f%t2*f%x2 - 3.D0*f%t1*f%x1 )

    !DO iz=1,nz
    !   DO iy=1,ny
    !      DO ix=1,nx
     !        tmp_sn_d2sn = sdens(ix,iy,iz,1,1)*d2s(ix,iy,iz,1,1) &
    !                     + sdens(ix,iy,iz,2,1)*d2s(ix,iy,iz,2,1) &
     !                    + sdens(ix,iy,iz,3,1)*d2s(ix,iy,iz,3,1)
    !         tmp_sp_d2sp = sdens(ix,iy,iz,1,2)*d2s(ix,iy,iz,1,2) &
      !                   + sdens(ix,iy,iz,2,2)*d2s(ix,iy,iz,2,2) &
      !                   + sdens(ix,iy,iz,3,2)*d2s(ix,iy,iz,3,2)
      !       tmp_sn_d2sp = sdens(ix,iy,iz,1,1)*d2s(ix,iy,iz,1,2) &
      !                   + sdens(ix,iy,iz,2,1)*d2s(ix,iy,iz,2,2) &
      !                   + sdens(ix,iy,iz,3,1)*d2s(ix,iy,iz,3,2)
     !        tmp_sp_d2sn = sdens(ix,iy,iz,1,2)*d2s(ix,iy,iz,1,1) &
     !                    + sdens(ix,iy,iz,2,2)*d2s(ix,iy,iz,2,1) &
     !                    + sdens(ix,iy,iz,3,2)*d2s(ix,iy,iz,3,1)
    !         efinodd = efinodd + wxyz * ( A_fin*(tmp_sn_d2sn + tmp_sp_d2sp) &
        !                                + B_fin*(tmp_sn_d2sp + tmp_sp_d2sn) )
    !      ENDDO
    !   ENDDO
    !ENDDO

    !=========================================================
    !  Step 4: E_sT^odd = A_st * s·T (isoscalar) + B_st * s_q·T_q (isovector-like)
    !=========================================================
    A_st = (1.D0/8.D0) * (f%t1*f%x1 + f%t2*f%x2)
    B_st = (1.D0/8.D0) * (f%t2 - f%t1)

    DO iz=1,nz
       DO iy=1,ny
          DO ix=1,nx
             Eeffodd = Eeffodd + wxyz * (&
                A_st * (& 
                   (sdens(ix,iy,iz,1,1)+sdens(ix,iy,iz,1,2))*(tdens(ix,iy,iz,1,1)+tdens(ix,iy,iz,1,2)) +& 
                   (sdens(ix,iy,iz,2,1)+sdens(ix,iy,iz,2,2))*(tdens(ix,iy,iz,2,1)+tdens(ix,iy,iz,2,2)) +& 
                   (sdens(ix,iy,iz,3,1)+sdens(ix,iy,iz,3,2))*(tdens(ix,iy,iz,3,1)+tdens(ix,iy,iz,3,2)) ) +& 
                B_st * (& 
                   sdens(ix,iy,iz,1,1)*tdens(ix,iy,iz,1,1) + sdens(ix,iy,iz,2,1)*tdens(ix,iy,iz,2,1) + &
                   sdens(ix,iy,iz,3,1)*tdens(ix,iy,iz,3,1) +& 
                   sdens(ix,iy,iz,1,2)*tdens(ix,iy,iz,1,2) + sdens(ix,iy,iz,2,2)*tdens(ix,iy,iz,2,2) +&
                   sdens(ix,iy,iz,3,2)*tdens(ix,iy,iz,3,2) ) )
          ENDDO
       ENDDO
    ENDDO

    !=========================================================
    !  Step 5: Edeltas
    !=========================================================

    !DO iq=1,2
    !   CALL rmulx(der1x,sdens(:,:,:,1,iq),worka(:,:,:,iq),0)
    !   CALL rmuly(der1y,sdens(:,:,:,2,iq),worka(:,:,:,iq),1)
    !   CALL rmulz(der1z,sdens(:,:,:,3,iq),worka(:,:,:,iq),1)
   ! ENDDO

   ! egradsodd = wxyz*SUM(C_0*(worka(:,:,:,1)+worka(:,:,:,2))**2 + C_1*(worka(:,:,:,1)-worka(:,:,:,2))**2)


    !=========================================================
   !  Step 6: total time-odd
    !=========================================================
   etodd = e0odd + e3odd  + Eeffodd! + efinodd+ egradsodd
   !=========================================================
   !  Step 7: lambda  contribution 
   !=========================================================
   ehf_lambda=0.0D0
   DO iz=1,nz
      DO iy=1,ny
         DO ix=1,nx
            ehf_lambda=ehf_lambda+wxyz*( &
                 h2ml*tau(ix,iy,iz,3) &
                 + a0*rho(ix,iy,iz,3)*(rho(ix,iy,iz,1)+rho(ix,iy,iz,2)) &
                 + a1*(rho(ix,iy,iz,3)*(tau(ix,iy,iz,1)+tau(ix,iy,iz,2)) &
                 + (rho(ix,iy,iz,1)+rho(ix,iy,iz,2))*tau(ix,iy,iz,3)) &
                 - 0.5D0*a2*(rho(ix,iy,iz,3)*(div2rho(ix,iy,iz,1)+div2rho(ix,iy,iz,2))&
                 + (rho(ix,iy,iz,1)+rho(ix,iy,iz,2))*div2rho(ix,iy,iz,3)) &
                 + a3*rho(ix,iy,iz,3)*(rho(ix,iy,iz,1)+rho(ix,iy,iz,2))**(1.0D0+alpha) &
                 - a4*(rho(ix,iy,iz,3)*(worka(ix,iy,iz,1)+worka(ix,iy,iz,2)) &
                 + (rho(ix,iy,iz,1)+rho(ix,iy,iz,2))*worka(ix,iy,iz,3)) )
         ENDDO
      ENDDO
   ENDDO
      ! printed  V_lambda-r
      !lambda potential
   if(((iter==maxiter).or.(sumflu/nstmax<serr.AND.(iter>1))).and.prho) then
      V_lambda(:,:,:) = a0*(rho(:,:,:,1)+rho(:,:,:,2)) &
                       + a1*(tau(:,:,:,1)+tau(:,:,:,2)) &
                       - a2*(div2rho(:,:,:,1)+div2rho(:,:,:,2)) &
                       + a3*(rho(:,:,:,1)+rho(:,:,:,2))**(1.0D0+alpha) &
                       - a4*(worka(:,:,:,1)+worka(:,:,:,2)) 
      ! 初始化
      V_r = 0.0
      r_count = 0.0

      ! 最大半径取 3D 网格范围
      r_max = max(max(x(nx), y(ny)), z(nz))
      dr = r_max / (n_r - 1)   ! 这里减 1 保证最后一个点是 r_max

      ! 第一个点 r=0，对应中心
      r_center(1) = 0.0
      V_r(1) = V_lambda(nx/2, ny/2, nz/2)
      r_count(1) = 1.0

      ! 遍历 3D 网格
      do iz = 1, nz
         do iy = 1, ny
            do ix = 1, nx
               r_bin = sqrt(x(ix)**2 + y(iy)**2 + z(iz)**2)
               ! 找到对应 bin
               bin = int(r_bin / dr) + 1
               if (bin < 2) bin = 2         ! 第一个 bin 已经处理
               if (bin > n_r) bin = n_r
               ! 累加密度
               V_r(bin) = V_r(bin) + V_lambda(ix, iy, iz)
               r_count(bin) = r_count(bin) + 1.0
            end do
         end do
      end do

      ! 求每个 bin 平均密度
      do bin = 2, n_r
         if (r_count(bin) > 0.0) V_r(bin) = V_r(bin) / r_count(bin)
         r_center(bin) = (bin-1)*dr
      end do

      ! 输出文件名
      if(b_d17 < 10) then
         WRITE(fname,'(A,f4.2,A)') 'V_lambda', b_d17, '.txt'
      else
         WRITE(fname,'(A,f5.2,A)') 'V_lambda', b_d17, '.txt'
      end if

      ! 输出两列 r, V_lambda
      OPEN(UNIT=99, FILE=fname, STATUS='REPLACE', FORM='FORMATTED')
      do bin = 1, n_r
         WRITE(99,'(F12.6,1X,F12.6)') r_center(bin), V_r(bin)
      end do
      CLOSE(99)
   end if
!=========================================================










   ehfint=ehft+ehf0+ehf1+ehf2+ehf3+ehfls+ehfc+ehf_lambda-epair(1)-epair(2)+&
          Zeeman_spin+Zeeman_orb+Zeeman_A2+etodd+Zeeman_lambda
  END SUBROUTINE integ_energy
!---------------------------------------------------------------------------  
! DESCRIPTION: sum_energy
!> @brief
!!This subroutine mainly computes the Koopman sum, but also
!!sums up a number of other single-particle properties. 
!>
!> @details
!!For systematics, the latter should be done in a different
!!place, but at present is left here.
!!
!!The summation of the total energy uses \c spenerg to compute 
!!\f[ \sum_k (\epsilon_k-\tfrac1{2}v_k)=\tfrac1{2}\sum_k(2t_k+v_k)=
!!\tfrac1{2}\sum_k(t_k+\epsilon_k). \f]
!!The last sum is calculated, the rearrangement corrections 
!!are added and the pairing energies subtracted.
!!
!!The subroutine then sums up the single-particle energy fluctuation
!!\c sp_efluct1 and \c sp_efluct2, dividing them by the nucleon
!!number.  Finally the orbital and spin angular
!!momentum components are summed to form the total ones. 
!--------------------------------------------------------------------------- 
  SUBROUTINE sum_energy
    USE Moment, ONLY: pnrtot
    INTEGER :: i
    ehf=SUM(wocc*(sp_kinetic+sp_energy))/2.D0+e3corr+ecorc &
         -epair(1)-epair(2)+Zeeman_spin+Zeeman_orb+Zeeman_A2+etodd
    tke=SUM(wocc*sp_kinetic)
    efluct1=SUM(wocc*sp_efluct1)/pnrtot
    efluct2=SUM(wocc*sp_efluct2)/pnrtot
    DO i=1,3
       orbital(i)=SUM(wocc*sp_orbital(i,:))
       spin(i)=SUM(wocc*sp_spin(i,:))
       total_angmom(i)=orbital(i)+spin(i)
    END DO

  END SUBROUTINE sum_energy
  !***************************************************************************
END MODULE Energies

