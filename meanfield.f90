  !------------------------------------------------------------------------------
! MODULE: Meanfield
!------------------------------------------------------------------------------
! DESCRIPTION: 
!> @brief
!!This module calculates all the ingredients needed for the
!!energy functional and for applying the single-particle Hamiltonian to
!!a wave function.
!>
!>@details
!!The work is done by two subroutines: \c skyrme for the calculation
!!of all the fields, which can be scalar or vector and
!!isospin-dependent. \c hpsi then is the routine applying the
!!single-particle Hamiltonian to one single-particle wave function.
!!
!!Note the division of labor between \c skyrme and
!!\c add_density of module \c Densities: everything that
!!constructs fields - densities and current densities - from the
!!single-particle wave functions is done in \c add_density, which is
!!called in a loop over the states, while \c skyrme does the further
!!manipulations to complete the fields entering the single-particle
!!Hamiltonian by combining the densities and their derivatives. It does
!!not need access to the wave functions.
!------------------------------------------------------------------------------
Module Meanfield
  USE Params, ONLY: db,tcoul,b_d17,C_0,C_1,a0,a1,a2,a3,alpha,a4,h2ml
  USE Densities
  USE Forces 
  USE Grids, ONLY: nx,ny,nz,der1x,der2x,der1y,der2y,der1z,der2z,dx,dy,x,y
  USE Coulomb, ONLY: poisson,wcoul
  USE Moment, ONLY:cmtot
  IMPLICIT NONE
  REAL(db),ALLOCATABLE,DIMENSION(:,:,:,:)   :: upot   !<this is the local part of the mean field 
  !!\f$ U_q \f$. It is a scalar field with isospin index.
  REAL(db),ALLOCATABLE,DIMENSION(:,:,:,:)   :: bmass  !<this is the effective mass \f$ B_q \f$.
  !!It is a scalar, isospin-dependent field.
  REAL(db),ALLOCATABLE,DIMENSION(:,:,:,:)   :: divaq  !<this is the divergence of \c aq,
  !!i.e., \f$ \nabla\cdot\vec A_q \f$. Its is a scalar, isospin-dependent field.
  REAL(db),ALLOCATABLE,DIMENSION(:,:,:,:,:) :: aq     !<This is the vector filed \f$ \vec A_q \f$. 
  !!It is a vector, isospin-dependent field.  
  REAL(db),ALLOCATABLE,DIMENSION(:,:,:,:,:) :: spot   !<the field \f$ \vec{S}_q \f$.
  !!It is a vector, isospin-dependent field.
  REAL(db),ALLOCATABLE,DIMENSION(:,:,:,:,:) :: wlspot !<the field \f$ \vec W_q \f$. 
  !!It is a vector, isospin-dependent field.
  REAL(db),ALLOCATABLE,DIMENSION(:,:,:,:,:) :: dbmass !<contains the gradient of \c bmass.
  !!It is a vector, isospin-dependent field.
  REAL(db),ALLOCATABLE,DIMENSION(:,:,:,:,:,:)   :: cq     !<the coefficient of the spin density in the S·σ term.
  PRIVATE :: divaq,aq,wlspot,dbmass
CONTAINS
!---------------------------------------------------------------------------  
! DESCRIPTION: alloc_fields
!> @brief
!!This subroutine has the simple task of allocating all the fields that
!!are local to the module \c Meanfield.
!--------------------------------------------------------------------------- 
  SUBROUTINE alloc_fields
    ALLOCATE(upot(nx,ny,nz,3),bmass(nx,ny,nz,3),divaq(nx,ny,nz,3), &
         aq(nx,ny,nz,3,3),spot(nx,ny,nz,3,3),wlspot(nx,ny,nz,3,3), &
         dbmass(nx,ny,nz,3,3),cq(nx,ny,nz,3,3,3))
    upot=0.D0
    bmass=0.D0
    wlspot=0.D0
    aq=0.D0
    divaq=0.D0
    dbmass=0.D0
  END SUBROUTINE alloc_fields
!---------------------------------------------------------------------------  
! DESCRIPTION: skyrme
!> @brief
!!In this subroutine the various fields are calculated from the
!!densities that were previously generated in module \c Densities.
!>
!> @details
!!The expressions divide up the
!!contributions into an isospin-summed part with \f$ b \f$-coefficients
!!followed by the isospin-dependent one with \f$ b' \f$-coefficients. As it
!!would be a waste of space to store the summed densities and currents,
!!the expressions are divided up more conveniently. If we denote the
!!isospin index \f$ q \f$ by \c 0 and the index for the opposite isospin
!!\f$ q' \f$ by \c ic (as \c iq can take the values 1 or 2, it can
!!conveniently be calculated as <tt> 3-iq </tt>), this can be written for
!!example as
!!
!!\f$ b_1\rho-b_1'\rho_q\longrightarrow b_1(\rho_q+\rho_{q'})-b_1'\rho_q \f$
!!\f$ \longrightarrow \f$ <tt>(b1-b1p)*rho(:,:,:,iq)+b1*rho(:,:,:,ic)</tt>
!!This decomposition is used in all applicable cases.
!!
!!For intermediate results the fields \c workden (scalar) and \c workvec 
!!(vector) are used.
!!
!!Now the subroutine proceeds in the following steps:
!!  -# all the parts in \c upot involving an \f$ \alpha \f$-dependent power of the
!!     density are collected. Note that in order to avoid having to
!!     calculate different powers, \f$ \rho^\alpha \f$ is factored out. The
!!     division by the total density uses the small number \c epsilon to
!!     avoid division by zero.
!!  -# the divergence of \f$ \vec J \f$ (\c sodens) is
!!     calculated for both isospins in \c workden and the contributions
!!     are added to \c upot.
!!  -# the Coulomb potential is calculated using
!!     subroutine \c poisson  (see module \c Coulomb). It and the
!!     Slater exchange correction (only if the \c ex parameter in the
!!     force is nonzero) are added to \c upot for protons, <tt> iq=2 </tt>.
!!  -# the Laplacian is applied to the densities and
!!     the result stored in \c workden. Then the remaining terms are constructed.  
!!     Note that the  \c iq -loop is combined with the following steps.
!!  -# the effective mass is calculated.
!!  -# the gradient of the density is calculated and the
!!     spin-orbit vector \f$ \vec W_q \f$ is constructed in \c wlspot.
!!  -# the curl of the spin density vector is calculated
!!     and stored in \c workvec.
!!  -# the vector \f$ \vec A_q \f$ is calculated from the current density 
!!     and the curl of the spin density.
!!  -# the curl of the current density is calculated and stored in \c spot.
!!  -# now the two isospin contributions in \c spot
!!     are combined in the proper way.
!!     This way of handling it avoids the introduction of an additional
!!     work vector for \f$ \nabla\times\vec\jmath_q \f$.
!!  -# the divergence of \f$ \vec A_q \f$ is calculated and stored in \c divaq.
!!  -# finally, the gradient of the effective mass term
!!     \f$ B_q \f$ is calculated and stored in the vector variable \c dbmass.
!!  .
!!This concludes the calculation of all scalar and vector fields needed
!!for the application of the Skyrme force. 
!--------------------------------------------------------------------------- 
  SUBROUTINE skyrme 
    USE Trivial, ONLY: rmulx,rmuly,rmulz
    REAL(db),PARAMETER :: epsilon=1.0d-25  
    REAL(db) :: rotspp,rotspn
    REAL(db),ALLOCATABLE :: workden(:,:,:,:),workvec(:,:,:,:,:),Laplace_sq(:,:,:,:,:)
    INTEGER :: ix,iy,iz,ic,iq,icomp,mu,nu
    ALLOCATE(workden(nx,ny,nz,3),workvec(nx,ny,nz,3,3),Laplace_sq(nx,ny,nz,3,2))
    !  Step 1: 3-body contribution to upot.
    DO iq=1,2  
       ic=3-iq  
       upot(:,:,:,iq)=(rho(:,:,:,1)+rho(:,:,:,2))**f%power * &
            ((b3*(f%power+2.D0)/3.D0-2.D0*b3p/3.D0)*rho(:,:,:,iq) &
            +b3*(f%power+2.D0)/3.D0*rho(:,:,:,ic) &
            -(b3p*f%power/3.D0)*(rho(:,:,:,1)**2+rho(:,:,:,2)**2)/ &
            (rho(:,:,:,1)+rho(:,:,:,2)+epsilon))
    ENDDO



    ! Step 2: add divergence of spin-orbit current to upot
    DO iq=1,3
       CALL rmulx(der1x,sodens(:,:,:,1,iq),workden(:,:,:,iq),0)
       CALL rmuly(der1y,sodens(:,:,:,2,iq),workden(:,:,:,iq),1)
       CALL rmulz(der1z,sodens(:,:,:,3,iq),workden(:,:,:,iq),1)
    ENDDO
    DO iq=1,2  
       ic=3-iq 
       upot(:,:,:,iq)=upot(:,:,:,iq) &
          -(b4+b4p)*workden(:,:,:,iq)-b4*workden(:,:,:,ic) &
            -a4*workden(:,:,:,3)  ! add ΛN interaction contribution part1 to upot
    ENDDO
   !setp 2.5: add hyperon interaction contribution part1 to upot(3)
    upot(:,:,:,3)=-a4*(workden(:,:,:,1)+workden(:,:,:,2))

    ! Step 3: Coulomb potential
    IF(tcoul) THEN
       CALL poisson(rho(:,:,:,2))
       upot(:,:,:,2)=upot(:,:,:,2)+wcoul
       IF(f%ex/=0) &
            upot(:,:,:,2)=upot(:,:,:,2)-slate*rho(:,:,:,2)**(1.0D0/3.0D0)
    ENDIF
    ! Step 4: remaining terms of upot
    DO iq=1,3
       CALL rmulx(der2x,rho(:,:,:,iq),workden(:,:,:,iq),0)  
       CALL rmuly(der2y,rho(:,:,:,iq),workden(:,:,:,iq),1)  
       CALL rmulz(der2z,rho(:,:,:,iq),workden(:,:,:,iq),1)
    ENDDO

  !step 4.25: add hyperon interaction contribution part2 to upot(3)
    upot(:,:,:,3)=upot(:,:,:,3)+a0*(rho(:,:,:,1) + rho(:,:,:,2)) &
                 + a1*(tau(:,:,:,1)+tau(:,:,:,2)) &
                 - a2*(workden(:,:,:,1)+workden(:,:,:,2))&
                 +a3*(rho(:,:,:,1)+rho(:,:,:,2))**(alpha+1.0D0)
                 !+2.0/3.0*a3*((rho(:,:,:,1)+rho(:,:,:,2))**2.0&
                 !+2.0*rho(:,:,:,1)*rho(:,:,:,2)) !the second one from ΛNN force

    !step 4.5: add ΛN interaction contributions part2 to upot(1) and upot(2)
    DO iq=1,2  
       ic=3-iq  
       upot(:,:,:,iq)=upot(:,:,:,iq)+(b0-b0p)*rho(:,:,:,iq)+b0*rho(:,:,:,ic) &
                                ! t1,t2, and tau-dependent part      !
            +(b1-b1p)*tau(:,:,:,iq)+b1*tau(:,:,:,ic) &
                                ! two-body laplacian*rho-dependent part
            -(b2-b2p)*workden(:,:,:,iq)-b2*workden(:,:,:,ic)&
            +a0*rho(:,:,:,3)+a1*tau(:,:,:,3)-a2*workden(:,:,:,3)&
            +a3*(1+alpha)*rho(:,:,:,3)*(rho(:,:,:,1)+rho(:,:,:,2))**alpha 
           ! +2.0/3.0*a3*(2.0*rho(:,:,:,iq)+4.0*rho(:,:,:,ic)) !the second one from ΛNN force
       ! Step 5: effective mass
       bmass(:,:,:,iq)=f%h2m(iq)+(b1-b1p)*rho(:,:,:,iq)+b1*rho(:,:,:,ic)&
       +a1*rho(:,:,:,3)  ! add ΛN interaction contribution to effective mass
    ENDDO
        ! Step 6: calculate grad(rho) and wlspot
    do iq=1,3
       CALL rmulx(der1x,rho(:,:,:,iq),workvec(:,:,:,1,iq),0)
       CALL rmuly(der1y,rho(:,:,:,iq),workvec(:,:,:,2,iq),0)
       CALL rmulz(der1z,rho(:,:,:,iq),workvec(:,:,:,3,iq),0)
    enddo
    DO iq=1,2
       ic=3-iq
       wlspot(:,:,:,:,iq)= &
            (b4+b4p)*workvec(:,:,:,:,iq)+b4*workvec(:,:,:,:,ic)&
            +a4*workvec(:,:,:,:,3)  ! add ΛN interaction contribution to spin-orbit potential
    END DO
      ! effective mass and spin-orbit potential for hyperon
    bmass(:,:,:,3)=h2ml +a1*(rho(:,:,:,1)+rho(:,:,:,2))  
    wlspot(:,:,:,:,3)=a4*(workvec(:,:,:,:,1)+workvec(:,:,:,:,2))  
    ! Step 7: calculate curl of spin density vector, store in workvec
    DO iq=1,2  
       CALL rmuly(der1y,sdens(:,:,:,3,iq),workvec(:,:,:,1,iq),0)
       CALL rmulz(der1z,sdens(:,:,:,2,iq),workvec(:,:,:,1,iq),-1)
       CALL rmulz(der1z,sdens(:,:,:,1,iq),workvec(:,:,:,2,iq),0)
       CALL rmulx(der1x,sdens(:,:,:,3,iq),workvec(:,:,:,2,iq),-1)
       CALL rmulx(der1x,sdens(:,:,:,2,iq),workvec(:,:,:,3,iq),0)
       CALL rmuly(der1y,sdens(:,:,:,1,iq),workvec(:,:,:,3,iq),-1)
    ENDDO
    ! Step 8: calculate A_q vector
    DO iq=1,2
       ic=3-iq
       aq(:,:,:,:,iq)=-2.0D0*(b1-b1p)*current(:,:,:,:,iq) &
            -2.0D0*b1*current(:,:,:,:,ic) &
            -(b4+b4p)*workvec(:,:,:,:,iq)-b4*workvec(:,:,:,:,ic)
    ENDDO
    ! Step 9: calculate the curl of the current density, stopr in spot
    DO iq=1,2  
       CALL rmuly(der1y,current(:,:,:,3,iq),spot(:,:,:,1,iq),0)
       CALL rmulz(der1z,current(:,:,:,2,iq),spot(:,:,:,1,iq),-1)
       CALL rmulz(der1z,current(:,:,:,1,iq),spot(:,:,:,2,iq),0)
       CALL rmulx(der1x,current(:,:,:,3,iq),spot(:,:,:,2,iq),-1)
       CALL rmulx(der1x,current(:,:,:,2,iq),spot(:,:,:,3,iq),0)
       CALL rmuly(der1y,current(:,:,:,1,iq),spot(:,:,:,3,iq),-1)
    ENDDO
    ! Step 10: combine isospin contributions
    DO icomp=1,3  
       DO iz=1,nz
          DO iy=1,ny
             DO ix=1,nx  
                rotspp=spot(ix,iy,iz,icomp,1)  
                rotspn=spot(ix,iy,iz,icomp,2)  
                spot(ix,iy,iz,icomp,1)=-(b4+b4p)*rotspp-b4*rotspn  
                spot(ix,iy,iz,icomp,2)=-(b4+b4p)*rotspn-b4*rotspp  
             ENDDO
          ENDDO
       ENDDO
    ENDDO
     ! Step 11: calculate divergence of aq in divaq 
     DO iq=1,2  
       CALL rmulx(der1x,aq(:,:,:,1,iq),divaq(:,:,:,iq),0)
       CALL rmuly(der1y,aq(:,:,:,2,iq),divaq(:,:,:,iq),1)
       CALL rmulz(der1z,aq(:,:,:,3,iq),divaq(:,:,:,iq),1)
     ENDDO
     ! Step 12: calculate the gradient of the effective mass in dbmass
     DO iq=1,3  
       CALL rmulx(der1x,bmass(:,:,:,iq),dbmass(:,:,:,1,iq),0)
       CALL rmuly(der1y,bmass(:,:,:,iq),dbmass(:,:,:,2,iq),0)
       CALL rmulz(der1z,bmass(:,:,:,iq),dbmass(:,:,:,3,iq),0)
     ENDDO
     ! Step 13: Laplacian of spin density components
     !DO iq=1,2
     !  DO icomp=1,3
     !    CALL rmulx(der2x,sdens(:,:,:,icomp,iq),Laplace_sq(:,:,:,icomp,iq),0)
     !    CALL rmuly(der2y,sdens(:,:,:,icomp,iq),Laplace_sq(:,:,:,icomp,iq),1)
     !    CALL rmulz(der2z,sdens(:,:,:,icomp,iq),Laplace_sq(:,:,:,icomp,iq),1)
     !  ENDDO
     !ENDDO
     ! Step 13.5: gradient of div(s_q) for C_0*(div s_0)^2 + C_1*(div s_1)^2
    ! DO iq=1,2
    !   CALL rmulx(der1x,sdens(:,:,:,1,iq),workden(:,:,:,iq),0)
    !   CALL rmuly(der1y,sdens(:,:,:,2,iq),workden(:,:,:,iq),1)
    !   CALL rmulz(der1z,sdens(:,:,:,3,iq),workden(:,:,:,iq),1)
    !   CALL rmulx(der1x,workden(:,:,:,iq),workvec(:,:,:,1,iq),0)
    !   CALL rmuly(der1y,workden(:,:,:,iq),workvec(:,:,:,2,iq),0)
    !   CALL rmulz(der1z,workden(:,:,:,iq),workvec(:,:,:,3,iq),0)
    ! ENDDO
     ! Step 14: add time-odd spin terms to spot 
     DO icomp=1,3
       DO iq=1,2
         ic=3-iq
         spot(:,:,:,icomp,iq)=spot(:,:,:,icomp,iq) &
            !+1.0D0/16.0D0*(3.0D0*f%t1*(1.0D0-f%x1)+f%t2*(1.0D0+f%x2))*Laplace_sq(:,:,:,icomp,iq) &
            !+1.0D0/16.0D0*(f%t2*f%x2-3.0D0*f%t1*f%x1)*Laplace_sq(:,:,:,icomp,ic) &
            +0.5D0*f%t0*(f%x0*(sdens(:,:,:,icomp,1)+sdens(:,:,:,icomp,2))-sdens(:,:,:,icomp,iq)) &
            +1.0D0/12.0D0*f%t3*(rho(:,:,:,1)+rho(:,:,:,2))**f%power* &
            (f%x3*(sdens(:,:,:,icomp,1)+sdens(:,:,:,icomp,2))-sdens(:,:,:,icomp,iq))&
            +1.0/8.0*(f%t1*f%x1+f%t2*f%x2)* (tdens(:,:,:,icomp,1)+tdens(:,:,:,icomp,2)) &
            +1.0/8.0*(f%t2-f%t1)* (tdens(:,:,:,icomp,iq)) !&
            !-2.0D0*((C_0+C_1)*workvec(:,:,:,icomp,iq)+(C_0-C_1)*workvec(:,:,:,icomp,ic))
       ENDDO
     ENDDO
    ! step 15: add time-odd spin terms to cq
    cq(:,:,:,:,:,:) = 0d0
    do iq = 1,2
      do mu = 1,3              ! 空间导数方向
        do nu = 1,3            ! 自旋分量
          cq(:,:,:,nu,mu,iq) = (1d0/8d0)*(f%t1*f%x1 + f%t2*f%x2) * ( sdens(:,:,:,nu,1) + sdens(:,:,:,nu,2) ) &
                        + (1d0/8d0)*(f%t2 - f%t1) * sdens(:,:,:,nu,iq)
        end do
      end do
    end do
    DEALLOCATE(workden,workvec,Laplace_sq)
  END SUBROUTINE skyrme
!---------------------------------------------------------------------------  
! DESCRIPTION: hpsi
!> @brief
!!This subroutine applies the single-particle Hamiltonian to a
!!single-particle wave function \c pinn to produce an output wave
!!function \c pout. The argument \c iq indicates the isospin for
!!the wave function and \c eshift is an energy shift which is zero in
!!the dynamic calculation but crucial to the static algorithm (see 
!!\c grstep in module \c Static).
!>
!> @details
!!For an understanding of this module the role of the following local
!!variables is crucial.  They are
!!  - \c is: this is used in the loops over spin to indicate the
!!    spin component: <tt> is=1 </tt> for spin up and <tt> is=2 </tt> for spin down.
!!  - \c ic: denotes the index for the opposite spin; it is
!!    calculated as <tt> ic=3-is </tt>. Note the similar handling of
!!    the two isospin projections using \c iq and \c icomp in
!!    subroutine \c skyrme.
!!  - \c sigis: this variable denotes the sign of the spin
!!    projection. It is calculated as <tt> sigis=3-2*is </tt> and thus is \c +1
!!    for spin up (<tt> is=1 </tt>) and \c - for spin down (<tt> is=2 </tt>).
!!
!!The general structure of the subroutine is as follows: first the part
!!of the Hamiltonian not involving derivatives is applied, followed by the
!!terms involving derivatives in order \f$ x \f$, \f$ y \f$, $z$.
!!Since the structure of the Hamiltonian involves only first or second
!!derivatives in one spatial direction in each term, the derivatives can
!!be calculated for one direction and then the working space can be
!!reused for the next one.
!!
!!The expressions for the different spatial derivatives are quite
!!analogous, so that only the $x$-direction will be discussed at length
!!below.
!!
!!The expressions is repeated here:
!!\f[
!!  \hat h=U_q(\vec r)-\nabla\cdot\left[B_q(\vec r)\nabla\right]
!!  +\I\vec W_q\cdot(\vec\sigma\times\nabla)
!!  +\vec S_q\cdot\vec\sigma
!!  -\frac{\I}{2} \left[(\nabla\cdot\vec A_q)+2\vec A_q\cdot\nabla\right].
!!\f]
!! 
!!  -# the non-derivative parts not involving spin. These
!!     arise from \f$ U_q \f$ and \f$ -\tfrac{\I}{2}\,\nabla\cdot\vec A_q \f$, which
!!     are combined into a complex expression. The energy shift \c eshift is also included.
!!  -# the spin current coupling is constructed by simply
!!     using the explicit definition of the Pauli matrices and multiplying
!!     the resulting matrix onto the spinor wave function.
!!  -# the first and second derivative in the \f$ x \f$-direction are evaluated 
!!     and stored in the arrays \c pswk and \c pswk2. The last term in the Hamiltonian 
!!     gives rise to the two contributions
!!     \f[ -\frac{\partial B_q}{\partial x}\frac{\partial}{\partial x}-B_q 
!!     \frac{\partial^2}{\partial x^2}, \f] 
!!     of which the second is evaluated
!!     straightforwardly, while the first one is combined with the spin-orbit
!!     contribution. The part of \f$ \I\vec W_q\cdot(\vec\sigma\times\nabla) \f$
!!     that contains an \f$ x \f$-derivative is
!!     \f[ (\I W_y\sigma_z-\I W_z\sigma_y)\frac{\partial}{\partial x}=
!!     \begin{pmatrix} \I W_y&-W_z\\ W_z&-\I W_y
!!     \end{pmatrix}\frac{\partial}{\partial x} \f]
!!     This is programmed employing the variable \c sigis to account for
!!     the different signs in the rows of the matrix.
!!  -# for the derivatives in the $y$-direction the
!!     procedure is similar; the spin-orbit part is now
!!     \f[ (\I W_z\sigma_x-\I W_x\sigma_z)\frac{\partial}{\partial y}=
!!     \begin{pmatrix} -\I W_x&\I W_z\\ \I W_z&\I W_x
!!     \end{pmatrix}\frac{\partial}{\partial y} \f]
!!  -# for the derivatives in the \f$ z \f$-direction the
!!     procedure is again similar; the spin-orbit part is now
!!     \f[ (\I W_x\sigma_y-\I W_y\sigma_x)\frac{\partial}{\partial z}=
!!     \begin{pmatrix} 0&W_x-\I W_y\\ -W_x-\I W_y & 0
!!     \end{pmatrix}\frac{\partial}{\partial z} \f]
!>
!> @param[in] iq
!> INTEGER, takes the isospin.
!> @param[in] eshift
!> REAL(db), takes the energy shift.
!> @param[in,out] pinn
!> COMPLEX(db), array, takes the wave function the Hamiltonian is supposed to be applied to.
!> @param[out] pout
!> COMPLEX(db), array, returns the Hamiltonian applied to the wave function.
!--------------------------------------------------------------------------- 
  SUBROUTINE hpsi(iq,eshift,pinn,pout)
    USE Trivial, ONLY: cmulx, cmuly, cmulz
    USE Levels, ONLY: cdervx,cdervy,cdervz
    INTEGER :: iq
    REAL(db) :: eshift
    COMPLEX(db),DIMENSION(:,:,:,:) :: pinn,pout
    INTENT(IN) :: iq,eshift
    INTENT(INOUT) :: pinn
    INTENT(OUT) :: pout
    INTEGER :: is,ic,i,j,k
    REAL(db) :: sigis,xx(nx),yy(ny)
    COMPLEX(db),ALLOCATABLE,DIMENSION(:,:,:,:) :: pswk,pswk2,psx,psy,Lzpsi,tmp, dxy, dyx
    ALLOCATE(pswk(nx,ny,nz,2),pswk2(nx,ny,nz,2))
    ALLOCATE(psx(nx,ny,nz,2),psy(nx,ny,nz,2))
    ALLOCATE(Lzpsi(nx,ny,nz,2))
    ALLOCATE(tmp(nx,ny,nz,2), dxy(nx,ny,nz,2), dyx(nx,ny,nz,2))
    ! Step 1: non-derivative parts not involving spin
    DO is=1,2
      pout(:,:,:,is)=CMPLX(upot(:,:,:,iq)-eshift, &
           -0D0,db)*pinn(:,:,:,is)
    ENDDO
    ! Step 2: the spin-current coupling
    pout(:,:,:,1)=pout(:,:,:,1)  &
         + CMPLX(spot(:,:,:,1,iq),-spot(:,:,:,2,iq),db) &
         *pinn(:,:,:,2)  + spot(:,:,:,3,iq)*pinn(:,:,:,1)
    pout(:,:,:,2)=pout(:,:,:,2) &
         + CMPLX(spot(:,:,:,1,iq),spot(:,:,:,2,iq),db) &
         *pinn(:,:,:,1) - spot(:,:,:,3,iq)*pinn(:,:,:,2)
    !add magnetic field
    if(iq==1) then
      pout(:,:,:,1) = pout(:,:,:,1) + CMPLX(0.60228d0, 0, db) * b_d17 * pinn(:,:,:,1)
      pout(:,:,:,2) = pout(:,:,:,2) - CMPLX(0.60228d0, 0, db) * b_d17 * pinn(:,:,:,2)
    else if(iq==2) then
      pout(:,:,:,1) = pout(:,:,:,1) - CMPLX(0.880416d0, 0, db) * b_d17 * pinn(:,:,:,1)
      pout(:,:,:,2) = pout(:,:,:,2) + CMPLX(0.880416d0, 0, db) * b_d17 * pinn(:,:,:,2)
    elseif(iq==3) then
      pout(:,:,:,1) = pout(:,:,:,1) + CMPLX(0.193245d0, 0.d0, db) * b_d17 * pinn(:,:,:,1)
      pout(:,:,:,2) = pout(:,:,:,2) - CMPLX(0.193245d0, 0.d0, db) * b_d17 * pinn(:,:,:,2)
    endif

    IF (iq==2) THEN

      xx = x - cmtot(1)
      yy = y - cmtot(2)
      ! ====== Lz term (Hermitian-symmetrized) ======
      ! First: psx = d/dx psi, psy = d/dy psi
      IF (TFFT) THEN
        CALL cdervx(pinn,psx)
        CALL cdervy(pinn,psy)
      ELSE
        CALL cmulx(der1x,pinn,psx,0)
        CALL cmuly(der1y,pinn,psy,0)
      ENDIF

      ! Compute d/dy ( x * psi )
      tmp = pinn
      DO is=1,2
        DO k=1,nz
          DO j=1,ny
            DO i=1,nx
              tmp(i,j,k,is) = xx(i) * pinn(i,j,k,is)
            END DO
          END DO
        END DO
      END DO
      IF (TFFT) THEN
        CALL cdervy(tmp,dxy)
      ELSE
        CALL cmuly(der1y,tmp,dxy,0)
      ENDIF

      ! Compute d/dx ( y * psi )
      tmp = pinn
      DO is=1,2
        DO k=1,nz
          DO j=1,ny
            DO i=1,nx
              tmp(i,j,k,is) = yy(j) * pinn(i,j,k,is)
            END DO
          END DO
        END DO
      END DO
      IF (TFFT) THEN
        CALL cdervx(tmp,dyx)
      ELSE
        CALL cmulx(der1x,tmp,dyx,0)
      ENDIF

      ! Lzpsi = -(i/2) * [ x*d_y psi + d_y(x psi) - y*d_x psi - d_x(y psi) ]
      DO is=1,2
        DO k=1,nz
          DO j=1,ny
            DO i=1,nx
              Lzpsi(i,j,k,is) = CMPLX(0D0,-0.5D0,db) * ( &
                   xx(i)*psy(i,j,k,is) + dxy(i,j,k,is) &
                 - yy(j)*psx(i,j,k,is) - dyx(i,j,k,is) )
            END DO
          END DO
        END DO
      END DO

      ! Add orbital Zeeman:  Δh = - (μ_N B) Lz
      DO is=1,2
        pout(:,:,:,is) = pout(:,:,:,is) - 0.31524512D0 * b_d17 * Lzpsi(:,:,:,is)
      END DO

      do is=1,2
        DO k=1,nz
          DO j=1,ny
            DO i=1,nx
              pout(i,j,k,is)=pout(i,j,k,is) + 0.31524512D0 * b_d17*0.31524512D0 &
              * b_d17*1.2049d-2*(xx(i)*xx(i)+yy(j)*yy(j))*pinn(i,j,k,is)
            ENDDO
          ENDDO
        ENDDO
      END DO
    ENDIF

    ! Step 3: derivative terms in x
    IF(TFFT) THEN
      CALL cdervx(pinn,pswk,d2psout=pswk2)
    ELSE
      CALL cmulx(der1x,pinn,pswk,0)
      CALL cmulx(der2x,pinn,pswk2,0)
    ENDIF
    ! ===== sigma·C term (Hermitian) :  -∂x [ (sigma·C_x) ∂x psi ] =====
    ! Here pswk = ∂x pinn is already available

    tmp(:,:,:,1) =  cq(:,:,:,3,1,iq) * pswk(:,:,:,1) &
                  + CMPLX( cq(:,:,:,1,1,iq), -cq(:,:,:,2,1,iq), db ) * pswk(:,:,:,2)

    tmp(:,:,:,2) =  CMPLX( cq(:,:,:,1,1,iq),  cq(:,:,:,2,1,iq), db ) * pswk(:,:,:,1) &
                  - cq(:,:,:,3,1,iq) * pswk(:,:,:,2)

    IF (TFFT) THEN
      CALL cdervx(tmp, dxy)
    ELSE
      CALL cmulx(der1x, tmp, dxy, 0)
    END IF

    pout(:,:,:,:) = pout(:,:,:,:) - dxy(:,:,:,:)

    DO is=1,2
      ic=3-is
      sigis=(3-2*is)*0.5D0
      pout(:,:,:,is)=pout(:,:,:,is) &
           -CMPLX(dbmass(:,:,:,1,iq),0.5D0*aq(:,:,:,1,iq) &
           -sigis*wlspot(:,:,:,2,iq),db)*pswk(:,:,:,is)  &
           -sigis*wlspot(:,:,:,3,iq)*pswk(:,:,:,ic) &
           -bmass(:,:,:,iq)*pswk2(:,:,:,is)
    ENDDO
    pswk2(:,:,:,1) = CMPLX(0D0,-0.5D0,db)*&
         (aq(:,:,:,1,iq)-wlspot(:,:,:,2,iq))*pinn(:,:,:,1)&
         -0.5D0*wlspot(:,:,:,3,iq)*pinn(:,:,:,2)
    pswk2(:,:,:,2) = CMPLX(0D0,-0.5D0,db)*&
         (aq(:,:,:,1,iq)+wlspot(:,:,:,2,iq))*pinn(:,:,:,2)&
         +0.5D0*wlspot(:,:,:,3,iq)*pinn(:,:,:,1)
    IF(TFFT) THEN
      CALL cdervx(pswk2,pswk)
    ELSE
      CALL cmulx(der1x,pswk2,pswk,0)
    ENDIF
    pout(:,:,:,:)=pout(:,:,:,:) + pswk(:,:,:,:)

    ! Step 4: derivative terms in y
    IF(TFFT) THEN
      CALL cdervy(pinn,pswk,d2psout=pswk2)
    ELSE
      CALL cmuly(der1y,pinn,pswk,0)
      CALL cmuly(der2y,pinn,pswk2,0)
    ENDIF

    ! ===== sigma·C term (Hermitian) :  -∂y [ (sigma·C_y) ∂y psi ] =====

    tmp(:,:,:,1) =  cq(:,:,:,3,2,iq) * pswk(:,:,:,1) &
                  + CMPLX( cq(:,:,:,1,2,iq), -cq(:,:,:,2,2,iq), db ) * pswk(:,:,:,2)

    tmp(:,:,:,2) =  CMPLX( cq(:,:,:,1,2,iq),  cq(:,:,:,2,2,iq), db ) * pswk(:,:,:,1) &
                  - cq(:,:,:,3,2,iq) * pswk(:,:,:,2)

    IF (TFFT) THEN
      CALL cdervy(tmp, dxy)
    ELSE
      CALL cmuly(der1y, tmp, dxy, 0)
    END IF

    pout(:,:,:,:) = pout(:,:,:,:) - dxy(:,:,:,:)

    DO is=1,2
      ic=3-is
      sigis=(3-2*is)*0.5D0
      pout(:,:,:,is)=pout(:,:,:,is) &
           -CMPLX(dbmass(:,:,:,2,iq),0.5D0*aq(:,:,:,2,iq) &
           +sigis*wlspot(:,:,:,1,iq),db)*pswk(:,:,:,is) &
           +CMPLX(0.D0,0.5D0*wlspot(:,:,:,3,iq),db)*pswk(:,:,:,ic) &
           -bmass(:,:,:,iq)*pswk2(:,:,:,is)
    ENDDO
    pswk2(:,:,:,1) = CMPLX(0D0,-0.5D0,db)*&
         (aq(:,:,:,2,iq)+wlspot(:,:,:,1,iq))*pinn(:,:,:,1)&
         +CMPLX(0D0,0.5D0,db)*wlspot(:,:,:,3,iq)*pinn(:,:,:,2)
    pswk2(:,:,:,2) = CMPLX(0D0,-0.5D0,db)*&
         (aq(:,:,:,2,iq)-wlspot(:,:,:,1,iq))*pinn(:,:,:,2)&
         +CMPLX(0D0,0.5D0*wlspot(:,:,:,3,iq),db)*pinn(:,:,:,1)
    IF(TFFT) THEN
      CALL cdervy(pswk2,pswk)
    ELSE
      CALL cmuly(der1y,pswk2,pswk,0)
    ENDIF
    pout(:,:,:,:)=pout(:,:,:,:) + pswk(:,:,:,:)
    ! Step 5: derivative terms in z
    IF(TFFT) THEN
      CALL cdervz(pinn,pswk,d2psout=pswk2)
    ELSE
      CALL cmulz(der1z,pinn,pswk,0)
      CALL cmulz(der2z,pinn,pswk2,0)
    ENDIF

    ! ===== sigma·C term (Hermitian) :  -∂z [ (sigma·C_z) ∂z psi ] =====

    tmp(:,:,:,1) =  cq(:,:,:,3,3,iq) * pswk(:,:,:,1) &
                  + CMPLX( cq(:,:,:,1,3,iq), -cq(:,:,:,2,3,iq), db ) * pswk(:,:,:,2)

    tmp(:,:,:,2) =  CMPLX( cq(:,:,:,1,3,iq),  cq(:,:,:,2,3,iq), db ) * pswk(:,:,:,1) &
                  - cq(:,:,:,3,3,iq) * pswk(:,:,:,2)

    IF (TFFT) THEN
      CALL cdervz(tmp, dxy)
    ELSE
      CALL cmulz(der1z, tmp, dxy, 0)
    END IF

    pout(:,:,:,:) = pout(:,:,:,:) - dxy(:,:,:,:)

    DO is=1,2
      ic=3-is
      sigis=(3-2*is)*0.5D0
      pout(:,:,:,is)=pout(:,:,:,is) &
           -CMPLX(dbmass(:,:,:,3,iq),0.5D0*aq(:,:,:,3,iq),db)*pswk(:,:,:,is) &
           +CMPLX(sigis*wlspot(:,:,:,1,iq),-0.5D0*wlspot(:,:,:,2,iq),db)* &
           pswk(:,:,:,ic)-bmass(:,:,:,iq)*pswk2(:,:,:,is)
    ENDDO
    pswk2(:,:,:,1) = CMPLX(0D0,-0.5D0,db)*aq(:,:,:,3,iq)*pinn(:,:,:,1)&
         +CMPLX(0.5D0*wlspot(:,:,:,1,iq),-0.5D0*wlspot(:,:,:,2,iq),db)*pinn(:,:,:,2)
    pswk2(:,:,:,2) = CMPLX(0D0,-0.5D0,db)*aq(:,:,:,3,iq)*pinn(:,:,:,2)&
         +CMPLX(-0.5D0*wlspot(:,:,:,1,iq),-0.5D0*wlspot(:,:,:,2,iq),db)*pinn(:,:,:,1)
    IF(TFFT) THEN
      CALL cdervz(pswk2,pswk)
    ELSE
      CALL cmulz(der1z,pswk2,pswk,0)
    ENDIF
    pout(:,:,:,:)=pout(:,:,:,:) + pswk(:,:,:,:)

    DEALLOCATE(pswk,pswk2,Lzpsi,psx,psy,tmp,dxy,dyx)
  END SUBROUTINE hpsi
  !***********************************************************************
END Module Meanfield
