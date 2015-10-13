SUBROUTINE evalFluxes(flx,quadVals,elemAvg,coeffs,uin,qWeights,nelem,nNodes,dt,dx)
  ! ------------------------------------------
  !
  ! ------------------------------------------
  USE testParameters
  IMPLICIT NONE
  ! Inputs
  INTEGER, INTENT(IN) :: nelem,nNodes
  DOUBLE PRECISION, INTENT(IN) :: dt,dx
  DOUBLE PRECISION, DIMENSION(0:nNodes,0:nelem), INTENT(IN) :: uin
  DOUBLE PRECISION, DIMENSION(1:nelem), INTENT(IN) :: elemAvg
  DOUBLE PRECISION, DIMENSION(0:nNodes), INTENT(IN) :: qWeights

  ! Outputs
  DOUBLE PRECISION, DIMENSION(0:nelem), INTENT(OUT) :: flx
  DOUBLE PRECISION, DIMENSION(0:nNodes,1:nelem), INTENT(OUT) :: quadVals
  DOUBLE PRECISION, DIMENSION(0:nNodes,1:nelem), INTENT(INOUT) :: coeffs

  ! Local variables
  INTEGER :: j,k
  DOUBLE PRECISION, DIMENSION(0:1,0:nelem+1) :: edgeVals
  DOUBLE PRECISION :: tmp

  DOUBLE PRECISION, DIMENSION(0:nNodes,1:nelem) :: coeffs_tmp

  INTERFACE
    SUBROUTINE numFlux(flx,edgeVals,uin,nelem,nNodes)
    	IMPLICIT NONE
    	! -- Inputs
    	INTEGER, INTENT(IN) :: nelem,nNodes
    	DOUBLE PRECISION, DIMENSION(0:1,0:nelem+1), INTENT(IN) :: edgeVals
    	DOUBLE PRECISION, DIMENSION(0:nNodes,0:nelem), INTENT(IN) :: uin
      ! -- Outputs
      DOUBLE PRECISION,DIMENSION(0:nelem), INTENT(OUT) :: flx
    END SUBROUTINE numFlux

    SUBROUTINE evalExpansion(quadVals,edgeVals,qIn,nelem,nNodes)
    ! Evaluate ansatz solution at quadrature nodes and edges
    ! Note:  Assumes basis is Lagrange interpolating polynomials at quad nodes -> phi_j (xi_k) = qIn(k,j)
        IMPLICIT NONE
        ! Inputs
        INTEGER, INTENT(IN) :: nelem,nNodes
        DOUBLE PRECISION, DIMENSION(0:nNodes,1:nelem), INTENT(INOUT) :: qIn

        ! Outputs
        DOUBLE PRECISION, DIMENSION(0:nNodes,1:nelem), INTENT(OUT) :: quadVals
        DOUBLE PRECISION, DIMENSION(0:1,0:nelem+1), INTENT(OUT) :: edgeVals
      END SUBROUTINE evalExpansion

      SUBROUTINE fctPolyMod(coeffs,flx,elemAvgs,uin,dx,dt,nelem,nNodes)
        ! Makes FCT-based correction to fluxes to ensure that element means remain non-negative
        ! Then modifies local polynomials so that FCT fluxes will be consistent
        ! NOTE:: THIS IS NOT CONSERVATIVE! Requires calling conservPolyMod() after to keep conservation
        IMPLICIT NONE
        ! Inputs
        INTEGER, INTENT(IN) :: nelem,nNodes
        DOUBLE PRECISION, INTENT(IN) :: dx,dt
        DOUBLE PRECISION, DIMENSION(1:nelem), INTENT(IN) :: elemAvgs
        DOUBLE PRECISION, DIMENSION(0:nNodes,0:nelem), INTENT(IN) :: uin
        ! Outputs
        DOUBLE PRECISION, DIMENSION(0:nelem), INTENT(INOUT) :: flx
        DOUBLE PRECISION, DIMENSION(0:nNodes,1:nelem), INTENT(INOUT) :: coeffs
      END SUBROUTINE fctPolyMod

      SUBROUTINE conservPolyMod(coeffs,prevEdgeVals,qWeights,nelem,nNodes)
        ! Makes additive conservative modification to internal (1:nNodes-1) coefficients
        ! corresponding to internal nodal quadrature locations within each element
        ! NOTE:: This only changes the values of INTERNAL points
        ! Inputs -
        !  . prevEdgeVals : array of left and right edge values of polynomial BEFORE modification
        !  . qWeights : array of GLL quadrature weights
        ! Outputs -
        !  . coeffs : array of nodal coefficients after desired edge node modification
        IMPLICIT NONE
        ! Inputs
        INTEGER, INTENT(IN) :: nelem,nNodes
        DOUBLE PRECISION, DIMENSION(0:nNodes), INTENT(IN) :: qWeights
        DOUBLE PRECISION, DIMENSION(0:1,0:nelem+1), INTENT(IN) :: prevEdgeVals
        ! Outputs
        DOUBLE PRECISION, DIMENSION(0:nNodes,1:nelem), INTENT(INOUT) :: coeffs
      END SUBROUTINE conservPolyMod
  END INTERFACE

  ! Get incoming edge and quadrature values
  CALL evalExpansion(quadVals,edgeVals,coeffs,nelem,nNodes)
  CALL numFlux(flx,edgeVals,uin,nelem,nNodes)

  ! For non-FCT methods, return unmodified fluxes
  IF(.not. doFCT) THEN
    RETURN
  ENDIF

  SELECT CASE(limitingMeth)
    CASE(4)
      ! modified FCT
      CALL fctPolyMod(coeffs,flx,elemAvg,uin,dx,dt,nelem,nNodes)
  END SELECT ! limiting meth
  CALL conservPolyMod(coeffs,edgeVals,qWeights,nelem,nNodes)

  ! Update quadurature and edge values using modified polynomial coefficients
  CALL evalExpansion(quadVals,edgeVals,coeffs,nelem,nNodes)
!  CALL numFlux(flx,edgeVals,uin,nelem,nNodes)

!  DO j=1,nelem
!    mean = elemAvg(j)
!    tmp = 0.5D0*SUM(qWeights(:)*coeffs(:,j))
!    write(*,*) 'j=',j
!    write(*,*) 'DM after polymod =',tmp-mean
!    IF(abs(tmp-mean) .gt. eps) THEN
!      write(*,*) 'j=',j
!      write(*,*) 'before=',mean
!      write(*,*) 'after=',tmp
!      write(*,*) '***'
!      STOP
!    ENDIF
!  ENDDO !j

END SUBROUTINE evalFluxes

SUBROUTINE fctPolyMod(coeffs,flx,elemAvg,uin,dx,dt,nelem,nNodes)
  ! Makes FCT-based correction to fluxes to ensure that element means remain non-negative
  ! Then modifies local polynomials so that FCT fluxes will be consistent
  ! NOTE:: THIS IS NOT CONSERVATIVE! Requires calling conservPolyMod() after to keep conservation
  IMPLICIT NONE
  ! Inputs
  INTEGER, INTENT(IN) :: nelem,nNodes
  DOUBLE PRECISION, INTENT(IN) :: dx,dt
  DOUBLE PRECISION, DIMENSION(1:nelem), INTENT(IN) :: elemAvg
  DOUBLE PRECISION, DIMENSION(0:nNodes,0:nelem), INTENT(IN) :: uin
  ! Outputs
  DOUBLE PRECISION, DIMENSION(0:nelem), INTENT(INOUT) :: flx
  DOUBLE PRECISION, DIMENSION(0:nNodes,1:nelem), INTENT(INOUT) :: coeffs

  ! Local variables
  INTEGER :: j
  DOUBLE PRECISION :: Qj,Pj,eps,mean
  DOUBLE PRECISION, DIMENSION(0:nelem+1) :: fluxRatio

  eps = 1D-10

  DO j=1,nelem
    mean = elemAvg(j)
    Qj = mean*dx/dt
    Pj = MAX(flx(j),0D0)-MIN(flx(j-1),0D0) + eps
    fluxRatio(j) = MIN(1D0,Qj/Pj)
  ENDDO !j
  fluxRatio(0) = fluxRatio(nelem)
  fluxRatio(nelem+1) = fluxRatio(1)

  DO j=0,nelem-1
    ! Pick appropriate flux modification ratio and adjust coefficient
    IF(flx(j) .ge. 0D0) THEN
      flx(j) = fluxRatio(j)*flx(j)
    ELSE
      flx(j) = fluxRatio(j+1)*flx(j)
    ENDIF
  ENDDO !j
  ! Maintain flux periodicity
  flx(nelem) = flx(0)

  ! TODO: This section implicitly refers to coeffs using periodic extension
  ! This reference should be made explicit (by adding ghost cells to coeffs earlier)

  ! Modify coefficient on upstream side of interfaces
  IF(uIn(nNodes,0) .gt. 0D0) THEN
    coeffs(nNodes,nelem) = flx(0)/uIn(nNodes,0)
  ELSE IF(uIn(nNodes,0) .lt. 0D0) THEN
    coeffs(0,1) = flx(0)/uIn(nNodes,0)
  ELSE
    ! Make no modifications (u==0)
  ENDIF

  DO j=1,nelem-1
    ! Modify coefficient on upstream side of interface
    IF(uIn(nNodes,j) .gt. 0D0) THEN
      coeffs(nNodes,j) = flx(j)/uIn(nNodes,j)
    ELSE IF(uIn(nNodes,j) .lt. 0D0) THEN
      coeffs(0,j+1) = flx(j)/uIn(nNodes,j)
    ELSE
      ! Make no modifications (u==0)
    ENDIF
  ENDDO !j

END SUBROUTINE fctPolyMod

SUBROUTINE conservPolyMod(coeffs,prevEdgeVals,qWeights,nelem,nNodes)
  ! Makes additive conservative modification to internal (1:nNodes-1) coefficients
  ! corresponding to internal nodal quadrature locations within each element
  ! NOTE:: This only changes the values of INTERNAL points
  ! Inputs -
  !  . prevEdgeVals : array of left and right edge values of polynomial BEFORE modification
  !  . qWeights : array of GLL quadrature weights
  ! Outputs -
  !  . coeffs : array of nodal coefficients after desired edge node modification
  IMPLICIT NONE
  ! Inputs
  INTEGER, INTENT(IN) :: nelem,nNodes
  DOUBLE PRECISION, DIMENSION(0:nNodes), INTENT(IN) :: qWeights
  DOUBLE PRECISION, DIMENSION(0:1,0:nelem+1), INTENT(IN) :: prevEdgeVals
  ! Outputs
  DOUBLE PRECISION, DIMENSION(0:nNodes,1:nelem), INTENT(INOUT) :: coeffs
  ! Local variables
  INTEGER :: j,k
  DOUBLE PRECISION :: delLeft,delRight,massCor

  ! Correct mass for modified polynomials to maintain conservation
  DO j=1,nelem
    ! Change at left edge
    delLeft = prevEdgeVals(0,j)-coeffs(0,j)
    ! Change at right edge
    delRight = prevEdgeVals(1,j)-coeffs(nNodes,j)
    ! Net mass change due to modification
    massCor = (qWeights(0)*delLeft+qWeights(nNodes)*delRight)/(nNodes-1)

    ! Make modification so that each (internal) node recieves an equal
    ! portion of the mass deficit
    DO k=1,nNodes-1
      coeffs(k,j) = coeffs(k,j)+massCor/qWeights(k)
    ENDDO !k
  ENDDO !j
END SUBROUTINE conservPolyMod
