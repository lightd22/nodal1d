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
  DOUBLE PRECISION :: mean,eps,delLeft,delRight,massCor,tmp
  DOUBLE PRECISION :: Qj,Pj,C
  DOUBLE PRECISION, DIMENSION(0:nelem) :: fluxLow,fluxCorrector,corrector
  DOUBLE PRECISION, DIMENSION(0:nelem+1) :: fluxRatio

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
  END INTERFACE

  eps = 1D-10

  ! Get incoming edge and quadrature values
  CALL evalExpansion(quadVals,edgeVals,coeffs,nelem,nNodes)
  CALL numFlux(flx,edgeVals,uin,nelem,nNodes)

  ! For non-FCT methods, return unmodified fluxes
  IF(.not. doFCT) THEN
    RETURN
  ENDIF

  DO j=1,nelem
    mean = elemAvg(j)
    Qj = mean*dx/dt
    Pj = MAX(flx(j),0D0)-MIN(flx(j-1),0D0) + eps
    fluxRatio(j) = MIN(1D0,Qj/Pj)
  ENDDO !j
  fluxRatio(0) = fluxRatio(nelem)
  fluxRatio(nelem+1) = fluxRatio(1)

  ! TODO: This section implicitly refers to coeffs using periodic extension
  ! This reference should be made explicit (by adding ghost cells to coeffs earlier)
  IF(flx(0) .ge. 0D0) THEN
    flx(0) = fluxRatio(0)*flx(0)
  ELSE
    flx(0) = fluxRatio(1)*flx(0)
  ENDIF
  ! Maintain flux periodicity
  flx(nelem) = flx(0)

  ! Modify coefficient on upstream side of interface
  IF(uIn(nNodes,0) .gt. 0D0) THEN
    coeffs(nNodes,nelem) = flx(0)/uIn(nNodes,0)
  ELSE IF(uIn(nNodes,0) .lt. 0D0) THEN
    coeffs(0,1) = flx(0)/uIn(nNodes,0)
  ELSE
    ! Make no modifications (u==0)
  ENDIF

  DO j=1,nelem-1
    ! Pick appropriate flux modification ratio and adjust coefficient
    IF(flx(j) .ge. 0D0) THEN
      flx(j) = fluxRatio(j)*flx(j)
    ELSE
      flx(j) = fluxRatio(j+1)*flx(j)
    ENDIF

    ! Modify coefficient on upstream side of interface
    IF(uIn(nNodes,j) .gt. 0D0) THEN
      coeffs(nNodes,j) = flx(j)/uIn(nNodes,j)
    ELSE IF(uIn(nNodes,j) .lt. 0D0) THEN
      coeffs(0,j+1) = flx(j)/uIn(nNodes,j)
    ELSE
      ! Make no modifications (u==0)
    ENDIF
  ENDDO !j

  ! Correct mass for modified polynomials to maintain conservation
  DO j=1,nelem
    ! Change at left edge
    delLeft = edgeVals(0,j)-coeffs(0,j)
    ! Change at right edge
    delRight = edgeVals(1,j)-coeffs(nNodes,j)
    ! Net mass change due to modification
    massCor = (qWeights(0)*delLeft+qWeights(nNodes)*delRight)/(nNodes-1)

    ! Make modification so that each (internal) node recieves an equal
    ! portion of the mass deficit
    DO k=1,nNodes-1
      coeffs(k,j) = coeffs(k,j)+massCor/qWeights(k)
    ENDDO !k
  ENDDO !j

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
