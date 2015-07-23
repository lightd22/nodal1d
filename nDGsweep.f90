!--------------------------------------------------------------------------------
! Nodal DG coefficient update for 1D advection
! By: Devin Light 4/21/14
! -------------------------------------------------------------------------------

SUBROUTINE nDGsweep(q,nelem,dxel,nNodes,qNodes,qWeights,u,lagDeriv,doposlimit,dt,&
                    nZSNodes,quadZSNodes,quadZSWeights,lagValsZS)
    IMPLICIT NONE
    INTEGER, PARAMETER :: DOUBLE = KIND(1D0) ! Specification of DOUBLE

    ! External Functions
	  REAL(KIND=DOUBLE), EXTERNAL :: dadt ! RHS function for evolution ODE for kth expansion coefficent

    ! Inputs
    INTEGER, INTENT(IN) :: nelem,nNodes
    INTEGER, INTENT(IN) :: nZSNodes
    LOGICAL, INTENT(IN) :: doposlimit
    REAL(KIND=DOUBLE), INTENT(IN) :: dxel,dt
    REAL(KIND=DOUBLE), DIMENSION(0:nNodes), INTENT(IN):: qNodes,qWeights
    REAL(KIND=DOUBLE), DIMENSION(0:nZSNodes), INTENT(IN) :: quadZSNodes,quadZSWeights
    REAL(KIND=DOUBLE), DIMENSION(0:nNodes,0:nZSNodes), INTENT(IN) :: lagValsZS
    REAL(KIND=DOUBLE), DIMENSION(0:nNodes,0:nNodes), INTENT(IN) :: lagDeriv
    REAL(KIND=DOUBLE), DIMENSION(0:nNodes,1:nelem), INTENT(IN) :: u

    ! Outputs
    REAL(KIND=DOUBLE), DIMENSION(0:nNodes,1:nelem), INTENT(INOUT) :: q

    ! Local Variables
    REAL(KIND=DOUBLE), DIMENSION(0:nNodes,1:nelem) :: qFwd,qStar,quadVals
    REAL(KIND=DOUBLE), DIMENSION(0:nNodes,0:nelem) :: utmp
    REAL(KIND=DOUBLE), DIMENSION(0:1,0:nelem+1) :: edgeVals
    REAL(KIND=DOUBLE), DIMENSION(0:nelem) :: flx
    REAL(KIND=DOUBLE), DIMENSION(0:nNodes) :: tmp

    REAL(KIND=DOUBLE), DIMENSION(1:nelem) :: avgVals
    REAL(KIND=DOUBLE), DIMENSION(0:nZSNodes) :: qVals

    INTEGER :: k,j,stage,l,jstar

    INTERFACE
      SUBROUTINE limitNodePositivity(qIn,avgVals,qWeights,nelem,nNodes,nQuad)
        ! Modifies approximating polynomial so that nodal values are non-negative for output
        USE testParameters
        IMPLICIT NONE
        ! Inputs
        INTEGER, INTENT(IN) :: nelem,nNodes,nQuad
        REAL(KIND=8), DIMENSION(0:nNodes,1:nelem), INTENT(INOUT) :: qIn
        REAL(KIND=8), DIMENSION(1:nelem), INTENT(IN) :: avgVals
        REAL(KIND=8), DIMENSION(0:nQuad), INTENT(IN) :: qWeights
      END SUBROUTINE limitNodePositivity

      SUBROUTINE limitMeanPositivity(qIn,avgVals,qWeights,nelem,nNodes,nQuad)
        ! Modifies approximating polynomial so that element mean value remains non-negative
        ! according to Zhang and Shu (2010) Thm 2.2
        USE testParameters
        IMPLICIT NONE
        ! Inputs
        INTEGER, INTENT(IN) :: nelem,nNodes,nQuad
        REAL(KIND=8), DIMENSION(0:nNodes,1:nelem), INTENT(INOUT) :: qIn
        REAL(KIND=8), DIMENSION(1:nelem), INTENT(IN) :: avgVals
        REAL(KIND=8), DIMENSION(0:nQuad), INTENT(IN) :: qWeights
      END SUBROUTINE limitMeanPositivity
    END INTERFACE

    jstar = -1

    utmp(:,1:nelem) = u
    utmp(:,0) = u(:,nelem)

    qStar = q
    ! Do SSPRK3 Update
    DO stage=1,3
      IF(doposlimit) THEN ! Rescale for element-mean positivity
        DO j=1,nelem
          avgVals(j) = 0.5D0*SUM( qStar(:,j)*qWeights )
          IF(avgVals(j) .lt. 0D0) THEN
            write(*,'(A,I2,A,I2,A,E10.4)') ' Average of element ', j,' is negative in stage ',stage, &
                                           ' avgVal = ',avgVals(j)
            STOP
          ENDIF
        ENDDO
        CALL limitMeanPositivity(qStar,avgVals,quadZSWeights,nelem,nNodes,nZSNodes)

!        IF(nZSnodes .eq. nNodes) THEN
!          CALL positivityLimit(qStar,qNodes,qWeights,nelem,nNodes,nNodes,lagValsZS)
!        ELSE
!          CALL positivityLimit(qStar,quadZSNodes,quadZSWeights,nelem,nNodes,nZSNodes,lagValsZS)
!        ENDIF
      ENDIF
        CALL evalExpansion(quadVals,edgeVals,qStar,nelem,nNodes)
        CALL numFlux(flx,edgeVals,utmp,nelem,nNodes)

        ! Forward Step
        DO j=1,nelem
            DO k=0,nNodes
                tmp = lagDeriv(k,:)
                qFwd(k,j) = qStar(k,j) + (dt/dxel)*dadt(quadVals(:,j),flx,utmp(:,j),qWeights,tmp,k,j,nelem,nNodes)
            ENDDO !k
        ENDDO ! j

  		SELECT CASE(stage)
  		CASE(1)
  			qStar = qFwd
  		CASE(2)
  			qStar = 0.75D0*q + 0.25D0*qFwd
  		CASE(3)
  			qStar = q/3d0 + 2D0*qFwd/3D0
  		END SELECT
    ENDDO !stage

    ! After RK-update, rescale polynomial to remove all non-negatives
    IF(doposlimit) THEN
      DO j=1,nelem
        avgVals(j) = 0.5D0*SUM( qStar(:,j)*qWeights )
      ENDDO!j
      CALL limitNodePositivity(qStar,avgVals,qWeights,nelem,nNodes,nNodes)
!      IF(nZSnodes .eq. nNodes) THEN
!        CALL positivityLimit(q,qNodes,qWeights,nelem,nNodes,nNodes,lagValsZS)
!      ELSE
!        CALL positivityLimit(q,quadZSNodes,quadZSWeights,nelem,nNodes,nZSNodes,lagValsZS)
!      ENDIF
    ENDIF

    q = qStar

END SUBROUTINE nDGsweep

SUBROUTINE evalExpansion(quadVals,edgeVals,qIn,nelem,nNodes)
! Evaluate ansatz solution at quadrature nodes and edges
! Note:  Assumes basis is Lagrange interpolating polynomials at quad nodes -> phi_j (xi_k) = qIn(k,j)
    IMPLICIT NONE
    ! Inputs
    INTEGER, INTENT(IN) :: nelem,nNodes
    REAL(KIND=8), DIMENSION(0:nNodes,1:nelem), INTENT(INOUT) :: qIn

    ! Outputs
    REAL(KIND=8), DIMENSION(0:nNodes,1:nelem), INTENT(OUT) :: quadVals
    REAL(KIND=8), DIMENSION(0:1,0:nelem+1), INTENT(OUT) :: edgeVals

    ! Local Variables

    ! Ansatz value at quad locations for the jth element is just coefficient value
    quadVals = qIn

    edgeVals(0,1:nelem) = qIn(0,:) ! left edge value is just left-most coefficient
    edgeVals(1,1:nelem) = qIn(nNodes,:) ! right edge value is just right-most coefficient

	! Extend edgeVals periodically
	edgeVals(:,0) = edgeVals(:,nelem)
	edgeVals(:,nelem+1) = edgeVals(:,1)

END SUBROUTINE evalExpansion

SUBROUTINE numFlux(flx,edgeVals,uin,nelem,nNodes)
	IMPLICIT NONE
	INTEGER, PARAMETER :: DOUBLE = KIND(1D0)
	! -- Inputs
	INTEGER, INTENT(IN) :: nelem,nNodes
	REAL(KIND=DOUBLE), DIMENSION(0:1,0:nelem+1), INTENT(IN) :: edgeVals
	REAL(KIND=DOUBLE), DIMENSION(0:nNodes,0:nelem), INTENT(IN) :: uin

	! -- Outputs
	REAL(KIND=DOUBLE),DIMENSION(0:nelem), INTENT(OUT) :: flx

	! -- Local variables
	INTEGER :: j

	DO j=0,nelem
		flx(j) = 0.5D0*edgeVals(1,j)*(uin(nNodes,j)+DABS(uin(nNodes,j)))+0.5D0*edgeVals(0,j+1)*(uin(nNodes,j)-DABS(uin(nNodes,j)))
	ENDDO
END SUBROUTINE numFlux

REAL(KIND=8) FUNCTION dadt(quadVals,flx,u,qWeights,lagDeriv,k,j,nelem,nNodes)
    IMPLICIT NONE
    ! Inputs
    INTEGER, INTENT(IN) :: nelem,nNodes,k,j
    REAL(KIND=8), DIMENSION(0:nNodes) :: quadVals,u,qWeights,lagDeriv
    REAL(KIND=8), DIMENSION(0:nelem) :: flx

    dadt = SUM( u(:)*quadVals(:)*lagDeriv(:)*qWeights(:) )

    IF( k .eq. 0) THEN
        dadt = dadt + flx(j-1)
    ENDIF
    IF( k .eq. nNodes) THEN
        dadt = dadt - flx(j)
    ENDIF
    dadt = 2D0*dadt/qWeights(k)

END FUNCTION dadt
