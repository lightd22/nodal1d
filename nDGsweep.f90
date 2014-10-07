!--------------------------------------------------------------------------------
! Nodal DG coefficient update for 1D advection
! By: Devin Light 4/21/14
! -------------------------------------------------------------------------------

SUBROUTINE nDGsweep(q,nelem,dxel,nNodes,qNodes,qWeights,u,lagDeriv,dozhangshu,dt,&
                    nZSNodes,quadZSNodes,quadZSWeights,lagValsZS)
    IMPLICIT NONE
    INTEGER, PARAMETER :: DOUBLE = KIND(1D0) ! Specification of DOUBLE

    ! External Functions
	REAL(KIND=DOUBLE), EXTERNAL :: dadt ! RHS function for evolution ODE for kth expansion coefficent

    ! Inputs
    INTEGER, INTENT(IN) :: nelem,nNodes
    INTEGER, INTENT(IN) :: nZSNodes
    LOGICAL, INTENT(IN) :: dozhangshu
    REAL(KIND=DOUBLE), INTENT(IN) :: dxel,dt
    REAL(KIND=DOUBLE), DIMENSION(0:nNodes), INTENT(IN):: qNodes,qWeights
    REAL(KIND=DOUBLE), DIMENSION(0:nZSNodes), INTENT(IN) :: quadZSNodes,quadZSWeights
    REAL(KIND=DOUBLE), DIMENSION(0:nNodes,1:nZSNodes-1), INTENT(IN) :: lagValsZS
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
    INTEGER :: k,j,stage

    utmp(:,1:nelem) = u
    utmp(:,0) = u(:,nelem)

    qStar = q
    ! Do SSPRK3 Update
    DO stage=1,3
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
            q = qStar
			qStar = qFwd
		CASE(2)
			qStar = 0.75D0*q + 0.25D0*qFwd
		CASE(3)
			qStar = q/3d0 + 2D0*qFwd/3D0
		END SELECT
        
        IF(dozhangshu) THEN ! Rescale for positivity for next step
            CALL polyMod(qStar,qNodes,qWeights,nelem,nNodes)
!            CALL polyMod(qStar,quadZSNodes,quadZSWeights,nelem,nZSNodes)
        ENDIF
    ENDDO !stage
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

SUBROUTINE polyMod(qIn,qNodes,qWeights,nelem,nNodes)
! Rescales DG polynomial around element averages for positivity, based on Zhang and Shu (2010)
    IMPLICIT NONE
    ! Inputs
    INTEGER, INTENT(IN) :: nelem,nNodes
    REAL(KIND=8), DIMENSION(0:nNodes,1:nelem), INTENT(INOUT) :: qIn
    REAL(KIND=8), DIMENSION(0:nNodes), INTENT(IN) :: qWeights,qNodes
    ! Local Variables
    INTEGER :: j
    REAL(KIND=8) :: avgVal,theta,valMin,valMax

    DO j=1,nelem
        avgVal = 0.5D0*SUM( qWeights(:)*qIn(:,j) )
        valMin = MINVAL(qIn(:,j))
!       valMax = MAXVAL(qIn(:,j))

        theta = MIN( abs(avgVal/(valMin-avgVal)),1D0 )
!       theta = MIN( abs(avgVal/(valMin-avgVal)),abs((1D0-avgVal)/(valMax-avgVal)),1D0 )

        qIn(:,j) = theta*(qIn(:,j)-avgVal) + avgVal
    ENDDO !j
    
END SUBROUTINE polyMod


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
