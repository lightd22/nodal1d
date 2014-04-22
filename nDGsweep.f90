!--------------------------------------------------------------------------------
! Nodal DG coefficient update for 1D advection
! By: Devin Light 4/21/14
! -------------------------------------------------------------------------------

SUBROUTINE nDGsweep(q,nelem,dxel,nNodes,qNodes,qWeights,u,lagDeriv,dozhangshu,dt)
    IMPLICIT NONE
    INTEGER, PARAMETER :: DOUBLE = KIND(1D0) ! Specification of DOUBLE

    ! External Functions
	REAL(KIND=DOUBLE), EXTERNAL :: dadt ! RHS function for evolution ODE for kth expansion coefficent

    ! Inputs
    INTEGER, INTENT(IN) :: nelem,nNodes
    LOGICAL, INTENT(IN) :: dozhangshu
    REAL(KIND=DOUBLE), INTENT(IN) :: dxel,dt
    REAL(KIND=DOUBLE), DIMENSION(0:nNodes), INTENT(IN):: qNodes,qWeights
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
        CALL evalExpansion(quadVals,edgeVals,qStar,nelem,nNodes,dozhangshu)
        CALL numFlux(flx,edgeVals,utmp,nelem,nNodes)

        ! Forward Step
        DO j=1,nelem
            DO k=0,nNodes
                tmp = lagDeriv(k,:)
                qFwd(k,j) = qStar(k,j) + (dt/dxel)*dadt(quadVals(:,j),flx,u(:,j),qWeights,tmp,k,j,nelem,nNodes)
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
999 continue
    q = qStar

END SUBROUTINE nDGsweep

SUBROUTINE evalExpansion(quadVals,edgeVals,qIn,nelem,nNodes,dozhangshu)
! Evaluate ansatz solution at quadrature nodes and edges
! Note: Currently assumes basis is Lagrange interpolating polynomials at quad nodes -> phi_j (xi_k) = qIn(k,j)
    IMPLICIT NONE
    ! Inputs
    INTEGER, INTENT(IN) :: nelem,nNodes
    LOGICAL, INTENT(IN) :: dozhangshu
    REAL(KIND=8), DIMENSION(0:nNodes,1:nelem), INTENT(IN) :: qIn

    ! Outputs
    REAL(KIND=8), DIMENSION(0:nNodes,1:nelem), INTENT(OUT) :: quadVals
    REAL(KIND=8), DIMENSION(0:1,0:nelem+1), INTENT(OUT) :: edgeVals

    ! Local Variables
    INTEGER :: j

    ! Ansatz value at quad locations for the jth element is just coefficient value
    quadVals = qIn
    
    edgeVals(0,1:nelem) = qIn(0,:) ! left edge value is just left-most coefficient
    edgeVals(1,1:nelem) = qIn(nNodes,:) ! right edge value is just right-most coefficient

    IF(dozhangshu) THEN
        ! To be filled in
        write(*,*) ' WARNING::: THIS ISNT IMPLEMENTED YET'
    ENDIF

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
