!--------------------------------------------------------------------------------
! Nodal DG coefficient update for 1D advection
! By: Devin Light 4/21/14
! -------------------------------------------------------------------------------

SUBROUTINE nDGsweep(q,nelem,dxel,nNodes,qNodes,qWeights,u,lagDeriv,dozhangshu,dt)
    IMPLICIT NONE
    INTEGER, PARAMETER :: DOUBLE = KIND(1D0) ! Specification of DOUBLE
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
    REAL(KIND=DOUBLE), DIMENSION(0:nNodes,1:nelem) :: qFwd, qStar,quadVals
    REAL(KIND=DOUBLE), DIMENSION(0:1,0:nelem+1) :: edgeVals
    INTEGER :: k,j,stage

    qStar = q
    ! Do SSPRK3 Update
    DO stage=1,3
        CALL evalExpansion(quadVals,edgeVals,qStar,nelem,nNodes,dozhangshu)
        CALL numFlux(flx,edgeVals,u,nelem,nNodes)

        ! Forward Step
        DO j=1,nelem
            DO k=0,nNodes
                qFwd(k,j) = qStar + (dt/dxel)*dadt(quadVals(:,j),flx,u(:,j),qWeights,lagDeriv(k,:),k,j,nelem,nNodes)
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
    edgeVals(1,1:nelem) = qIn(1,:) ! right edge value is just right-most coefficient

    IF(dozhangshu) THEN
        ! To be filled in
        write(*,*) ' WARNING::: THIS ISNT IMPLEMENTED YET'
    ENDIF

	! Extend edgeVals periodically
	edgeVals(:,0) = edgeVals(:,nelem)
	edgeVals(:,nelem+1) = edgeVals(:,1)

END SUBROUTINE evalExpansion

SUBROUTINE numFlux(flx,edgeVals,u,nelem,nNodes)
	IMPLICIT NONE
	INTEGER, PARAMETER :: DOUBLE = KIND(1D0)
	! -- Inputs
	INTEGER, INTENT(IN) :: nelem,nNodes
	REAL(KIND=DOUBLE), DIMENSION(0:1,0:nelem+1), INTENT(IN) :: edgeVals
	REAL(KIND=DOUBLE), DIMENSION(0:nNodes,1:nelem), INTENT(IN) :: u

	! -- Outputs	
	REAL(KIND=DOUBLE),DIMENSION(0:nelem), INTENT(OUT) :: flx

	! -- Local variables
	INTEGER :: j

	DO j=0,nelem
		flx(j) = 0.5D0*edgeVals(1,j)*(u(nNodes,j)+DABS(u(nNodes,j)))+0.5D0*edgeVals(0,j+1)*(u(nNodes,j)-DABS(u(nNodes,j)))
	ENDDO
END SUBROUTINE numFlux
