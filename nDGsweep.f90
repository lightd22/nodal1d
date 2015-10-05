!--------------------------------------------------------------------------------
! Nodal DG coefficient update for 1D advection
! By: Devin Light 4/21/14
! -------------------------------------------------------------------------------

SUBROUTINE nDGsweep(q,nelem,dxel,nNodes,qNodes,qWeights,u,lagDeriv,doposlimit,dt,&
                    nZSNodes,quadZSNodes,quadZSWeights,lagValsZS)
    USE testParameters, ONLY: dofct
    IMPLICIT NONE
    ! External Functions
	  DOUBLE PRECISION, EXTERNAL :: dadt ! RHS function for evolution ODE for kth expansion coefficent

    ! Inputs
    INTEGER, INTENT(IN) :: nelem,nNodes
    INTEGER, INTENT(IN) :: nZSNodes
    LOGICAL, INTENT(IN) :: doposlimit
    DOUBLE PRECISION, INTENT(IN) :: dxel,dt
    DOUBLE PRECISION, DIMENSION(0:nNodes), INTENT(IN):: qNodes,qWeights
    DOUBLE PRECISION, DIMENSION(0:nZSNodes), INTENT(IN) :: quadZSNodes,quadZSWeights
    DOUBLE PRECISION, DIMENSION(0:nNodes,0:nZSNodes), INTENT(IN) :: lagValsZS
    DOUBLE PRECISION, DIMENSION(0:nNodes,0:nNodes), INTENT(IN) :: lagDeriv
    DOUBLE PRECISION, DIMENSION(0:nNodes,1:nelem), INTENT(IN) :: u

    ! Outputs
    DOUBLE PRECISION, DIMENSION(0:nNodes,1:nelem), INTENT(INOUT) :: q

    ! Local Variables
    DOUBLE PRECISION, DIMENSION(0:nNodes,1:nelem) :: qFwd,qStar,quadVals
    DOUBLE PRECISION, DIMENSION(0:nNodes,0:nelem) :: utmp
    DOUBLE PRECISION, DIMENSION(0:1,0:nelem+1) :: edgeVals
    DOUBLE PRECISION, DIMENSION(0:nelem) :: flx
    DOUBLE PRECISION, DIMENSION(0:nNodes) :: tmp

    DOUBLE PRECISION, DIMENSION(1:nelem) :: avgVals,prevAvg
    DOUBLE PRECISION, DIMENSION(0:nZSNodes) :: qVals

    INTEGER :: k,j,stage,l
    LOGICAL :: stopFlag = .FALSE.
    DOUBLE PRECISION :: projAvg,actAvg

    INTERFACE
      SUBROUTINE limitNodePositivity(qIn,avgVals,qWeights,nelem,nNodes,nQuad)
        ! Modifies approximating polynomial so that nodal values are non-negative for output
        USE testParameters
        IMPLICIT NONE
        ! Inputs
        INTEGER, INTENT(IN) :: nelem,nNodes,nQuad
        DOUBLE PRECISION, DIMENSION(0:nNodes,1:nelem), INTENT(INOUT) :: qIn
        DOUBLE PRECISION, DIMENSION(1:nelem), INTENT(IN) :: avgVals
        DOUBLE PRECISION, DIMENSION(0:nQuad), INTENT(IN) :: qWeights
      END SUBROUTINE limitNodePositivity

      SUBROUTINE limitMeanPositivity(qIn,avgVals,qWeights,nelem,nNodes,nQuad)
        ! Modifies approximating polynomial so that element mean value remains non-negative
        ! according to Zhang and Shu (2010) Thm 2.2
        USE testParameters
        IMPLICIT NONE
        ! Inputs
        INTEGER, INTENT(IN) :: nelem,nNodes,nQuad
        DOUBLE PRECISION, DIMENSION(0:nNodes,1:nelem), INTENT(INOUT) :: qIn
        DOUBLE PRECISION, DIMENSION(1:nelem), INTENT(IN) :: avgVals
        DOUBLE PRECISION, DIMENSION(0:nQuad), INTENT(IN) :: qWeights
      END SUBROUTINE limitMeanPositivity

      SUBROUTINE evalFluxes(flx,elemAvg,edgeVals,uin,nelem,nNodes,dt,dx,coeffs)
        ! ------------------------------------------
        !
        ! ------------------------------------------
        USE testParameters
        IMPLICIT NONE
        ! Inputs
        INTEGER, INTENT(IN) :: nelem,nNodes
        DOUBLE PRECISION, INTENT(IN) :: dt,dx
        DOUBLE PRECISION, DIMENSION(0:1,0:nelem+1), INTENT(IN) :: edgeVals
        DOUBLE PRECISION, DIMENSION(0:nNodes,0:nelem), INTENT(IN) :: uin
        DOUBLE PRECISION, DIMENSION(1:nelem), INTENT(IN) :: elemAvg
        ! Outputs
        DOUBLE PRECISION, DIMENSION(0:nelem), INTENT(OUT) :: flx
        DOUBLE PRECISION, DIMENSION(0:nNodes,1:nelem), INTENT(INOUT) :: coeffs

      END SUBROUTINE evalFluxes

      SUBROUTINE computeAverages(avgs,coeffs,quadWeights,nelem,nNodes)
        IMPLICIT NONE
        ! Inputs
        INTEGER, INTENT(IN) :: nelem,nNodes
        DOUBLE PRECISION, DIMENSION(0:nNodes,1:nelem), INTENT(IN) :: coeffs
        DOUBLE PRECISION, DIMENSION(0:nNodes), INTENT(IN) :: quadWeights
        ! Outputs
        DOUBLE PRECISION, DIMENSION(1:nelem), INTENT(OUT) :: avgs
      END SUBROUTINE computeAverages
    END INTERFACE

    utmp(:,1:nelem) = u
    utmp(:,0) = u(:,nelem)

    qStar = q
    CALL computeAverages(prevAvg,qStar,qWeights,nelem,nNodes)
    avgVals = prevAvg
    ! SSPRK3 Update
    DO stage=1,3

      ! Rescale for element-mean positivity (if necessary)
      IF(doposlimit) THEN
        DO j=1,nelem
          avgVals(j) = 0.5D0*SUM( qStar(:,j)*qWeights )
          IF(avgVals(j) .lt. 0D0) THEN
            write(*,'(A,I2,A,I2,A,D10.4)') ' Average of element ', j,&
                                           ' is negative in stage ',stage, &
                                           ' avgVal = ',avgVals(j)
            stopFlag = .TRUE.
          ENDIF
        ENDDO !j
        CALL limitMeanPositivity(qStar,avgVals,quadZSWeights,nelem,nNodes,nZSNodes)
      ENDIF

      CALL evalExpansion(quadVals,edgeVals,qStar,nelem,nNodes)
      CALL evalFluxes(flx,avgVals,edgeVals,utmp,nelem,nNodes,dt,dxel,qStar)
      !CALL numFlux(flx,edgeVals,utmp,nelem,nNodes)

      ! Forward Step
      DO j=1,nelem
        DO k=0,nNodes
          tmp = lagDeriv(k,:)
          qFwd(k,j) = qStar(k,j) + dt*dadt(quadVals(:,j),flx,utmp(:,j),qWeights,tmp,k,j,nelem,nNodes)/dxel
        ENDDO !k
        qFwd(0,j) = qFwd(0,j)+2D0*dt*flx(j-1)/qWeights(0)/dxel
        qFwd(nNodes,j) = qFwd(nNodes,j)-2D0*dt*flx(j)/qWeights(nNodes)/dxel
        ! Track change to averages
        avgVals(j) = avgVals(j)-(dt/dxel)*(flx(j)-flx(j-1))
      ENDDO ! j

  		SELECT CASE(stage)
  		CASE(1)
  			qStar = qFwd
  		CASE(2)
  			qStar = 0.75D0*q + 0.25D0*qFwd
        avgVals = 0.75D0*prevAvg+0.25D0*avgVals
  		CASE(3)
  			qStar = q/3D0 + 2D0*qFwd/3D0
        avgVals = prevAvg/3D0+2D0*avgVals/3D0
  		END SELECT

      IF(doposlimit) THEN
        DO j=1,nelem
!          avgVals(j) = 0.5D0*SUM(qStar(:,j)*qWeights)
          IF(avgVals(j) .lt. 0d0) THEN
            write(*,'(A,I2,A,I2)') 'Average of element',j,' is negative after stage',stage
            write(*,'(A,D10.4)') 'Problem average=',avgVals(j)
            stopFlag = .TRUE.
          ENDIF
        ENDDO !j
      ENDIF

      IF(stopFlag) THEN
        write(*,*) '********* CRITICAL ERROR *********'
        write(*,*) 'S T O P P I N G...'
        write(*,*) '********* CRITICAL ERROR *********'
        STOP
      ENDIF
    ENDDO !stage

    ! After RK-update, rescale polynomial to remove all non-negatives
    IF(doposlimit) THEN
      DO j=1,nelem
        avgVals(j) = 0.5D0*SUM( qStar(:,j)*qWeights )
      ENDDO!j
      CALL limitNodePositivity(qStar,avgVals,qWeights,nelem,nNodes,nNodes)
    ENDIF
    q = qStar

END SUBROUTINE nDGsweep

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
	! -- Inputs
	INTEGER, INTENT(IN) :: nelem,nNodes
	DOUBLE PRECISION, DIMENSION(0:1,0:nelem+1), INTENT(IN) :: edgeVals
	DOUBLE PRECISION, DIMENSION(0:nNodes,0:nelem), INTENT(IN) :: uin

	! -- Outputs
	DOUBLE PRECISION,DIMENSION(0:nelem), INTENT(OUT) :: flx

	! -- Local variables
	INTEGER :: j

	DO j=0,nelem
		flx(j) = 0.5D0*edgeVals(1,j)*(uin(nNodes,j)+DABS(uin(nNodes,j)))+0.5D0*edgeVals(0,j+1)*(uin(nNodes,j)-DABS(uin(nNodes,j)))
	ENDDO
END SUBROUTINE numFlux

DOUBLE PRECISION FUNCTION dadt(quadVals,flx,u,qWeights,lagDeriv,k,j,nelem,nNodes)
    IMPLICIT NONE
    ! Inputs
    INTEGER, INTENT(IN) :: nelem,nNodes,k,j
    DOUBLE PRECISION, DIMENSION(0:nNodes) :: quadVals,u,qWeights,lagDeriv
    DOUBLE PRECISION, DIMENSION(0:nelem) :: flx

    dadt = SUM( u(:)*quadVals(:)*lagDeriv(:)*qWeights(:) )

!    IF( k .eq. 0) THEN
!        dadt = dadt + flx(j-1)
!    ENDIF
!    IF( k .eq. nNodes) THEN
!        dadt = dadt - flx(j)
!    ENDIF
    dadt = 2D0*dadt/qWeights(k)

END FUNCTION dadt
