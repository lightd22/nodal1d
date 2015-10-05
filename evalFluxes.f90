SUBROUTINE evalFluxes(flx,elemAvg,edgeVals,uin,nelem,nNodes,dt,dx)
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

  ! Local variables
  INTEGER :: j
  DOUBLE PRECISION :: meanLow,permissableNetFlux,actualNetFlux,eps,foo
  DOUBLE PRECISION, DIMENSION(0:nelem) :: fluxLow,fluxCor,corrector
  DOUBLE PRECISION, DIMENSION(0:nelem+1) :: fluxRatio

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
  END INTERFACE

  eps = epsilon(1D0)

  ! Compute high order fluxes
  CALL numFlux(flx,edgeVals,uin,nelem,nNodes)

  ! For non-FCT methods, return high order fluxes
  IF(.not. doFCT) THEN
    RETURN
  ENDIF

  ! Evalutate corrected flux
  fluxLow = 0D0
  fluxCor = flx-fluxLow

  DO j=1,nelem
    meanLow = elemAvg(j)
    permissableNetFlux = meanLow*dx/dt
    actualNetFlux = MAX(0D0,fluxCor(j))-MIN(0D0,fluxCor(j-1))+eps
    fluxRatio(j) = MIN(1D0,permissableNetFlux/actualNetFlux)
  ENDDO !j
  ! Periodic extension
  fluxRatio(0) = fluxRatio(nelem)
  fluxRatio(nelem+1) = fluxRatio(1)

  ! If flux at edge is negative, use limiting ratio from right element
  ! since that is where we are pulling mass from
  DO j=0,nelem
    IF(fluxCor(j) .ge. 0D0) THEN
      corrector(j) = fluxRatio(j)
    ELSE
      corrector(j) = fluxRatio(j+1)
    ENDIF
  ENDDO !j

  ! Update FCT fluxes
  flx = fluxLow + corrector*fluxCor
  !flx = corrector*flx

END SUBROUTINE evalFluxes
