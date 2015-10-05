SUBROUTINE computeAverages(avgs,coeffs,quadWeights,nelem,nNodes)
  IMPLICIT NONE
  ! Inputs
  INTEGER, INTENT(IN) :: nelem,nNodes
  DOUBLE PRECISION, DIMENSION(0:nNodes,1:nelem), INTENT(IN) :: coeffs
  DOUBLE PRECISION, DIMENSION(0:nNodes), INTENT(IN) :: quadWeights
  ! Outputs
  DOUBLE PRECISION, DIMENSION(1:nelem), INTENT(OUT) :: avgs
  ! Local variables
  INTEGER :: j

  DO j=1,nelem
    avgs(j) = 0.5D0*SUM(quadWeights(:)*coeffs(:,j))
  ENDDO
  
END SUBROUTINE computeAverages
