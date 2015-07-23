SUBROUTINE positivityLimit(qIn,avgVals,qNodes,qWeights,nelem,nNodes,nQuad,lagVals,stat)
  ! Rescales polynomial so that nodal values are non-negative
  USE testParameters
  IMPLICIT NONE
  ! Inputs
  INTEGER, INTENT(IN) :: nelem,nNodes,nQuad,stat
  REAL(KIND=8), DIMENSION(0:nNodes,1:nelem), INTENT(INOUT) :: qIn
  REAL(KIND=8), DIMENSION(1:nelem), INTENT(IN) :: avgVals
  REAL(KIND=8), DIMENSION(0:nQuad), INTENT(IN) :: qWeights,qNodes
  REAL(KIND=8), DIMENSION(0:nNodes,0:nQuad), INTENT(IN) :: lagVals
  ! Local Variables
  INTEGER :: j,l,flag
  INTEGER, DIMENSION(0:nNodes) :: msk
  REAL(KIND=8), DIMENSION(0:nQuad) :: qVals
  REAL(KIND=8) :: avgVal,theta,valMin,valMax,eps,Mt,Mp,lowrBnd,uprBnd,qStar

  eps = epsilon(1d0)
  ! First check incoming element averages
  IF(limitingMeth == 1) THEN
    DO j=1,nelem
      ! -- Evaluate DG polynomial at quad locations
      ! Note: This is not necessary if the quad nodes used are same nodes that the basis is interpolating --
      !       In this case the polynomial coefficients are the nodal values. Two exceptions are the edge
      !       values (-1 and 1) which are always part of any GLL integration. These may be read from coefficients directly
      qVals(0) = qIn(0,j)
      qVals(nQuad) = qIn(nNodes,j)
      DO l=1,nQuad-1
        qVals(l) = SUM(qIn(:,j)*lagVals(:,l))
      ENDDO !l

      avgVal = 0.5D0*SUM( qWeights(:)*qVals(:) )
!     valMin = MINVAL(qVals(:))-eps
      valMin = MIN( MINVAL(qIn(1:nNodes-1,j)) , MINVAL(qVals(:)) ) - eps

      ! -- Compute rescaling factor
      theta = MIN( abs(avgVal/(valMin-avgVal)),1D0 )

      qIn(:,j) = theta*(qIn(:,j)-avgVal) + avgVal
    ENDDO !j
  ELSE IF(limitingMeth==2) THEN
    DO j=1,nelem
      ! Use nodes themselves to do limiting, rather than evaluating polynomial at multiple locations
      avgVal = 0.5D0*SUM( qWeights(:)*qIn(:,j) )
      valMin = MINVAL(qIn(:,j))-eps

      IF(avgVal .lt. 0D0) THEN
          write(*,*) 'Element AVG in Z&S Limiter is negative!! Negative average value = ',avgVal
      ENDIF

      ! -- Compute rescale factor
      theta = MIN( abs(avgVal/(valMin-avgVal)),1D0 )
!     theta = MIN( abs(avgVal/(valMin-avgVal)),abs((1D0-avgVal)/(valMax-avgVal)),1D0 )

      ! -- Rescale polynomial
      qIn(:,j) = theta*(qIn(:,j)-avgVal) + avgVal
    ENDDO !j
  ELSE IF(limitingMeth==3) THEN
    ! Use "mass aware truncation" for rescaling
    DO j=1,nelem
      Mp = 0D0
      Mt = 0D0

      DO l=0,nNodes
          Mt = Mt + qWeights(l)*qIn(l,j)
          qIn(l,j) = MAX(0D0,qIn(l,j)) ! Zero out negative nodes
          Mp = Mp + qWeights(l)*qIn(l,j)
      ENDDO !l
      theta = MAX(Mt,0D0)/MAX(Mp,TINY(1D0))
      qIn(:,j) = theta*qIn(:,j) ! Reduce remaining positive nodes by reduction factor
    ENDDO !j
  ELSE IF(limitingMeth==4) THEN
    ! Use "mass aware" truncation bounds preserving limiter
    lowrBnd = 0D0
    uprBnd = 1D0
    msk = 1
    DO j=1,nelem
      Mp = 0D0
      Mt = 0D0
      DO l=0,nNodes
        Mt = Mt + qWeights(l)*qIn(l,j)
        Mp = Mp + qWeights(l)*qIn(l,j)
        qStar = MAX(lowrBnd,MIN(uprBnd,qIn(l,j)))

        flag = 0.5*(uprBnd-lowrBnd-(abs(uprBnd-qIn(l,j))+abs(qIn(l,j)-lowrBnd)) )
        qIn(l,j) = qStar
        msk(l) = flag
      ENDDO!l
      theta = Mt/MAX(Mp,eps)
      WHERE(msk.eq.0)
        qIn(:,j) = theta*qIn(:,j)
      ENDWHERE
      IF( MAXVAL(qIn(:,j)) .gt. uprBnd) THEN
        write(*,*) 'After limiting, maximal value is:',maxval(qIn(:,j))
        write(*,*) 'Flagging array:',msk
      ENDIF

    ENDDO !j
  ENDIF !limitingMeth

END SUBROUTINE positivityLimit
