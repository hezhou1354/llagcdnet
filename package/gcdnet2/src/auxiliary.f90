! DESCRIPTION: 
!
!    These functions are minor modifications from the glmnet package:
!
!    Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010). 
!    Regularization Paths for Generalized Linear Models via Coordinate Descent. 
!    Journal of Statistical Software, 33(1), 1-22. 
!    URL http://www.jstatsoft.org/v33/i01/.
!
! --------------------------------------------------------------------------
! standard: An auxiliary function for standardize x matrix.
! --------------------------------------------------------------------------
!
! USAGE:
! 
! call standard (nobs,nvars,x,ju,isd,xmean,xnorm,maj)   
! 
! INPUT ARGUMENTS:
! 
!    nobs = number of observations
!    nvars = number of predictor variables
!    x(nobs, nvars) = matrix of predictors, of dimension N * p; each row is an observation vector.
!    ju(nvars) = flag of predictor variables
!                ju(j) = 0 => this predictor has zero variance
!                ju(j) = 1 => this predictor does not have zero variance
!    isd = standarization flag:
!          isd = 0 => do not standardize predictor variables
!          isd = 1 => standardize predictor variables
!          NOTE: no matter isd is 1 or 0, matrix x is always centered by column. That is, col.mean(x) = 0.
!    
! OUTPUT:
!
!    x(nobs, nvars) = standarized matrix x
!    xmean(nvars) = column mean of x matrix
!    xnorm(nvars) = column standard deviation of x matrix
!    maj(nvars) = column variance of x matrix
!
! --------------------------------------------------------------------------
! chkvars: An auxiliary function for variable check.
! --------------------------------------------------------------------------
!
! USAGE:
! 
! call chkvars (nobs, nvars, x, ju)
! 
! INPUT ARGUMENTS:
! 
!    nobs = number of observations
!    nvars = number of predictor variables
!    x(nobs, nvars) = matrix of predictors, of dimension N * p; each row is an observation vector.
!    y(nobs) = response variable. This argument should be a two-level factor {-1, 1} 
!            for classification.
!    
! OUTPUT:
!
!    ju(nvars) = flag of predictor variables
!                ju(j) = 0 => this predictor has zero variance
!                ju(j) = 1 => this predictor does not have zero variance
!

! --------------------------------------------------
SUBROUTINE standard(nobs,nvars,x,ju,isd,xmean,xnorm,maj)     
! --------------------------------------------------
    IMPLICIT NONE
    ! - - - arg types - - -
    INTEGER::nobs
    INTEGER::nvars
    INTEGER::isd
    INTEGER::ju(nvars)
    DOUBLE PRECISION::x(nobs,nvars)
    DOUBLE PRECISION::xmean(nvars)
    DOUBLE PRECISION::xnorm(nvars)
    DOUBLE PRECISION::maj(nvars)
    ! - - - local declarations - - -
    INTEGER:: j
! - - - begin - - -                                
    DO j=1,nvars                                  
        IF(ju(j)==1) THEN                         
            xmean(j)=sum(x(:,j))/nobs     !mean                        
            x(:,j)=x(:,j)-xmean(j)    
            maj(j)=dot_product(x(:,j),x(:,j))/nobs                                              
              IF(isd==1) THEN
                xnorm(j)=sqrt(maj(j))    !standard deviation               
                x(:,j)=x(:,j)/xnorm(j)
                maj(j)=1.0D0
            ENDIF                                                        
        ENDIF                                     
    ENDDO                             
END SUBROUTINE standard


! --------------------------------------------------
SUBROUTINE chkvars (nobs, nvars, x, ju)
! --------------------------------------------------
      IMPLICIT NONE
    ! - - - arg types - - -
      INTEGER :: nobs
      INTEGER :: nvars
      INTEGER :: ju (nvars)
      DOUBLE PRECISION :: x (nobs, nvars)
    ! - - - local declarations - - -
      INTEGER :: i
      INTEGER :: j
      DOUBLE PRECISION :: t
! - - - begin - - -
      DO j = 1, nvars
         ju (j) = 0
         t = x (1, j)
         DO i = 2, nobs
            IF (x(i, j) /= t) THEN
               ju (j) = 1
               EXIT
            END IF
         END DO
      END DO
END SUBROUTINE chkvars


subroutine normp ( z, p, pdf )

!*****************************************************************************80
!
!! NORMP computes the cumulative density of the standard normal distribution.
!
!  Discussion:
!
!    This is algorithm 5666 from Hart, et al.
!
!  Modified:
!
!    13 January 2008
!
!  Author:
!
!    Original FORTRAN77 version by Alan Miller.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
!    Charles Mesztenyi, John Rice, Henry Thacher,
!    Christoph Witzgall,
!    Computer Approximations,
!    Wiley, 1968,
!    LC: QA297.C64.
!
!  Parameters:
!
!    Input, DOUBLE PRECISION Z, divides the real line into two
!    semi-infinite intervals, over each of which the standard normal
!    distribution is to be integrated.
!
!    Output, DOUBLE PRECISION P, Q, the integrals of the standard normal
!    distribution over the intervals ( - Infinity, Z] and
!    [Z, + Infinity ), respectively.
!
!    Output, DOUBLE PRECISION PDF, the value of the standard normal
!    distribution at Z.
!
  implicit none

  DOUBLE PRECISION :: cutoff = 7.071D+00
  DOUBLE PRECISION expntl
  DOUBLE PRECISION p
  DOUBLE PRECISION :: p0 = 220.2068679123761D+00
  DOUBLE PRECISION :: p1 = 221.2135961699311D+00
  DOUBLE PRECISION :: p2 = 112.0792914978709D+00
  DOUBLE PRECISION :: p3 = 33.91286607838300D+00
  DOUBLE PRECISION :: p4 = 6.373962203531650D+00
  DOUBLE PRECISION :: p5 = 0.7003830644436881D+00
  DOUBLE PRECISION :: p6 = 0.03526249659989109D+00
  DOUBLE PRECISION pdf
  DOUBLE PRECISION q
  DOUBLE PRECISION :: q0 = 440.4137358247522D+00
  DOUBLE PRECISION :: q1 = 793.8265125199484D+00
  DOUBLE PRECISION :: q2 = 637.3336333788311D+00
  DOUBLE PRECISION :: q3 = 296.5642487796737D+00
  DOUBLE PRECISION :: q4 = 86.78073220294608D+00
  DOUBLE PRECISION :: q5 = 16.06417757920695D+00
  DOUBLE PRECISION :: q6 = 1.755667163182642D+00
  DOUBLE PRECISION :: q7 = 0.08838834764831844D+00
  DOUBLE PRECISION :: root2pi = 2.506628274631001D+00
  DOUBLE PRECISION z
  DOUBLE PRECISION zabs

  zabs = abs ( z )
!
!  37 < |Z|.
!
  if ( 37.0D+00 < zabs ) then

    pdf = 0.0D+00
    p = 0.0D+00
!
!  |Z| <= 37.
!
  else

    expntl = exp ( - 0.5D+00 * zabs * zabs )
    pdf = expntl / root2pi
!
!  |Z| < CUTOFF = 10 / sqrt(2).
!
    if ( zabs < cutoff ) then

      p = expntl * (((((( &
          p6   * zabs &
        + p5 ) * zabs &
        + p4 ) * zabs &
        + p3 ) * zabs &
        + p2 ) * zabs &
        + p1 ) * zabs &
        + p0 ) / ((((((( &
          q7   * zabs &
        + q6 ) * zabs &
        + q5 ) * zabs &
        + q4 ) * zabs &
        + q3 ) * zabs &
        + q2 ) * zabs &
        + q1 ) * zabs &
      + q0 )
!
!  CUTOFF <= |Z|.
!
    else

      p = pdf / ( &
        zabs + 1.0D+00 / ( &
        zabs + 2.0D+00 / ( &
        zabs + 3.0D+00 / ( &
        zabs + 4.0D+00 / ( &
        zabs + 0.65D+00 )))))

    end if

  end if

  if ( z < 0.0D+00 ) then
    q = 1.0D+00 - p
  else
    q = p
    p = 1.0D+00 - q
  end if

  return
end
