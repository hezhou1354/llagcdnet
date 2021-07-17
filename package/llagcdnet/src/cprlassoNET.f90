! --------------------------------------------------------------------------
! cprlassoNET.f90: the GCD algorithm for composite probit regression.
! --------------------------------------------------------------------------
! 
! USAGE:
! 
! call cprlassoNET (lam2, nobs, nvars, x, cpry, jd, pf, pf2, dfmax, &
! & pmax, nlam, flmin, ulam, eps, isd, maxit, &
! ! begin He
! & nthrs, wt, &
! ! end He
! & nalam, b0, beta, ibeta, &
! & nbeta, alam, npass, jerr)
! 
! INPUT ARGUMENTS:
! 
!    lam2 = regularization parameter for the quadratic penalty of the coefficients
!    nobs = number of observations
!    nvars = number of predictor variables
!    x(nobs, nvars) = matrix of predictors, of dimension N * p; each row is an observation vector.
!    ! begin He
!    cpry(nobs, nvars) = response matrix, of dimension N * p; each row is an observation vector.
!    ! end He
!    jd(jd(1)+1) = predictor variable deletion flag
!                  jd(1) = 0  => use all variables
!                  jd(1) != 0 => do not use variables jd(2)...jd(jd(1)+1)
!    pf(nvars) = relative L1 penalties for each predictor variable
!                pf(j) = 0 => jth variable unpenalized
!    pf2(nvars) = relative L2 penalties for each predictor variable
!                pf2(j) = 0 => jth variable unpenalized
!    dfmax = limit the maximum number of variables in the model.
!            (one of the stopping criterion)
!    pmax = limit the maximum number of variables ever to be nonzero. 
!           For example once beta enters the model, no matter how many 
!           times it exits or re-enters model through the path, it will 
!           be counted only once. 
!    nlam = the number of lambda values
!    flmin = user control of lambda values (>=0)
!            flmin < 1.0 => minimum lambda = flmin*(largest lambda value)
!            flmin >= 1.0 => use supplied lambda values (see below)
!    ulam(nlam) = user supplied lambda values (ignored if flmin < 1.0)
!    eps = convergence threshold for coordinate majorization descent. 
!          Each inner coordinate majorization descent loop continues 
!          until the relative change in any coefficient is less than eps.
!    isd = standarization flag:
!          isd = 0 => regression on original predictor variables
!          isd = 1 => regression on standardized predictor variables
!          Note: output solutions always reference original
!                variables locations and scales.
!    maxit = maximum number of outer-loop iterations allowed at fixed lambda value. 
!            (suggested values, maxit = 100000)
!    ! begin He
!    nthrs = number of thresholds for composite probit regression
!    wt(nthrs) = weights for each threshold
!                wt(j) = 0 => jth composite not included
!    ! end He
! 
! OUTPUT:
! 
!    nalam = actual number of lambda values (solutions)
!    b0(nalam) = intercept values for each solution
!    beta(pmax, nalam) = compressed coefficient values for each solution
!    ibeta(pmax) = pointers to compressed coefficients
!    nbeta(nalam) = number of compressed coefficients for each solution
!    alam(nalam) = lambda values corresponding to each solution
!    npass = actual number of passes over the data for all lambda values
!    jerr = error flag:
!           jerr  = 0 => no error
!           jerr > 0 => fatal error - no output returned
!                    jerr < 7777 => memory allocation error
!                    jerr = 7777 => all used predictors have zero variance
!                    jerr = 10000 => maxval(vp) <= 0.0
!                    ! beging He
!                    jerr = 10001 => maxval(wt) <= 0.0
!                    ! end He
!           jerr < 0 => non fatal error - partial output:
!                    Solutions for larger lambdas (1:(k-1)) returned.
!                    jerr = -k => convergence for kth lambda value not reached
!                           after maxit (see above) iterations.
!                    jerr = -10000-k => number of non zero coefficients along path
!                           exceeds pmax (see above) at kth lambda value.
! 
! LICENSE: GNU GPL (version 2 or later)
! 
! AUTHORS:
!    Yi Yang (yi.yang6@mcgill.ca)  
!    and Hui Zou (hzou@stat.umn.edu).
!    School of Statistics, University of Minnesota.
! 
! REFERENCES:
!    Yang, Y. and Zou, H. (2012). An Efficient Algorithm for Computing 
!    The HHSVM and Its Generalizations.
!    Journal of Computational and Graphical Statistics, 22, 396-415.


! --------------------------------------------------
SUBROUTINE cprlassoNET (lam2, nobs, nvars, x, cpry, jd, pf, pf2, dfmax, pmax, &
& nlam, flmin, ulam, eps, isd, maxit, &
! begin He
& nthrs, wt, &
! end He
& nalam, b0, beta, ibeta, nbeta, alam, npass, jerr)
! --------------------------------------------------
      IMPLICIT NONE
    ! - - - arg types - - -
      INTEGER :: nobs
      INTEGER :: nvars
      INTEGER :: dfmax
      INTEGER :: pmax
      INTEGER :: nlam
      INTEGER :: isd
      ! begin He
      INTEGER :: nthrs
      ! end He
      INTEGER :: nalam
      INTEGER :: npass
      INTEGER :: jerr
      INTEGER :: maxit
      INTEGER :: jd (*)
      INTEGER :: ibeta (pmax)
      INTEGER :: nbeta (nlam)
      DOUBLE PRECISION :: lam2
      DOUBLE PRECISION :: flmin
      DOUBLE PRECISION :: eps
      DOUBLE PRECISION :: x (nobs, nvars)
      DOUBLE PRECISION :: cpry (nobs, nthrs)
      DOUBLE PRECISION :: pf (nvars)
      DOUBLE PRECISION :: pf2 (nvars)
      DOUBLE PRECISION :: ulam (nlam)
      DOUBLE PRECISION :: beta (pmax, nlam)
      DOUBLE PRECISION :: b0 (nthrs, nlam)
      DOUBLE PRECISION :: alam (nlam)
      ! begin He
      DOUBLE PRECISION :: wt (nthrs)
      ! end He
    ! - - - local declarations - - -
      ! begin He
      INTEGER :: t
      ! end He
      INTEGER :: j
      INTEGER :: l
      INTEGER :: nk
      INTEGER :: ierr
      INTEGER, DIMENSION (:), ALLOCATABLE :: ju
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: xmean
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: xnorm
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: maj
! - - - begin - - -
! - - - allocate variables - - -
      jerr = 0
      ALLOCATE (ju(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (xmean(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (maj(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (xnorm(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      IF (jerr /= 0) RETURN
      CALL chkvars (nobs, nvars, x, ju)
      IF (jd(1) > 0) ju (jd(2:(jd(1)+1))) = 0
      IF (maxval(ju) <= 0) THEN
         jerr = 7777
         RETURN
      END IF
      IF (maxval(pf) <= 0.0D0) THEN
         jerr = 10000
         RETURN
      END IF
      IF (maxval(pf2) <= 0.0D0) THEN
         jerr = 10000
         RETURN
      END IF
      ! begin He
      IF (maxval(wt) <= 0.0D0) THEN
         jerr = 10000
         RETURN
      END IF
      wt = Max (0.0D0, wt)
      ! end He
      pf = Max (0.0D0, pf)
      pf2 = Max (0.0D0, pf2)
      CALL standard (nobs, nvars, x, ju, isd, xmean, xnorm, maj)
      CALL cprlassoNETpath (lam2, maj, nobs, nvars, x, cpry, ju, pf, pf2, dfmax, &
     & pmax, nlam, flmin, ulam, eps, maxit, nalam, b0, beta, ibeta, &
     & nbeta, alam, npass, jerr, &
     ! begin He
     & nthrs, wt)
     ! end He
      IF (jerr > 0) RETURN! check error after calling function
! - - - organize beta afterward - - -
      DO l = 1, nalam
         nk = nbeta (l)
         IF (isd == 1) THEN
            DO j = 1, nk
               beta (j, l) = beta (j, l) / xnorm (ibeta(j))
            END DO
         END IF
         ! begin He
         DO t = 1, nthrs
            b0 (t, l) = b0 (t, l) + dot_product (beta(1:nk, l), &
            & xmean(ibeta(1:nk)))
         END DO
         ! end He
      END DO
      DEALLOCATE (ju, xmean, xnorm, maj)
      RETURN
END SUBROUTINE cprlassoNET
! --------------------------------------------------
SUBROUTINE cprlassoNETpath (lam2, maj, nobs, nvars, x, cpry, ju, pf, pf2, dfmax, &
& pmax, nlam, flmin, ulam, eps, maxit, nalam, b0, beta, m, nbeta, alam, &
& npass, jerr, &
! begin He
& nthrs, wt)
! end He
! --------------------------------------------------
      IMPLICIT NONE
        ! - - - arg types - - -
      DOUBLE PRECISION, PARAMETER :: big = 9.9E30
      DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
      ! begin He: 1/sqrt(2*pi)
      DOUBLE PRECISION, PARAMETER :: sqrtpi = 3.9894228040143267794D-1
      ! end He
      INTEGER, PARAMETER :: mnlam = 6
      INTEGER :: mnl
      INTEGER :: nobs
      INTEGER :: nvars
      INTEGER :: dfmax
      INTEGER :: pmax
      INTEGER :: nlam
      INTEGER :: maxit
      ! begin He
      INTEGER :: nthrs
      ! end He
      INTEGER :: nalam
      INTEGER :: npass
      INTEGER :: jerr
      INTEGER :: ju (nvars)
      INTEGER :: m (pmax)
      INTEGER :: nbeta (nlam)
      DOUBLE PRECISION :: lam2
      DOUBLE PRECISION :: eps
      DOUBLE PRECISION :: x (nobs, nvars)
      DOUBLE PRECISION :: cpry (nobs, nthrs)
      DOUBLE PRECISION :: pf (nvars)
      DOUBLE PRECISION :: pf2 (nvars)
      ! begin He
      DOUBLE PRECISION :: wt (nthrs)
      ! end He
      DOUBLE PRECISION :: beta (pmax, nlam)
      DOUBLE PRECISION :: ulam (nlam)
      DOUBLE PRECISION :: b0 (nthrs, nlam)
      DOUBLE PRECISION :: alam (nlam)
      DOUBLE PRECISION :: maj (nvars)
    ! - - - local declarations - - -
      ! begin He
      DOUBLE PRECISION :: d
      ! end He
      DOUBLE PRECISION :: dif
      DOUBLE PRECISION :: oldb
      DOUBLE PRECISION :: u
      DOUBLE PRECISION :: v
      DOUBLE PRECISION :: al
      DOUBLE PRECISION :: alf
      DOUBLE PRECISION :: flmin
      ! begin He
      DOUBLE PRECISION :: dl (nobs, nthrs)
      DOUBLE PRECISION :: pdf
      DOUBLE PRECISION :: cdf
      DOUBLE PRECISION :: rr
      DOUBLE PRECISION :: ul (nthrs)
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: bzero
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: oldbzero
      ! end He
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: oldbeta
      
      DOUBLE PRECISION, DIMENSION (:, :), ALLOCATABLE :: r
      ! begin He
      INTEGER :: i
      INTEGER :: t
      ! end He
      INTEGER :: k
      INTEGER :: j
      INTEGER :: l
      INTEGER :: vrg
      INTEGER :: ierr
      INTEGER :: ni
      INTEGER :: me
      INTEGER, DIMENSION (:), ALLOCATABLE :: mm
! - - - begin - - -
! - - - allocate variables - - -
      ! begin He
      ALLOCATE (b(1:nvars), STAT=jerr)
      ALLOCATE (bzero(1:nthrs), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (oldbeta(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (oldbzero(1:nthrs), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (mm(1:nvars), STAT=ierr)
      jerr = jerr + ierr
      ALLOCATE (r(1:nobs, 1:nthrs), STAT=ierr)
      jerr = jerr + ierr
      ! end He
      IF (jerr /= 0) RETURN
! - - - some initial setup - - -
      ! begin He
      r = 0.0D0
      b = 0.0D0
      bzero = 0.0D0
      oldbeta = 0.0D0
      oldbzero = 0.0D0
      ! end He
      m = 0
      mm = 0
      npass = 0
      ni = npass
      mnl = Min (mnlam, nlam)
      ! begin He
      maj = 1.0D0 * maj
      alf = 1.0D0
      ! end He
      IF (flmin < 1.0D0) THEN
         flmin = Max (mfl, flmin)
         alf = flmin ** (1.0D0/(DBLE(nlam)-1.0D0))
      END IF
! --------- lambda loop ----------------------------
      DO l = 1, nlam
! --------- computing lambda ----------------------------
         IF (flmin >= 1.0D0) THEN
            al = ulam (l)
         ELSE
            IF (l > 2) THEN
               al = al * alf
            ELSE IF (l == 1) THEN
               al = big
            ELSE IF (l == 2) THEN
               al = 0.0D0
               ! begin He:
               DO t = 1, nthrs
                  DO i = 1, nobs
                     rr = r(i, t)
                     CALL normp(rr, cdf, pdf)
                     dl(i, t) = pdf / cdf
                  END DO
               END DO
               DO j = 1, nvars
                  IF (ju(j) /= 0) THEN
                     IF (pf(j) > 0.0D0) THEN
                        DO t = 1, nthrs
                           ul(t) = dot_product (cpry(:, t) * dl(:, t), x(:, j))
                        END DO
                        u = sum (wt * ul)
                        al = Max (al, Abs(u)/pf(j))
                     END IF
                  END IF
               END DO
               al = al * alf / nobs
               ! end He
            END IF
         END IF
        ! --------- outer loop ----------------------------
         DO
            ! begin He
            oldbzero = bzero
            IF (ni > 0) oldbeta (m(1:ni)) = b (m(1:ni))
            ! end He
        ! --middle loop-------------------------------------
            DO
               npass = npass + 1
               dif = 0.0D0
               DO k = 1, nvars
                  IF (ju(k) /= 0) THEN
                     ! begin He
                     oldb = b (k)
                     DO t = 1, nthrs
                        DO i = 1, nobs
                           rr = r(i, t)
                           CALL normp(rr, cdf, pdf)
                           dl(i, t) = pdf / cdf
                        END DO
                        ul(t) = dot_product (cpry(:, t) * dl(:, t), x(:, k))
                     END DO
                     u = sum (wt * ul)
                     u = maj (k) * b (k) + u / nobs
                     ! end He
                     v = al * pf (k)
                     v = Abs (u) - v
                     IF (v > 0.0D0) THEN
                        b (k) = sign (v, u) / (maj(k) + pf2(k) * lam2)
                     ELSE
                        b (k) = 0.0D0
                     END IF
                     d = b (k) - oldb
                     IF (Abs(d) > 0.0D0) THEN
                        ! begin He
                        dif = Max (dif, d**2)
                        !dif = Max (dif, d)
                        DO t = 1, nthrs
                           r(:, t) = r(:, t) + cpry(:, t) * x(:, k) * d
                        END DO
                        ! end He
                        IF (mm(k) == 0) THEN
                           ni = ni + 1
                           IF (ni > pmax) EXIT
                           mm (k) = ni
                           m (ni) = k !indicate which one is non-zero
                        END IF
                     END IF
                  END IF
               END DO
               ! begin He
               DO t = 1, nthrs
                  DO i = 1, nobs
                     rr = r(i, t)
                     CALL normp(rr, cdf, pdf)
                     dl(i, t) = pdf / cdf
                  END DO
                  d = dot_product (dl(:, t) , cpry(:, t)) / nobs
                  IF (Abs(d) > 0.0D0) THEN
                     bzero(t) = bzero(t) - d
                     r(:, t) = r(:, t) + cpry(:, t) * d
                     dif = Max (dif, d**2)
                  END IF
               END DO
               ! end He
               IF (ni > pmax) EXIT
               IF (dif < eps) EXIT
               IF (npass > maxit) THEN
                  jerr = -l
                  RETURN
               END IF
        ! --inner loop----------------------
               DO
                  npass = npass + 1
                  dif = 0.0D0
                  DO j = 1, ni
                     k = m (j)
                     oldb = b (k)
                     ! begin He
                     DO t = 1, nthrs
                        DO i = 1, nobs
                           rr = r(i, t)
                           CALL normp(rr, cdf, pdf)
                           dl(i, t) = pdf / cdf
                        END DO
                        ul(t) = dot_product (cpry(:, t) * dl(:, t), x(:, k))
                     END DO
                     u = sum (wt * ul)
                     u = maj (k) * b (k) + u / nobs
                     ! end He
                     v = al * pf (k)
                     v = Abs (u) - v
                     IF (v > 0.0D0) THEN
                        b (k) = sign (v, u) / (maj(k) + pf2(k) * lam2)
                     ELSE
                        b (k) = 0.0D0
                     END IF
                     d = b (k) - oldb
                     IF (Abs(d) > 0.0D0) THEN
                        ! begin He
                        dif = Max (dif, d**2)
                        !dif = Max (dif, d)
                        DO t = 1, nthrs
                           r(:, t) = r(:, t) + cpry(:, t) * x (:, k) * d
                        END DO
                        ! end He
                     END IF
                  END DO
                  ! begin He
                  DO t = 1, nthrs
                     DO i = 1, nobs
                        rr = r(i, t)
                        CALL normp(rr, cdf, pdf)
                        dl(i, t) = pdf / cdf
                     END DO
                     d = dot_product (dl(:, t) , cpry(:, t)) / nobs
                     IF (Abs(d) > 0.0D0) THEN
                        bzero (t) = bzero (t) - d
                        r(:, t) = r(:, t) + cpry(:, t) * d
                        dif = Max (dif, d**2)
                        !dif = Max (dif, Abs(d))
                     END IF
                  END DO
                  ! end He
                  IF (dif < eps) EXIT
                  IF (npass > maxit) THEN
                     jerr=-l
                     RETURN
                  ENDIF
               END DO
            END DO
            IF (ni > pmax) EXIT
        !--- this is the final check ------------------------
            vrg = 1
            ! begin He
            DO t = 1, nthrs
               IF ((bzero(t)-oldbzero(t))**2 >= eps) THEN
                  vrg = 0
                  EXIT
               END IF
            END DO
            DO j = 1, ni
               IF ((b(m(j))-oldbeta(m(j)))**2 >= eps) THEN
                  vrg = 0
                  EXIT
               END IF
            END DO
            IF (vrg == 1) EXIT
         END DO
    ! final update variable save results------------
         IF (ni > pmax) THEN
            jerr = - 10000 - l
            EXIT
         END IF
         IF (ni > 0) beta (1:ni, l) = b (m(1:ni))
         nbeta (l) = ni
         b0 (:, l) = bzero
         alam (l) = al
         nalam = l
         IF (l < mnl) CYCLE
         IF (flmin >= 1.0D0) CYCLE
         me = count (beta(1:ni, l) /= 0.0D0)
         IF (me > dfmax) EXIT
      END DO
      DEALLOCATE (b, bzero, oldbeta, oldbzero, r, mm)
      RETURN
! begin He
END SUBROUTINE cprlassoNETpath
! end He

