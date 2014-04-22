PROGRAM EXECUTE
    USE nDGmod, netcdf
    IMPLICIT NONE
    	INTEGER, PARAMETER :: DOUBLE=KIND(1D0)
	INTEGER :: start_res,ntest
	REAL(KIND=DOUBLE) :: mu
	LOGICAL :: dodghybrid, doposlimit, debug, dozshulimit
	INTEGER :: nout

	start_res = 8
	mu = 0.1D0
	nout = 100

	debug = .FALSE.
	dozshulimit = .FALSE.
	dodghybrid = .FALSE.
	doposlimit = .FALSE.


	write(*,*) '======'
	write(*,*) 'TEST 1: Square Wave'
	write(*,*) '======'
	CALL test1d_nodal(1,start_res,2,3,nout,mu,debug)

	write(*,*) '======'
	write(*,*) 'TEST 2: cos'
	write(*,*) '======'
	CALL test1d_nodal(2,start_res,2,3,nout,mu,debug)

	write(*,*) '======'
	write(*,*) 'TEST 3: cos**2'
	write(*,*) '======'
	CALL test1d_nodal(3,start_res,2,3,nout,mu,debug)

	write(*,*) '======'
	write(*,*) 'TEST 4: cos**4'
	write(*,*) '======'
	CALL test1d_nodal(4,start_res,2,3,nout,mu,debug)

!	write(*,*) '======'
!	write(*,*) 'TEST 5: 2 sine waves'
!	write(*,*) '======'
!	CALL test1d_modal(5,start_res,2,3,nout,mu,debug)

	WRITE(*,*)
	WRITE(*,*) 'PROGRAM COMPLETE!'

	CONTAINS
	SUBROUTINE test1d_nodal(ntest,nex0,nscale,nlvl,noutput,maxcfl,debug)
		IMPLICIT NONE
		! -- Inputs
		INTEGER, INTENT(IN) :: nex0,ntest,nscale,nlvl,noutput
		REAL(KIND=DOUBLE), INTENT(IN) :: maxcfl
		LOGICAL, INTENT(IN) :: debug

		! -- Local variables
	    REAL(KIND=8), DIMENSION(nlvl) :: e1, e2, ei
		REAL(KIND=8) :: cnvg1, cnvg2, cnvgi,tfinal, cons
		INTEGER :: nsteps, nxfv, nmethod_final, nmethod, imethod, nelem,nout, N, ierr
		LOGICAL :: dozshulimit
		INTEGER, DIMENSION(10) :: tmp_method
		CHARACTER(len=9) :: outdir
		CHARACTER(len=40) :: cdf_out

		REAL(KIND=8) :: xLeft,xRight,domWidth

		REAL(KIND=DOUBLE), DIMENSION(:), ALLOCATABLE :: qNodes,qWeights,ecent,nodeSpacing,tmpErr
		REAL(KIND=DOUBLE), DIMENSION(:,:), ALLOCATABLE :: lagDeriv
		REAL(KIND=DOUBLE), ALLOCATABLE, DIMENSION(:,:) :: q0,q,u
		REAL(KIND=DOUBLE), DIMENSION(:), ALLOCATABLE :: M0,Mf

		REAL(KIND=DOUBLE) :: dxel, PI, dt,t
		REAL(KIND=8) :: tmp_qmax,tmp_qmin
		REAL(KIND=4) :: t0,tf
		REAL(KIND=4), DIMENSION(2) :: tstart,tend

		INTEGER :: i,j,k,l,p ! Looping variables

        PI = DACOS(-1D0)

		if(nlvl.lt.1) STOP 'nlev should be at least 1 in test1d_modal'

		nmethod_final = 1
		tmp_method = 0
		tmp_method(1) = 1 

		! -----
		! Set domain
		! -----

		xLeft = 0D0
		xRight = 1D0
		domWidth = xRight - xLeft

		DO nmethod=1,nmethod_final
			dozshulimit = .FALSE.

			imethod = tmp_method(nmethod)

			SELECT CASE(imethod)
				CASE(1)
				  outdir = 'ndgunlim/'
          	      dozshulimit = .FALSE.
				  N = 4
				  WRITE(*,*) '1D Nodal, No limiting, N=', N
				CASE(2)
				  outdir = 'ndgzhshu/'
				  N = 4
				  dozshulimit = .TRUE.
				  WRITE(*,*) '1D Nodal, Zhang and Shu Positivity limiting, N=', N
			END SELECT
		
		ALLOCATE(qnodes(0:N),qweights(0:N),lagDeriv(0:N,0:N),nodeSpacing(0:N-1))

		! --  Compute Gauss-Lobatto quadrature nodes and weights
		CALL gllquad_nodes(N,qnodes)
		CALL gllquad_weights(N,qnodes,qweights)

        nodeSpacing = qnodes(1:N)-qnodes(0:N-1)

        ! -- Fill in matrix of first derivative of Lagrange node polynomials at quadrature nodes
        CALL Dmat(N,qnodes,lagDeriv)

		DO p=1,nlvl
			t0 = etime(tstart)

			nelem = nex0*nscale**(p-1)

			ALLOCATE(ecent(1:nelem),M0(1:nelem),Mf(1:nelem),q0(0:N,1:nelem),q(0:N,1:nelem,u(0:N,1:nelem))
			! -- Set up x-grid via elements
			dxel = domWidth/DBLE(nelem)
			ecent(1) = xLeft + dxel/2D0
			DO j=2,nelem
				ecent(j) = ecent(j-1)+dxel
			END DO

			! Initialize rhoq, u, and filenames
			CALL inits(ntest,q,u,nelem,dxel,N,ecent,qnodes,tfinal,cdf_out)
			cdf_out = outdir // cdf_out

            q0 = q

			! Set up timestep
			dxm = dxel*MINVAL(nodeSpacing)/2D0

			IF(noutput .eq. -1) THEN
				nsteps = CEILING( (tfinal/maxcfl)*(MAXVAL(DABS(u))/dxm) )
				nout = nsteps
			ELSE
				nsteps = noutput*CEILING( (tfinal/maxcfl)*(MAXVAL(DABS(u))/dxm)/DBLE(noutput) )
				nout = noutput
			ENDIF

			dt = tfinal/DBLE(nsteps)
            mu = maxval(dabs(u))*dt/dxm

			IF(p .eq. 1) THEN ! Set up netCDF file
				CALL output1d(q0,xQuad,qweights,qnodes,N,nelem,tfinal,mu,cdf_out,nout,-1)
			ENDIF
			CALL output1d(q0,xQuad,qweights,qnodes,N,nelem,0D0,mu,cdf_out,p,0) ! Set up variables for this value of p ; Write x and initial conditions

            ! Compute initial mass from ICs for conservation testing
            DO j=1,nelem
                M0(j) = 0.5D0*SUM(qWeights(:)*q0(:,j))
            ENDDO !j

			! Time integration loop
			t = 0d0
			tmp_qmax = MAXVAL(q0)
			tmp_qmin = MINVAL(q0)

            DO l=1,nsteps
                CALL nDGsweep(q,nelem,dxel,N,qnodes,qweights,u,lagDeriv,dozshulimit,dt)
                t = t + dt
            
                IF((MOD(l,nsteps/nout).eq.0).OR.(l.eq.nsteps)) THEN
                    output1d(q,xQuad,qweights,qnodes,N,nelem,t,mu,cdf_out,p,2)
                ENDIF ! output check

                tmp_qmax = MAX(MAXVAL(q),tmp_qmax)
                tmp_qmin = MIN(MINVAL(q),tmp_qmin)
            ENDDO !l
            tf = etime(tend) - t0

            ALLOCATE(tmpErr(1:nelem),STAT=ierr)

            ! Compute errors
            DO j=1,nelem
                tmpErr(j) = SUM(qweights(:)*ABS(q(:,j)-q0(:,j)))
            ENDDO !j
			e1(p) = dxel*SUM(tmpErr)/2D0

            DO j=1,nelem
                tmpErr(j) = SUM(qweights(:)*(q(:,j)-q0(:,j))**2)
            ENDDO !j
			e2(p) = SQRT(dxel*SUM(tmpErr)/2D0)

			ei(p) = MAXVAL(ABS(q(:,:)-q0(:,:)))

            ! Compute final mass for conservation testing
            DO j=1,nelem
                Mf(j) = 0.5D0*SUM(qWeights(:)*q0(:,j))
            ENDDO !j
            cons = SUM(Mf-M0)/DBLE(nelem)

			if (p.eq.1) then
			write(UNIT=6,FMT='(A107)') &
'nex      E1        E2       Einf      convergence rate  overshoot  undershoot   cons cputime time step'
			cnvg1 = 0.d0
			cnvg2 = 0.d0
			cnvgi = 0.d0
		    else
        		cnvg1 = -log(e1(p)/e1(p-1))/log(dble(nscale))
        		cnvg2 = -log(e2(p)/e2(p-1))/log(dble(nscale))
        		cnvgi = -log(ei(p)/ei(p-1))/log(dble(nscale))
		    end if
       write(*,990) nelem, e1(p), e2(p), ei(p), &
            cnvg1, cnvg2, cnvgi, &
            tmp_qmax-MAXVAL(q0), &
            MINVAL(q0)-tmp_qmin, &
            cons, tf, nsteps


			IF(p .eq. nlvl) THEN
				CALL output1d(q,xQuad,qweights,qnodes,N,nelem,t,mu,cdf_out,p,1) ! Close netCDF files
			ENDIF

			DEALLOCATE(ecent,nodeSpacing,q,u,q0,M0,Mf,tmpErr)
        ENDDO !p
    ENDDO ! nmethod

    END SUBROUTINE test1d_nodal
    
    SUBROUTINE inits(ntest,q,u,nelem,dxel,N,ecent,qnodes,tfinal,cdf_out)
        IMPLICIT NONE
        ! Inputs
        INTEGER, INTENT(IN) :: N,nelem,ntest
        REAL(KIND=8), INTENT(IN) :: dxel
        REAL(KIND=8), DIMENSION(1:nelem), INTENT(IN) :: ecent
        REAL(KIND=8), DIMENSION(0:N), INTENT(IN) :: qnodes

        ! Outputs
        REAL(KIND=8), INTENT(OUT) :: tfinal
        REAL(KIND=8), DIMENSION(0:N,1:nelem) :: q,u
		CHARACTER(LEN=40), INTENT(OUT) :: cdf_out

        ! Local Variables
        REAL(KIND=8), DIMENSION(0:N,1:nelem) :: xQuad
        REAL(KIND=8) :: PI,r
        INTEGER :: j

		PI = DACOS(-1D0)

		u = 1D0
		tfinal = 10D0

		DO j=1,nelem
				xQuad(:,j) = ecent(j)+qnodes(:)*dxel/2D0
		ENDDO !j

        SELECT CASE(ntest)
			CASE(1) ! Square wave
				cdf_out = 'dg1d_sqwave.nc'
                q = 0D0
				WHERE( xQuad .lt. 0.6D0 .and. xQuad .gt. 0.4D0)
					q = 1D0
				END WHERE

			CASE(2) ! cos
				cdf_out = 'dg1d_cos.nc'
                q = 0D0
                r = 4D0*ABS(xQuad-0.25D0)
                WHERE(r .lt. 1D0)
                    q = 0.5D0*(1D0+DCOS(PI*r))
                END WHERE 

            CASE(3) ! cos**2
				cdf_out = 'dg1d_cos2.nc'
                q = 0D0
                r = 4D0*ABS(xQuad-0.25D0)
                WHERE(r .lt. 1D0)
                    q = (0.5D0*(1D0+DCOS(PI*r)))**2
                END WHERE 

            CASE(4) ! cos**4
				cdf_out = 'dg1d_cos2.nc'
                q = 0D0
                r = 4D0*ABS(xQuad-0.25D0)
                WHERE(r .lt. 1D0)
                    q = (0.5D0*(1D0+DCOS(PI*r)))**4
                END WHERE 
            
            CASE(5) ! 2 sine waves
                cdf_out = 'dg1d_2wave.nc'
                q = 0D0
                q = 2D0 + DSIN(6D0*PI*xQuad) + DSIN(8D0*PI*xQuad)

        END SELECT !ntest
    END SUBROUTINE inits

	SUBROUTINE output1d(q,x,wghts,nodes,N,nelem,tval_in,mu,cdf_out,ilvl,stat)
		IMPLICIT NONE
		! inputs
		INTEGER, INTENT(IN) :: N,nelem,ilvl,stat
		CHARACTER(len=40), INTENT(IN) :: cdf_out
		REAL(KIND=8), INTENT(IN) :: tval_in,mu
		REAL(KIND=8), DIMENSION(0:N,1:nelem), INTENT(IN) :: x,q
		REAL(KIND=8), DIMENSION(0:N), INTENT(IN) :: wghts,nodes
		
		! outputs

		! local variables
		INTEGER :: cdfid ! ID for netCDF file
		INTEGER, PARAMETER :: NDIMS = 2
		INTEGER :: ierr
	    INTEGER :: idq,idt,idx,idweight,idnode,dimids(NDIMS),idmu
	    INTEGER :: x_dimid, t_dimid,node_dimid
		INTEGER, DIMENSION(1:NDIMS) :: start, count
		CHARACTER(len=8) :: nxname,xname,qname,muname

		REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: tmp,temp
		INTEGER :: i,j,nxout

	    SAVE cdfid, idq, t_dimid, start, count
        
        nxout = (N+1)*nelem

		IF(stat .eq. -1) THEN
			! Create netCDF file and time variables
			ierr = NF90_CREATE(TRIM(cdf_out),NF90_CLOBBER,cdfid)

			ierr = NF90_REDEF(cdfid)
			ierr = NF90_DEF_DIM(cdfid, "nt", ilvl+1, t_dimid)
			ierr = NF90_DEF_DIM(cdfid, "nnodes", N+1, node_dimid)

			ierr = NF90_DEF_VAR(cdfid, "qweights",NF90_FLOAT, node_dimid, idweight)
			ierr = NF90_DEF_VAR(cdfid, "qnodes",NF90_FLOAT, node_dimid, idnode)
			ierr = NF90_DEF_VAR(cdfid, "time", NF90_FLOAT, t_dimid,idt)

			ierr = NF90_ENDDEF(cdfid)

			! Calculate time at output levels (note ilvl=noutput)
			ALLOCATE(tmp(1:ilvl+1), STAT=ierr)
			DO i=0,ilvl
				tmp(i+1) = DBLE(i)*tval_in/DBLE(ilvl)
			ENDDO

			! Write t, nodes, and weights
			ierr = NF90_PUT_VAR(cdfid,idt,tmp)
			ierr = NF90_PUT_VAR(cdfid,idweight,wghts)
            ierr = NF90_PUT_VAR(cdfid,idnode,nodes)

			DEALLOCATE(tmp, STAT=ierr)

			RETURN

		ELSEIF(stat .eq. 0) THEN
			! Create dimensions and variables for this level of runs (ilvl = p)
			start = 1
			count = 1
           
			! Define names of variables
			WRITE(nxname,'(a2,i1)') 'nx',ilvl
			WRITE(xname, '(a1,i1)') 'x',ilvl
			WRITE(qname, '(a1,i1)') 'q',ilvl
            WRITE(muname,'(a2,i1)') 'mu',ilvl

			ierr = NF90_REDEF(cdfid)
			ierr = NF90_DEF_DIM(cdfid, TRIM(nxname), nxout, x_dimid)

			dimids(1) = x_dimid
			dimids(2) = t_dimid

			ierr = NF90_DEF_VAR(cdfid, TRIM(qname),NF90_FLOAT,dimids,idq)
			ierr = NF90_DEF_VAR(cdfid, TRIM(xname),NF90_FLOAT,x_dimid,idx)
            ierr = NF90_DEF_VAR(cdfid, TRIM(muname), NF90_FLOAT,0,idmu)

			ierr = NF90_enddef(cdfid)

			! Write x values
			ierr = NF90_PUT_VAR(cdfid, idx, x)
            ierr = NF90_PUT_VAR(cdfid,idmu,mu)

			start(2) = 1

		ELSEIF(stat .eq. 1) THEN
			ierr = NF90_CLOSE(cdfid)
			RETURN
		ENDIF

		! Write out concentration field rhoq
		count(1) = nxout
    
		ALLOCATE(temp(1:nxout), STAT=ierr)
        DO j=1,nelem
        		temp(1+(j-1)*(N+1):j*(N+1)) = q(:,j)
        ENDDO !j

		ierr = NF90_PUT_VAR(cdfid,idq,temp,start,count)

		! Increment t level 
		start(2) = start(2) + 1 

		DEALLOCATE(temp, STAT=ierr)

	END SUBROUTINE output1d
END PROGRAM EXECUTE
