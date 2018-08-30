      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
C      REAL*16
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
C
      DIMENSION SIGMA(3,3),DSIG(3,3),SIGMATRI(3,3),YD1(3,3),
	1 DEPS(3,3),DOMEGA(3,3),OMEGA(3,3),EPSID(3,3),
	2 EPSE(3,3),EPSEL(3,3),DEPSE(3,3),DEPSID(3,3),
     3 SIGENER(3,3),TEMP1(3,3),DEPSEL(3,3),EPS(3,3),
     4 OMEGA_V(6),EPSEL_V(6),EPSE_V(6),EPSID_V(6),
	5 YD1_V(6),ENERGY(7),OMEGA_O(3,3),EPSEL_O(3,3),EPSE_O(3,3),
	6 EPSID_O(3,3),YD1_O(3,3),EPS_O(3,3),SIGMA_O(3,3),
	7 ZERO_T(3,3),AMATDO(3,3,3,3),AMATDO_2(6,6),SO(3,3,3,3),
	8 AMATDOM(3,3,3,3),AMATDOM_2(6,6),SOM(3,3,3,3),
	9 SIGMA_i(3,3),OMEGA_i(3,3),EPSID_i(3,3)
	
	REAL*8  E0,ENU0,A1,A2,A3,A4,C0,C1,ALPHA
      COMMON/PARAMETERS/E0,ENU0,A1,A2,A3,A4,C0,C1,ALPHA
C
      DATA ZERO,HALF,ONE,TWO,THREE /0.0D0,0.5D0,1.0D0,2.0D0,3.0D0/
      DATA THETA /0.5D0/
C
      PARAMETER (FTOL=1.D-6,PI=3.1415926526D0,ITER=800,nwrite=1,
     1    FTOL2=1.D-8)
C
C     MATERIAL PARAMETERS-(obtain its value from material input)
C
      E0 = PROPS(1)
      ENU0 = PROPS(2)
      A1 = PROPS(3)
      A2 = PROPS(4)
      A3 = PROPS(5)
      A4 = PROPS(6)
      C0 = PROPS(7)
      C1 = PROPS(8)
      ALPHA = PROPS(9)
      OMEGA_11 = PROPS(10)
      OMEGA_22 = PROPS(11)
	  OMEGA_33 = PROPS(12)
C=============================================================== 
C      DEFAULT NTENS = 6 FOR 3D
C      NDI = 3 NUMBER OF DIRECT STRESS COMPONENS
C      NSHR= 3 NUMBER OF ENGINEERING SHEAR STRESS COMPONENTS
C ---  Get state variables from STATEV(vector) 
      DO I = 1, NTENS 
        OMEGA_V(I) = STATEV(I)
        EPSEL_V(I) = STATEV(I+NTENS)
        EPSE_V(I)  = STATEV(I+2*NTENS)
        EPSID_V(I) = STATEV(I+3*NTENS)
        YD1_V(I)   = STATEV(I+4*NTENS)
        ENERGY(I)  = STATEV(I+5*NTENS)
      END DO
	  ENERGY(7) = STATEV(1+6*NTENS)
      FD = STATEV(2+6*NTENS)
C
C --- TRANSFORM VECTOR VALUE INTO 2D TENSOR
C
      CALL MAT1_MAT2(NTENS,OMEGA_V,OMEGA_O,ONE)
      IF (TIME(2).EQ.0) THEN
C --- check the value to determine 2D/3D simulation 
C	   WRITE(7,*),'NTENS=', NTENS
C      WRITE(7,*),'NSHR=', NSHR
        OMEGA_O(1,1) = OMEGA_11
        OMEGA_O(2,2) = OMEGA_22
		OMEGA_O(3,3) = OMEGA_33
      END IF
      CALL MAT1_MAT2(NTENS,EPSEL_V,EPSEL_O,HALF)
      CALL MAT1_MAT2(NTENS,EPSE_V,EPSE_O,HALF)
      CALL MAT1_MAT2(NTENS,EPSID_V,EPSID_O,HALF)
      CALL MAT1_MAT2(NTENS,YD1_V,YD1_O,ONE)
      CALL MAT1_MAT2(NTENS,DSTRAN,DEPS,HALF)
      CALL MAT1_MAT2(NTENS,STRAN,EPS_O,HALF)
      CALL MAT1_MAT2(NTENS,STRESS,SIGMA_O,ONE)
C
C --- CALCULATION OF THE EFFECTIVE FOURTH-ORDER ELASTIC STIFFNESS TENSOR
C
      DO I=1,3
        DO J=1,3
          ZERO_T(I,J)=0.0D0
        END DO
      END DO
C
      CALL DMAT(AMATDO,SO,ZERO_T)
      CALL MAT4_MAT2(AMATDO,AMATDO_2,1)
      CALL DMAT(AMATDOM,SOM,OMEGA_O)
      CALL MAT4_MAT2(AMATDOM,AMATDOM_2,1)
C
C     ELASTICAL TRIAL
C
      CALL Aijkl_Bkl(AMATDOM,DEPS,DSIG)
C
C      WRITE(7,*),'DEPS=',DEPS
C      WRITE(7,*),'DSTRAN=',STRAN
C      WRITE(7,*),'DSIG=',DSIG
C
      CALL Aij_PLUS_Bij(SIGMA_O,DSIG,SIGMATRI)      
      CALL FDDP(FD1,SIGMATRI,OMEGA_O,YD1)
C       WRITE(7,*),'FD=',FD1
C      WRITE(7,*),'SIGMATRI=',SIGMATRI
C
C     DAMAGE OR NOT?
C
      IF (FD1.LE.FTOL) THEN
C
C       NO DAMAGE GENERATED, ELASTIC INCREMENT
C        WRITE(7,*)'sigmal=1’,
C
        DO II=1,3
          DO JJ=1,3
            SIGMA(II,JJ)=SIGMATRI(II,JJ)
            DOMEGA(II,JJ)=0.0D0
            OMEGA(II,JJ)=OMEGA_O(II,JJ)
            EPSID(II,JJ)=EPSID_O(II,JJ)
          END DO
        END DO
        CALL DMAT(AMATDOM,SOM,OMEGA)
        CALL Aijkl_Bkl(SOM,SIGMA,EPSE)
        CALL Aijkl_Bkl(SO, SIGMA,EPSEL)
        CALL FDDP(FD2,SIGMA,OMEGA,YD1)
C
        DO II=1,3
          DO JJ=1,3
            DEPSID(II,JJ)=0.0D0
            SIGENER(II,JJ)=SIGMA_O(II,JJ)+THETA*DSIG(II,JJ)
            DEPSE(II,JJ)=DEPS(II,JJ)
            DEPSEL(II,JJ)=EPSEL(II,JJ)-EPSEL_O(II,JJ)
          END DO
        END DO
C - - - ENERGY CALCULATION
        CALL Aij_Bij(SIGENER,DEPS,ENERGY(1))
        CALL Aij_Bij(SIGENER,DEPSE,ENERGY(2))
        CALL Aij_Bij(SIGENER,DEPSEL,ENERGY(3))
        CALL Aij_Bij(SIGENER,DEPSID,ENERGY(4))
        CALL Aij_Bij(YD1,DOMEGA,ENERGY(5))
        TRSIGENER=SIGENER(1,1)+SIGENER(2,2)+SIGENER(3,3)
        CALL Aik_Bkj(SIGENER,SIGENER,TEMP1)
        TRTEMP1=TEMP1(1,1)+TEMP1(2,2)+TEMP1(3,3)
        DG_DTRDO=(A1+A3)*TRSIGENER**2.+(A2+A4)*TRTEMP1
        ENERGY(6)=DG_DTRDO*TRDOMEGA
        CALL Aij_Bij(SIGMA,EPSE,TEMP2)
        ENERGY(7)=0.5*TEMP2
C
      ELSEIF (FD1.GT.FTOL) THEN
C
C     DMAGE OCCURS,IRREVERSIBLE INCREMENT
C      WRITE(7,*)'signal=2’,
C
        INC=1
        FD2=FD1
        CALL Aij_PLUS_Bij(SIGMATRI,ZERO_T,SIGMA_i)
        CALL Aij_PLUS_Bij(OMEGA_O, ZERO_T,OMEGA_i)
        CALL Aij_PLUS_Bij(EPSID_O, ZERO_T,EPSID_i)
        DO WHILE (ABS(FD2).GT. FTOL2 .AND. INC.LT.ITER)
          CALL DLAMBDA(OMEGA_i,SIGMA_i,DEPSID,DOMEGA,DSIG)
C
          CALL Aij_MINUS_Bij(SIGMA_i,DSIG,SIGMA)
          CALL Aij_PLUS_Bij(OMEGA_i,DOMEGA,OMEGA)
          CALL Aij_PLUS_Bij(EPSID_i,DEPSID,EPSID)
C
          CALL FDDP(FD2,SIGMA,OMEGA,YD1)
          IF (ABS(FD2).GT.FTOL2) THEN
            CALL Aij_PLUS_Bij(SIGMA,ZERO_T,SIGMA_i)
            CALL Aij_PLUS_Bij(OMEGA,ZERO_T,OMEGA_i)
            CALL Aij_PLUS_Bij(EPSID,ZERO_T,EPSID_i)
            INC = INC+1
c           WRITE(7,*),'SIGMA_i=',SIGMA_i
c           WRITE(7,*),'OMEGA_i=',OMEGA_i
c           WRITE(7,*),'EPSID_i=',EPSID_i
          END IF
        END DO
C
        CALL DMAT(AMATDOM,SOM,OMEGA)
        CALL Aijkl_Bkl(SOM,SIGMA,EPSE)
        CALL Aijkl_Bkl(SO, SIGMA,EPSEL)
        CALL FDDP(FD2,SIGMA,OMEGA,YD1)
C
        DO II=1,3
          DO JJ=1,3
            SIGENER(II,JJ)=SIGMA_O(II,JJ)
     1 +THETA*(SIGMA(II,JJ)-SIGMA_O(II,JJ))
            DEPSE(II,JJ)=EPSE(II,JJ)-EPSE_O(II,JJ)
            DEPSEL(II,JJ)=EPSEL(II,JJ)-EPSEL_O(II,JJ)
            DEPSID(II,JJ)=EPSID(II,JJ)-EPSID_O(II,JJ)
            DOMEGA(II,JJ)=OMEGA(II,JJ)-OMEGA_O(II,JJ)
          END DO
        END DO
        TRDOMEGA=DOMEGA(1,1)+DOMEGA(2,2)+DOMEGA(3,3)
C - - - ENERGY CALCULATION
        CALL Aij_Bij(SIGENER,DEPS,ENERGY(1))
        CALL Aij_Bij(SIGENER,DEPSE,ENERGY(2))
        CALL Aij_Bij(SIGENER,DEPSEL,ENERGY(3))
        CALL Aij_Bij(SIGENER,DEPSID,ENERGY(4))
        CALL Aij_Bij(YD1,DOMEGA,ENERGY(5))
        TRSIGENER=SIGENER(1,1)+SIGENER(2,2)+SIGENER(3,3)
        CALL Aik_Bkj(SIGENER,SIGENER,TEMP1)
        TRTEMP1=TEMP1(1,1)+TEMP1(2,2)+TEMP1(3,3)
        DG_DTRDO=(A1+A3)*TRSIGENER**2.+(A2+A4)*TRTEMP1
        ENERGY(6)=DG_DTRDO*TRDOMEGA
        CALL Aij_Bij(SIGMA,EPSE,TEMP2)
        ENERGY(7)=0.5*TEMP2  
C
C       write(7,*) 'OMEGA=',OMEGA
C	    write(7,*) 'SIGMA=',SIGMA
C	    write(7,*) 'EPSID=',EPSID
C	    write(7,*) 'DEPSE=',DEPSE
      END IF
C
C      write(7,*) 'OMEGA=',OMEGA
C	 write(7,*) 'SIGMA=',SIGMA
C	 write(7,*) 'EPSID=',EPSID
C	 write(7,*) 'DEPSE=',DEPSE
C 
C 
C --- UPDATING STATE VARIABLES
C
      CALL MAT2_MAT1(NTENS,OMEGA,OMEGA_V,ONE)
      CALL MAT2_MAT1(NTENS,EPSEL,EPSEL_V,TWO)
      CALL MAT2_MAT1(NTENS,EPSE,EPSE_V,TWO)
      CALL MAT2_MAT1(NTENS,EPSID,EPSID_V,TWO)
      CALL MAT2_MAT1(NTENS,SIGMA,STRESS,ONE)
      CALL MAT2_MAT1(NTENS,YD1,YD1_V,ONE)
C
C  STORE STATE VARIABLES
C
      DO I = 1,NTENS
        STATEV(I) = OMEGA_V(I)
        STATEV(I+NTENS) = EPSEL_V(I)
        STATEV(I+2*NTENS) = EPSE_V(I)
        STATEV(I+3*NTENS) = EPSID_V(I)
        STATEV(I+4*NTENS) = YD1_V(I)
	    STATEV(I+5*NTENS) = ENERGY(I)+STATEV(I+5*NTENS)
      END DO
	  STATEV(1+6*NTENS) = ENERGY(7)
      CALL FDDP(FD2,SIGMA,OMEGA,YD1)
      IF (FD2.GT.1.E-2) THEN
        STATEV(2+6*NTENS) = 1.E-2
      ELSE
        STATEV(2+6*NTENS) = FD2
      END IF
C
C  CREATE NEW JACOBIAN
C --- CONVERSION OF THE FORTH ORDER STIFFNESS TENSOR TO A SECOND
C     ORDER TENSOR
C
      CALL MAT4_MAT2(AMATDOM,AMATDOM_2,1)
C 	  write(7,*) 'DDSDDE=',AMATDOM_2
C
C --- UPDATING THE JACOBIAN MATRIX FOR ACCELERATING CONVERGENCE
C
      DO I = 1 , NTENS
        DO J = 1 , NTENS
          DDSDDE (I , J) = AMATDOM_2 (I , J)
        ENDDO
      ENDDO
C
C
      RETURN
      END
C
C
C ====================================================================
C                              SUBROUTINES 
C
C
      SUBROUTINE DMAT(AMATDOM,AMATS,OMEGA)
CC ====================================================================
C
C                              D M A T 
C
C ====================================================================
      INCLUDE 'ABA_PARAM.INC'
	  REAL*8  E0,ENU0,A1,A2,A3,A4,C0,C1,ALPHA
C
      DIMENSION OMEGA(3,3),AMATDOM(3,3,3,3),AMATDOM_2(6,6),
     1 E(3,3),AMATS(3,3,3,3),AMATS_2(6,6)
      COMMON/PARAMETERS/E0,ENU0,A1,A2,A3,A4,C0,C1,ALPHA  
C
      DATA ONE,TWO,HALF / 1.0D0,2.0D0,0.5D0 /
C
C      INITIALIZING TENSORS AND MATRIX
C
       DO I=1,6
         DO J=1,6
           AMATS_2(I,J)=0.0D0
           AMATDOM_2(I,J)=0.0D0
         END DO
       END DO
C
       DO I=1,3
         DO J=1,3
           DO K=1,3
             DO L=1,3
               AMATS(I,J,K,L)=0.0D0
               AMATDOM(I,J,K,L)=0.0D0
             END DO
           END DO 
         END DO
       END DO
C
       CALL DIDENTITY_2(E)
C
       TROMEGA = OMEGA(1,1) + OMEGA(2,2) + OMEGA(3,3)
       B1 = (1.0+ENU0)/E0/2.
       B2 = ENU0/E0
C
       DO I=1,3
         DO J=1,3
           DO K=1,3
             DO L=1,3
                AMATS(I,J,K,L) = B1*(E(I,K)*E(J,L)+E(I,L)*E(J,K))-
     1 B2*E(I,J)*E(K,L)+2.0*A1*TROMEGA*E(I,J)*E(K,L)+
     2 HALF*A2*(E(I,K)*OMEGA(J,L)+E(I,L)*OMEGA(J,K)+
     3 OMEGA(I,K)*E(J,L)+OMEGA(I,L)*E(J,K))+
     4 A3*(E(I,J)*OMEGA(K,L)+OMEGA(I,J)*E(K,L))+
     5 A4*TROMEGA*(E(I,K)*E(J,L)+E(I,L)*E(J,K))
             END DO
           END DO
         END DO
       END DO
C
       CALL MAT4_MAT2(AMATS,AMATS_2,2)
       CALL INVERSE(AMATS_2,6,6,AMATDOM_2)
       CALL MAT2_MAT4(AMATDOM_2,AMATDOM,1)
C
C	  
      RETURN
      END
C
      SUBROUTINE MAT1_MAT2(NTENS,VECTOR,TENSOR,FACT)
C ====================================================================
C 
C                MAT1_MAT2 : VECTOR TO TENSOR  
C
C ====================================================================
      INCLUDE 'ABA_PARAM.INC'
	  REAL*8  E0,ENU0,A1,A2,A3,A4,C0,C1,ALPHA
C
      DIMENSION VECTOR(NTENS),TENSOR(3,3)
C
      DO I=1,3
        DO J=1,3
          TENSOR(I,J)=0.0D0
        END DO
      END DO
C
      TENSOR(1 , 1) = VECTOR( 1 )
      TENSOR(2 , 2) = VECTOR( 2 )
	  TENSOR(3 , 3) = VECTOR( 3 )
      TENSOR(1 , 2) = VECTOR( 4 )*FACT
      TENSOR(2 , 1) = VECTOR( 4 )*FACT
      TENSOR(1 , 3) = VECTOR( 5 )*FACT
      TENSOR(3 , 1) = VECTOR( 5 )*FACT
      TENSOR(2 , 3) = VECTOR( 6 )*FACT
      TENSOR(3 , 2) = VECTOR( 6 )*FACT
C
      RETURN
      END
C
C
      SUBROUTINE MAT2_MAT1(NTENS,TENSOR,VECTOR,FACT)
C
C ====================================================================
C
C =================== MAT2_MAT1: TENSOR TO VECTOR=====================
C
C ====================================================================
      INCLUDE 'ABA_PARAM.INC'
	  REAL*8  E0,ENU0,A1,A2,A3,A4,C0,C1,ALPHA
C
      DIMENSION VECTOR(NTENS),TENSOR(3,3)
C
      DO I=1,NTENS
          VECTOR(I)=0.0D0
      END DO
C
      VECTOR( 1 ) = TENSOR(1 , 1)
      VECTOR( 2 ) = TENSOR(2 , 2)
      VECTOR( 3 ) = TENSOR(3 , 3)
      VECTOR( 4 ) = TENSOR(1 , 2)*FACT
	  VECTOR( 5 ) = TENSOR(1 , 3)*FACT
	  VECTOR( 6 ) = TENSOR(2 , 3)*FACT
C
      RETURN
      END
C
C
C
      SUBROUTINE MAT2_MAT4(DMATRIX,TENSOR,ICOE)
C======================================================================================
C
C                             MAT2_MAT4
C
C======================================================================================
      INCLUDE 'ABA_PARAM.INC'
	  REAL*8  E0,ENU0,A1,A2,A3,A4,C0,C1,ALPHA
C  
      DIMENSION DMATRIX(6,6),TENSOR(3,3,3,3)
C
C     INITALIZATION
C
       DO I=1,3
         DO J=1,3
           DO K=1,3
             DO L=1,3 
               TENSOR(I,J,K,L)=0.0D0
             END DO
           END DO
         END DO
       END DO
C
      IF (ICOE.EQ.1) THEN
         COE1=1.
         COE2=1.
      ELSEIF(ICOE.EQ.2) THEN
         COE1=2.
         COE2=4.
      END IF
C      
      TENSOR(1,1,1,1) = DMATRIX(1,1)
      TENSOR(1,1,2,2) = DMATRIX(1,2)
      TENSOR(1,1,3,3) = DMATRIX(1,3)
      TENSOR(1,1,1,2) = DMATRIX(1,4)/COE1
      TENSOR(1,1,2,1) = DMATRIX(1,4)/COE1
      TENSOR(1,1,2,3) = DMATRIX(1,5)/COE1
      TENSOR(1,1,3,2) = DMATRIX(1,5)/COE1
      TENSOR(1,1,1,3) = DMATRIX(1,6)/COE1
      TENSOR(1,1,3,1) = DMATRIX(1,6)/COE1
C
      TENSOR(2,2,1,1) = DMATRIX(2,1)
      TENSOR(2,2,2,2) = DMATRIX(2,2)
      TENSOR(2,2,3,3) = DMATRIX(2,3)
      TENSOR(2,2,1,2) = DMATRIX(2,4)/COE1
      TENSOR(2,2,2,1) = DMATRIX(2,4)/COE1
      TENSOR(2,2,2,3) = DMATRIX(2,5)/COE1
      TENSOR(2,2,3,2) = DMATRIX(2,5)/COE1
      TENSOR(2,2,1,3) = DMATRIX(2,6)/COE1
      TENSOR(2,2,3,1) = DMATRIX(2,6)/COE1
C
      TENSOR(3,3,1,1) = DMATRIX(3,1)
      TENSOR(3,3,2,2) = DMATRIX(3,2)
      TENSOR(3,3,3,3) = DMATRIX(3,3)
      TENSOR(3,3,1,2) = DMATRIX(3,4)/COE1
      TENSOR(3,3,2,1) = DMATRIX(3,4)/COE1
      TENSOR(3,3,2,3) = DMATRIX(3,5)/COE1
      TENSOR(3,3,3,2) = DMATRIX(3,5)/COE1
      TENSOR(3,3,1,3) = DMATRIX(3,6)/COE1
      TENSOR(3,3,3,1) = DMATRIX(3,6)/COE1
C
      TENSOR(1,2,1,1) = DMATRIX(4,1)/COE1
      TENSOR(1,2,2,2) = DMATRIX(4,2)/COE1
      TENSOR(1,2,3,3) = DMATRIX(4,3)/COE1
      TENSOR(1,2,1,2) = DMATRIX(4,4)/COE2
      TENSOR(1,2,2,1) = DMATRIX(4,4)/COE2
      TENSOR(1,2,2,3) = DMATRIX(4,5)/COE2
      TENSOR(1,2,3,2) = DMATRIX(4,5)/COE2
      TENSOR(1,2,1,3) = DMATRIX(4,6)/COE2
      TENSOR(1,2,3,1) = DMATRIX(4,6)/COE2
C
      TENSOR(2,3,1,1) = DMATRIX(5,1)/COE1
      TENSOR(2,3,2,2) = DMATRIX(5,2)/COE1
      TENSOR(2,3,3,3) = DMATRIX(5,3)/COE1
      TENSOR(2,3,1,2) = DMATRIX(5,4)/COE2
      TENSOR(2,3,2,1) = DMATRIX(5,4)/COE2
      TENSOR(2,3,2,3) = DMATRIX(5,5)/COE2
      TENSOR(2,3,3,2) = DMATRIX(5,5)/COE2
      TENSOR(2,3,1,3) = DMATRIX(5,6)/COE2
      TENSOR(2,3,3,1) = DMATRIX(5,6)/COE2
C
      TENSOR(1,3,1,1) = DMATRIX(6,1)/COE1
      TENSOR(1,3,2,2) = DMATRIX(6,2)/COE1
      TENSOR(1,3,3,3) = DMATRIX(6,3)/COE1
      TENSOR(1,3,1,2) = DMATRIX(6,4)/COE2
      TENSOR(1,3,2,1) = DMATRIX(6,4)/COE2
      TENSOR(1,3,2,3) = DMATRIX(6,5)/COE2
      TENSOR(1,3,3,2) = DMATRIX(6,5)/COE2
      TENSOR(1,3,1,3) = DMATRIX(6,6)/COE2
      TENSOR(1,3,3,1) = DMATRIX(6,6)/COE2
C      
      TENSOR(2,1,1,1) = DMATRIX(4,1)/COE1
      TENSOR(2,1,2,2) = DMATRIX(4,2)/COE1
      TENSOR(2,1,3,3) = DMATRIX(4,3)/COE1
      TENSOR(2,1,1,2) = DMATRIX(4,4)/COE2
      TENSOR(2,1,2,1) = DMATRIX(4,4)/COE2
      TENSOR(2,1,2,3) = DMATRIX(4,5)/COE2
      TENSOR(2,1,3,2) = DMATRIX(4,5)/COE2
      TENSOR(2,1,1,3) = DMATRIX(4,6)/COE2
      TENSOR(2,1,3,1) = DMATRIX(4,6)/COE2
C
      TENSOR(3,2,1,1) = DMATRIX(5,1)/COE1
      TENSOR(3,2,2,2) = DMATRIX(5,2)/COE1
      TENSOR(3,2,3,3) = DMATRIX(5,3)/COE1
      TENSOR(3,2,1,2) = DMATRIX(5,4)/COE2
      TENSOR(3,2,2,1) = DMATRIX(5,4)/COE2
      TENSOR(3,2,2,3) = DMATRIX(5,5)/COE2
      TENSOR(3,2,3,2) = DMATRIX(5,5)/COE2
      TENSOR(3,2,1,3) = DMATRIX(5,6)/COE2
      TENSOR(3,2,3,1) = DMATRIX(5,6)/COE2
C
      TENSOR(3,1,1,1) = DMATRIX(6,1)/COE1
      TENSOR(3,1,2,2) = DMATRIX(6,2)/COE1
      TENSOR(3,1,3,3) = DMATRIX(6,3)/COE1
      TENSOR(3,1,1,2) = DMATRIX(6,4)/COE2
      TENSOR(3,1,2,1) = DMATRIX(6,4)/COE2
      TENSOR(3,1,2,3) = DMATRIX(6,5)/COE2
      TENSOR(3,1,3,2) = DMATRIX(6,5)/COE2
      TENSOR(3,1,1,3) = DMATRIX(6,6)/COE2
      TENSOR(3,1,3,1) = DMATRIX(6,6)/COE2
C
C
      RETURN
      END
C
C
      SUBROUTINE MAT4_MAT2(TENSOR,DMATRIX,ICOE)
C
C ====================================================================
C                        MAT4_MAT2                                                                  
C I        THIS PROGRAM TRANSFORMS THE FOURTH ORDER COMPLIANCE       I
C I        TENSOR TO A SECOND ORDER MATRIX                           I
C I                                                                  I
C ====================================================================
C
      INCLUDE 'ABA_PARAM.INC'
	  REAL*8  E0,ENU0,A1,A2,A3,A4,C0,C1,ALPHA
C
      DIMENSION TENSOR(3,3,3,3),DMATRIX(6,6)
C
      DATA ZERO,TWO /0.0D0,2.0D0/
C
C     D2 = THE SECOND ORDER STIFFNESS MATRIX
C
      DO I=1,6
        DO J=1,6
          DMATRIX(I,J)=0.0D0
        END DO
      END DO

      IF (ICOE.EQ.1) THEN
         COE1=1.
         COE2=1.
      ELSEIF(ICOE.EQ.2) THEN
         COE1=2.
         COE2=4.
      END IF
C
      DMATRIX(1,1)=TENSOR(1,1,1,1)
      DMATRIX(1,2)=TENSOR(1,1,2,2)
      DMATRIX(1,3)=TENSOR(1,1,3,3)
      DMATRIX(1,4)=TENSOR(1,1,1,2)*COE1
      DMATRIX(1,5)=TENSOR(1,1,2,3)*COE1
      DMATRIX(1,6)=TENSOR(1,1,1,3)*COE1
C
      DMATRIX(2,1)=TENSOR(2,2,1,1)
      DMATRIX(2,2)=TENSOR(2,2,2,2)
      DMATRIX(2,3)=TENSOR(2,2,3,3)
      DMATRIX(2,4)=TENSOR(2,2,1,2)*COE1
      DMATRIX(2,5)=TENSOR(2,2,2,3)*COE1
      DMATRIX(2,6)=TENSOR(2,2,1,3)*COE1
C
      DMATRIX(3,1)=TENSOR(3,3,1,1)
      DMATRIX(3,2)=TENSOR(3,3,2,2)
      DMATRIX(3,3)=TENSOR(3,3,3,3)
      DMATRIX(3,4)=TENSOR(3,3,1,2)*COE1
      DMATRIX(3,5)=TENSOR(3,3,2,3)*COE1
      DMATRIX(3,6)=TENSOR(3,3,1,3)*COE1
C 
C - — - here, engineering shear strain is used
C
      DMATRIX(4,1)=TENSOR(1,2,1,1)*COE1
      DMATRIX(4,2)=TENSOR(1,2,2,2)*COE1
      DMATRIX(4,3)=TENSOR(1,2,3,3)*COE1
      DMATRIX(4,4)=TENSOR(1,2,1,2)*COE2
      DMATRIX(4,5)=TENSOR(1,2,2,3)*COE2
      DMATRIX(4,6)=TENSOR(1,2,1,3)*COE2
C
      DMATRIX(5,1)=TENSOR(2,3,1,1)*COE1
      DMATRIX(5,2)=TENSOR(2,3,2,2)*COE1
      DMATRIX(5,3)=TENSOR(2,3,3,3)*COE1
      DMATRIX(5,4)=TENSOR(2,3,1,2)*COE2
      DMATRIX(5,5)=TENSOR(2,3,2,3)*COE2
      DMATRIX(5,6)=TENSOR(2,3,1,3)*COE2
C
      DMATRIX(6,1)=TENSOR(1,3,1,1)*COE1
      DMATRIX(6,2)=TENSOR(1,3,2,2)*COE1
      DMATRIX(6,3)=TENSOR(1,3,3,3)*COE1
      DMATRIX(6,4)=TENSOR(1,3,1,2)*COE2
      DMATRIX(6,5)=TENSOR(1,3,2,3)*COE2
      DMATRIX(6,6)=TENSOR(1,3,1,3)*COE2
C
C
      RETURN
      END
C
C
C
      SUBROUTINE INVERSE(A,N,NP,AINV)
C========================================================================
C
C    CALCULATE THE SECOND ORDER TENSOR A'S INVERSE, AINV
C    A^{-1} = AINV    
C    this subroutine inverses a (n x n) A matrix
C	 following a Gauss-Jordan elimination process
C
C========================================================================
      INCLUDE 'ABA_PARAM.INC'
	  REAL*8  E0,ENU0,A1,A2,A3,A4,C0,C1,ALPHA
C  
      DIMENSION A(NP,NP),IPIV(NP),INDXR(NP),INDXC(NP),
     1 A0(NP,NP),AINV(NP,NP)
C
      DO J=1,N
        IPIV(J)=0
      END DO
C
C     storage of the original A matrix
C
      DO I=1,N
        DO J=1,N
          A0(I,J)=A(I,J)
        END DO
      END DO
C
C	find a pivot (largest absolute value) among the rows of A that have not already been reduced
C
      DO I=1,N
        BIG=0.0D0
        DO J=1,N
          IF(IPIV(J).NE.1)THEN
            DO K=1,N
                IF(IPIV(K).EQ.0)THEN
                  IF(ABS(A(J,K)).GE.BIG)THEN
                    BIG=ABS(A(J,K))
                    IROW=J
                    ICOL=K
                    PIV=A(J,K)
                  END IF
                ELSEIF(IPIV(K).GT.1)THEN
                  write (7,*) 'Singular Matrix'
                END IF
            END DO
          END IF
        END DO
C
        IPIV(ICOL)=IPIV(ICOL)+1
        INDXR(I)=IROW
        INDXC(I)=ICOL
C	  
C     interchange the rows to put the pivot on the diagonal
C
        IF(IROW.NE.ICOL)THEN
          DO L=1,N
            DUM=A(IROW,L)
            A(IROW,L)=A(ICOL,L)
            A(ICOL,L)=DUM
          END DO
        END IF
C
C     reduction of the row of the pivot
C       
        IF(PIV.EQ.0) write (7,*) 'Singular Matrix2'
C       
        PIVIN=1./PIV          ! numerical stabilization
C
        A(ICOL,ICOL)=1.       ! numerical stabilization
        DO L=1,N
           A(ICOL,L)=A(ICOL,L)*PIVIN
        END DO
C
C     reduction of the column of the pivot
C
        DO LL=1,N
           IF(LL.NE.ICOL)THEN
             DUM=A(LL,ICOL)
             A(LL,ICOL)=0.    ! numerical stabilization
             DO L=1,N
                A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
             END DO
           END IF
        END DO
      END DO
C
C     unscramble the columns to get A^{-1}
C		
      DO J=N,1,-1   ! reverse DO loop
        DO K=1,N
          DUM=A(K,INDXR(J))
          A(K,INDXR(J))=A(K,INDXC(J))
          A(K,INDXC(J))=DUM
        END DO
      END DO
C
C	restitution process of A and Ainv
C
      DO I=1,N
        DO J=1,N
          AINV(I,J)=A(I,J)
        END DO
      END DO
C
      DO I=1,N
        DO J=1,N
          A(I,J)=A0(I,J)
        END DO
      END DO
C     
      RETURN
      END
C
C
      SUBROUTINE DIDENTITY_2(DELTA)
C========================================================================
C                                                                       =
C                               DIDENTITY_2                             =
C                                                                       =
C========================================================================
      INCLUDE 'ABA_PARAM.INC'
	  REAL*8  E0,ENU0,A1,A2,A3,A4,C0,C1,ALPHA
C
      DIMENSION DELTA(3,3)
      DATA ZERO,ONE /0.0D0,1.0D0/
C         
      DO I=1,3
        DO J=1,3
          DELTA(I,J)=ZERO
        END DO
		DELTA(I,I)=ONE
      END DO
C
      RETURN
      END
C
      SUBROUTINE Aijkl_Bkl(A,B,C)
C========================================================================
C                                                                       =
C                              Aijkl_Bkl                                =
C                                                                       =
C========================================================================
      INCLUDE 'ABA_PARAM.INC'
	  REAL*8  E0,ENU0,A1,A2,A3,A4,C0,C1,ALPHA
C
      DIMENSION A(3,3,3,3),B(3,3),C(3,3)
      DATA ZERO /0.0D0/
C
      DO I=1,3
        DO J=1,3
          C(I,J)=ZERO
          DO K=1,3
            DO L=1,3
              C(I,J)=C(I,J)+A(I,J,K,L)*B(K,L)
            END DO
          END DO
        END DO
      END DO
C
      RETURN
      END
C
      SUBROUTINE FDDP(FD,SIGMA,OMEGA,YD1)
C=========================================================================
C
C       Damage Yield Function (Pseudo-Drucker-Pager)
C
C=========================================================================
      INCLUDE 'ABA_PARAM.INC'
	  REAL*8  E0,ENU0,A1,A2,A3,A4,C0,C1,ALPHA
C
      DIMENSION OMEGA(3,3),SIGMA(3,3),YD1(3,3),E(3,3),
     1 SIGSIG(3,3),P_1(3,3,3,3),P_1Y(3,3),SIJ(3,3)
C
      COMMON/PARAMETERS/E0,ENU0,A1,A2,A3,A4,C0,C1,ALPHA 
C
      DATA ONE,TWO,HALF / 1.0D0,2.0D0,0.5D0 /
C
      CALL DIDENTITY_2(E)
C
      TRSIGMA = SIGMA(1,1)+SIGMA(2,2)+SIGMA(3,3)
      TROMEGA = OMEGA(1,1)+OMEGA(2,2)+OMEGA(3,3)
C
      CALL Aik_Bkj(SIGMA,SIGMA,SIGSIG)
C
      TRSIGSIG =SIGSIG(1,1)+SIGSIG(2,2)+SIGSIG(3,3)
C
      DO I=1,3
        DO J=1,3
          YD1(I,J)= A1*TRSIGMA**2.*E(I,J)+A2*SIGSIG(I,J)+
     1 A3*TRSIGMA*SIGMA(I,J)+A4*TRSIGSIG*E(I,J)
        END DO
      END DO
C  
C
      CALL MATP_1(P_1,SIGMA)
C
      CALL Aijkl_Bkl(P_1,YD1,P_1Y)
C
      TRY = P_1Y(1,1)+P_1Y(2,2)+P_1Y(3,3)
C
      DO I= 1,3
        DO J=1,3
          SIJ(I,J) = P_1Y(I,J)-1.D0/3.D0*TRY*E(I,J)
        END DO
      END DO
C
      CALL Aij_Bij(SIJ,SIJ,SS)
      FD = SQRT(HALF*SS)+ALPHA*TRY-C0-C1*TROMEGA   
C— - -THE SIGN BEFORE ALPHA IS "+" DUE TO MECHANICAL CONVENTION
C
      RETURN
      END
C
C
      SUBROUTINE MATP_1(P_1,SIGMA)
C====================================================================================================
C
C                          MATP_1	  
C
C=====================================================================================================
      INCLUDE 'ABA_PARAM.INC'
	  REAL*8  E0,ENU0,A1,A2,A3,A4,C0,C1,ALPHA
C
      DIMENSION P_1(3,3,3,3),SIGMA(3,3),AN1(3),AN2(3),AN3(3),S(3),IV1(3)
      REAL*8 EIGVALR(3),EIGVALI(3),EIGVEC(3,3),FV1(3),SIGMA0(3,3)
C
      PARAMETER (FTOL=1.D-6,PI=3.1415926526D0)
C
       DO I=1,3
         AN1(I)=0.
         AN2(I)=0.
         AN3(I)=0. 
         S(I)=0.
         EIGVALR(I)=0.
         EIGVALI(I)=0.
         FV1(I)=0.
         IV1(I)=0.
         DO J=1,3
           EIGVEC(I,J)=0.
           SIGMA0(I,J)=SIGMA(I,J)
           DO K=1,3
             DO L=1,3 
               P_1(I,J,K,L)=0.
             END DO
           END DO 
         END DO
       END DO 
C
C
      CALL GEIGEN(3,SIGMA0,EIGVALR,EIGVALI,EIGVEC,FV1,IV1)
C 
C
      DO I=1,3
        IF (ABS(EIGVALR(I)).LT.FTOL) EIGVALR(I)=0.D0
      END DO
C
      DO I=1,3
        AN1(I)=EIGVEC(I,1)
        AN2(I)=EIGVEC(I,2)
        AN3(I)=EIGVEC(I,3)
      END DO
C	  
C	  write(7,*) 'STRESS=',SIGMA
C	  write(7,*) 'EIGVALUE=',EIGVALR
C	  write(7,*) 'AN1=',AN1
C	  write(7,*) 'AN2=',AN2
C	  write(7,*) 'AN3=',AN3
C
      DO I=1,3 
        IF (EIGVALR(I).GE.0.) THEN
          S(I)=1.
        ELSE
          S(I)=-1.
        END IF
      END DO
C	  
      
       DO I=1,3
         DO J=1,3
           DO K=1,3
             DO L=1,3 
               P_1(I,J,K,L)=S(1)*AN1(I)*AN1(J)*AN1(K)*AN1(L)+
     1 S(2)*AN2(I)*AN2(J)*AN2(K)*AN2(L)+
     2 S(3)*AN3(I)*AN3(J)*AN3(K)*AN3(L)
             END DO
           END DO 
         END DO
       END DO
      
      RETURN
      END
  
C      SUBROUTINE MATP_1(P_1,SIGMA)
C====================================================================================================
C
C                          MATP_1	(using built in Fortran code)  
C
C=====================================================================================================
C     INCLUDE 'ABA_PARAM.INC'
C
C      DIMENSION P_1(3,3,3,3),SIGMA(3,3),AN1(3),AN2(3),AN3(3)
C      REAL*8 EIGVALR(3),EIGVEC(3,3)
C
C      PARAMETER (FTOL=1.D-6,PI=3.1415926526D0,LSTR=1,NDI=3,NSHR=1)
C*************
C      CALL SPRIND(SIGMA,EIGVALR,EIGVEC,LSTR,NDI,NSHR)
C************* 
C
C      write(7,*) 'EIGVALR=',EIGVALR
C	  write(7,*) 'EIGVEC=',EIGVEC
C	  
C      DO I=1,3
C        IF (ABS(EIGVALR(I)).LT.FTOL) EIGVALR(I)=0.D0
C      END DO
C
C      DO I=1,3
C        AN1(I)=EIGVEC(I,1)
C        AN2(I)=EIGVEC(I,2)
C        AN3(I)=EIGVEC(I,3)
C      END DO
C
C      DO I=1,3 
C        IF (EIGVALR(I).GE.0.) THEN
C          S(I)=1.
C        ELSE
C          S(I)=-1.
C        END IF
C      END DO
C	  
C       DO I=1,3
C         DO J=1,3
C           DO K=1,3
C             DO L=1,3 
C               P_1(I,J,K,L)=S(1)*AN1(I)*AN1(J)*AN1(K)*AN1(L)+
C     1    S(2)*AN2(I)*AN2(J)*AN2(K)*AN2(L)+
C     2    S(3)*AN3(I)*AN3(J)*AN3(K)*AN3(L)
C             END DO
C           END DO 
C         END DO
C       END DO
C      
C     RETURN
C      END
C
      SUBROUTINE Aij_Bij(A,B,C)
C====================================================================================
C                                                                                   =
C                                   Aij_Bij                                         =
C                                                                                   =
C====================================================================================
       INCLUDE 'ABA_PARAM.INC'
	   REAL*8  E0,ENU0,A1,A2,A3,A4,C0,C1,ALPHA
C
      DIMENSION A(3,3),B(3,3)
      DATA ZERO /0.0D0/
C
      C=ZERO
      DO I=1,3
        DO J=1,3
          C=C+A(I,J)*B(I,J)
        END DO
      END DO
C
      RETURN
      END
C
      SUBROUTINE Aij_MINUS_Bij(A,B,C)
C====================================================================================
C                                                                      *
C                                  Aij_MINUS_Bij                        *
C                                                                      *
C====================================================================================
      INCLUDE 'ABA_PARAM.INC'
	  REAL*8  E0,ENU0,A1,A2,A3,A4,C0,C1,ALPHA
C
      DIMENSION A(3,3),B(3,3),C(3,3)
C
      DO I=1,3
        DO J=1,3
          C(I,J)=0
          C(I,J)=A(I,J)-B(I,J)
        END DO
      END DO
C
      RETURN
      END
C
      SUBROUTINE Aij_PLUS_Bij(A,B,C)
C====================================================================================
C                                                                      *
C                                  Aij_PLUS_Bij                        *
C                                                                      *
C====================================================================================
      INCLUDE 'ABA_PARAM.INC'
	  REAL*8  E0,ENU0,A1,A2,A3,A4,C0,C1,ALPHA
C
      DIMENSION A(3,3),B(3,3),C(3,3)
C
      DO I=1,3
        DO J=1,3
          C(I,J)=0
          C(I,J)=A(I,J)+B(I,J)
        END DO
      END DO
C
      RETURN
      END
C
      SUBROUTINE Aik_Bkj(A,B,C)
CC====================================================================================================
C
C                          Aik_Bkj=Cij	  
C
C=====================================================================================================
C
      INCLUDE 'ABA_PARAM.INC'
	  REAL*8  E0,ENU0,A1,A2,A3,A4,C0,C1,ALPHA
C
      DIMENSION A(3,3),B(3,3),C(3,3)
C
C
      DO I=1,3
        DO J=1,3
          C(I,J)=0.0D0
          DO K=1,3
            C(I,J)=C(I,J)+A(I,K)*B(K,J)
          END DO
        END DO
      END DO
C
      RETURN
      END
C
      SUBROUTINE Aijmn_Bmnkl(A,B,C)
CC====================================================================================================
C
C                          Aijmnk_Bmnkl=Cijkl	  
C
C=====================================================================================================
C
      INCLUDE 'ABA_PARAM.INC'
	  REAL*8  E0,ENU0,A1,A2,A3,A4,C0,C1,ALPHA
C
      DIMENSION A(3,3,3,3),B(3,3,3,3),C(3,3,3,3)
C
C
      DO I=1,3
        DO J=1,3
          DO K=1,3
            DO L=1,3
              C(I,J,K,L)=0.0D0
              DO M=1,3
                DO N=1,3
                  C(I,J,K,L)=C(I,J,K,L)+A(I,J,M,N)*B(M,N,K,L)
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
C
      RETURN
      END
C
      SUBROUTINE Aijklmn_Bmn(A,B,C)
C==============================================================================
C
C                      Aijklmn_Bmn=Cijkl
C
C==============================================================================
C
      INCLUDE 'ABA_PARAM.INC'
	  REAL*8  E0,ENU0,A1,A2,A3,A4,C0,C1,ALPHA
C
      DIMENSION A(3,3,3,3,3,3),B(3,3),C(3,3,3,3)
C
C
      DO I=1,3
        DO J=1,3
          DO K=1,3
            DO L=1,3
              C(I,J,K,L)=0.0D0
              DO M=1,3
                DO N=1,3
                  C(I,J,K,L)=C(I,J,K,L)+A(I,J,K,L,M,N)*B(M,N)
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO

      RETURN
      END
C
      SUBROUTINE INVMAT4(A,AINV)
C====================================================================================================
C
C                          INVERSE OF 4TH ORDER TENSOR	  
C
C=====================================================================================================
C
      INCLUDE 'ABA_PARAM.INC'
	  REAL*8  E0,ENU0,A1,A2,A3,A4,C0,C1,ALPHA
C
      DIMENSION E(6,6),A(3,3,3,3),AINV(3,3,3,3),C(3,3,3,3),A_2(6,6),
     1 AINV_2(6,6)
C
      DO I=1,6
        DO J=1,6
           IF (I.EQ.J) THEN
              E(I,J)=1.
           ELSE
              E(I,J)=0.
           END IF
        END DO
      END DO
      CALL MAT4_MAT2(A,A_2,2)
      CALL INVERSE(A_2,6,6,AINV_2)
      CALL MAT2_MAT4(AINV_2,AINV,1)
      RETURN
      END
C
      SUBROUTINE MATP_2(P_2,SIGMA)
C====================================================================================================
C
C                          MATP_2	  
C
C=====================================================================================================
C
      INCLUDE 'ABA_PARAM.INC'
	  REAL*8  E0,ENU0,A1,A2,A3,A4,C0,C1,ALPHA
C
      DIMENSION P_2(3,3,3,3),SIGMA(3,3),AN1(3),AN2(3),AN3(3),S(3),IV1(3)
      REAL*8 EIGVALR(3),EIGVALI(3),EIGVEC(3,3),FV1(3),SIGMA0(3,3)
C
      PARAMETER (FTOL=1.D-6,PI=3.1415926526D0)
       DO I=1,3
         AN1(I)=0.
         AN2(I)=0.
         AN3(I)=0. 
         S(I)=0.
         EIGVALR(I)=0.
         EIGVALI(I)=0.
         FV1(I)=0.
         IV1(I)=0.
         DO J=1,3
           EIGVEC(I,J)=0.
           SIGMA0(I,J)=SIGMA(I,J)
           DO K=1,3
             DO L=1,3 
               P_2(I,J,K,L)=0.
             END DO
           END DO 
         END DO
       END DO 
C
C
      CALL GEIGEN(3,SIGMA0,EIGVALR,EIGVALI,EIGVEC,FV1,IV1)
C
C      DO I=1,3
C        IF (ABS(EIGVALR(I)).LT.(1.D-10)) EIGVALR(I)=0.D0
C      END DO
      DO I=1,3
        AN1(I)=EIGVEC(I,1)
        AN2(I)=EIGVEC(I,2)
        AN3(I)=EIGVEC(I,3)
      END DO
      DO I=1,3 
        RES=EIGVALR(I)-MIN(EIGVALR(1),EIGVALR(2),EIGVALR(3))
        IF (RES.GT.0.0D0) THEN     !CHECK
          S(I)=1.0D0       
        ELSE
          S(I)=0.0D0
        END IF
      END DO
C	   
       DO I=1,3
         DO J=1,3
           DO K=1,3
             DO L=1,3 
               P_2(I,J,K,L)=S(1)*AN1(I)*AN1(J)*AN1(K)*AN1(L)+
     1 S(2)*AN2(I)*AN2(J)*AN2(K)*AN2(L)+
     2 S(3)*AN3(I)*AN3(J)*AN3(K)*AN3(L)
             END DO
           END DO 
         END DO
       END DO
C    
      RETURN
      END
C
      SUBROUTINE Aijkl_Bij(A,B,C)
C=======================================================================
C                                                                      *
C                          Aijkl_Bij                                   *
C                                                                      *
C=======================================================================
      INCLUDE 'ABA_PARAM.INC'
	  REAL*8  E0,ENU0,A1,A2,A3,A4,C0,C1,ALPHA
C
      DIMENSION A(3,3,3,3),B(3,3),C(3,3)
      DATA ZERO /0.0D0/
C
      DO K=1,3
        DO L=1,3
          C(K,L)=ZERO
          DO I=1,3
            DO J=1,3
              C(K,L)=C(K,L)+A(I,J,K,L)*B(I,J)
            END DO
          END DO
        END DO
      END DO
C
      RETURN
      END
C
      SUBROUTINE Aij_Bkl(A,B,C)
C=======================================================================
C                                                                      *
C                           Aij_Bkl                                    *
C                                                                      *
C=======================================================================
      INCLUDE 'ABA_PARAM.INC'
	  REAL*8  E0,ENU0,A1,A2,A3,A4,C0,C1,ALPHA
C
      DIMENSION A(3,3),B(3,3),C(3,3,3,3)
C
      DO I=1,3
        DO J=1,3
          DO K=1,3
            DO L=1,3
              C(I,J,K,L)=A(I,J)*B(K,L)
            END DO
          END DO
        END DO
      END DO
C
      RETURN
      END
C
      SUBROUTINE DLAMBDA(OMEGA,SIGMA,DEPSID,DOMEGA,DSIG)
C=======================================================================
C                                                                      *
C                           DLAMBDA                                    *
C                                                                      *
C=======================================================================
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION SIGMA(3,3), OMEGA(3,3), DS_DO(3,3,3,3,3,3), DG_DY(3,3),
     1 DF_DY(3,3), DY_DSIG(3,3,3,3),DF_DSIG(3,3),DSO_DGY(3,3,3,3),
     2 SIG_DSO_DGY(3,3),BRACKET(3,3), AMATDOM(3,3,3,3), SOM(3,3,3,3),
     3 D_BRACKET(3,3),DF_DO(3,3),DEPSID(3,3),DOMEGA(3,3),DSIG(3,3),
     4 YD1(3,3),E(3,3),SIGSIG(3,3), EEP1(3,3,3,3),      
     5 P_1(3,3,3,3), P1YD1(3,3), EP1(3,3), 
     6 P_3(3,3,3,3), F1IJ(3,3), F1P3(3,3), 
     7 P_2(3,3,3,3), F2IJ(3,3), F2P2(3,3)      
C
      REAL*8  E0,ENU0,A1,A2,A3,A4,C0,C1,ALPHA
      COMMON/PARAMETERS/E0,ENU0,A1,A2,A3,A4,C0,C1,ALPHA 
C
      CALL FDDP(FD1,SIGMA,OMEGA,YD1)
C      write(7,*) 'FD1=',FD1
C
      CALL DIDENTITY_2(E)
C========================================================
C      calculate DS_DOMEGA
C========================================================
      DO I=1,3
        DO J=1,3
          DO K=1,3
            DO L=1,3
              DO M=1,3
                DO N=1,3
                  DS_DO(I,J,K,L,M,N)=2.*A1*E(I,J)*E(K,L)*E(M,N)+
     1 0.25*A2*(E(I,K)*(E(J,M)*E(L,N)+E(J,N)*E(L,M))+
     2 E(I,L)*(E(J,M)*E(K,N)+E(J,N)*E(K,M))+
     3 E(J,L)*(E(I,M)*E(K,N)+E(I,N)*E(K,M))+
     4 E(J,K)*(E(I,M)*E(L,N)+E(I,N)*E(L,M)))+
     5 0.5*A3*(E(I,J)*(E(K,M)*E(L,N)+E(K,N)*E(L,M))+
     6 E(K,L)*(E(I,M)*E(J,N)+E(I,N)*E(J,M)))+
     7 A4*(E(I,K)*E(J,L)+E(I,L)*E(J,K))*E(M,N)
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
C========================================================
C      calculate DY_DSIGMA
C========================================================     
      TRSIGMA=SIGMA(1,1)+SIGMA(2,2)+SIGMA(3,3)
      DO I=1,3
        DO J=1,3
          DO K=1,3
            DO L=1,3
                DY_DSIG(I,J,K,L)=2.*A1*TRSIGMA*E(I,J)*E(K,L)+
     1 0.5*A2*(E(I,K)*SIGMA(L,J)+E(I,L)*SIGMA(K,J)+SIGMA(I,K)*E(J,L)+
     2 SIGMA(I,L)*E(J,K))+A3*(E(K,L)*SIGMA(I,J)
     3 +0.5*TRSIGMA*(E(I,K)*E(J,L)+E(I,L)*E(J,K)))
     4 +2.*A4*SIGMA(K,L)*E(I,J)
            END DO
          END DO
        END DO
      END DO
C========================================================
C      calculate DF_DY,DG_DY,DF_DOMEGA
C========================================================  
      CALL Aik_Bkj(SIGMA,SIGMA,SIGSIG)
      TRSIGSIG =SIGSIG(1,1)+SIGSIG(2,2)+SIGSIG(3,3)
      DO I=1,3
        DO J=1,3
          YD1(I,J)= A1*TRSIGMA**2.*E(I,J)+A2*SIGSIG(I,J)+
     1 A3*TRSIGMA*SIGMA(I,J)+A4*TRSIGSIG*E(I,J)
        END DO
      END DO	  
      CALL MATP_1(P_1,SIGMA)
      CALL Aijkl_Bkl(P_1,YD1,P1YD1)
      CALL Aij_Bij(P1YD1,E,P1YD1E)
      CALL Aijkl_Bij(P_1,E,EP1)
      CALL Aij_Bkl(E,EP1,EEP1)
      DO I=1,3
        DO J=1,3
          DO K=1,3
            DO L=1,3
              P_3(I,J,K,L)=P_1(I,J,K,L)-1./3.*EEP1(I,J,K,L)
            END DO
          END DO
        END DO
      END DO
      DO I=1,3
        DO J=1,3
          F1IJ(I,J)=P1YD1(I,J)-1./3.*P1YD1E*E(I,J)
        END DO
      END DO    
      CALL Aijkl_Bij(P_3,F1IJ,F1P3)
      CALL Aij_Bij(F1IJ,F1IJ,F1F1)
      DO I=1,3
        DO J=1,3
          DF_DY(I,J)=F1P3(I,J)/SQRT(2.*F1F1)+ALPHA*EP1(I,J)
        END DO
      END DO
C
C
      CALL MATP_2(P_2,SIGMA)
      CALL Aijkl_Bkl(P_2,YD1,F2IJ)
      CALL Aijkl_Bij(P_2,F2IJ,F2P2)
      CALL Aij_Bij(F2IJ,F2IJ,F2F2)
      DO I=1,3
        DO J=1,3
          DF_DO(I,J)=-C1*E(I,J)
          IF(F2F2.EQ.0.) THEN
             DG_DY(I,J)=0.
          ELSE
             DG_DY(I,J)=F2P2(I,J)/SQRT(2.0*F2F2)
          END IF
        END DO
      END DO
      CALL Aijkl_Bij(DY_DSIG,DF_DY,DF_DSIG)
C	  
C      write(7,*) 'DG_DY=',DG_DY
C      write(7,*) 'DF_DSIG=',DF_DSIG
C      write(7,*) 'DS_DO=',DS_DO
C      write(7,*) 'DF_DO=',DF_DO
C
C========================================================
C      calculate MULTIPLIER
C======================================================== 
      CALL Aijklmn_Bmn(DS_DO,DG_DY,DSO_DGY)
C
      CALL Aijkl_Bij(DSO_DGY,SIGMA,SIG_DSO_DGY)
C
      CALL Aij_PLUS_Bij(DF_DSIG,SIG_DSO_DGY,BRACKET)
C
      CALL DMAT(AMATDOM,SOM,OMEGA)
C
      CALL Aijkl_Bkl(AMATDOM,BRACKET,D_BRACKET)
C     
C      WRITE(7,*),'D_BRACKET=',D_BRACKET
C
      DENOMINATOR=0.0
      DO I=1,3
         DO J=1,3
           DENOMINATOR=DENOMINATOR+DF_DSIG(I,J)*D_BRACKET(I,J)-
     1 DF_DO(I,J)*DG_DY(I,J)
        END DO
      END DO
C       write(7,*) 'DENOMINATOR=',DENOMINATOR   
C
      DO I=1,3
        DO J=1,3
          DEPSID(I,J)=0
          DOMEGA(I,J)=0
          DSIG(I,J)=0
        END DO
      END DO
C
      ALAM=0
      ALAM=ALAM+FD1/DENOMINATOR
C       write(7,*) 'ALAM=',ALAM
C      
      DO I=1,3
        DO J=1,3
          DEPSID(I,J)=DEPSID(I,J)+ALAM*DF_DSIG(I,J)
          DOMEGA(I,J)=DOMEGA(I,J)+ALAM*DG_DY(I,J)
          DSIG(I,J)=DSIG(I,J)+ALAM*D_BRACKET(I,J)
        END DO
      END DO
C      
C      WRITE(7,*),'DEPSID=',DEPSID
C 	  WRITE(7,*),'DOMEGA=',DOMEGA
C      WRITE(7,*) 'DSIG=',DSIG      
C
      RETURN
      END
C
      SUBROUTINE GEIGEN( N,A,WR,WI,Z,FV1,IV1)
C=======================================================================
C                                                                      *
C                           Aij_Bkl                                    *
C                                                                      *
C=======================================================================
C   25/01/77            MEMBER NAME  GEIGEN   (SO)     FORTRAN
C   Input
C	N - dimension
C	A - input array
C   Output
C	A   - array of eigenvectors
C    WR  - array of eigenvalues (real)
C	WI  - array of eigenvalues (imagin)
C	Z   - working array
C	FV1 - working array
C	IV1 - working array
C======================================================================
CCCCCCCOMMON/ KONTR/KTRW(80)
C
      INCLUDE 'ABA_PARAM.INC'
      INTEGER*4 IV1(N)
      REAL*8 A(N,N),WR(N),WI(N),Z(N,N),FV1(N)
CCC
CCC   SUBROUTINES FROM EISPECK-PACKAGE
      CALL BALANC(N,N,A, IS1,IS2,FV1)
      CALL ELMHES(N,N,IS1,IS2,A,IV1)
      CALL GEIELTRAN(N,N,IS1,IS2,A,IV1,Z)
      CALL HQR2(N,N,IS1,IS2,A,WR,WI,Z,IERR)
CG    IF(IERR.NE.0)  CALL WRITEI('IERR',1,IERR,1)
      CALL BALBAK(N,N,IS1,IS2,FV1,N,Z)
C
CCC---
CG    IF(KTRW(4).EQ.1)   CALL WRITED('WR  ',1,WR,N)
CG    IF(KTRW(4).EQ.1)   CALL WRITED('WI  ',1,WI,N)
      RETURN
      END
C
      SUBROUTINE BALANC(NM,N,A,LOW,IGH,SCALE)
C===============================================================================
C===============================================================================
      INCLUDE 'ABA_PARAM.INC'
      INTEGER I,J,K,L,M,N,JJ,NM,IGH,LOW,IEXC
      REAL*8 A(NM,N),SCALE(N)
      REAL*8 C,F,G,R,S,B2,RADIX
      REAL*8 DABS
      LOGICAL NOCONV
C
C     ***** RADIX IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C          THE BASE OF THE MACHINE FLOATING POINT REPRESENTATION
C     ******
      RADIX=16
C
      B2=RADIX*RADIX
      K=1
      L=N
      GO TO 100
C     ***** IN-LINE PROCEDURE FOR ROW AND
C          COLUMN EXCHANGE*****
20    SCALE(M)=J
      IF(J.EQ.M) GO TO 50
C
       DO 30 I=1,L
      F=A(I,J)
      A(I,J)=A(I,M)
      A(I,M)=F
30    CONTINUE
C
       DO  40 I=K,N
      F=A(J,I)
      A(J,I)=A(M,I)
      A(M,I)=F
40    CONTINUE
C
50    GO TO (80,130),IEXC
C     *****SEARCH FOR ROWS ISOLATING AN EIGENVALUE
C          AND PUSH THEM DOWN*****
80    IF(L.EQ.1) GO TO 280
      L=L-1
C     *****FOR J=L STEP -1 UNTIL 1 DO -- *****
100   DO 120 JJ=1,L
      J=L+1-JJ
C
       DO 110 I=1,L
      IF(I.EQ.J) GO TO 110
      IF(A(J,I).NE.0.0) GO TO 120
110   CONTINUE
C
      M=L
      IEXC=1
      GO TO 20
120   CONTINUE
C
      GO TO 140
C     *****SEARCH FOR COLUMS ISOLATING AN EIGENVALUE
C          AND PUSH THEM LEFT*****
130   K=K+1
C
140   DO 170 J=K,L
C
       DO 150 I=K,L
      IF(I.EQ.J) GO TO 150
      IF(A(I,J).NE.0.0) GO TO 170
150   CONTINUE
C
      M=K
      IEXC=2
      GO TO 20
170   CONTINUE
C     *****NOW BALANCE THE SUBMATRIX IN ROWS K TO L*****
       DO 180 I=K,L
180   SCALE(I)=1.0
C     *****ITERATIVE LOOP FOR NORM REDUCTION*****
190   NOCONV=.FALSE.
C
       DO 270 I=K,L
      C=0.0
      R=0.0
C
       DO 200 J=K,L
      IF(J.EQ.I) GO TO 200
      C=C+DABS(A(J,I))
      R=R+DABS(A(I,J))
200   CONTINUE
C
      G=R/RADIX
      F=1.0
      S=C+R
210   IF(C.GE.G) GO TO 220
      F=F*RADIX
      C=C*B2
      GO TO 210
220   G=R*RADIX
230   IF(C.LT.G) GO TO 240
      F=F/RADIX
      C=C/B2
      GO TO 230
C     *****NOW BALANCE*****
240   IF((C+R)/F.GE.0.95*S) GO TO 270
      G=1.0/F
      SCALE(I)=SCALE(I)*F
      NOCONV=.TRUE.
C
       DO 250 J=K,N
250   A(I,J)=A(I,J)*G
C
       DO 260 J=1,L
260   A(J,I)=A(J,I)*F
C
270   CONTINUE
C
      IF(NOCONV) GO TO 190
C
280   LOW=K
      IGH=L
      RETURN
      END
C
C   02/03/76            MEMBER NAME  ELMHES   (QUELLE)      FORTRAN
      SUBROUTINE ELMHES(NM,N,LOW,IGH,A,INT)
C=================================================================================
C=================================================================================
      INCLUDE 'ABA_PARAM.INC'
      INTEGER I,J,M,N,LA,NM,IGH,KP1,LOW,NM1,MP1
      REAL*8 A(NM,N)
      REAL*8 X,Y
      REAL*8 DABS
      INTEGER INT(IGH)
C
      LA=IGH-1
      KP1=LOW+1
      IF(LA.LT.KP1) GO TO 200
C
       DO 180 M=KP1,LA
      MM1=M-1
      X=0.0
      I=M
C
       DO 100 J=M,IGH
      IF(DABS(A(J,MM1)).LE.DABS(X)) GO TO 100
      X=A(J,MM1)
      I=J
100   CONTINUE
C
      INT(M)=I
      IF(I.EQ.M) GO TO 130
C     *****INTERCHANGE ROWS AND COLUMNS OF A*****
       DO 110 J=MM1,N
      Y=A(I,J)
      A(I,J)=A(M,J)
      A(M,J)=Y
110   CONTINUE
C
       DO 120 J=1,IGH
      Y=A(J,I)
      A(J,I)=A(J,M)
      A(J,M)=Y
120   CONTINUE
C     *****END INTERCHANGE*****
130   IF(X.EQ.0.0) GO TO 180
      MP1=M+1
C
       DO 160 I=MP1,IGH
      Y=A(I,MM1)
      IF(Y.EQ.0.0) GO TO 160
      Y=Y/X
      A(I,MM1)=Y
C
       DO 140 J=M,N
140   A(I,J)=A(I,J)-Y*A(M,J)
C
       DO 150 J=1,IGH
150   A(J,M)=A(J,M)+Y*A(J,I)
C
160   CONTINUE
C
180   CONTINUE
C
200   RETURN
      END
C   02/03/76            MEMBER NAME  ELTRAN   (QUELLE)      FORTRAN
      SUBROUTINE GEIELTRAN(NM,N,LOW,IGH,A,INT,Z)
C======================================================================================
C-======================================================================================
      INCLUDE 'ABA_PARAM.INC'
      INTEGER I,J,N,KL,NM,MP,NM1,IGH,LOW,MP1
      REAL*8 A(NM,IGH),Z(NM,N)
      INTEGER INT(IGH)
C
C     *****INITIALIZE Z TO IDENTITY MATRIX*****
       DO 80 I=1,N
C
       DO 60 J=1,N
60    Z(I,J)=0.0
C
      Z(I,I)=1.0
80    CONTINUE
C
      KL=IGH-LOW-1
      IF(KL.LT.1) GO TO 200
C     *****FOR MP=IGH-1 STEP -1 UNTIL LOW+1 DO -- *****
       DO 140 MM=1,KL
      MP=IGH-MM
      MP1=MP+1
C
       DO 100 I=MP1,IGH
100   Z(I,MP)=A(I,MP-1)
C
      I=INT(MP)
      IF(I.EQ.MP) GO TO 140
C
       DO 130 J=MP,IGH
      Z(MP,J)=Z(I,J)
      Z(I,J)=0.0
130   CONTINUE
C
      Z(I,MP)=1.0
140   CONTINUE
C
200   RETURN
      END
C   02/03/76            MEMBER NAME  HQR2     (QUELLE)      FORTRAN
      SUBROUTINE HQR2(NM,N,LOW,IGH,H,WR,WI,Z,IERR)
C=====================================================================================
C=====================================================================================
      INCLUDE 'ABA_PARAM.INC'
      INTEGER I,J,K,L,M,N,EN,II,JJ,LL,MM,NA,NM,NN,
     X        IGH,ITS,LOW,MP2,ENM2,IERR
      REAL*8 H(NM,N),WR(N),WI(N),Z(NM,N)
      REAL*8 P,Q,R,S,T,W,X,Y,RA,SA,VI,VR,ZZ,NORM,MACHEP
      REAL*8 DSQRT,DABS,DSIGN
      INTEGER MIN0
      LOGICAL NOTLAS
      COMPLEX*16 Z3
      COMPLEX*16 DCMPLX
      REAL*8 T3(2)
      EQUIVALENCE(Z3,T3(1))
C
C     *****MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C          THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
C     *****
      MACHEP= 1.0D-16
C
      IERR=0
C     *****STORE ROOTS ISOLATED BY BALANC*****
       DO 50 I=1,N
      IF(I.GE.LOW.AND.I.LE.IGH) GO TO 50
      WR(I)=H(I,I)
      WI(I)=0.0
50    CONTINUE
C
      EN=IGH
      T=0.0
C     *****SEARCH FOR NEXT EIGENVALUES*****
60    IF(EN.LT.LOW) GO TO 340
      ITS=0
      NA=EN-1
      ENM2=NA-1
C     *****LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
C          FOR L=EN STEP -1 UNTIL LOW DO -- *****
70    DO 80 LL=LOW,EN
      L=EN+LOW-LL
      IF(L.EQ.LOW) GO TO 100
      IF(DABS(H(L,L-1)).LE.MACHEP*(DABS(H(L-1,L-1))
     X    +DABS(H(L,L)))) GO TO 100
80    CONTINUE
C     *****FORM SHIFT*****
100   X=H(EN,EN)
      IF(L.EQ.EN) GO TO 270
      Y=H(NA,NA)
      W=H(EN,NA)*H(NA,EN)
      IF(L.EQ.NA) GO TO 280
      IF(ITS.EQ.30) GO TO 1000
      IF(ITS.NE.10.AND.ITS.NE.20) GO TO 130
C     *****FORM EXCEPTIONAL SHIFT*****
      T=T+X
C
       DO 120 I=LOW,EN
120   H(I,I)=H(I,I)-X
C
      S=DABS(H(EN,NA))+DABS(H(NA,ENM2))
      X=0.75*S
      Y=X
      W=-0.4375*S*S
130   ITS=ITS+1
C     LOOK FOR TWO CONSECUTIVE SMALL
C          SUB-DIAGONAL ELEMENTS.
C          FOR M=EN-2 STEP -1 UNTIL L DO -- *****
       DO 140 MM=L,ENM2
      M=ENM2+L-MM
      ZZ=H(M,M)
      R=X-ZZ
      S=Y-ZZ
      P=(R*S-W)/H(M+1,M)+H(M,M+1)
      Q=H(M+1,M+1)-ZZ-R-S
      R=H(M+2,M+1)
      S=DABS(P)+DABS(Q)+DABS(R)
      P=P/S
      Q=Q/S
      R=R/S
      IF(M.EQ.L) GO TO 150
      IF(DABS(H(M,M-1))*(DABS(Q)+DABS(R)).LE.MACHEP*DABS(P)
     X  *(DABS(H(M-1,M-1))+DABS(ZZ)+DABS(H(M+1,M+1)))) GO TO 150
140   CONTINUE
C
150   MP2=M+2
C
       DO 160 I=MP2,EN
      H(I,I-2)=0.0
      IF(I.EQ.MP2) GO TO 160
      H(I,I-3)=0.0
160   CONTINUE
C     *****DOUBLE QR STEP INVOLVING ROWS L TO EN AND
C          COLUMNS M TO EN*****
       DO 260 K=M,NA
      NOTLAS=K.NE.NA
      IF(K.EQ.M) GO TO 170
      P=H(K,K-1)
      Q=H(K+1,K-1)
      R=0.0
      IF(NOTLAS) R=H(K+2,K-1)
      X=DABS(P)+DABS(Q)+DABS(R)
      IF(X.EQ.0.0) GO TO 260
      P=P/X
      Q=Q/X
      R=R/X
170   S=DSIGN(DSQRT(P*P+Q*Q+R*R),P)
      IF(K.EQ.M) GO TO 180
      H(K,K-1)=-S*X
      GO TO 190
180   IF(L.NE.M) H(K,K-1)=-H(K,K-1)
190   P=P+S
      X=P/S
      Y=Q/S
      ZZ=R/S
      Q=Q/P
      R=R/P
C     *****ROW MODIFICATION*****
       DO 210 J=K,N
      P=H(K,J)+Q*H(K+1,J)
      IF(.NOT.NOTLAS) GO TO 200
      P=P+R*H(K+2,J)
      H(K+2,J)=H(K+2,J)-P*ZZ
200   H(K+1,J)=H(K+1,J)-P*Y
      H(K,J)=H(K,J)-P*X
210   CONTINUE
C
      J=MIN0(EN,K+3)
C     *****COLUMN MODIFICATION*****
       DO 230 I=1,J
      P=X*H(I,K)+Y*H(I,K+1)
      IF(.NOT.NOTLAS) GO TO 220
      P=P+ZZ*H(I,K+2)
      H(I,K+2)=H(I,K+2)-P*R
220   H(I,K+1)=H(I,K+1)-P*Q
      H(I,K)=H(I,K)-P
230   CONTINUE
C     *****ACCUMULATE TRANSFORMATIONS*****
       DO 250 I=LOW,IGH
      P=X*Z(I,K)+Y*Z(I,K+1)
      IF(.NOT.NOTLAS) GO TO 240
      P=P+ZZ*Z(I,K+2)
      Z(I,K+2)=Z(I,K+2)-P*R
240   Z(I,K+1)=Z(I,K+1)-P*Q
      Z(I,K)=Z(I,K)-P
250   CONTINUE
C
260   CONTINUE
C
      GO TO 70
C     *****ONE ROOT FOUND*****
270   H(EN,EN)=X+T
      WR(EN)=H(EN,EN)
      WI(EN)=0.0
      EN=NA
      GO TO 60
C     *****TWO ROOTS FOUND*****
280   P=(Y-X)/2.0
      Q=P*P+W
      ZZ=DSQRT(DABS(Q))
      H(EN,EN)=X+T
      X=H(EN,EN)
      H(NA,NA)=Y+T
      IF(Q.LT.0.0) GO TO 320
C     *****REAL PAIR*****
      ZZ=P+DSIGN(ZZ,P)
      WR(NA)=X+ZZ
      WR(EN)=WR(NA)
      IF(ZZ.NE.0.0) WR(EN)=X-W/ZZ
      WI(NA)=0.0
      WI(EN)=0.0
      X=H(EN,NA)
      R=DSQRT(X*X+ZZ*ZZ)
      P=X/R
      Q=ZZ/R
C     *****ROW MODIFICATION*****
       DO 290 J=NA,N
      ZZ=H(NA,J)
      H(NA,J)=Q*ZZ+P*H(EN,J)
      H(EN,J)=Q*H(EN,J)-P*ZZ
290   CONTINUE
C     *****COLUMN MODIFICATION*****
       DO 300 I=1,EN
      ZZ=H(I,NA)
      H(I,NA)=Q*ZZ+P*H(I,EN)
      H(I,EN)=Q*H(I,EN)-P*ZZ
300   CONTINUE
C     *****ACCUMULATE TRANSFORMATIONS*****
       DO 310 I=LOW,IGH
      ZZ=Z(I,NA)
      Z(I,NA)=Q*ZZ+P*Z(I,EN)
      Z(I,EN)=Q*Z(I,EN)-P*ZZ
310   CONTINUE
C
      GO TO 330
C     *****COMPLEX PAIR*****
320   WR(NA)=X+P
      WR(EN)=X+P
      WI(NA)=ZZ
      WI(EN)=-ZZ
330   EN=ENM2
      GO TO 60
C     *****ALL ROOTS FOUND.BACKSUBSTITUTE TO FIND
C          VECTORS OF UPPER TRIANGULAR FORM*****
340   NORM=0.0
      K=1
C
       DO 360 I=1,N
C
       DO 350 J=K,N
350   NORM=NORM+DABS(H(I,J))
C
      K=I
360   CONTINUE
C
      IF(NORM.EQ.0.0) GO TO 1001
C     *****FOR EN=N STEP -1 UNTIL 1 DO -- *****
       DO 800 NN=1,N
      EN=N+1-NN
      P=WR(EN)
      Q=WI(EN)
      NA=EN-1
      IF(Q) 710,600,800
C     *****REAL VECTOR*****
600   M=EN
      H(EN,EN)=1.0
      IF(NA.EQ.0) GO TO 800
C     *****FOR I=EN-1 STEP -1 UNTIL 1 DO --*****
       DO 700 II=1,NA
      I=EN-II
      W=H(I,I)-P
      R=H(I,EN)
      IF(M.GT.NA) GO TO 620
C
       DO 610 J=M,NA
610   R=R+H(I,J)*H(J,EN)
C
620   IF(WI(I).GE.0.0) GO TO 630
      ZZ=W
      S=R
      GO TO 700
630   M=I
      IF(WI(I).NE.0.0) GO TO 640
      T=W
      IF(W.EQ.0.0) T=MACHEP*NORM
      H(I,EN)=-R/T
      GO TO 700
C************* SOLVE REAL EQUATIONS *******************
640   X=H(I,I+1)
      Y=H(I+1,I)
      Q=(WR(I)-P)*(WR(I)-P)+WI(I)*WI(I)
      T=(X*S-ZZ*R)/Q
      H(I,EN)=T
      IF(DABS(X).LE.DABS(ZZ)) GOTO 650
      H(I+1,EN)=(-R-W*T)/X
      GOTO 700
650   H(I+1,EN)=(-S-Y*T)/ZZ
700   CONTINUE
C     *****END REAL VECTOR*****
      GO TO 800
C     *****COMPLEX VECTOR*****
710   M=NA
C     *****LAST VECTOR COMPONENT CHOSEN IMAGINARY SO THAT
C          EIGENVECTOR MATRIX IS TRIANGULAR*****
      IF(DABS(H(EN,NA)).LE.DABS(H(NA,EN))) GO TO 720
      H(NA,NA)=Q/H(EN,NA)
      H(NA,EN)=-(H(EN,EN)-P)/H(EN,NA)
      GO TO 730
720   Z3=DCMPLX(0.D+0,-H(NA,EN))/DCMPLX(H(NA,NA)-P,Q)
      H(NA,NA)=T3(1)
      H(NA,EN)=T3(2)
730   H(EN,NA)=0.0
      H(EN,EN)=1.0
      ENM2=NA-1
      IF(ENM2.EQ.0) GO TO 800
C
       DO 790 II=1,ENM2
      I=NA-II
      W=H(I,I)-P
      RA=0.0
      SA=H(I,EN)
C
       DO 760 J=M,NA
      RA=RA+H(I,J)*H(J,NA)
      SA=SA+H(I,J)*H(J,EN)
760   CONTINUE
C
      IF(WI(I).GE.0.0) GO TO 770
      ZZ=W
      R=RA
      S=SA
      GO TO 790
770   M=I
      IF(WI(I).NE.0.0) GO TO 780
      Z3=DCMPLX(-RA,-SA)/DCMPLX(W,Q)
      H(I,NA)=T3(1)
      H(I,EN)=T3(2)
      GO TO 790
C     *****SOLVE COMPLEX EQUATIONS*****
780   X=H(I,I+1)
      Y=H(I+1,I)
      VR=(WR(I)-P)*(WR(I)-P)+WI(I)*WI(I)-Q*Q
      VI=(WR(I)-P)*2.0*Q
      IF(VR.EQ.0.0.AND.VI.EQ.0.0) VR=MACHEP*NORM
     X  *(DABS(W)+DABS(Q)+DABS(X)+DABS(Y)+DABS(ZZ))
      Z3=DCMPLX(X*R-ZZ*RA+Q*SA,X*S-ZZ*SA-Q*RA)/DCMPLX(VR,VI)
      H(I,NA)=T3(1)
      H(I,EN)=T3(2)
      IF(DABS(X).LE.DABS(ZZ)+DABS(Q)) GO TO 785
      H(I+1,NA)=(-RA-W*H(I,NA)+Q*H(I,EN))/X
      H(I+1,EN)=(-SA-W*H(I,EN)-Q*H(I,NA))/X
      GO TO 790
785   Z3=DCMPLX(-R-Y*H(I,NA),-S-Y*H(I,EN))/DCMPLX(ZZ,Q)
      H(I+1,NA)=T3(1)
      H(I+1,EN)=T3(2)
790   CONTINUE
C     *****END COMPLEX VECTOR*****
800   CONTINUE
C     *****END BACK SUBSTITUTION.
C          VECTORS OF ISOLATED ROOTS*****
       DO 840 I=1,N
      IF(I.GE.LOW.AND.I.LE.IGH) GO TO 840
C
       DO 820 J=1,N
820   Z(I,J)=H(I,J)
C
840   CONTINUE
C     *****MULTIPLY BY TRANSFORMATION MATRIX TO GIVE
C          VECTORS OF ORIGINAL FULL MATRIX.
C          FOR J=N STEP -1 UNTIL LOW DO --*****
       DO 880 JJ=LOW,N
      J=N+LOW-JJ
      M=MIN0(J,IGH)
C
       DO 880 I=LOW,IGH
      ZZ=0.0
C
       DO 860 K=LOW,M
860   ZZ=ZZ+Z(I,K)*H(K,J)
C
      Z(I,J)=ZZ
880   CONTINUE
C
      GO TO 1001
C     *****SET ERROR -- NO CONVERGENCE TO AN
C          EIGENVALUE AFTER 30 ITERATIONS*****
1000  IERR=EN
1001  RETURN
      END
C   02/03/76            MEMBER NAME  BALBAK   (QUELLE)      FORTRAN
      SUBROUTINE BALBAK(NM,N,LOW,IGH,SCALE,M,Z)
C========================================================================
C========================================================================
      INCLUDE 'ABA_PARAM.INC'
      INTEGER I,J,K,M,N,II,NM,IGH,LOW
      REAL*8 SCALE(N),Z(NM,M)
      REAL*8 S,SUM
C
      IF(IGH.EQ.LOW) GO TO 120
C
       DO 110 I=LOW,IGH
      S=SCALE(I)
C     *****LEFT HAND EIGENVECTORS ARE BACK TRANSFORMED
C          IF THE FOREGOING STATEMENT IS REPLACED BY
C          S=1.0/SCALE(I).*****
       DO 100 J=1,M
100   Z(I,J)=Z(I,J)*S
C
110   CONTINUE
C     *****FOR I=LOW-1 STEP -1 UNTIL 1,
C          IGH+1 STEP 1 UNTIL N DO -- *****
120   DO 140 II=1,N
      I=II
      IF(I.GE.LOW.AND.I.LE.IGH) GO TO 140
      IF(I.LT.LOW) I=LOW-II
      K=SCALE(I)
      IF(K.EQ.I) GO TO 140
C
       DO 130 J=1,M
      S=Z(I,J)
      Z(I,J)=Z(K,J)
      Z(K,J)=S
130   CONTINUE
C
140   CONTINUE
C
       DO  170 J=1,M
      SUM=0.0
       DO  150 I=1,N
      SUM=SUM+Z(I,J)*Z(I,J)
150   CONTINUE
      SUM=DSQRT(SUM)
       DO  160 I=1,N
      Z(I,J)=Z(I,J)/SUM
160   CONTINUE
170   CONTINUE
C
      RETURN
      END
C