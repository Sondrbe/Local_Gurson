      subroutine vumat(
! Read only -
     .     nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     .     stepTime, totalTime, dt, cmname, coordMp, charLength,
     .     props, density, strainInc, relSpinInc,
     .     tempOld, stretchOld, defgradOld, fieldOld,
     .     stressOld, stateOld, enerInternOld, enerInelasOld,
     .     tempNew, stretchNew, defgradNew, fieldNew,
! Write only -
     .     stressNew, stateNew, enerInternNew, enerInelasNew)
     
      include 'vaba_param.inc'
      
      DIMENSION PROPS(NPROPS), DENSITY(NBLOCK), COORDMP(NBLOCK),
     . CHARLENGTH(NBLOCK), STRAININC(NBLOCK, NDIR+NSHR),
     . RELSPININC(NBLOCK, NSHR), TEMPOLD(NBLOCK),
     . STRETCHOLD(NBLOCK, NDIR+NSHR),DEFGRADOLD(NBLOCK,NDIR+NSHR+NSHR),
     . FIELDOLD(NBLOCK, NFIELDV), STRESSOLD(NBLOCK, NDIR+NSHR),
     . STATEOLD(NBLOCK, NSTATEV), ENERINTERNOLD(NBLOCK),
     . ENERINELASOLD(NBLOCK), TEMPNEW(NBLOCK),
     . STRETCHNEW(NBLOCK, NDIR+NSHR),DEFGRADNEW(NBLOCK,NDIR+NSHR+NSHR),
     . FIELDNEW(NBLOCK, NFIELDV), STRESSNEW(NBLOCK,NDIR+NSHR),
     . STATENEW(NBLOCK, NSTATEV), ENERINTERNNEW(NBLOCK),
     . ENERINELASNEW(NBLOCK)
!
      CHARACTER*8 CMNAME
!
!-----------------------------------------------------------------------
!     Internal VUMAT variables
!-----------------------------------------------------------------------
      REAL*8, PARAMETER :: PI = 3.1415927
      parameter ( zero = 0.d0, one = 1.d0, two = 2.d0,
     .   third = 1.d0 / 3.d0, half = 0.5d0, op5 = 1.5d0)
      REAL*8 SigDevOld(6)           ! Old deviatoric stress
      REAL*8 SigDevNew(6)           ! New deviatoric stress
      REAL*8 kappa, Eff_Void
      REAL*8 stresstr(6)
!
!
!-----------------------------------------------------------------------
!     Read parameters from ABAQUS material card
!-----------------------------------------------------------------------
      YOUNG  = props(1)      ! Elastic/Youngs Modulus (e.g: E=210MPa for steel)
      POISS  = props(2)      ! Poisson's ratio (often: \mu=0.33)
      SIGMA0 = props(3)      ! Initial Yield stress (at EpsPEqv=0)
      BETA1  = props(4)      ! Modified Gursons 1st constant: q_1
      BETA2  = props(5)      ! Modified Gursons 2nd constant: q_2
      BETA3  = props(6)      ! Modified Gursons 3rd constant: q_3
      TOL    = props(8)      ! TOLerance on the plastic return mapping
      MXITER = props(9)      ! Maximum number of iterations in the return map
      T_1     = props(10)    ! 
      Q_1     = props(11)    ! 3 term Voce hardening rule...
      T_2     = props(12)    ! T_i * (1 - EXP(-Q_i*EpsPEqv))
      Q_2     = props(13)    ! ... "" ...
      T_3     = props(14)
      Q_3     = props(15)
      FN     = props(16)     ! f_N, the volume fraction of void nucleating particles.
      SN     = props(17)     ! S_N, i.e the standard deviation in the nucleation distribution.
      EPSN   = props(18)     ! \epsilon_N, i.e mean plastic strain for the distribution.
      f_c    = PROPS(19)     ! Void volume fraction upon which coalescence occurs.
      f_f    = PROPS(20)     ! Critical void volume fraction upon which failure occurs.
      k_s    = PROPS(21)     ! Shear damage term, k_s=0 means no contribution
!      
!      
!-----------------------------------------------------------------------
!     Compute elasticity matrix
!-----------------------------------------------------------------------       
      G      = 0.5*YOUNG/(1.0+POISS)
      K      = YOUNG/(3.0*(1.0-2.0*POISS))


! Elastic constants:
      twomu  = YOUNG / ( one + POISS )
      alamda = POISS * twomu / ( one - two * POISS )
      thremu = op5 * twomu
!
!---- Calculate the preprocessor step --------------------------------------------                  
      IF ((totalTime.eq.0.0).and.(stepTime.eq.0.0)) THEN
        DO km = 1, NBLOCK
          trace = strainInc(km,1) + strainInc(km,2) + strainInc(km,3)
          stressNew(km,1) = stressOld(km,1)
     .          + twomu * strainInc(km,1) + alamda * trace
          stressNew(km,2) = stressOld(km,2)
     .          + twomu * strainInc(km,2) + alamda * trace
          stressNew(km,3) = stressOld(km,3)
     .          + twomu * strainInc(km,3) + alamda * trace
          stressNew(km,4)=stressOld(km,4) + twomu * strainInc(km,4)
          IF ( nshr .gt. 1 ) THEN
            stressNew(km,5)=stressOld(km,5) + twomu * strainInc(km,5)
            stressNew(km,6)=stressOld(km,6) + twomu * strainInc(km,6)
          END IF
        END DO
      ELSE
!      
!--------------------------------------------------------------------------------------
!     Continue with the return mapping algorithm:
!--------------------------------------------------------------------------------------      
      DO km = 1, nblock
        trace = strainInc(km,1) + strainInc(km,2) + strainInc(km,3)
        stresstr(1) = stressOld(km,1)
     .        + twomu * strainInc(km,1) + alamda * trace
        stresstr(2) = stressOld(km,2)
     .        + twomu * strainInc(km,2) + alamda * trace
        stresstr(3) = stressOld(km,3)
     .        + twomu * strainInc(km,3) + alamda * trace
        stresstr(4)=stressOld(km,4) + twomu * strainInc(km,4)
        IF ( nshr .gt. 1 ) THEN
          stresstr(5)=stressOld(km,5) + twomu * strainInc(km,5)
          stresstr(6)=stressOld(km,6) + twomu * strainInc(km,6)
        END IF      
!      
!-----------------------------------------------------------------------
!       Read the void content either from material card or memory:
!-----------------------------------------------------------------------
        IF (totalTime.eq.dt) THEN !read the void content from material properties, as it is the first increment in the simulation
          VoidTr = PROPS(7)
        ELSE
          VoidTr = STATEOLD(km, 1)
        END IF
!        
!-----------------------------------------------------------------------
!       Calculate pTr
!-----------------------------------------------------------------------
        SigHNew = (stresstr(1)+stresstr(2)+stresstr(3))/3.0 
        pTr = -SigHNew  
!    
!-----------------------------------------------------------------------
!       Calculate qTr
!-----------------------------------------------------------------------    
        SigDevNew(1) = stresstr(1)-SigHNew
        SigDevNew(2) = stresstr(2)-SigHNew
        SigDevNew(3) = stresstr(3)-SigHNew
        SigDevNew(4) = stresstr(4)
        IF ( nshr .gt. 1 ) THEN
          SigDevNew(5) = stresstr(5)
          SigDevNew(6) = stresstr(6)
        END IF      
        IF ( nshr .eq. 1 ) THEN
          SigEqv = sqrt((3.0/2.0)*(SigDevNew(1)*SigDevNew(1)
     .           +SigDevNew(2)*SigDevNew(2)
     .           +SigDevNew(3)*SigDevNew(3))
     .       +3.0*(SigDevNew(4)*SigDevNew(4)))
        ELSE 
          SigEqv = sqrt((3.0/2.0)*(SigDevNew(1)*SigDevNew(1)
     .           +SigDevNew(2)*SigDevNew(2)
     .           +SigDevNew(3)*SigDevNew(3))
     .      +3.0*(SigDevNew(4)*SigDevNew(4)
     .           +SigDevNew(5)*SigDevNew(5)
     .           +SigDevNew(6)*SigDevNew(6)))        
        END IF
        qTr = SigEqv  
!
!-----------------------------------------------------------------------
!       Third deviatoric stress invariant
!-----------------------------------------------------------------------        
        IF ( nshr .eq. 1 ) THEN
          J_3 = SigDevNew(1)*SigDevNew(2)*SigDevNew(3) 
     .        - SigDevNew(4)*SigDevNew(4)*SigDevNew(3)
        ELSE IF (nshr .gt. 1) THEN 
          J_3 = SigDevNew(1)*SigDevNew(2)*SigDevNew(3) 
     .        - SigDevNew(1)*SigDevNew(5)*SigDevNew(5) 
     .        - SigDevNew(4)*SigDevNew(4)*SigDevNew(3)
     .        + 2*SigDevNew(4)*SigDevNew(5)*SigDevNew(6) 
     .        - SigDevNew(2)*SigDevNew(6)*SigDevNew(6)
        END IF        
!        
!-----------------------------------------------------------------------
!       State variable from previous time step
!-----------------------------------------------------------------------
        EpsPEqvTr = STATEOLD(km,2)
!        
!-----------------------------------------------------------------------
!       Update yield stress
!-----------------------------------------------------------------------
        CALL VOCE(SigM, SIGMA0, EpsPEqvTr, T_1, Q_1, T_2, 
     .          Q_2, T_3, Q_3)
!        
!-----------------------------------------------------------------------
!       Calculate the effective void volume fraction f^*
!-----------------------------------------------------------------------
        IF ((VoidTr .ge. f_c).and.(f_f.ne.0)) THEN
          C_Factor = ((1.0/BETA1) - f_c) / (f_f - f_c)
          Eff_Void = f_c + C_Factor * (VoidTr - f_c)
        ELSE 
          Eff_Void = VoidTr
        END IF         
!
!---------- In order to avoid numerical errors -------------------------              
        IF (Eff_Void .ge. 0.99*(1.0/BETA1)) THEN
          Eff_Void = 0.99*(1.0/BETA1)
        END IF   
!        
!-----------------------------------------------------------------------
!       Compute Bracket Term
!-----------------------------------------------------------------------
        TH=(-3.0*BETA2*pTr)/(2.0*SigM)
!        
!-----------------------------------------------------------------------
!     Compute Yield Function
!-----------------------------------------------------------------------
        F = (qTr**2.0/SigM**2.0)+(2.0*Eff_Void*BETA1*cosh(TH))
     .    -1.0-BETA3*Eff_Void**2.0
!
!-----------------------------------------------------------------------
!       Check for plasticity
!-----------------------------------------------------------------------
        IF(F.GT.0.0)THEN ! Plastic flow
!-----------------------------------------------------------------------
!         Initialize Variables
!-----------------------------------------------------------------------
          p = pTr
          q = qTr
          Void = VoidTr
          EpsPEqv = EpsPEqvTr
          DeltaEpsP = 0.
          DeltaEpsQ = 0.
          omega_shear  = 1.0 - (27*J_3 / (2.0*qTr**3))**2
          IF(q.LE.1.0E-20)THEN
            q=1.0E-20
          ENDIF
!
!-----------------------------------------------------------------------
!         Start iteration
!-----------------------------------------------------------------------
          DO ITER=1,MXITER
            IF(q.LE.1.0E-20)THEN
              q=1.0E-20
            ENDIF
!-----------------------------------------------------------------------
!           Calculate Derivatives
!-----------------------------------------------------------------------
            SigMDerivX2 = T_1*Q_1*EXP(-Q_1*EpsPEqv)                     ! dSigM / deps_eq_p
     .                  + T_2*Q_2*EXP(-Q_2*EpsPEqv)
     .                  + T_3*Q_3*EXP(-Q_3*EpsPEqv)
            GDerivQ  = 2.0*q/SigM**2.0                                  ! dG/dq
            GDerivP  = -3.0*Eff_Void*BETA1*BETA2/SigM*SINH(TH)          ! dG/dp
            GDerivX1 = 2.0*BETA1*COSH(TH)-2.0*BETA3*Eff_Void            ! dG/df^*
            GDerivX2 = -2.0*q**2.0 * SigMDerivX2 /SigM**3.0             ! dG/deps_eq_p
     .                + 3.0*Eff_Void*BETA1*BETA2*p*SigMDerivX2
     .                / SigM**2.0*SINH(TH)

            G2DerivPQ = 0.0
            G2DerivP2 = (9.0*Eff_Void*BETA1*BETA2**2.0)/(2.0*SigM**2.0)
     .                   *COSH(TH)
            G2DerivQ2 = 2.0/SigM**2.0
            G2DerivQX1 = 0.0
            G2DerivQX2 = -4.0*q*SigMDerivX2/SigM**3.0
            G2DerivPX1 = -3.0*BETA1*BETA2/SigM*SINH(TH)
            G2DerivPX2 = (3.0*SigMDerivX2*Eff_Void*BETA1*BETA2/SigM**2.0)
     .                 *(SINH(TH)+TH*COSH(TH))
            h1DerivDeltaEpsP = 1.0-Void
            h1DerivP = 0.0
            h2DerivDeltaEpsP = (-p)/((1-Void)*SigM)
            h2DerivP = -DeltaEpsP/((1-Void)*SigM)
            h1DerivDeltaEpsQ = k_s*Void*omega_shear     ! 0 Initially
            h1DerivQ = 0.0
            h2DerivDeltaEpsQ = q/((1-Void)*SigM)
            h2DerivQ = (DeltaEpsQ) /((1-Void)*SigM)

            h1DerivX1 = -DeltaEpsP
            IF((p.LE.0.0).and.(FN.ne.0))THEN
                h1DerivX2 = -DeltaEpsPEqv*FN/(SN*sqrt(2.0*PI))
     .    *(EpsPEqv-EPSN)/(SN**2.0)*exp(-0.5*((EpsPEqv-EPSN)/SN)**2.0)
            ELSE
                h1DerivX2 = 0.0
            ENDIF
            h2DerivX1 = (-p*DeltaEpsP+q*DeltaEpsQ)
     .                   /(SigM*(1.0-Void)**2.0)
            h2DerivX2 = SigMDerivX2*(p*DeltaEpsP-q*DeltaEpsQ)
     .                 /((1.0-Void)*SigM**2.0)

            C11Inv = 1.0-h1DerivX1
            C12Inv = h1DerivX2
            C21Inv = h2DerivX1
            C22Inv = 1.0-h2DerivX2
            DetC = C11Inv*C22Inv-C12Inv*C21Inv
            C11 = C22Inv/DetC
            C12 = -C12Inv/DetC
            C21 = -C21Inv/DetC
            C22 = C11Inv/DetC
            X1DerivDeltaEpsP =
     .         +C11*((h1DerivDeltaEpsP)+K*(h1DerivP))
     .         +C12*((h2DerivDeltaEpsP)+K*(h2DerivP))
            X2DerivDeltaEpsP =
     .         +C22*((h2DerivDeltaEpsP)+K*(h2DerivP))
     .         +C21*((h1DerivDeltaEpsP)+K*(h1DerivP))
            X1DerivDeltaEpsQ =
     .         +C11*((h1DerivDeltaEpsQ)-3.0*G*(h1DerivQ))
     .         +C12*((h2DerivDeltaEpsQ)-3.0*G*(h2DerivQ))
            X2DerivDeltaEpsQ =
     .         +C22*((h2DerivDeltaEpsQ)-3.0*G*(h2DerivQ))
     .         +C21*((h1DerivDeltaEpsQ)-3.0*G*(h1DerivQ))
!-----------------------------------------------------------------------
!           Calculate Constants For Update
!-----------------------------------------------------------------------
            A11 = GDerivQ+DeltaEpsP*(K*G2DerivPQ
     .               +G2DerivQX1*X1DerivDeltaEpsP
     .               +G2DerivQX2*X2DerivDeltaEpsP)
     .               +DeltaEpsQ*(K*G2DerivP2
     .               +G2DerivPX1*X1DerivDeltaEpsP
     .               +G2DerivPX2*X2DerivDeltaEpsP)
            A12 = GDerivP+DeltaEpsP*(-3.0*G*G2DerivQ2
     .               +G2DerivQX1*X1DerivDeltaEpsQ
     .               +G2DerivQX2*X2DerivDeltaEpsQ)
     .               +DeltaEpsQ*(-3.0*G*G2DerivPQ
     .               +G2DerivPX1*X1DerivDeltaEpsQ
     .               +G2DerivPX2*X2DerivDeltaEpsQ)
            A21 = K*GDerivP+GDerivX1*X1DerivDeltaEpsP
     .                       +GDerivX2*X2DerivDeltaEpsP
            A22 = -3.0*G*GDerivQ+GDerivX1*X1DerivDeltaEpsQ
     .                            +GDerivX2*X2DerivDeltaEpsQ
            b1 = -DeltaEpsP*GDerivQ-DeltaEpsQ*GDerivP
            b2 = -F
            T1 = b2-A21*b1/A11
            T2 = A22-A21*A12/A11
            Cq = T1/T2
            Cp = (b1-A12*Cq)/A11
!
!-----------------------------------------------------------------------
!           Update Values
!-----------------------------------------------------------------------
            DeltaEpsP = DeltaEpsP+Cp
            DeltaEpsQ = DeltaEpsQ+Cq
!
!---------- Update pressure and eqv stress ----------------------------------              
            p=pTr+K*DeltaEpsP
            q=qTr-3.0*G*DeltaEpsQ
!
!---------- Update equivalent plastic strain ----------------------------------             
            DeltaEpsPEqv = (-p*DeltaEpsP+q*DeltaEpsQ)/((1.0-Void)*SigM)
!     
!---------- Update Void Growth term ----------------------------------   
            DeltaVoidG = DeltaEpsP*(1.0-Void)
!              
!---------- Update Void Nucleation term ---------------------------------- 
            IF ((p.LE.0.0).and.(FN.ne.0)) THEN
              DeltaVoidN = DeltaEpsPEqv*FN/(SN*sqrt(2.0*PI))
     .                   * exp(-0.5*((EpsPEqv-EPSN)/SN)**2)
            ELSE
              DeltaVoidN = 0.0
            END IF
!
!---------- Update Void Shear term ---------------------------------- 
            IF (k_s .ne. 0) THEN
              DeltaVoidS = k_s*VoidTr*omega_shear*DeltaEpsQ
            ELSE
              DeltaVoidS = 0.      
            END IF 
!
!---------- Update the equivalent plastic strain ---------------------------------------
            EpsPEqv = EpsPEqvTr + DeltaEpsPEqv
!              
!---------- Update yield strength ------------------------------------------------ 
            CALL VOCE(SigM, SIGMA0, EpsPEqv, T_1, Q_1, T_2, Q_2, T_3, Q_3)

!---------- Update bracket term ------------------------------------------------ 
            TH = (-3.0*BETA2*p)/(2.0*SigM)
!
!----------- Sum all void contributions ------------------------------------------------              
            DeltaVoid = DeltaVoidG + DeltaVoidN + DeltaVoidS
            Void = VoidTr+DeltaVoid
!
!            
!           VOID COALESCENCE!              
!----------- Calculate the effective void volume fraction f^* --------------------------   
            IF ((Void .ge. f_c).and.(f_f.ne.0)) THEN
              C_Factor = ((1/BETA1) - f_c) / (f_f - f_c)
              Eff_Void = f_c + C_Factor * (Void - f_c)
            ELSE 
              Eff_Void = Void
            END IF         
!
!---------- In order to avoid numerical errors -------------------------              
            IF (Eff_Void .ge. 0.99*(1/BETA1)) THEN
              Eff_Void = 0.99*(1/BETA1)
            END IF    
!            
!---------- Check the modified yield function --------------------------
            F = (q**2.0/SigM**2.0)+
     .             (2.0*Eff_Void*BETA1*COSH(TH))          
     .             -1.0-BETA3*Eff_Void**2.0  
!     
!-----------------------------------------------------------------------
!           Compute convergence criterion
!-----------------------------------------------------------------------
            RESNOR = ABS(F/SigM)
!-----------------------------------------------------------------------
!           Check for convergence
!-----------------------------------------------------------------------
            IF(RESNOR.LE.TOL)THEN ! RMAP has converged
!-----------------------------------------------------------------------
! Now, introducing the \kappa parameter in order to speed up the void deformations.
! This code is only executed if there are plastic deformations in the element, which is correct.
! The adjusted update in the void volume fraction is local, but the rate is adjusted dependent on its 
! neighbors plastics increments...
!-----------------------------------------------------------------------              
! Based on the plastic increment ratio NL/L, we will speed up / slow down the void growth:
!
                kappa = 1.      !If \kappa=1 is used, the local solution is obtained
                Void  =  VoidTr + kappa * DeltaVoid
!                
!-------------- In order to avoid numerical errors -----------------------   
                IF (f_f.ne.0) THEN  ! Then we have coalescence, and Void should be less than f_f
                  IF (Void .ge. 0.99*f_f) THEN
                    Void = 0.99*f_f
                  END IF  
                ELSE 
                  IF (Void .ge. 0.99*(1.0/BETA1)) THEN
                    Void = 0.99*(1.0/BETA1)
                  END IF  
                END IF
!
!-----------------------------------------------------------------------
!               Update the stress tensor
!-----------------------------------------------------------------------
                STRESSNEW(km,1) = (stresstr(1)+pTr)*(q/qTr) - p
                STRESSNEW(km,2) = (stresstr(2)+pTr)*(q/qTr) - p
                STRESSNEW(km,3) = (stresstr(3)+pTr)*(q/qTr) - p
                STRESSNEW(km,4) = stresstr(4)*(q/qTr)
                IF ( nshr .gt. 1 ) THEN
                  STRESSNEW(km,5) = stresstr(5)*(q/qTr)
                  STRESSNEW(km,6) = stresstr(6)*(q/qTr)
                END IF
!-----------------------------------------------------------------------
!              Update the history variables
!-----------------------------------------------------------------------
                IF(abs(Void).LE.1.0E-10)THEN
                  Void=0.0                 
                ENDIF
                STATENEW(km,1) = Void      ! VoidFrac
                STATENEW(km,2) = EpsPEqv   ! EpsPEqv
                STATENEW(km,3) = q         ! SigEqv
                STATENEW(km,4) = p         ! SigH
                STATENEW(km,5) = SigM      ! SigM
                STATENEW(km,6) = F         ! F
                STATENEW(km,8) = EpsPEqv   ! Local plastic strain
!---------------------------------------------------------------------------
! Delete the element if the void content has reached the critical value:  
                IF (STATEOLD(km,7) .eq. 1) THEN         ! If the element hasn't been deleted yet
                  IF (f_f .ne. 0) THEN 
                    IF (Void .ge. 0.98*f_f) THEN 
                      STATENEW(km,7) = 0
                    ELSE 
                      STATENEW(km,7) = 1
                    END IF 
                  ELSE 
                    IF (Void.ge.(0.98*(1/BETA1))) THEN
                      STATENEW(km,7) = 0                  ! If critical void content: delete the element
                    ELSE 
                      STATENEW(km,7) = 1
                    END IF 
                  END IF 
                ELSE    
                  STATENEW(km,7) = 0                    ! The element was deleted in an earlier increment..
                END IF          
!---------------------------------------------------------------------------                
                GOTO 90
              ELSE ! RMAP has not converged yet
                IF(ITER.eq.(MXITER-1))THEN
                  PRINT*, "RMAP has not converged!"
                  PRINT*, "Convergence: ", RESNOR
                  PRINT*, Void
                  STOP
                ENDIF
              ENDIF
            ENDDO  !The plastic convergence loop in the vumat
          ELSE
!
!-----------------------------------------------------------------------
! If this case, the element has undergone purely elastic deformation:  
!-----------------------------------------------------------------------
!        Update the stress tensor
!-----------------------------------------------------------------------
            STRESSNEW(km,1)   = stresstr(1)
            STRESSNEW(km,2)   = stresstr(2)
            STRESSNEW(km,3)   = stresstr(3)
            STRESSNEW(km,4)   = stresstr(4)
            IF ( nshr .gt. 1 ) THEN
              STRESSNEW(km,5) = stresstr(5)
              STRESSNEW(km,6) = stresstr(6)
            END IF
!-----------------------------------------------------------------------
!        Update the history variables
!-----------------------------------------------------------------------
            STATENEW(km,1) = VoidTr         ! VoidFrac
            STATENEW(km,2) = EpsPEqvTr      ! EpsPEqv
            STATENEW(km,3) = SigEqv         ! SigEqv
            STATENEW(km,4) = pTr            ! SigH
            STATENEW(km,5) = SigM           ! SigM
            STATENEW(km,6) = F              ! F
            STATENEW(km,7) = 1              ! Delete the element or not..
            STATENEW(km,9) = EpsPEqvTr      ! EpsPEqv
          ENDIF
  90  CONTINUE

!-----------------------------------------------------------------------
!      Store the specific internal energy -
!-----------------------------------------------------------------------  
          if ( nshr .eq. 1 ) then
            stressPower = half * (
     .      ( stressOld(km,1) + stressNew(km,1) ) * strainInc(km,1) +
     .      ( stressOld(km,2) + stressNew(km,2) ) * strainInc(km,2) +
     .      ( stressOld(km,3) + stressNew(km,3) ) * strainInc(km,3) ) +
     .      ( stressOld(km,4) + stressNew(km,4) ) * strainInc(km,4)
          else
            stressPower = half * (
     .      ( stressOld(km,1) + stressNew(km,1) ) * strainInc(km,1) +
     .      ( stressOld(km,2) + stressNew(km,2) ) * strainInc(km,2) +
     .      ( stressOld(km,3) + stressNew(km,3) ) * strainInc(km,3) ) +
     .      ( stressOld(km,4) + stressNew(km,4) ) * strainInc(km,4) +
     .      ( stressOld(km,5) + stressNew(km,5) ) * strainInc(km,5) +
     .      ( stressOld(km,6) + stressNew(km,6) ) * strainInc(km,6)
          end if
          enerInternNew(km) = enerInternOld(km) + stressPower / density(km)  
!          
!-----------------------------------------------------------------------
!     END the material point loop and vumat branch:
!-----------------------------------------------------------------------
      END DO                ! End of the material point loop
      END IF                ! End of the preprocessor / vumat branch
!
!-----------------------------------------------------------------------
!     END SUBROUTINE
!-----------------------------------------------------------------------
      RETURN
      END
!      











!      
!
!
! 3 term Voce hardening rule:
      SUBROUTINE VOCE(sigma_new, sigma0, peeqOld, T_1, Q_1, T_2, 
     .Q_2, T_3, Q_3)
     
        include 'vaba_param.inc'
!
        term1 = T_1 * (1 - EXP(-Q_1 * peeqOld))
        term2 = T_2 * (1 - EXP(-Q_2 * peeqOld))
        term3 = T_3 * (1 - EXP(-Q_3 * peeqOld))
      
        sigma_new = sigma0 + term1 + term2 + term3
      
      RETURN 
      END
!
!