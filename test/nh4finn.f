C       Subroutine POT to compute (H2) + (NH2-) potential energy 
C       for a given geometry
C       Atom order : H--H--H--H--N
C       Units for  : Cartesian in Angstrom and energy in eV
C       Total 45932 CCSD(T)/F12a/aug-cc-pVTZ points fitted by FI-NN method
C       Fitting RMSE= 5.0307meV, which is average of 2 PESs; 
C       Update by Mengyi Pan
C       Subroutine "PREPOT" must be called only once before "POT"
C       The data file 'PARA.txt' and 'W0**.txt' should be in the same  
C       directory as the program
C       Date: 1-June-2021 
C       Corresponding author: Hongwei Song
C       Email:   hwsong@wipm.ac.cn

C
C       NNPES PARAMETERS INFO.
C       MULTIDIMENSION ARRAYS ARE ARRANGED BY ... NEUROS, LAYERS, NPES, DISTRICTS
C
        MODULE NNDATA
        IMPLICIT NONE
        INTEGER :: NA,DIN,NLAYS,NPES
        INTEGER,ALLOCATABLE :: NEU(:)
        REAL*8, ALLOCATABLE :: XRANGE(:,:,:),YRANGE(:,:),W(:,:,:,:),
     &  B(:,:,:),RLI(:,:,:)
        CHARACTER*100 :: PFNAME
        END MODULE

C
C       READ ALL LMNNPARA
C
        SUBROUTINE PREPOT
        USE NNDATA
        IMPLICIT NONE

        INTEGER :: IPES,ILAYS,jneuron,jp,i,j

        OPEN(10,FILE='PARA.txt',STATUS='OLD')
        READ(10,*) NA,DIN,NLAYS,NPES

        ALLOCATE(NEU(NLAYS+1))
        NEU(1)=DIN
        READ(10,*) NEU(2:NLAYS+1)

        ALLOCATE(XRANGE(2,DIN,NPES))
        ALLOCATE(YRANGE(2,NPES))
        ALLOCATE(W(MAXVAL(NEU),MAXVAL(NEU),NLAYS,NPES))
        ALLOCATE(B(MAXVAL(NEU),NLAYS,NPES))
        ALLOCATE(RLI(MAXVAL(NEU),NLAYS+1,NPES))
        
        XRANGE=0.D0
        YRANGE=0.D0
        W=0.D0
        B=0.D0
        RLI=0.D0
        IPES=0;ILAYS=0

        DO IPES=1,NPES
            READ(10,*) PFNAME
            OPEN(11,FILE=PFNAME,STATUS='OLD')
            DO ILAYS=1,NLAYS
              READ(11,*)W(1:NEU(ILAYS),1:NEU(ILAYS+1),ILAYS,IPES)
              READ(11,*) 
              READ(11,*) B(1:NEU(ILAYS+1),ILAYS,IPES)
              READ(11,*)
            ENDDO
            READ(11,*) XRANGE(:,:,IPES)
            READ(11,*)
            READ(11,*) YRANGE(:,IPES)
            CLOSE(11)
        ENDDO
        CLOSE(10)

        RETURN
        END
!---------------------------------------------------------------------------
C       NH4NNPOT, R(:)/ANGSTROM, V(NPES)/EV      
!---------------------------------------------------------------------------
        SUBROUTINE NNPOT(R,VOUT)        
        USE NNDATA
        IMPLICIT NONE
        REAL*8,INTENT(IN)  :: R(DIN)
        REAL*8 V(NPES)
        REAL*8 VOUT,RMSE
        INTEGER :: IP,I,J,K,jjj
        REAL*8 :: VFIT,TEST
        
        DO IP=1,NPES
C
C       INPUT LAYER
C
        DO I=1,DIN
        RLI(I,1,IP)=2.D0*(R(I)-XRANGE(1,I,IP))/(XRANGE(2,I,IP)
     &              -XRANGE(1,I,IP))-1.D0
        ENDDO
C
C       layer 1-3
C
        DO I=1,NLAYS
           RLI(1:NEU(I+1),I+1,IP)=B(1:NEU(I+1),I,IP)
           DO J=1,NEU(I+1)
                IF (I==NLAYS)THEN
                   DO K=1,NEU(I)
                      RLI(J,I+1,IP)=RLI(J,I+1,IP)+W(K,J,I,IP)
     &               *RLI(K,I,IP)
                   ENDDO
                ELSE
                   DO K=1,NEU(I)
                      RLI(J,I+1,IP)=RLI(J,I+1,IP)+W(K,J,I,IP)
     &               *RLI(K,I,IP)
                   ENDDO
                   TEST=-2.D0*RLI(J,I+1,IP)
                   IF(TEST<-100.D0) THEN
                       RLI(J,I+1,IP)=1.D0
                   ELSE IF(TEST>-100.D0 .AND. TEST<100.D0) THEN
                       RLI(J,I+1,IP)=2.D0/(1.D0+EXP(TEST
     &                ))-1.D0
                   ELSE
                       RLI(J,I+1,IP)=-1.D0
                   ENDIF
                ENDIF
           ENDDO
        ENDDO
        VFIT=RLI(1,NLAYS+1,IP)
        V(IP)=0.5D0*(VFIT+1.D0)*(YRANGE(2,IP)-YRANGE(1,IP))
     &  +YRANGE(1,IP)
        ENDDO

        VOUT=0.D0

        DO I=1,NPES                        
           VOUT=VOUT+V(I)            ! V average, uint in eV
        ENDDO
        VOUT=VOUT/(NPES*1.D0)

        RETURN
        END

!**************************************************************************
C       POTENTIAL ENERGY SUBROUTINE INTERFACE
C       INPUT: C/ANGSTROM, V/EV
        SUBROUTINE POT(CI,VOUT)
        USE NNDATA, ONLY:NA,NPES
        IMPLICIT NONE
        INTEGER :: I,K
        REAL*8  :: V(NPES),VOUT,t
        REAL*8  :: CI(3,NA),C(3,NA),R(NA*(NA-1)/2)
        REAL*8  :: R1,V1,R2(3),V2
        REAL*8  :: RT(NA*(NA-1)/2)                            
        REAL*8  :: COST,VT,RMSE,RR,VMIN
        REAL,PARAMETER::alpha=1.d0               
        real*8  ::RINT(31)                       
        C=CI
        CALL DISTANC(NA,C,R)

        RT(:)=dexp(-R(:)/alpha)    
       
        call fi_a4b(RT,RINT)       
        CALL NNPOT(RINT,VOUT,RMSE)   

        RETURN
        END
!**************************************************************	
        SUBROUTINE DISTANC(NA,C,R)
        IMPLICIT NONE
        INTEGER :: I,J,IJ,K,NA
        REAL*8  :: C(3,NA),R(NA*(NA-1)/2),t

        R=0.D0
        K=0
        DO I=1,NA-1
           DO J=I+1,NA
              K=K+1
              DO IJ=1,3
                 R(K)=R(K)+(C(IJ,I)-C(IJ,J))**2
              ENDDO
              R(K)=DSQRT(R(K))
           ENDDO
        ENDDO

        RETURN
        END

!****************************************************************
        SUBROUTINE DESTRUCT
        USE NNDATA

        DEALLOCATE(W)
        DEALLOCATE(RLI)
        DEALLOCATE(B)
        DEALLOCATE(NEU)
        DEALLOCATE(XRANGE)

        RETURN
        END
!***************************************************************
!---------------------------------------------------------------------

! compute the fundamental invariants
! atom order : A A A A B
! r (input) real(kind=8) the internuclear distances
! p (output) real(kind=8) the fundamental invariants
      subroutine fi_a4b(r,p)
      implicit none
      real(kind=8),intent(in) :: r(10)
      real(kind=8),intent(out) :: p(31)
      p(1)=r(1)+r(3)+r(2)+r(8)+r(5)+r(6)
      p(2)=r(4)+r(10)+r(9)+r(7)
      p(3)=r(1)**2+r(3)**2+r(2)**2+r(8)**2+r(5)**2+r(6)**2
      p(4)=r(4)**2+r(10)**2+r(9)**2+r(7)**2
      p(5)=r(1)*r(2)+r(1)*r(3)+r(3)*r(6)+r(3)*r(8)+r(2)*r(3)
     &     +r(2)*r(5)+r(2)*r(8)+r(6)*r(8)+r(5)*r(8)+r(1)*r(5)
     &     +r(5)*r(6)+r(1)*r(6)
      p(6)=r(1)*r(4)+r(3)*r(4)+r(3)*r(10)+r(2)*r(4)+r(2)*r(9)
     &     +r(8)*r(9)+r(8)*r(10)+r(5)*r(9)+r(5)*r(7)+r(6)*r(7)
     &     +r(6)*r(10)+r(1)*r(7)
      p(7)=r(1)**3+r(3)**3+r(2)**3+r(8)**3+r(5)**3+r(6)**3
      p(8)=r(4)**3+r(10)**3+r(9)**3+r(7)**3
      p(9)=r(1)**2*r(2)+r(1)**2*r(3)+r(1)*r(3)**2+r(3)**2*r(6)
     &     +r(3)**2*r(8)+r(2)*r(3)**2+r(2)**2*r(3)+r(1)*r(2)**2
     &     +r(2)**2*r(5)+r(2)**2*r(8)+r(2)*r(8)**2+r(3)*r(8)**2
     &     +r(6)*r(8)**2+r(5)*r(8)**2+r(5)**2*r(8)+r(2)*r(5)**2
     &     +r(1)*r(5)**2+r(5)**2*r(6)+r(5)*r(6)**2+r(6)**2*r(8)
     &     +r(3)*r(6)**2+r(1)*r(6)**2+r(1)**2*r(6)+r(1)**2*r(5)
      p(10)=r(1)**2*r(4)+r(3)**2*r(4)+r(3)**2*r(10)+r(2)**2*r(4)
     &      +r(2)**2*r(9)+r(8)**2*r(9)+r(8)**2*r(10)+r(5)**2*r(9)
     &      +r(5)**2*r(7)+r(6)**2*r(7)+r(6)**2*r(10)+r(1)**2*r(7)
      p(11)=r(4)**2*r(5)+r(4)**2*r(6)+r(1)*r(10)**2+r(2)*r(10)**2
     &      +r(4)**2*r(8)+r(1)*r(9)**2+r(3)*r(9)**2+r(5)*r(10)**2
     &      +r(6)*r(9)**2+r(2)*r(7)**2+r(7)**2*r(8)+r(3)*r(7)**2
      p(12)=r(1)*r(2)*r(3)+r(3)*r(6)*r(8)+r(2)*r(5)*r(8)+r(1)*r(5)*r(6)
      p(13)=r(1)*r(2)*r(4)+r(1)*r(3)*r(4)+r(3)*r(6)*r(10)+
     &      r(3)*r(8)*r(10)+r(2)*r(3)*r(4)+r(2)*r(5)*r(9)+r(2)*r(8)*r(9)
     &      +r(6)*r(8)*r(10)+r(5)*r(8)*r(9)+r(1)*r(5)*r(7)+
     &      r(5)*r(6)*r(7)+r(1)*r(6)*r(7)
      p(14)=r(1)*r(4)*r(7)+r(3)*r(4)*r(10)+r(2)*r(4)*r(9)
     &      +r(8)*r(9)*r(10)+r(5)*r(7)*r(9)+r(6)*r(7)*r(10)
      p(15)=r(1)**4+r(3)**4+r(2)**4+r(8)**4+r(5)**4+r(6)**4
      p(16)=r(4)**4+r(10)**4+r(9)**4+r(7)**4
      p(17)=r(1)**3*r(2)+r(1)**3*r(3)+r(1)*r(3)**3+r(3)**3*r(6)+
     &      r(3)**3*r(8)+r(2)*r(3)**3+r(2)**3*r(3)+r(1)*r(2)**3+
     &      r(2)**3*r(5)+r(2)**3*r(8)+r(2)*r(8)**3+r(3)*r(8)**3+
     &      r(6)*r(8)**3+r(5)*r(8)**3+r(5)**3*r(8)+r(2)*r(5)**3+
     &      r(1)*r(5)**3+r(5)**3*r(6)+r(5)*r(6)**3+r(6)**3*r(8)+
     &      r(3)*r(6)**3+r(1)*r(6)**3+r(1)**3*r(6)+r(1)**3*r(5)
      p(18)=r(1)**3*r(4)+r(3)**3*r(4)+r(3)**3*r(10)+r(2)**3*r(4)+
     &      r(2)**3*r(9)+r(8)**3*r(9)+r(8)**3*r(10)+r(5)**3*r(9)+
     &      r(5)**3*r(7)+r(6)**3*r(7)+r(6)**3*r(10)+r(1)**3*r(7)
      p(19)=r(4)**3*r(5)+r(4)**3*r(6)+r(1)*r(10)**3+r(2)*r(10)**3+
     &      r(4)**3*r(8)+r(1)*r(9)**3+r(3)*r(9)**3+r(5)*r(10)**3+
     &      r(6)*r(9)**3+r(2)*r(7)**3+r(7)**3*r(8)+r(3)*r(7)**3
      p(20)=r(1)**2*r(4)**2+r(3)**2*r(4)**2+r(3)**2*r(10)**2+
     &      r(2)**2*r(4)**2+r(2)**2*r(9)**2+r(8)**2*r(9)**2+
     &      r(8)**2*r(10)**2+r(5)**2*r(9)**2+r(5)**2*r(7)**2+
     &      r(6)**2*r(7)**2+r(6)**2*r(10)**2+r(1)**2*r(7)**2
      p(21)=r(1)**2*r(2)*r(4)+r(1)**2*r(3)*r(4)+r(1)*r(3)**2*r(4)+
     &      r(3)**2*r(6)*r(10)+r(3)**2*r(8)*r(10)+r(2)*r(3)**2*r(4)+
     &      r(2)**2*r(3)*r(4)+r(1)*r(2)**2*r(4)+r(2)**2*r(5)*r(9)+
     &      r(2)**2*r(8)*r(9)+r(2)*r(8)**2*r(9)+r(3)*r(8)**2*r(10)+
     &      r(6)*r(8)**2*r(10)+r(5)*r(8)**2*r(9)+r(5)**2*r(8)*r(9)+
     &      r(2)*r(5)**2*r(9)+r(1)*r(5)**2*r(7)+r(5)**2*r(6)*r(7)+
     &      r(5)*r(6)**2*r(7)+r(6)**2*r(8)*r(10)+r(3)*r(6)**2*r(10)+
     &      r(1)*r(6)**2*r(7)+r(1)**2*r(6)*r(7)+r(1)**2*r(5)*r(7)
      p(22)=r(1)**2*r(2)*r(7)+r(1)**2*r(3)*r(7)+r(1)*r(3)**2*r(10)+
     &      r(3)**2*r(4)*r(6)+r(3)**2*r(4)*r(8)+r(2)*r(3)**2*r(10)+
     &      r(2)**2*r(3)*r(9)+r(1)*r(2)**2*r(9)+r(2)**2*r(4)*r(5)+
     &      r(2)**2*r(4)*r(8)+r(2)*r(8)**2*r(10)+r(3)*r(8)**2*r(9)+
     &      r(6)*r(8)**2*r(9)+r(5)*r(8)**2*r(10)+r(5)**2*r(7)*r(8)+
     &      r(2)*r(5)**2*r(7)+r(1)*r(5)**2*r(9)+r(5)**2*r(6)*r(9)+
     &      r(5)*r(6)**2*r(10)+r(6)**2*r(7)*r(8)+r(3)*r(6)**2*r(7)+
     &      r(1)*r(6)**2*r(10)+r(1)**2*r(4)*r(6)+r(1)**2*r(4)*r(5)
      p(23)=r(1)**2*r(4)*r(7)+r(3)**2*r(4)*r(10)+r(2)**2*r(4)*r(9)+
     &      r(8)**2*r(9)*r(10)+r(5)**2*r(7)*r(9)+r(6)**2*r(7)*r(10)
      p(24)=r(4)**2*r(5)*r(6)+r(4)**2*r(6)*r(8)+r(1)*r(2)*r(10)**2+
     &      r(4)**2*r(5)*r(8)+r(1)*r(3)*r(9)**2+r(3)*r(6)*r(9)**2+
     &      r(2)*r(5)*r(10)**2+r(1)*r(6)*r(9)**2+r(2)*r(7)**2*r(8)+
     &      r(3)*r(7)**2*r(8)+r(1)*r(5)*r(10)**2+r(2)*r(3)*r(7)**2
      p(25)=r(1)**5+r(3)**5+r(2)**5+r(8)**5+r(5)**5+r(6)**5
      p(26)=r(1)**4*r(4)+r(3)**4*r(4)+r(3)**4*r(10)+r(2)**4*r(4)+
     &      r(2)**4*r(9)+r(8)**4*r(9)+r(8)**4*r(10)+r(5)**4*r(9)+
     &      r(5)**4*r(7)+r(6)**4*r(7)+r(6)**4*r(10)+r(1)**4*r(7)
      p(27)=r(1)**3*r(4)**2+r(3)**3*r(4)**2+r(3)**3*r(10)**2+
     &      r(2)**3*r(4)**2+r(2)**3*r(9)**2+r(8)**3*r(9)**2+
     &      r(8)**3*r(10)**2+r(5)**3*r(9)**2+r(5)**3*r(7)**2+
     &      r(6)**3*r(7)**2+r(6)**3*r(10)**2+r(1)**3*r(7)**2
      p(28)=r(4)**3*r(5)**2+r(4)**3*r(6)**2+r(1)**2*r(10)**3+
     &      r(2)**2*r(10)**3+r(4)**3*r(8)**2+r(1)**2*r(9)**3+
     &      r(3)**2*r(9)**3+r(5)**2*r(10)**3+r(6)**2*r(9)**3+
     &      r(2)**2*r(7)**3+r(7)**3*r(8)**2+r(3)**2*r(7)**3
      p(29)=r(1)**3*r(2)*r(4)+r(1)**3*r(3)*r(4)+r(1)*r(3)**3*r(4)+
     &      r(3)**3*r(6)*r(10)+r(3)**3*r(8)*r(10)+r(2)*r(3)**3*r(4)+
     &      r(2)**3*r(3)*r(4)+r(1)*r(2)**3*r(4)+r(2)**3*r(5)*r(9)+
     &      r(2)**3*r(8)*r(9)+r(2)*r(8)**3*r(9)+r(3)*r(8)**3*r(10)+
     &      r(6)*r(8)**3*r(10)+r(5)*r(8)**3*r(9)+r(5)**3*r(8)*r(9)+
     &      r(2)*r(5)**3*r(9)+r(1)*r(5)**3*r(7)+r(5)**3*r(6)*r(7)+
     &      r(5)*r(6)**3*r(7)+r(6)**3*r(8)*r(10)+r(3)*r(6)**3*r(10)+
     &      r(1)*r(6)**3*r(7)+r(1)**3*r(6)*r(7)+r(1)**3*r(5)*r(7)
      p(30)=r(1)**3*r(4)*r(7)+r(3)**3*r(4)*r(10)+r(2)**3*r(4)*r(9)+
     &      r(8)**3*r(9)*r(10)+r(5)**3*r(7)*r(9)+r(6)**3*r(7)*r(10)
      p(31)=r(4)**3*r(5)*r(7)+r(4)**3*r(6)*r(7)+r(4)**3*r(6)*r(10)+
     &      r(1)*r(4)*r(10)**3+r(2)*r(4)*r(10)**3+r(4)**3*r(8)*r(10)+
     &      r(4)**3*r(8)*r(9)+r(4)**3*r(5)*r(9)+r(1)*r(4)*r(9)**3+
     &      r(3)*r(4)*r(9)**3+r(3)*r(9)**3*r(10)+r(2)*r(9)*r(10)**3+
     &      r(5)*r(9)*r(10)**3+r(6)*r(9)**3*r(10)+r(6)*r(7)*r(9)**3+
     &      r(1)*r(7)*r(9)**3+r(2)*r(7)**3*r(9)+r(7)**3*r(8)*r(9)+
     &      r(7)**3*r(8)*r(10)+r(5)*r(7)*r(10)**3+r(1)*r(7)*r(10)**3+
     &      r(3)*r(7)**3*r(10)+r(3)*r(4)*r(7)**3+r(2)*r(4)*r(7)**3
      end subroutine fi_a4b
!------------------------------------------------------------------------

