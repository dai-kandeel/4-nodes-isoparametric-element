c Copyright ANSYS.  All Rights Reserved.
c User defined Isoparametric Element
c MMS @ UDG
      
            module AllinOne
                implicit none
                contains

            subroutine Bmatrixx(xi,eta,x1,x2,x3,x4,y1,y2,y3,y4,
     &                 BMatrix,det)
                  implicit none
                  real xi,eta, par_1,par_2
                  DOUBLE PRECISION, dimension(2,2) :: jac
                  DOUBLE PRECISION, dimension(2,2) :: inv_jac
                  INTEGER  , dimension (2,2) :: ipiv
                  DOUBLE PRECISION , dimension (2,2) :: work
                  DOUBLE PRECISION , dimension (3,8) :: BMatrix
                  real :: x1,x2,x3,x4,y1,y2,y3,y4,N1_xi,N2_xi,N3_xi,
     &                    N4_xi,N1_eta,N2_eta,N3_eta,N4_eta,det

                  !Define the Jacobian matrix

                  par_1=xi/4.0D0 - 1.0D0/4.0D0
                  par_2=(xi/4.0D0) + 1.0D0/4.0D0
                  jac(1,1)=(x1*(eta-1.0D0)/4.0D0)-(x2*(eta-1.0D0)/4.0D0)
     &                    +(x3*(eta+1.0D0)/4.0D0)-(x4*(eta+1.0D0)/4.0D0)
                  jac(1,2)=(y1*(eta-1.0D0)/4.0D0)-
     &                     (y2*(eta-1.0D0)/4.0D0)+(y3*(eta+1.0D0)/4.0D0)
     &                      -(y4*(eta+1.0D0)/4.0D0)
                  jac(2,1)=(x1*par_1)-(x2*par_2)+(x3*par_2)-(x4*par_1)
                  jac(2,2)=(y1*par_1)-(y2*par_2)+(y3*par_2)-(y4*par_1)

                  !compute the determinant of the Jacobian

                  det=(jac(1,1)*jac(2,2)) - (jac(1,2)*jac(2,1))

                  !Get the inverse of the Jacobian
                  inv_jac(1,1)=jac(2,2)
                  inv_jac(2,2)=jac(1,1)
                  inv_jac(1,2)=-jac(1,2)
                  inv_jac(2,1)=-jac(2,1)
                  inv_jac=(1/det)*inv_jac

                  !Define derivatives
                  N1_xi=0.25*eta-0.25
                  N2_xi=0.25-0.25*eta
                  N3_xi=0.25*eta+0.25
                  N4_xi=-0.25*eta-0.25

                  N1_eta=0.25*xi-0.25
                  N2_eta=-0.25*xi-0.25
                  N3_eta=0.25*xi+0.25
                  N4_eta=0.25-0.25*xi

                  !Get the B-matrix
                  BMatrix(1,1)=(N1_xi*inv_jac(1,1))+((N1_eta)
     &            *inv_jac(1,2))
                  BMatrix(1,2)=0.0D0
                  BMatrix(1,3)=(N2_xi*inv_jac(1,1))+((N2_eta)
     &            *inv_jac(1,2))
                  BMatrix(1,4)=0.0D0
                  BMatrix(1,5)=(N3_xi*inv_jac(1,1))+((N3_eta)
     &            *inv_jac(1,2))
                  BMatrix(1,6)=0.0D0
                  BMatrix(1,7)=(N4_xi*inv_jac(1,1))+((N4_eta)
     &            *inv_jac(1,2))
                  BMatrix(1,8)=0.0D0

                  BMatrix(2,1)=0.0D0
                  BMatrix(2,2)=(N1_xi*inv_jac(2,1))+((N1_eta)
     &            *inv_jac(2,2))
                  BMatrix(2,3)=0.0D0
                  BMatrix(2,4)=(N2_xi*inv_jac(2,1))+((N2_eta)
     &            *inv_jac(2,2))
                  BMatrix(2,5)=0.0D0
                  BMatrix(2,6)=(N3_xi*inv_jac(2,1))+((N3_eta)
     &            *inv_jac(2,2))
                  BMatrix(2,7)=0.0D0
                  BMatrix(2,8)=(N4_xi*inv_jac(2,1))+((N4_eta)
     &            *inv_jac(2,2))

                  BMatrix(3,1)=(N1_xi*inv_jac(2,1))+((N1_eta)
     &            *inv_jac(2,2))
                  BMatrix(3,2)=(N1_xi*inv_jac(1,1))+((N1_eta)
     &            *inv_jac(1,2))
                  BMatrix(3,3)=(N2_xi*inv_jac(2,1))+((N2_eta)
     &            *inv_jac(2,2))
                  BMatrix(3,4)=(N2_xi*inv_jac(1,1))+((N2_eta)
     &            *inv_jac(1,2))
                  BMatrix(3,5)=(N3_xi*inv_jac(2,1))+((N3_eta)
     &            *inv_jac(2,2))
                  BMatrix(3,6)=(N3_xi*inv_jac(1,1))+((N3_eta)
     &            *inv_jac(1,2))
                  BMatrix(3,7)=(N4_xi*inv_jac(2,1))+((N4_eta)
     &            *inv_jac(2,2))
                  BMatrix(3,8)=(N4_xi*inv_jac(1,1))+((N4_eta)
     &            *inv_jac(1,2))


            end subroutine

            subroutine stiffer(x1,x2,x3,x4,y1,y2,y3,y4,t,gau_n,E,v,
     &                eStiff,KeyOP)
                implicit none
                integer  KeyOp, i,j, gau_n !gau_n is number of gauss points in 1 direction
                real :: det,x1,x2,x3,x4,y1,y2,y3,y4,xi,eta
                DOUBLE PRECISION , dimension (3,8) :: BMatrix
                DOUBLE PRECISION , dimension (8,8) :: eStiff
                real , dimension (3,3) :: D  !material matrix
                real , dimension (8,3) :: Gf1
                real , dimension (8,8) :: Gf
                real , dimension (8,8) :: eq
                real , dimension(2) :: g
                real , dimension(2) :: w
                DOUBLE PRECISION :: v,E,t
                
                g(1)=-1/sqrt(3.0D0)  !gauss points coordinates in natural space
                g(2)=1/sqrt(3.0D0)

                w(1)=1           !weights
                w(2)=1


                !Plane stress
                IF (KeyOp.EQ.0) THEN

                    D(1,1)=1.0D0
                    D(1,2)=v
                    D(1,3)=0.0D0

                    D(2,1)=v
                    D(2,2)=1.0D0
                    D(2,3)=0.0D0

                    D(3,1)=0.0D0
                    D(3,2)=0.0D0
                    D(3,3)=(1-v)/2
                    D= E/(1-(v*v)) * D

                  !plane strain
                Else IF (KeyOp.EQ.1) THEN
                    
                    D(1,1)=1.0D0-v
                    D(1,2)=v
                    D(1,3)=0.0D0

                    D(2,1)=v
                    D(2,2)=1.0D0-v
                    D(2,3)=0.0D0

                    D(3,1)=0.0D0
                    D(3,2)=0.0D0
                    D(3,3)=(1.0D0-2.0D0*v)/2.0D0
                    D=(E/((1.0D0+v)*(1.0D0-2.0D0*v))) * D

                END IF
                eStiff=0
                Do i =1,gau_n                   
                    Do j=1,gau_n
                        xi=g(i)
                        eta=g(j)
                        call Bmatrixx(xi,eta,x1,x2,x3,x4,y1,y2,y3,
     &                                y4,BMatrix,det)
                        Gf1=matmul(Transpose(BMatrix),D)!Gauss function
                        Gf1=(t*det)*Gf1

                        Gf=matmul(Gf1,BMatrix)
                        eq=w(i)*w(j)*Gf  !Gauss function multiplied by weights
                        eStiff=eStiff+eq
                    End Do
                End Do

            end subroutine


            end module

      subroutine UserElem (elId, matId, keyMtx, lumpm, nDim, nNodes,
     &                     Nodes, nIntPnts, nUsrDof, kEStress, 
     &                     keyAnsMat, keySym, nKeyOpt, KeyOpt,
     &                     temper, temperB, tRef, kTherm, 
     &                     nPress, Press, kPress, nReal, RealConst, 
     &                     nSaveVars, saveVars, xRef, xCur, 
     &                     TotValDofs, IncValDofs, ItrValDofs,
     &                     VelValDofs, AccValDofs,
     &                     kfstps, nlgeom, nrkey, outkey, elPrint, iott,
     &                     keyHisUpd, ldstep, isubst, ieqitr, timval, 
     &                     keyEleErr, keyEleCnv,
     &                     eStiff, eMass, eDamp, eSStiff,
     &                     fExt, fInt, elVol, elMass, elCG, 
     &                     nRsltBsc, RsltBsc, nRsltVar, RsltVar, 
     &                     nElEng, elEnergy) 
                   


      use AllinOne
#include "impcom.inc"

      EXTERNAL         ElemGetMat  

      INTEGER          elId, matId, keyMtx(10), lumpm, nDim, nNodes,
     &                 Nodes(nNodes), nIntPnts, nUsrDof, kEStress, 
     &                 keyAnsMat, keySym, nKeyOpt, KeyOpt(nKeyOpt),
     *                 kTherm, nPress, kPress, nReal, nSaveVars, 
     &                 kfstps, nlgeom, nrkey, outkey, 
     &                 elPrint, iott, keyHisUpd, 
     &                 ldstep, isubst, ieqitr, keyEleErr, keyEleCnv,
     &                 nRsltBsc, nRsltVar, nElEng


      DOUBLE PRECISION temper(nNodes), temperB(nNodes), tRef, 
     &                 Press(nPress), RealConst(nReal),
     &                 saveVars(nSaveVars), 
     &                 xRef(nDim,nNodes), xCur(nDim,nNodes),
     &                 TotValDofs(nUsrDof), IncValDofs(nUsrDof), 
     &                 ItrValDofs(nUsrDof), VelValDofs(nUsrDof),
     &                 AccValDofs(nUsrDof),      timval,
     &                 eStiff(nUsrDof,nUsrDof), eMass(nUsrDof,nUsrDof), 
     &                 eDamp(nUsrDof,nUsrDof), eSStiff(nUsrDof,nUsrDof), 
     &                 fExt(nUsrDof), fInt(nUsrDof), 
     &                 elVol, elMass, elCG(3),
     &                 RsltBsc(nRsltBsc), RsltVar(nRsltVar), 
     &                 elEnergy(nElEng),stifoo(8,8)


#include "locknm.inc"
c Subroutine UElMatx (Accessing Element Matrices and Load Vectors)
c Function GetRForce (Getting Nodal Reaction Force values) (204)
c Function GetStackDisp (Getting Current Displacement Values) (205)
c Function disget (Getting a cinstraint from the database) (236)
c Function disput (Storing a constraint at node) (236)
c Function forget (Getting a constraint from the database) (237)
c Function forput ( Storing a nodal load at node) (237)
c FUNCTION elmget (Getting an element's characterestics and nodes)

      EXTERNAL         vzero, vmove, vmult, vdot, vidot,
     &                 maxv, matxb, matba, maat, matsym, getMatProp,
     &                 erhandler, equivStrain, ElemJac, ElemMass,
     &                 ElemRsltNode, ElemShpFn, pplock, ppunlock,
     &                 disget, forget

      DOUBLE PRECISION vdot, vidot

      INTEGER          nUsrDof2, intPnt, iNode, nTens, flgSingular,dis,
     &                 k1, k2, k3, nComp, iDim, iDim1, iComp, yourLabel,
     &                 nNodesCorner, nDirect, kThermIP, yourNode, disget
     &                 ,real_node, forget, forg, dspget
                      
      DOUBLE PRECISION BMat(nDim*2,nUsrDof), Ex, nu, density, G, workDb,
     &                 con1, con2, cMat(nDim*2,nDim*2), shIsoC(nNodes),
     &                 shIso(nNodes), shDerIso(nDim,nNodes), wtIP(1),
     &                 workArr(360), elJac(nDim*nDim), detJac, dperr(2),
     &                 shDerEl(nDim,nNodes), dVol, Strain(nDim*2), 
     &                 Stress(nDim*2), wStrain(48), wStress(48),
     &                 nStrain(28), nStress(28), sigm, tem, prop(3),
     &                 IncStrain(nDim*2),  defG(3,3),
     &                 defG0(3,3), xCurIP(nDim), TemperIP, 
     &                 TemperIPB, StressTh(nDim*2), MatProp(5),
     &                 StrainPl(nDim*2), StrainCr(nDim*2), 
     &                 StrainTh(nDim*2), StrainSw, StressBk(nDim*2),
     &                 MatRotGlb(3,3), wStrainTh(48), wStrainPl(48),
     &                 wStrainCr(48), eMassb(nNodes,nNodes), EnergyD(3)
     &                 ,yourLvalue(4),yourDvalue(4)
      CHARACTER*4      label(3)

c --- temporary debug key
      INTEGER debug, ix


      
      INTEGER          keyop , ext_i, ext_j, indd, fi
      real             x1,x2,x3,x4,y1,y2,y3,y4,E,v,det
      DOUBLE PRECISION, dimension(8) :: TheDOFS
      DOUBLE PRECISION, dimension(8) :: THELOADS
      DOUBLE PRECISION ::nooo,yesss,t
      real, dimension(4) ::g_xi
      real, dimension(4) ::g_eta
      DOUBLE PRECISION , dimension (3,8) :: BMatrix
      real :: N1,N2,N3,N4
      real :: x_g,y_g !gauss points physical coordinates
      real :: g_xi_scalar,g_eta_scalar

      DOUBLE PRECISION :: p
      DOUBLE PRECISION, dimension(4,4) :: Extrapolation_matrix

      DOUBLE PRECISION, dimension(4) :: n1_strains
      DOUBLE PRECISION, dimension(4) :: n2_strains
      DOUBLE PRECISION, dimension(4) :: n3_strains
      DOUBLE PRECISION, dimension(4) :: n4_strains
      
      DOUBLE PRECISION, dimension(4) :: n1_stresses
      DOUBLE PRECISION, dimension(4) :: n2_stresses
      DOUBLE PRECISION, dimension(4) :: n3_stresses
      DOUBLE PRECISION, dimension(4) :: n4_stresses
      DOUBLE PRECISION xi_ext, eta_ext,  xiii, etaaa
      
      DOUBLE PRECISION, dimension(3,3) :: D
      DOUBLE PRECISION, dimension(3) :: Strain_mod
      DOUBLE PRECISION, dimension(3) :: Stress_mod

      DOUBLE PRECISION, dimension(8) ::  realdofs
      DOUBLE PRECISION, dimension(8) ::  THELOADS1
      DOUBLE PRECISION, dimension(8) ::  THELOADS2
      DOUBLE PRECISION, dimension(3) :: YourStrain
      DOUBLE PRECISION, dimension(4,4) :: M
      DOUBLE PRECISION, dimension(4,3) :: g_stress_tensor
      DOUBLE PRECISION, dimension(4,3) :: g_strain_tensor
      DOUBLE PRECISION, dimension(4,3) :: n_stress_tensor
      DOUBLE PRECISION, dimension(4,3) :: n_strain_tensor
      
      DOUBLE PRECISION, dimension(4) :: copy_n2_stresses
      DOUBLE PRECISION, dimension(4) :: copy_n4_stresses
      DOUBLE PRECISION, dimension(4) :: copy_n2_strains
      DOUBLE PRECISION, dimension(4) :: copy_n4_strains
      
c --- B E G I N   E X E C U T A B L E  C O D I N G
c - - -                                                                                       P    R    E           P    R    O    C    E    S    S    I    N    G
c --- initialization
      nlgeom=0
      nTens = nDim*2
      nComp = nDim*nDim
      nDirect = 3
      nUsrDof2 = nUsrDof*nUsrDof
      CALL vzero (BMat(1,1),nUsrDof*nTens)
      IF (keyMtx(1).EQ.1) CALL vzero (eStiff(1,1),nUsrDof2)
      IF (keyMtx(2).EQ.1) CALL vzero (eMass(1,1) ,nUsrDof2)
      IF (keyMtx(5).EQ.1) CALL vzero (fExt(1)    ,nUsrDof)
      IF (keyMtx(6).EQ.1) CALL vzero (fInt(1)    ,nUsrDof)
      IF (nlgeom.EQ.0) THEN
         DO iDim = 1, 3
            DO iDim1 = 1, 3
               defG0(iDim, iDim1) = 0.0D0
            END DO
            defG0(iDim, iDim) = 1.0D0
         END DO
         CALL vmove (defG0(1,1),defG(1,1),9)
      ELSE
c        Nonlinear logic not defined here
      END IF
      elVol = 0.d0
      elMass= 0.d0
      IF (keyMtx(2).EQ.1) CALL vzero (eMassb(1,1),nNodes*nNodes)
      CALL vzero (elEnergy(1), nElEng)

c --- debug
      debug = 0
      IF (debug.EQ.1) THEN
          write (*,*)
          write (*,*)'elId =',elId, 
     &    ' keyAnsMat=',keyAnsMat,' nDim=',nDim, ' kEStress=',
     &    kEStress, ' keySym=',keySym, ' nIntPnts=',nIntPnts
      END IF

c --- Define constitutive tensors for plane stress and plane strain

      IF (keyAnsMat.EQ.1) THEN
c        --- from the MP/MPDATA/TB/TBDATA command information
      ELSE

c        ---- from the R/RMORE command information (METHOD(3))
c        ---- can also be done at the integration points
         IF (nReal.GE.1 ) Ex     = RealConst(1)
         IF (nReal.GE.2) nu      = RealConst(2)
         IF (nReal.GE.3 ) t      = RealConst(3)
c --- Plane stress
         IF (KeyOpt(1).EQ.0) THEN
c .AND. Ex     .EQ.0.0d0    
             D(1,1)=1.0D0
             D(1,2)=nu
             D(1,3)=0.0D0

             D(2,1)=nu
             D(2,2)=1.0D0
             D(2,3)=0.0D0
 
             D(3,1)=0.0D0
             D(3,2)=0.0D0
             D(3,3)=(1-nu)/2.0D0
             D= Ex/(1-(nu*nu)) * D

c --- plane strain
         Else IF (KeyOpt(1).EQ.1) THEN
             D(1,1)=1.0D0-nu
             D(1,2)=nu
             D(1,3)=0.0D0

             D(2,1)=nu
             D(2,2)=1.0D0-nu
             D(2,3)=0.0D0

             D(3,1)=0.0D0
             D(3,2)=0.0D0
             D(3,3)=(1.0D0-2.0D0*nu)/2.0D0
             D=(Ex/((1.0D0+nu)*(1.0D0-2.0D0*nu))) * D
         END IF
 
      END IF
      kTherm = 0

c --- start integration loop >> gauss points coordinates in natural space




      g_xi(1)=-1.0D0/sqrt(3.0D0)  
      g_xi(2)=-1.0D0/sqrt(3.0D0)
      g_xi(3)=1.0D0/sqrt(3.0D0)
      g_xi(4)=1.0D0/sqrt(3.0D0)
      g_eta(1)=-1.0D0/sqrt(3.0D0)
      g_eta(2)=1.0D0/sqrt(3.0D0)
      g_eta(3)=-1.0D0/sqrt(3.0D0)
      g_eta(4)=1.0D0/sqrt(3.0D0)

c Get physical coordinates
      x1=xCur(1,1)
      y1=xCur(2,1)
      x2=xCur(1,2)
      y2=xCur(2,2)
      x3=xCur(1,3)
      y3=xCur(2,3)
      x4=xCur(1,4)
      y4=xCur(2,4)

c --- create stiffness matrix


C                                                                        I     M     P     O     R     T     A     N     T                                        N     O     T     E




C  THE LOAD VECTOR IS CALCULATED AT THE END BECAUSE IF IT IS CALCULATED HERE THE TOTVALDOFS MAY CHANGE THROUGH SOULTOION ITERATIONS, SO YOU MAY ENTER THE SAME INPUT FILE AND EACH TIME YOU RUN YOU GET A DIFFERENT SOLUTIONS
C FOR THE NODAL FORCE VECTOR


    


c --- Get nodal force and dispalcements vectors
C --- Get element's nodes numberes and remap them to 1 2 3 4 to be able to get thier constraint displacements and loads      
      
      
c --- for debugding -Get predefined displacements and loads and REACTION FORCES
      


      Do yourNode=1,nNodes
          DO yourLabel=1,nDim
c ---         Ux=1, Uy=2    Fx=1,   Fy=2
              real_node=Nodes(yourNode)
              dis=disget(real_node,yourLabel,yourDvalue)
              forg=forget(real_node,yourLabel,yourLvalue)
              nooo=dspget(real_node,1,yourLabel,yesss)
              realdofs((yourNode-1)*2+yourLabel)=yesss
c             IF dis =4 OR -4 , THERE IS A CONSTRAINT , 0 no constraint
c             print*, forg
              IF (dis.NE.0) Then
                  THEDOFS((yourNode-1)*2+yourLabel) = yourDvalue(1)
              END IF
c             IF forg =4 , THERE IS A CONSTRAINT , 0 no constrint
              IF (forg.NE.0) Then
                  THELOADS((yourNode-1)*2+yourLabel) = yourLvalue(1)
              END IF           
          END DO
      END Do
 
      print*, 'PRE DEFINED LOADS'
      print*, THEDOFS
    
      
      print*, 'PRE DEFINED DOFS'
      print*, THEDOFS

      
      print*, 'coordinatesare'
      print*, x1,x2,x3,x4,y1,y2,y3,y4     

    
    
c--- Nodal force vector 

c    NOTE: the TotalDofs here is a PURE INPUT by ansys , and the load and STRAIN/STRESS vectors are calculated using it with the defined stiffness and B matrix



    
      
c - - -                                                                                       P    O    S    T                 P    R    O    C    E    S    S    I    N    G
      DO 100 intPnt = 1, nIntPnts


c --- create B matrix
          g_xi_scalar=g_xi(intPnt)
          g_eta_scalar=g_eta(intPnt)
          Call Bmatrixx(g_xi_scalar,g_eta_scalar,x1,x2,x3,x4,y1,y2,
     &         y3,y4,BMatrix,det)

c --- calculate strains and stress

          CALL maxv(BMatrix, TOTValDofs, Strain_mod, 3,8)

          IF (keyAnsMat.EQ.1) THEN

             CALL ElemGetMat (elId, matId, nDim, nTens, nDirect,
     &                         intPnt, xCurIP(1), TemperIP,
     &                         TemperIPB, kThermIP, IncStrain(1),
     &                         defG0(1,1), defG(1,1),
     &                         cMat(1,1), MatProp(1), Stress(1), 
     &                         Strain(1), StressTh(1), StrainTh(1),
     &                         StrainPl(1), StrainCr(1), 
     &                         StressBk(1), StrainSw, EnergyD(1),
     &                         MatRotGlb(1,1))
             if (kThermIP .eq. 1) kTherm = 1
             density = MatProp(4)
             nu = MatProp(5)
          ELSE


             Strain_mod(3)=Strain_mod(3)/2
             Strain(1)=Strain_mod(1)
             Strain(2)=Strain_mod(2)
             Strain(4)=Strain_mod(3)
c--- Strains at gauss points              
             g_strain_tensor(intpnt,1)=Strain_mod(1)
             g_strain_tensor(intpnt,2)=Strain_mod(2)
             g_strain_tensor(intpnt,3)=Strain_mod(3)

  
             CALL maxv(D, Strain_mod, Stress_mod, 3, 
     &              3)
             
             Stress(1)=Stress_mod(1)
             Stress(2)=Stress_mod(2)
             Stress(4)=Stress_mod(3)
c--- Stresses at gauss points              
             g_stress_tensor(intpnt,1)=Stress_mod(1)
             g_stress_tensor(intpnt,2)=Stress_mod(2)
             g_stress_tensor(intpnt,3)=Stress_mod(3)
             
             
c           --- don't deal with thermal strains or stresses in this example
            kThermIP = 0
            call vzero (StressTh(1),nTens)
         END IF

c --- prepare to create mass matrix (avoided cause density =0)

         IF (keyMtx(2).EQ.1) THEN
            IF (density.NE.0.0d0) THEN
               workDb = density*dVol
               CALL maat (shIso(1),eMassb(1,1),nNodes,nNodes,workDb)
            ENDIF
         ENDIF

c --- create external force vector (avoided cause kthermal ip = 0)

         IF (keyMtx(5).EQ.1 .AND. kThermIP.EQ.1 .AND. outkey.EQ.0) THEN
            CALL vmult (StressTh(1), workArr(1), nTens, dVol)
            CALL matxb (BMatrix, workArr(1), fExt(1), nTens, nTens,
     &                  nUsrDof, nUsrDof, 1, -nTens)
         END IF

c --- calculate other element quantities

         elVol = elVol+dVol
         elMass = elMass+dVol*density
         IF (keyAnsMat.EQ.0) elEnergy(1) = elEnergy(1)
     &                 + 0.5d0*dVol*vdot(Strain(1), Stress(1), nTens)
         k1 = (intPnt-1)*nTens+1
c         CALL vmove (Stress(1), saveVars(k1), nTens)
         IF (outkey.EQ.1) THEN
            CALL vmove (Strain(1), wStrain(k1), nTens)
            CALL vmove (Stress(1), wStress(k1), nTens)
            IF (keyAnsMat.EQ.1) THEN
               CALL vmove (StrainTh(1), wStrainTh(k1), nTens)
               CALL vmove (StrainPl(1), wStrainPl(k1), nTens)
               CALL vmove (StrainCr(1), wStrainCr(k1), nTens)
            END IF
            IF (debug.EQ.1) THEN
               write (*,3010) intPnt, (Strain(ix),ix=1,nTens)
               write (*,3020) (Stress(ix),ix=1,nTens)
               write (*,3030) (StrainPl(ix),ix=1,nTens)
 3010          FORMAT (/1x, 'intPnt=',i2, 'Strain=',6(e15.8,2x))
 3020          FORMAT (1x, 8x, 'Stress=',6(e15.8,2x))
 3030          FORMAT (1x, 8x, 'StrainPl=',6(e15.8,2x))
            END IF
         END IF
 100  CONTINUE
      
      
                      
      CALL stiffer(x1,x2,x3,x4,y1,y2,y3,y4,t,2,
     &                  Ex,nu,eStiff,KeyOpt(1))

      
      fInt= matmul(eStiff, TOTValDofs) 
      
      print*, 'Nodal loads vector'
      print*, fInt


c --- calculate strains and stresses at the elemnt's nodes by extrapolating and output
c       to result files

      IF (outkey.EQ.1) THEN

          nNodesCorner = nNodes
          
c --- Create Extrapolation matrix 
          M(1,1) = (sqrt(3.0D0)/4.0D0+1.0D0/4.0D0)*(sqrt(3.0D0)+1.0D0)
          M(1,2) = -(sqrt(3.0D0)/4.0D0-1.0D0/4.0D0)*(sqrt(3.0D0)+1.0D0)
          M(1,3) = (sqrt(3.0D0)/4.0D0-1.0D0/4.0D0)*(sqrt(3.0D0)-1.0D0)
          M(1,4) = -(sqrt(3.0D0)/4.0D0+1.0D0/4.0D0)*(sqrt(3.0D0)-1.0D0)
          M(2,1) = -(sqrt(3.0D0)/4.0D0-1.0D0/4.0D0)*(sqrt(3.0D0)+1.0D0)
          M(2,2) = (sqrt(3.0D0)/4.0D0+1.0D0/4.0D0)*(sqrt(3.0D0)+1.0D0)
          M(2,3) = -(sqrt(3.0D0)/4.0D0+1.0D0/4.0D0)*(sqrt(3.0D0)-1.0D0)
          M(2,4) = (sqrt(3.0D0)/4.0D0-1.0D0/4.0D0)*(sqrt(3.0D0)-1.0D0)
          M(3,1) = (sqrt(3.0D0)/4.0D0-1.0D0/4.0D0)*(sqrt(3.0D0)-1.0D0)
          M(3,2) = -(sqrt(3.0D0)/4.0D0+1.0D0/4.0D0)*(sqrt(3.0D0)-1.0D0)
          M(3,3) = (sqrt(3.0D0)/4.0D0+1.0D0/4.0D0)*(sqrt(3.0D0)+1.0D0)
          M(3,4) = -(sqrt(3.0D0)/4.0D0-1.0D0/4.0D0)*(sqrt(3.0D0)+1.0D0)
          M(4,1) = -(sqrt(3.0D0)/4.0D0+1.0D0/4.0D0)*(sqrt(3.0D0)-1.0D0)
          M(4,2) = (sqrt(3.0D0)/4.0D0-1.0D0/4.0D0)*(sqrt(3.0D0)-1.0D0)
          M(4,3) = -(sqrt(3.0D0)/4.0D0-1.0D0/4.0D0)*(sqrt(3.0D0)+1.0D0)
          M(4,4) = (sqrt(3.0D0)/4.0D0+1.0D0/4.0D0)*(sqrt(3.0D0)+1.0D0)


c --- Extrapolate 
          
          
          n_strain_tensor=matmul(M,g_strain_tensor)
          n_stress_tensor=matmul(M,g_stress_tensor)

        
C get vectors for each node that contains SX,SY,SXY
          Do fi=1,3
              
              n1_stresses(fi)=n_stress_tensor(1,fi)
              n2_stresses(fi)=n_stress_tensor(2,fi)
              n3_stresses(fi)=n_stress_tensor(3,fi)
              n4_stresses(fi)=n_stress_tensor(4,fi)
          
              n1_strains(fi)=n_strain_tensor(1,fi)
              n2_strains(fi)=n_strain_tensor(2,fi)
              n3_strains(fi)=n_strain_tensor(3,fi)
              n4_strains(fi)=n_strain_tensor(4,fi)
          END DO
c move the 3rd component (SXY) and (EXY) to the 4th position          
          
          n1_stresses(4)=n1_stresses(3)
          n2_stresses(4)=n2_stresses(3)
          n3_stresses(4)=n3_stresses(3)
          n4_stresses(4)=n4_stresses(3)
          
          n1_strains(4)=n1_strains(3)
          n2_strains(4)=n2_strains(3)
          n3_strains(4)=n3_strains(3)
          n4_strains(4)=n4_strains(3)
          
c Calculate the 3rd components (SZ) and (EZ) according to the constitutive tensor (Plane stress/ Plane strain)
          
          IF (KeyOpt(1).EQ.0) THEN !plane stress
              
              n1_strains(3)=(-nu/(1-nu))*(n1_strains(1)+n1_strains(2))  
              n1_stresses(3)=0
              
              n2_strains(3)=(-nu/(1-nu))*(n2_strains(1)+n2_strains(2))  
              n2_stresses(3)=0
              
              n3_strains(3)=(-nu/(1-nu))*(n3_strains(1)+n3_strains(2))  
              n3_stresses(3)=0
              
              n4_strains(3)=(-nu/(1-nu))*(n4_strains(1)+n4_strains(2))  
              n4_stresses(3)=0
              
          Else if (KeyOpt(1).EQ.1) THEN !plane strain
              n1_strains(3)=0
              n1_stresses(3)=nu*(n1_stresses(1)+n1_stresses(2))
              
              n2_strains(3)=0
              n2_stresses(3)=nu*(n2_stresses(1)+n2_stresses(2))
              
              n3_strains(3)=0
              n3_stresses(3)=nu*(n3_stresses(1)+n3_stresses(2))
              
              n4_strains(3)=0
              n4_stresses(3)=nu*(n4_stresses(1)+n4_stresses(2))
          End if
          
c         Convert to a COUNTERCLOCKWISE ORDER TO FOLLOW ANSYS
          copy_n2_stresses=n2_stresses
          copy_n4_stresses=n4_stresses
          
          copy_n2_strains=n2_strains
          copy_n4_strains=n4_strains
          
          n2_stresses=copy_n4_stresses
          n4_stresses=copy_n2_stresses
          
          n2_strains=copy_n4_strains
          n4_strains=copy_n2_strains
          
          
          
          CALL vmove (n1_stresses, saveVars(1), 4)
          CALL vmove (n2_stresses, saveVars(5), 4)
          CALL vmove (n3_stresses, saveVars(9), 4)
          CALL vmove (n4_stresses, saveVars(13), 4)

        
          nStress([1:4])=n1_stresses
          nStress([5:8])=n2_stresses
          nStress([9:12])=n3_stresses
          nStress([13:16])=n4_stresses
   
          nStrain([1:4])=n1_strains
          nStrain([5:8])=n2_strains
          nStrain([9:12])=n3_strains
          nStrain([13:16])=n4_strains
          
          
          

          


c --- only calculate basic result variables when it is allowed
c        and necessary

          IF(nRsltBsc.GT.0) THEN
                           
              DO iNode = 1, nNodesCorner

                  k1 = (iNode-1)*nTens + 1
                  k2 = (iNode-1)*7 + 1
                  CALL vmove (nStress(k1), RsltBsc(k2), nTens)
                  sigm = (nStress(k1)+nStress(k1+1)+nStress(k1+2))/3.0d0

                  RsltBsc(k2+4) = 0.0d0
                  RsltBsc(k2+5) = 0.0d0
               
c --- Equivalent stress               
                  RsltBsc(k2+6) = SQRT(1.5d0*
     &                 ( (nStress(k1)-sigm)*(nStress(k1)-sigm)
     &                 + (nStress(k1+1)-sigm)*(nStress(k1+1)-sigm)
     &                 + (nStress(k1+2)-sigm)*(nStress(k1+2)-sigm)
     &                 + 2.0d0*nStress(k1+3)*nStress(k1+3)))

                  k2 = (nNodesCorner+iNode-1)*7 + 1
                  CALL vmove (nStrain(k1), RsltBsc(k2), nTens)
               
                  RsltBsc(k2+4) = 0.0d0
                  RsltBsc(k2+5) = 0.0d0

c --- Equivalent strain                
                  CALL equivStrain (nu, nStrain(k1), nTens, 
     &                            RsltBsc(k2+6))
              END DO    
          END IF
         
c          k1 = nNodesCorner*nTens
c          CALL vmove (nStrain(1), RsltVar(1), k1)
c          CALL vmove (nStress(1), RsltVar(k1+1), k1)
         
         

         
         IF (elPrint .EQ. 1) THEN
c --- print out the results in OUT file (requested by the OUTPR command)
            CALL pplock (LOCKOT)
            WRITE (iott,2000) elId
 2000       FORMAT (/1x, 'Material Point output for element',I8)
            WRITE (iott,2100)
 2100       FORMAT(/4x, 'Intg.Pt. "S"     Stresses')
            DO intPnt = 1, nIntPnts
               k1 = (intPnt-1)*nTens
               WRITE (iott, 2110) intPnt, (wStress(k1+k2),k2=1,nTens)
 2110          FORMAT (4x,I4, 4x, 6(E12.5,1x))
            END DO
            WRITE (iott,2200)
 2200       FORMAT(/4x, 'Intg.Pt. "EPTO"     Strains')
            DO intPnt = 1, nIntPnts
               k1 = (intPnt-1)*nTens
               WRITE (iott, 2110) intPnt, (wStrain(k1+k2),k2=1,nTens)
            END DO
            
            
            
            IF (keyAnsMat.EQ.1) THEN
               workDb = 0.0d0
               DO k1 = 1, nIntPnts*nTens
                  workDb = workDb + ABS(wStrainPl(k1))
               END DO
               IF (workDb.GT.0.0d0) THEN
                  WRITE (iott,2300)
 2300             FORMAT(/4x, 'Intg.Pt. "EPPL"     Strains')
                  DO intPnt = 1, nIntPnts
                     k1 = (intPnt-1)*nTens
                     WRITE (iott, 2110) intPnt, (wStrainPl(k1+k2),
     &                                  k2=1,nTens)
                  END DO
               END IF
               workDb = 0.0d0
               DO k1 = 1, nIntPnts*nTens
                  workDb = workDb + ABS(wStrainCr(k1))
               END DO
               IF (workDb.GT.0.0d0) THEN
                  WRITE (iott,2400)
 2400             FORMAT(/4x, 'Intg.Pt. "EPCR"     Strains')
                  DO intPnt = 1, nIntPnts
                     k1 = (intPnt-1)*nTens
                     WRITE (iott, 2110) intPnt, (wStrainCr(k1+k2),
     &                                  k2=1,nTens)
                  END DO
               END IF
               workDb = 0.0d0
               DO k1 = 1, nIntPnts*nTens
                  workDb = workDb + ABS(wStrainTh(k1))
               END DO
               IF (workDb.GT.1.0d-12) THEN
                  WRITE (iott,2500)
 2500             FORMAT(/4x, 'Intg.Pt. "EPTH"     Strains')
                  DO intPnt = 1, nIntPnts
                     k1 = (intPnt-1)*nTens
                     WRITE (iott, 2110) intPnt, (wStrainTh(k1+k2),
     &                                  k2=1,nTens)
                  END DO
               END IF
               write (iott,3000)
 3000          format(2/)
            END IF
            CALL ppunlock (LOCKOT)
         END IF
      END IF
     
 990  CONTINUE
 
      RETURN
      END

