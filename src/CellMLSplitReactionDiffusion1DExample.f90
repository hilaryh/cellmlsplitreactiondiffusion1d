!> \file
!> $Id: DiffusionExample.f90 1528 2010-09-21 01:32:29Z chrispbradley $
!> \author Chris Bradley
!> \brief This is an example program to solve a reaction diffusion equation using OpenCMISS calls.
!>
!> \section LICENSE
!>
!> Version: MPL 1.1/GPL 2.0/LGPL 2.1
!>
!> The contents of this file are subject to the Mozilla Public License
!> Version 1.1 (the "License"); you may not use this file except in
!> compliance with the License. You may obtain a copy of the License at
!> http://www.mozilla.org/MPL/
!>
!> Software distributed under the License is distributed on an "AS IS"
!> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
!> License for the specific language governing rights and limitations
!> under the License.
!>
!> The Original Code is OpenCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand and University of Oxford, Oxford, United
!> Kingdom. Portions created by the University of Auckland and University
!> of Oxford are Copyright (C) 2007 by the University of Auckland and
!> the University of Oxford. All Rights Reserved.
!>
!> Contributor(s):
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. If you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. If you do not delete
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>

!> \example ClassicalField/ReactionDiffusion/CellMLSplitReactionDiffusion1D/src/CellMLSplitReactionDiffusion1DExample.f90
!! Example program to solve a diffusion equation using OpenCMISS calls.
!! \htmlinclude ClassicalField/ReactionDiffusion/CellMLSplitReactionDiffusion1D/history.html
!<

!> Main program
PROGRAM CELLMLSPLITREACTIONDIFFUSION1DEXAMPLE

  USE OpenCMISS
  USE OpenCMISS_Iron
#ifndef NOMPIMOD
  USE MPI
#endif

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif


  !Test program parameters
  !you can change this parameter to change the length of the line in microns
  REAL(CMISSRP), PARAMETER :: LENGTH=10.0_CMISSRP
  
  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: DependentFieldUserNumber=9
  INTEGER(CMISSIntg), PARAMETER :: MaterialsFieldUserNumber=10
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetUserNumber=11
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=12
  INTEGER(CMISSIntg), PARAMETER :: SourceFieldUserNumber=13
  INTEGER(CMISSIntg), PARAMETER :: CellMLUserNumber=14
  INTEGER(CMISSIntg), PARAMETER :: CellMLModelsFieldUserNumber=15
  INTEGER(CMISSIntg), PARAMETER :: CellMLStateFieldUserNumber=16
  INTEGER(CMISSIntg), PARAMETER :: CellMLIntermediateFieldUserNumber=17
  INTEGER(CMISSIntg), PARAMETER :: CellMLParametersFieldUserNumber=18
  INTEGER(CMISSIntg), PARAMETER :: CompartmentTypeFieldUserNumber=19
  INTEGER(CMISSIntg), PARAMETER :: NaiSourceFieldUserNumber=20

INTEGER(CMISSIntg), PARAMETER :: FluoFieldUserNumber=21
INTEGER(CMISSIntg), PARAMETER :: FluoMaterialsFieldUserNumber=22
INTEGER(CMISSIntg), PARAMETER :: FluoEquationsSetUserNumber=23
INTEGER(CMISSIntg), PARAMETER :: FluoEquationsSetFieldUserNumber=24
INTEGER(CMISSIntg), PARAMETER :: FluoSourceFieldUserNumber=25

INTEGER(CMISSIntg), PARAMETER :: CaFluoFieldUserNumber=26
INTEGER(CMISSIntg), PARAMETER :: CaFluoMaterialsFieldUserNumber=27
INTEGER(CMISSIntg), PARAMETER :: CaFluoEquationsSetUserNumber=28
INTEGER(CMISSIntg), PARAMETER :: CaFluoEquationsSetFieldUserNumber=29
INTEGER(CMISSIntg), PARAMETER :: CaFluoSourceFieldUserNumber=30

INTEGER(CMISSIntg), PARAMETER :: mBufFieldUserNumber=31
INTEGER(CMISSIntg), PARAMETER :: mBufMaterialsFieldUserNumber=32
INTEGER(CMISSIntg), PARAMETER :: mBufEquationsSetUserNumber=33
INTEGER(CMISSIntg), PARAMETER :: mBufEquationsSetFieldUserNumber=34
INTEGER(CMISSIntg), PARAMETER :: mBufSourceFieldUserNumber=35

INTEGER(CMISSIntg), PARAMETER :: CamBufFieldUserNumber=36
INTEGER(CMISSIntg), PARAMETER :: CamBufMaterialsFieldUserNumber=37
INTEGER(CMISSIntg), PARAMETER :: CamBufEquationsSetUserNumber=38
INTEGER(CMISSIntg), PARAMETER :: CamBufEquationsSetFieldUserNumber=39
INTEGER(CMISSIntg), PARAMETER :: CamBufSourceFieldUserNumber=40

  !Program types
  
  !Program variables

  INTEGER(CMISSIntg) :: NUMBER_GLOBAL_X_ELEMENTS,CONDITION,ElemMap
  INTEGER(CMISSIntg) :: NUMBER_OF_DOMAINS,NODE_NUMBER,NodeDomain
  INTEGER(CMISSIntg),DIMENSION(2) :: BCNODES
  INTEGER(CMISSIntg) :: MPI_IERROR
  INTEGER :: node, GeometricMeshComponent
  REAL(CMISSRP) :: VALUE, VALUEa, VALUEb, VALUEc, VALUEd, VALUEe
  INTEGER(CMISSIntg) :: constantModelIndex

  !INTEGER(INTG) :: first_global_dof,first_local_dof,first_local_rank,last_global_dof,last_local_dof,last_local_rank,rank_idx
  !INTEGER(INTG) :: EQUATIONS_SET_INDEX
  !TYPE(DOMAIN_MAPPING_TYPE), POINTER :: DEPENDENT_DOF_MAPPING
  
  !CMISS variables

  TYPE(cmfe_BasisType) :: Basis
  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(cmfe_DecompositionType) :: Decomposition
  TYPE(cmfe_EquationsType) :: Equations
  TYPE(cmfe_EquationsSetType) :: EquationsSet
  TYPE(cmfe_FieldType) :: GeometricField,EquationsSetField,DependentField,MaterialsField,SourceField,NaiSourceField
  TYPE(cmfe_FieldsType) :: Fields
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditions
  TYPE(cmfe_GeneratedMeshType) :: GeneratedMesh  
  TYPE(cmfe_MeshType) :: Mesh
  TYPE(cmfe_ProblemType) :: Problem
  TYPE(cmfe_ControlLoopType) :: ControlLoop
  TYPE(cmfe_RegionType) :: Region,WorldRegion
  TYPE(cmfe_SolverType) :: Solver, LinearSolver
  TYPE(cmfe_SolverEquationsType) :: SolverEquations
  TYPE(cmfe_CellMLType) :: CellML
  TYPE(cmfe_CellMLEquationsType) :: CellMLEquations
  TYPE(cmfe_FieldType) :: CellMLModelsField,CellMLStateField,CellMLIntermediateField,CellMLParametersField
  TYPE(cmfe_FieldType) :: CompartmentTypeField

  ! Define buffer variables
  TYPE(cmfe_FieldType) :: FluoField,FluoMaterialsField,FluoEquationsSetField,FluoSourceField
  TYPE(cmfe_FieldType) :: CaFluoField,CaFluoMaterialsField,CaFluoEquationsSetField,CaFluoSourceField
  TYPE(cmfe_FieldType) :: mBufField,mBufMaterialsField,mBufEquationsSetField,mBufSourceField
  TYPE(cmfe_FieldType) :: CamBufField,CamBufMaterialsField,CamBufEquationsSetField,CamBufSourceField

  TYPE(cmfe_EquationsType) :: FluoEquations, CaFluoEquations,mBufEquations, CamBufEquations
  TYPE(cmfe_EquationsSetType) :: FluoEquationsSet, CaFluoEquationsSet,mBufEquationsSet, CamBufEquationsSet


  LOGICAL :: EXPORT_FIELD

#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Generic CMISS variables
  
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber
  INTEGER(CMISSIntg) :: EquationsSetIndex,CellMLIndex
  INTEGER(CMISSIntg) :: Err

  
#ifdef WIN32
  !Initialise QuickWin
  QUICKWIN_WINDOW_CONFIG%TITLE="General Output" !Window title
  QUICKWIN_WINDOW_CONFIG%NUMTEXTROWS=-1 !Max possible number of rows
  QUICKWIN_WINDOW_CONFIG%MODE=QWIN$SCROLLDOWN
  !Set the window parameters
  QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
  !If attempt fails set with system estimated values
  IF(.NOT.QUICKWIN_STATUS) QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
#endif

  !Intialise OpenCMISS
  CALL cmfe_Initialise(WorldCoordinateSystem,WorldRegion,Err)
  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,Err)
  !Get the computational nodes information
  CALL cmfe_ComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL cmfe_ComputationalNodeNumberGet(ComputationalNodeNumber,Err)

  !*******************************************************************************
  ! YOU CAN CHANGE THIS to change the refinement of your 1D line.
  ! Greg's data has a resolution of ~0.5 microns
  NUMBER_GLOBAL_X_ELEMENTS=20
  !*******************************************************************************

  NUMBER_OF_DOMAINS=NumberOfComputationalNodes

  !Set all diganostic levels on for testing

  CALL MPI_BCAST(NUMBER_GLOBAL_X_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)
  CALL MPI_BCAST(NUMBER_OF_DOMAINS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)

  !Start the creation of a new RC coordinate system
  CALL cmfe_CoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  !Set the coordinate system to be 1D
  CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystem,1,Err)
  !Finish the creation of the coordinate system
  CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystem,Err)
  

  !Start the creation of the region
  CALL cmfe_Region_Initialise(Region,Err)
  CALL cmfe_Region_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  !Set the regions coordinate system to the 1D RC coordinate system that we have created
  CALL cmfe_Region_CoordinateSystemSet(Region,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL cmfe_Region_CreateFinish(Region,Err)
  
  !Start the creation of a basis (default is linear lagrange)
  CALL cmfe_Basis_Initialise(Basis,Err)
  CALL cmfe_Basis_CreateStart(BasisUserNumber,Basis,Err)
  !Set the basis to be a linear Lagrange basis
  CALL cmfe_Basis_NumberOfXiSet(Basis,1,Err)
  !Finish the creation of the basis
  CALL cmfe_Basis_CreateFinish(Basis,Err)

  !Start the creation of a generated mesh in the region
  CALL cmfe_GeneratedMesh_Initialise(GeneratedMesh,Err)
  CALL cmfe_GeneratedMesh_CreateStart(GeneratedMeshUserNumber,Region,GeneratedMesh,Err)
  !Set up a regular 1D mesh
  CALL cmfe_GeneratedMesh_TypeSet(GeneratedMesh,CMFE_GENERATED_MESH_REGULAR_MESH_TYPE,Err)
  !Set the default basis
  CALL cmfe_GeneratedMesh_BasisSet(GeneratedMesh,Basis,Err)   
  !Define the mesh on the region
  CALL cmfe_GeneratedMesh_ExtentSet(GeneratedMesh,[LENGTH],Err)
  CALL cmfe_GeneratedMesh_NumberOfElementsSet(GeneratedMesh,[NUMBER_GLOBAL_X_ELEMENTS],Err)
  !Finish the creation of a generated mesh in the region
  CALL cmfe_Mesh_Initialise(Mesh,Err)
  CALL cmfe_GeneratedMesh_CreateFinish(GeneratedMesh,MeshUserNumber,Mesh,Err)
  
  !Create a decomposition
  CALL cmfe_Decomposition_Initialise(Decomposition,Err)
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL cmfe_Decomposition_TypeSet(Decomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,NUMBER_OF_DOMAINS,Err)
  !Finish the decomposition
  CALL cmfe_Decomposition_CreateFinish(Decomposition,Err)
  
  !Start to create a default (geometric) field on the region
  CALL cmfe_Field_Initialise(GeometricField,Err)
  CALL cmfe_Field_CreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  !Set the decomposition to use
  CALL cmfe_Field_MeshDecompositionSet(GeometricField,Decomposition,Err)
  !Set the domain to be used by the field components.
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,Err)
  !Finish creating the field
  CALL cmfe_Field_CreateFinish(GeometricField,Err)
  
  
  !Update the geometric field parameters
  CALL cmfe_GeneratedMesh_GeometricParametersCalculate(GeneratedMesh,GeometricField,Err)




!*************Dependent Field*******************
  !Create the cellml reaction with split reaction diffusion equations_set
  CALL cmfe_EquationsSet_Initialise(EquationsSet,Err)
  CALL cmfe_Field_Initialise(EquationsSetField,Err)
  CALL cmfe_EquationsSet_CreateStart(EquationsSetUserNumber,Region,GeometricField,[CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMFE_EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE,CMFE_EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE], &
    & EquationsSetFieldUserNumber,EquationsSetField,EquationsSet,Err)
  !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(EquationsSet,Err)

  !Create the equations set dependent field variables
  CALL cmfe_Field_Initialise(DependentField,Err)
  CALL cmfe_EquationsSet_DependentCreateStart(EquationsSet,DependentFieldUserNumber,DependentField,Err)
  CALL cmfe_Field_VariableLabelSet(DependentField,CMFE_FIELD_U_VARIABLE_TYPE,"Dependent Field",Err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateFinish(EquationsSet,Err)

  !Create the equations set material field variables
  !by default 2 comps for reac diff i.e. diff coeff in 1 direction set constant spatially = 1, and storage coeff set to 1
  CALL cmfe_Field_Initialise(MaterialsField,Err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(EquationsSet,MaterialsFieldUserNumber,MaterialsField,Err)
  CALL cmfe_Field_VariableLabelSet(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,"Materials Field",Err)
!Specify the interpolation
CALL cmfe_Field_ComponentInterpolationSet(MaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,1, &
& cmfe_FIELD_ELEMENT_BASED_INTERPOLATION,ERR)
  !Finish the equations set materials field variables
  CALL cmfe_EquationsSet_MaterialsCreateFinish(EquationsSet,Err)

! Buffer Fields

!################( Fluo )#############################################################--->

!Create Fluo  reaction diffusion with constant source equations_set
CALL cmfe_EquationsSet_Initialise(FluoEquationsSet,Err)
CALL cmfe_Field_Initialise(FluoEquationsSetField,Err)
CALL cmfe_EquationsSet_CreateStart(FluoEquationsSetUserNumber,Region,GeometricField,[CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
& CMFE_EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE,CMFE_EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE], &
& FluoEquationsSetFieldUserNumber,FluoEquationsSetField,FluoEquationsSet,Err)
CALL cmfe_EquationsSet_CreateFinish(FluoEquationsSet,Err)

CALL cmfe_Field_Initialise(FluoField,Err)
CALL cmfe_EquationsSet_DependentCreateStart(FluoEquationsSet,FluoFieldUserNumber,FluoField,Err)
CALL cmfe_Field_VariableLabelSet(FluoField,cmfe_FIELD_U_VARIABLE_TYPE,"Fluo Field",Err)
CALL cmfe_EquationsSet_DependentCreateFinish(FluoEquationsSet,Err)

CALL cmfe_Field_Initialise(FluoMaterialsField,Err)
CALL cmfe_EquationsSet_MaterialsCreateStart(FluoEquationsSet,FluoMaterialsFieldUserNumber,FluoMaterialsField,Err)
CALL cmfe_Field_VariableLabelSet(FluoMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,"Fluo Materials Field",Err)
CALL cmfe_EquationsSet_MaterialsCreateFinish(FluoEquationsSet,Err)

! Change second last variable to change diffusivity
CALL cmfe_Field_ComponentValuesInitialise(FluoMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
& cmfe_FIELD_VALUES_SET_TYPE,1,0.030_CMISSRP,Err)
CALL cmfe_Field_ComponentValuesInitialise(FluoMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
& 2,1.0_CMISSRP,Err) ! storage coefficient

!CALL cmfe_Field_ParameterSetUpdateStart(FluoMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
!CALL cmfe_Field_ParameterSetUpdateFinish(FluoMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)

!Set up source field for reaction diffusion equation set. Note that for the split problem subtype, the source field is not useall.
CALL cmfe_Field_Initialise(FluoSourceField,Err)
CALL cmfe_EquationsSet_SourceCreateStart(FluoEquationsSet,FluoSourceFieldUserNumber,FluoSourceField,Err)
CALL cmfe_Field_VariableLabelSet(FluoSourceField,cmfe_FIELD_U_VARIABLE_TYPE,"Fluo Source Field",Err)
!Finish the equations set source field variables
CALL cmfe_EquationsSet_SourceCreateFinish(FluoEquationsSet,Err)
CALL cmfe_Field_ComponentValuesInitialise(FluoSourceField,cmfe_FIELD_U_VARIABLE_TYPE, &
& cmfe_FIELD_VALUES_SET_TYPE,1,0.0_CMISSDP,Err)  !Set the equations set to be a standard Diffusion no source problem
!Finish creating the equations set

!################( CaFluo )#############################################################--->

!Create CaFluo reaction diffusion with split reaction diffusion equations_set
CALL cmfe_EquationsSet_Initialise(CaFluoEquationsSet,Err)
CALL cmfe_Field_Initialise(CaFluoEquationsSetField,Err)
CALL cmfe_EquationsSet_CreateStart(CaFluoEquationsSetUserNumber,Region,GeometricField,[CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
& CMFE_EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE,CMFE_EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE], &
& CaFluoEquationsSetFieldUserNumber,CaFluoEquationsSetField,CaFluoEquationsSet,Err)
! Finish creating the equations set
CALL cmfe_EquationsSet_CreateFinish(CaFluoEquationsSet,Err)
! Create the equations set dependent field variables
CALL cmfe_Field_Initialise(CaFluoField,Err)
CALL cmfe_EquationsSet_DependentCreateStart(CaFluoEquationsSet,CaFluoFieldUserNumber,CaFluoField,Err)
CALL cmfe_Field_VariableLabelSet(CaFluoField,cmfe_FIELD_U_VARIABLE_TYPE,"CaFluo Field",Err)
CALL cmfe_EquationsSet_DependentCreateFinish(CaFluoEquationsSet,Err)
! Create the equations set material field variables
CALL cmfe_Field_Initialise(CaFluoMaterialsField,Err)
CALL cmfe_EquationsSet_MaterialsCreateStart(CaFluoEquationsSet,CaFluoMaterialsFieldUserNumber,CaFluoMaterialsField,Err)
CALL cmfe_Field_VariableLabelSet(CaFluoMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,"CaFluo Materials Field",Err)
CALL cmfe_EquationsSet_MaterialsCreateFinish(CaFluoEquationsSet,Err)

! Change second last variable to change diffusivity
CALL cmfe_Field_ComponentValuesInitialise(CaFluoMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
& cmfe_FIELD_VALUES_SET_TYPE, &
& 1,0.030_CMISSRP,Err)
! Change the storage coefficient
CALL cmfe_Field_ComponentValuesInitialise(CaFluoMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
& 2,1.0_CMISSRP,Err) ! storage coefficient

!!!! Not in earlier version. Nec?
!CALL cmfe_Field_ParameterSetUpdateStart(CaFluoMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
!CALL cmfe_Field_ParameterSetUpdateFinish(CaFluoMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)

!Set up source field for reaction diffusion equation set. Note that for the split problem subtype, the source field is not useall.
CALL cmfe_Field_Initialise(CaFluoSourceField,Err)
CALL cmfe_EquationsSet_SourceCreateStart(CaFluoEquationsSet,CaFluoSourceFieldUserNumber,CaFluoSourceField,Err)
CALL cmfe_Field_VariableLabelSet(CaFluoSourceField,cmfe_FIELD_U_VARIABLE_TYPE,"CaFluo Source Field",Err)
!Finish the equations set source field variables
CALL cmfe_EquationsSet_SourceCreateFinish(CaFluoEquationsSet,Err)
CALL cmfe_Field_ComponentValuesInitialise(CaFluoSourceField,cmfe_FIELD_U_VARIABLE_TYPE, &
& cmfe_FIELD_VALUES_SET_TYPE,1,0.0_CMISSDP,Err)  !Set the equations set to be a standard Diffusion no source problem
!Finish creating the equations set

!################( mBuf )#############################################################--->

!Create mBuf  reaction diffusion with constant source equations_set
CALL cmfe_EquationsSet_Initialise(mBufEquationsSet,Err)
CALL cmfe_Field_Initialise(mBufEquationsSetField,Err)
CALL cmfe_EquationsSet_CreateStart(mBufEquationsSetUserNumber,Region,GeometricField,[CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
& CMFE_EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE,CMFE_EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE], &
& mBufEquationsSetFieldUserNumber,mBufEquationsSetField,mBufEquationsSet,Err)
CALL cmfe_EquationsSet_CreateFinish(mBufEquationsSet,Err)

CALL cmfe_Field_Initialise(mBufField,Err)
CALL cmfe_EquationsSet_DependentCreateStart(mBufEquationsSet,mBufFieldUserNumber,mBufField,Err)
CALL cmfe_Field_VariableLabelSet(mBufField,cmfe_FIELD_U_VARIABLE_TYPE,"mBuf Field",Err)
CALL cmfe_EquationsSet_DependentCreateFinish(mBufEquationsSet,Err)

CALL cmfe_Field_Initialise(mBufMaterialsField,Err)
CALL cmfe_EquationsSet_MaterialsCreateStart(mBufEquationsSet,mBufMaterialsFieldUserNumber,mBufMaterialsField,Err)
CALL cmfe_Field_VariableLabelSet(mBufMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,"mBuf Materials Field",Err)
CALL cmfe_EquationsSet_MaterialsCreateFinish(mBufEquationsSet,Err)

! Change second last variable to change diffusivity
CALL cmfe_Field_ComponentValuesInitialise(mBufMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
& cmfe_FIELD_VALUES_SET_TYPE,1,0.030_CMISSRP,Err)
CALL cmfe_Field_ComponentValuesInitialise(mBufMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
& 2,1.0_CMISSRP,Err) ! storage coefficient

!CALL cmfe_Field_ParameterSetUpdateStart(mBufMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
!CALL cmfe_Field_ParameterSetUpdateFinish(mBufMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)

!Set up source field for reaction diffusion equation set. Note that for the split problem subtype, the source field is not useall.
CALL cmfe_Field_Initialise(mBufSourceField,Err)
CALL cmfe_EquationsSet_SourceCreateStart(mBufEquationsSet,mBufSourceFieldUserNumber,mBufSourceField,Err)
CALL cmfe_Field_VariableLabelSet(mBufSourceField,cmfe_FIELD_U_VARIABLE_TYPE,"mBuf Source Field",Err)
!Finish the equations set source field variables
CALL cmfe_EquationsSet_SourceCreateFinish(mBufEquationsSet,Err)
CALL cmfe_Field_ComponentValuesInitialise(mBufSourceField,cmfe_FIELD_U_VARIABLE_TYPE, &
& cmfe_FIELD_VALUES_SET_TYPE,1,0.0_CMISSDP,Err)  !Set the equations set to be a standard Diffusion no source problem
!Finish creating the equations set

!################( CamBuf )#############################################################--->

!Create CamBuf  reaction diffusion with constant source equations_set
CALL cmfe_EquationsSet_Initialise(CamBufEquationsSet,Err)
CALL cmfe_Field_Initialise(CamBufEquationsSetField,Err)
CALL cmfe_EquationsSet_CreateStart(CamBufEquationsSetUserNumber,Region,GeometricField,[CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
& CMFE_EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE,CMFE_EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE], &
& CamBufEquationsSetFieldUserNumber,CamBufEquationsSetField,CamBufEquationsSet,Err)
CALL cmfe_EquationsSet_CreateFinish(CamBufEquationsSet,Err)

CALL cmfe_Field_Initialise(CamBufField,Err)
CALL cmfe_EquationsSet_DependentCreateStart(CamBufEquationsSet,CamBufFieldUserNumber,CamBufField,Err)
CALL cmfe_Field_VariableLabelSet(CamBufField,cmfe_FIELD_U_VARIABLE_TYPE,"CamBuf Field",Err)
CALL cmfe_EquationsSet_DependentCreateFinish(CamBufEquationsSet,Err)

CALL cmfe_Field_Initialise(CamBufMaterialsField,Err)
CALL cmfe_EquationsSet_MaterialsCreateStart(CamBufEquationsSet,CamBufMaterialsFieldUserNumber,CamBufMaterialsField,Err)
CALL cmfe_Field_VariableLabelSet(CamBufMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,"CamBuf Materials Field",Err)
CALL cmfe_EquationsSet_MaterialsCreateFinish(CamBufEquationsSet,Err)

! Change second last variable to change diffusivity
CALL cmfe_Field_ComponentValuesInitialise(CamBufMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE, &
& cmfe_FIELD_VALUES_SET_TYPE,1,0.030_CMISSRP,Err)
CALL cmfe_Field_ComponentValuesInitialise(CamBufMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
& 2,1.0_CMISSRP,Err) ! storage coefficient

!CALL cmfe_Field_ParameterSetUpdateStart(CamBufMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
!CALL cmfe_Field_ParameterSetUpdateFinish(CamBufMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)

!Set up source field for reaction diffusion equation set. Note that for the split problem subtype, the source field is not useall.
CALL cmfe_Field_Initialise(CamBufSourceField,Err)
CALL cmfe_EquationsSet_SourceCreateStart(CamBufEquationsSet,CamBufSourceFieldUserNumber,CamBufSourceField,Err)
CALL cmfe_Field_VariableLabelSet(CamBufSourceField,cmfe_FIELD_U_VARIABLE_TYPE,"CamBuf Source Field",Err)
!Finish the equations set source field variables
CALL cmfe_EquationsSet_SourceCreateFinish(CamBufEquationsSet,Err)
CALL cmfe_Field_ComponentValuesInitialise(CamBufSourceField,cmfe_FIELD_U_VARIABLE_TYPE, &
& cmfe_FIELD_VALUES_SET_TYPE,1,0.0_CMISSDP,Err)  !Set the equations set to be a standard Diffusion no source problem
!Finish creating the equations set


!*******************************************************************************
! YOU CAN CHANGE THIS to change the diffusivity

  CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
   & 1,0.53_CMISSRP,Err) !diff coeff in x
  ! Set diff at specific element i.e across nuclear envelope. nucpos
    ! But the nuclear envelope is ~45nm thick. 
    ! Currently, the nodes are spaced 500 nm apart.
    ElemMap = 13.0_CMISSDP
  CALL cmfe_Field_ParameterSetUpdateElement(MaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE, &
    & ElemMap,1,0.5_CMISSRP,Err)

! Prevent nuclear buffer diffusion across the nuclear envelope
  CALL cmfe_Field_ParameterSetUpdateElement(mBufMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE, &
    & ElemMap,1,0.0_CMISSRP,Err)
  CALL cmfe_Field_ParameterSetUpdateElement(CamBufMaterialsField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE, &
    & ElemMap,1,0.0_CMISSRP,Err)

!0.00067
!*******************************************************************************
  CALL cmfe_Field_ComponentValuesInitialise(MaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
   & 2,1.0_CMISSRP,Err) ! storage coefficient

!!! Added stuff
CALL cmfe_Field_Initialise(CompartmentTypeField,Err)
CALL cmfe_Field_CreateStart(CompartmentTypeFieldUserNumber,Region,CompartmentTypeField,Err)
CALL cmfe_Field_TypeSet(CompartmentTypeField,CMFE_FIELD_GENERAL_TYPE,Err)
CALL cmfe_Field_MeshDecompositionSet(CompartmentTypeField,Decomposition,Err)
CALL cmfe_Field_GeometricFieldSet(CompartmentTypeField,GeometricField,Err)
CALL cmfe_Field_NumberOfVariablesSet(CompartmentTypeField,1,Err)
CALL cmfe_Field_VariableTypesSet(CompartmentTypeField,[CMFE_FIELD_U_VARIABLE_TYPE],Err)
CALL cmfe_Field_DataTypeSet(CompartmentTypeField,CMFE_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_DP_TYPE,Err)
CALL cmfe_Field_DimensionSet(CompartmentTypeField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_SCALAR_DIMENSION_TYPE,Err)
CALL cmfe_Field_NumberOfComponentsSet(CompartmentTypeField,CMFE_FIELD_U_VARIABLE_TYPE,1,Err)
CALL cmfe_Field_VariableLabelSet(CompartmentTypeField,CMFE_FIELD_U_VARIABLE_TYPE,"CompartmentTypeField",Err)
CALL cmfe_Field_ComponentMeshComponentGet(GeometricField,cmfe_FIELD_U_VARIABLE_TYPE, &
& 1,GeometricMeshComponent,ERR)
!Default to the geometric interpolation setup
CALL cmfe_Field_ComponentMeshComponentSet(CompartmentTypeField,cmfe_FIELD_U_VARIABLE_TYPE,1, &
& GeometricMeshComponent,ERR)
!Specify the interpolation to be same as geometric interpolation
CALL cmfe_Field_ComponentInterpolationSet(CompartmentTypeField,cmfe_FIELD_U_VARIABLE_TYPE,1, &
& cmfe_FIELD_NODE_BASED_INTERPOLATION,ERR)

CALL cmfe_Field_CreateFinish(CompartmentTypeField,Err)
!Set all nodes to have value 1
CALL cmfe_Field_ComponentValuesInitialise(CompartmentTypeField,cmfe_FIELD_U_VARIABLE_TYPE, &
  & cmfe_FIELD_VALUES_SET_TYPE,1,1.0_CMISSDP,Err)

  !Set up source field for reaction diffusion equation set. Note that for the split problem subtype, the source field is not used at all.
  CALL cmfe_Field_Initialise(SourceField,Err)
  CALL cmfe_EquationsSet_SourceCreateStart(EquationsSet,SourceFieldUserNumber,SourceField,Err)
  CALL cmfe_Field_VariableLabelSet(SourceField,CMFE_FIELD_U_VARIABLE_TYPE,"Source Field",Err)
  !Finish the equations set source field variables
  CALL cmfe_EquationsSet_SourceCreateFinish(EquationsSet,Err)
  CALL cmfe_Field_ComponentValuesInitialise(SourceField,CMFE_FIELD_U_VARIABLE_TYPE, &
    & CMFE_FIELD_VALUES_SET_TYPE,1,0.0_CMISSRP,Err)

!CALL cmfe_Field_Initialise(NaiSourceField,Err)
!CALL cmfe_Field_CreateStart(NaiSourceFieldUserNumber,Region,NaiSourceField,Err)
!CALL cmfe_Field_TypeSet(NaiSourceField,CMFE_FIELD_GENERAL_TYPE,Err)
!CALL cmfe_Field_MeshDecompositionSet(NaiSourceField,Decomposition,Err)
!CALL cmfe_Field_GeometricFieldSet(NaiSourceField,GeometricField,Err)
!CALL cmfe_Field_NumberOfVariablesSet(NaiSourceField,1,Err)
!CALL cmfe_Field_VariableTypesSet(NaiSourceField,[CMFE_FIELD_U_VARIABLE_TYPE],Err)
!CALL cmfe_Field_DataTypeSet(NaiSourceField,CMFE_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_DP_TYPE,Err)
!CALL cmfe_Field_DimensionSet(NaiSourceField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_SCALAR_DIMENSION_TYPE,Err)
!CALL cmfe_Field_NumberOfComponentsSet(NaiSourceField,CMFE_FIELD_U_VARIABLE_TYPE,1,Err)
!CALL cmfe_Field_VariableLabelSet(NaiSourceField,CMFE_FIELD_U_VARIABLE_TYPE,"NaiField",Err)
!CALL cmfe_Field_ComponentMeshComponentGet(GeometricField,cmfe_FIELD_U_VARIABLE_TYPE, &
!& 1,GeometricMeshComponent,ERR)
!!Default to the geometric interpolation setup
!CALL cmfe_Field_ComponentMeshComponentSet(NaiSourceField,cmfe_FIELD_U_VARIABLE_TYPE,1, &
!& GeometricMeshComponent,ERR)
!!Specify the interpolation to be same as geometric interpolation
!CALL cmfe_Field_ComponentInterpolationSet(NaiSourceField,cmfe_FIELD_U_VARIABLE_TYPE,1, &
!& cmfe_FIELD_NODE_BASED_INTERPOLATION,ERR)
!
!CALL cmfe_Field_CreateFinish(NaiSourceField,Err)
!!Set all nodes to have value 10
!CALL cmfe_Field_ComponentValuesInitialise(NaiSourceField,cmfe_FIELD_U_VARIABLE_TYPE, &
!  & cmfe_FIELD_VALUES_SET_TYPE,1,10.0_CMISSDP,Err)

  !Start to set up CellML Fields

  !Create the CellML environment
  CALL cmfe_CellML_Initialise(CellML,Err)
  CALL cmfe_CellML_CreateStart(CellMLUserNumber,Region,CellML,Err)
  !Import a toy constant source model from a file
  CALL cmfe_CellML_ModelImport(CellML,"HinchBuf.cellml",constantModelIndex,Err)

  ! Now we have imported all the models we are able to specify which variables from the model we want:
  !   - to set from this side
  !These are effectively parameters that you know won't change in the course of the ode solving for one time step. i.e. fixed before running cellml, known in opencmiss and 
  !changed only in opencmiss - components of the parameters field
!  CALL cmfe_CellML_VariableSetAsKnown(CellML,constantModelIndex,"intracell/Na_i",Err)
  CALL cmfe_CellML_VariableSetAsKnown(CellML,constantModelIndex,"cellGeo/compartmentType",Err)

  !   - to get from the CellML side. variables in cellml model that are not state variables, but are dependent on independent and state variables. - components of intermediate field
!  CALL cmfe_CellML_VariableSetAsWanted(CellML,constantModelIndex,"fluer/CFluo",Err)

  !Finish the CellML environment
  CALL cmfe_CellML_CreateFinish(CellML,Err)


  !Start the creation of CellML <--> OpenCMISS field maps
  CALL cmfe_CellML_FieldMapsCreateStart(CellML,Err)
  !Now we can set up the field variable component <--> CellML model variable mappings.
  !here we map opencmiss fields to appropriate cellml field - either parameters/intermediates/state.
  !in monodomain, people typically want Vm, a state variable. In my case, Calcium is my state variable.
  !In my case, I will want to get an injection of calcium. Now, I have order splitting, which means that 
  !I don't actually use the source field (the problem is solved as an ode, and then as a pde separately.
  !I need to get calcium from my dependent field and solve the DAE to give new values of calcium in the dependent field again.
  ! this is then used as the initial set up for the dynamic solver. so map dependent field to appropriate component variable names.
  !cellml/opencmiss will look up the appropriate field.

  !If I didn't have order splitting, i.e. a proper reaction diffusion equation to solve, then I need to get the current ca conc. from the
  !dependent field, solve the dae, and then put the result of the dae into the source field. 
  CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
    & constantModelIndex,"intracell/Ca_i",CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,constantModelIndex,"intracell/Ca_i",CMFE_FIELD_VALUES_SET_TYPE, &
    & DependentField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,Err)

! Add buffer maps
CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,FluoField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE, &
& constantModelIndex,"intracell/Fluo",cmfe_FIELD_VALUES_SET_TYPE,Err)
CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,constantModelIndex,"intracell/Fluo",cmfe_FIELD_VALUES_SET_TYPE, &
& FluoField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE,Err)

CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,CaFluoField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE, &
& constantModelIndex,"intracell/CaFluo",cmfe_FIELD_VALUES_SET_TYPE,Err)
CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,constantModelIndex,"intracell/CaFluo",cmfe_FIELD_VALUES_SET_TYPE, &
& CaFluoField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE,Err)

CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,mBufField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE, &
& constantModelIndex,"intracell/B_nucm",cmfe_FIELD_VALUES_SET_TYPE,Err)
CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,constantModelIndex,"intracell/B_nucm",cmfe_FIELD_VALUES_SET_TYPE, &
& mBufField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE,Err)

CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,CamBufField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE, &
& constantModelIndex,"intracell/CaB_nucm",cmfe_FIELD_VALUES_SET_TYPE,Err)
CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,constantModelIndex,"intracell/CaB_nucm",cmfe_FIELD_VALUES_SET_TYPE, &
& CamBufField,cmfe_FIELD_U_VARIABLE_TYPE,1,cmfe_FIELD_VALUES_SET_TYPE,Err)

!  CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,NaiSourceField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
!    & constantModelIndex,"intracell/Na_i",CMFE_FIELD_VALUES_SET_TYPE,Err)
!  CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,constantModelIndex,"intracell/Na_i",CMFE_FIELD_VALUES_SET_TYPE, &
!    & NaiSourceField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,Err)

!  !!! Added stuff
  CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,CompartmentTypeField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
    & constantModelIndex,"cellGeo/compartmentType",CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,constantModelIndex,"cellGeo/compartmentType",CMFE_FIELD_VALUES_SET_TYPE, &
    & CompartmentTypeField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,Err)

  !Finish the creation of CellML <--> OpenCMISS field maps
  CALL cmfe_CellML_FieldMapsCreateFinish(CellML,Err)

  !set initial value of the dependent field/state variable, ca.
!*******************************************************************************
! YOU CAN CHANGE THIS to change the initial value of the dependent field at all the nodes
  !parametersetupdatenode changes a single node's value
 !componentvaluesinitialise changes the value of all nodes of a field

  CALL cmfe_Field_ComponentValuesInitialise(DependentField,CMFE_FIELD_U_VARIABLE_TYPE, &
    & CMFE_FIELD_VALUES_SET_TYPE,1,0.0001_CMISSRP,Err)

! Set initial values of buffers
CALL cmfe_Field_ComponentValuesInitialise(FluoField,cmfe_FIELD_U_VARIABLE_TYPE, &
& cmfe_FIELD_VALUES_SET_TYPE,1,0.002_CMISSRP,Err)
CALL cmfe_Field_ComponentValuesInitialise(CaFluoField,cmfe_FIELD_U_VARIABLE_TYPE, &
& cmfe_FIELD_VALUES_SET_TYPE,1,0.0_CMISSRP,Err)
CALL cmfe_Field_ComponentValuesInitialise(mBufField,cmfe_FIELD_U_VARIABLE_TYPE, &
& cmfe_FIELD_VALUES_SET_TYPE,1,0.0_CMISSRP,Err)
CALL cmfe_Field_ComponentValuesInitialise(CamBufField,cmfe_FIELD_U_VARIABLE_TYPE, &
& cmfe_FIELD_VALUES_SET_TYPE,1,0.0_CMISSRP,Err)


!*******************************************************************************

!Start the creation of the CellML models field. This field is an integer field that stores which nodes have which cellml model
  CALL cmfe_Field_Initialise(CellMLModelsField,Err)
  CALL cmfe_CellML_ModelsFieldCreateStart(CellML, CellMLModelsFieldUserNumber, &
    & CellMLModelsField,Err)
  !Finish the creation of the CellML models field
  CALL cmfe_CellML_ModelsFieldCreateFinish(CellML,Err)
  !The CellMLModelsField is an integer field that stores which model is being used by which node.
  !By default all field parameters have default model value of 1, i.e. the first model. But, this command below is for example purposes
  CALL cmfe_Field_ComponentValuesInitialise(CellMLModelsField,CMFE_FIELD_U_VARIABLE_TYPE, &
    & CMFE_FIELD_VALUES_SET_TYPE,1,1_CMISSIntg,Err)

  !Set up the models field
  !DO N=1,(NUMBER_GLOBAL_X_ELEMENTS+1)*(NUMBER_GLOBAL_Y_ELEMENTS+1)*(NUMBER_GLOBAL_Z_ELEMENTS+1)
  !  IF(N < 5) THEN
  !    CELL_TYPE = 1
  !  ELSE
  !    CELL_TYPE = 2
  !  ENDIF
  !  CALL cmfe_FieldParameterSetUpdateNode(CellMLModelsField, CMFE_FIELD_U_VARIABLE_TYPE, CMFE_FIELD_VALUES_SET_TYPE,1,N,1,CELL_TYPE,Err)
  !END DO
  !CALL cmfe_FieldParameterSetUpdateStart(CellMLModelsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  !CALL cmfe_FieldParameterSetUpdateFinish(CellMLModelsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)


  !Start the creation of the CellML state field
  CALL cmfe_Field_Initialise(CellMLStateField,Err)
  CALL cmfe_CellML_StateFieldCreateStart(CellML, & 
    & CellMLStateFieldUserNumber,CellMLStateField,Err)
  !Finish the creation of the CellML state field
  CALL cmfe_CellML_StateFieldCreateFinish(CellML,Err)

!  !Start the creation of the CellML intermediate field
!  CALL cmfe_Field_Initialise(CellMLIntermediateField,Err)
!  CALL cmfe_CellML_IntermediateFieldCreateStart(CellML, &
!    & CellMLIntermediateFieldUserNumber,CellMLIntermediateField,Err)
!  !Finish the creation of the CellML intermediate field
!  CALL cmfe_CellML_IntermediateFieldCreateFinish(CellML,Err)

  !Start the creation of CellML parameters field
  CALL cmfe_Field_Initialise(CellMLParametersField,Err)
  CALL cmfe_CellML_ParametersFieldCreateStart(CellML, &
    & CellMLParametersFieldUserNumber,CellMLParametersField,Err)
  !Finish the creation of CellML parameters
  CALL cmfe_CellML_ParametersFieldCreateFinish(CellML,Err)

 
  !Create the equations set equations
  CALL cmfe_Equations_Initialise(Equations,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(EquationsSet,Equations,Err)
  !Set the equations matrices sparsity type
  CALL cmfe_Equations_SparsityTypeSet(Equations,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  CALL cmfe_Equations_OutputTypeSet(Equations,CMFE_EQUATIONS_NO_OUTPUT,Err)
  !CALL cmfe_EquationsOutputTypeSet(Equations,cmfe_EquationsTimingOutput,Err)
  !CALL cmfe_EquationsOutputTypeSet(Equations,cmfe_EquationsMatrixOutput,Err)
  !CALL cmfe_EquationsOutputTypeSet(Equations,cmfe_EquationsElementMatrixOutput,Err)
  !Finish the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateFinish(EquationsSet,Err)

!!!! Buffer Equations
!Create the equations set equations for Fluo
CALL cmfe_Equations_Initialise(FluoEquations,Err)
CALL cmfe_EquationsSet_EquationsCreateStart(FluoEquationsSet,FluoEquations,Err)
CALL cmfe_Equations_SparsityTypeSet(FluoEquations,cmfe_EQUATIONS_SPARSE_MATRICES,Err)
CALL cmfe_Equations_OutputTypeSet(FluoEquations,cmfe_EQUATIONS_NO_OUTPUT,Err)
CALL cmfe_EquationsSet_EquationsCreateFinish(FluoEquationsSet,Err)

!Create the equations set equations for CaFluo
CALL cmfe_Equations_Initialise(CaFluoEquations,Err)
CALL cmfe_EquationsSet_EquationsCreateStart(CaFluoEquationsSet,CaFluoEquations,Err)
CALL cmfe_Equations_SparsityTypeSet(CaFluoEquations,cmfe_EQUATIONS_SPARSE_MATRICES,Err)
CALL cmfe_Equations_OutputTypeSet(CaFluoEquations,cmfe_EQUATIONS_NO_OUTPUT,Err)
CALL cmfe_EquationsSet_EquationsCreateFinish(CaFluoEquationsSet,Err)

!Create the equations set equations for mBuf
CALL cmfe_Equations_Initialise(mBufEquations,Err)
CALL cmfe_EquationsSet_EquationsCreateStart(mBufEquationsSet,mBufEquations,Err)
CALL cmfe_Equations_SparsityTypeSet(mBufEquations,cmfe_EQUATIONS_SPARSE_MATRICES,Err)
CALL cmfe_Equations_OutputTypeSet(mBufEquations,cmfe_EQUATIONS_NO_OUTPUT,Err)
CALL cmfe_EquationsSet_EquationsCreateFinish(mBufEquationsSet,Err)

!Create the equations set equations for CamBuf
CALL cmfe_Equations_Initialise(CamBufEquations,Err)
CALL cmfe_EquationsSet_EquationsCreateStart(CamBufEquationsSet,CamBufEquations,Err)
CALL cmfe_Equations_SparsityTypeSet(CamBufEquations,cmfe_EQUATIONS_SPARSE_MATRICES,Err)
CALL cmfe_Equations_OutputTypeSet(CamBufEquations,cmfe_EQUATIONS_NO_OUTPUT,Err)
CALL cmfe_EquationsSet_EquationsCreateFinish(CamBufEquationsSet,Err)


  !Create the problem
  CALL cmfe_Problem_Initialise(Problem,Err)
  CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_CLASSICAL_FIELD_CLASS, &
    & CMFE_PROBLEM_REACTION_DIFFUSION_EQUATION_TYPE,CMFE_PROBLEM_CELLML_REAC_INTEG_REAC_DIFF_STRANG_SPLIT_SUBTYPE],Problem,Err)
  !Finish the creation of a problem.
  CALL cmfe_Problem_CreateFinish(Problem,Err)

  !Create the problem control
  CALL cmfe_Problem_ControlLoopCreateStart(Problem,Err)
  !Get the control loop
  CALL cmfe_ControlLoop_Initialise(ControlLoop,Err)
  CALL cmfe_Problem_ControlLoopGet(Problem,CMFE_CONTROL_LOOP_NODE,ControlLoop,Err)
  !Set the times
!*******************************************************************************
! YOU CAN CHANGE THIS - you can change the initial time, final time and time step for your RD sim.
  CALL cmfe_ControlLoop_TimesSet(ControlLoop,0.0_CMISSRP,100.0_CMISSRP,2.0_CMISSRP,Err)
!*******************************************************************************
  CALL cmfe_ControlLoop_TimeOutputSet(ControlLoop,1,Err)
  CALL cmfe_ControlLoop_OutputTypeSet(ControlLoop,CMFE_CONTROL_LOOP_PROGRESS_OUTPUT,Err)
  !Finish creating the problem control loop
  CALL cmfe_Problem_ControlLoopCreateFinish(Problem,Err)


  !Set up the problem solvers for Strang splitting
  CALL cmfe_Problem_SolversCreateStart(Problem,Err)
  !First solver is a DAE solver
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL cmfe_Solver_DAESolverTypeSet(Solver,CMFE_SOLVER_DAE_EULER,Err)
  CALL cmfe_Solver_DAETimeStepSet(Solver,0.001_CMISSRP,Err)
  CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_MATRIX_OUTPUT,Err)

  !Second solver is the dynamic solver for solving the parabolic equation
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Solver_Initialise(LinearSolver,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,2,Solver,Err)
  !set theta - backward vs forward time step parameter
  CALL cmfe_Solver_DynamicThetaSet(Solver,1.0_CMISSRP,Err)
  !CALL cmfe_SolverOutputTypeSet(Solver,CMFE_SOLVER_NO_OUTPUT,Err)
  CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_TIMING_OUTPUT,Err)
  !get the dynamic linear solver from the solver
  CALL cmfe_Solver_DynamicLinearSolverGet(Solver,LinearSolver,Err)
  CALL cmfe_Solver_LibraryTypeSet(LinearSolver,CMFE_SOLVER_LAPACK_LIBRARY,Err)
  CALL cmfe_Solver_LinearTypeSet(LinearSolver,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
  CALL cmfe_Solver_LinearDirectTypeSet(LinearSolver,CMFE_SOLVER_DIRECT_LU,Err)
  !CALL cmfe_SolverLibraryTypeSet(LinearSolver,cmfe_SolverCMISSLibrary,Err)
  !CALL cmfe_SolverLinearTypeSet(LinearSolver,cmfe_SolverLinearDirectSolveType,Err)
  !CALL cmfe_SolverLibraryTypeSet(LinearSolver,cmfe_SolverMUMPSLibrary,Err)
  !CALL cmfe_Solver_LinearIterativeMaximumIterationsSet(LinearSolver,10000,Err)


  !Third solver is another DAE solver
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,3,Solver,Err)
  CALL cmfe_Solver_DAESolverTypeSet(Solver,CMFE_SOLVER_DAE_EULER,Err)
  CALL cmfe_Solver_DAETimeStepSet(Solver,0.001_CMISSRP,Err)
  CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_MATRIX_OUTPUT,Err)
  !CALL cmfe_SolverOutputTypeSet(Solver,cmfe_SolverTimingOutput,Err)
  !CALL cmfe_SolverOutputTypeSet(Solver,cmfe_SolverSolverOutput,Err)
  !CALL cmfe_SolverOutputTypeSet(Solver,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
  !Finish the creation of the problem solver
  CALL cmfe_Problem_SolversCreateFinish(Problem,Err)


  !Start the creation of the problem solver CellML equations
  CALL cmfe_Problem_CellMLEquationsCreateStart(Problem,Err)
  !Get the first solver  
  !Get the CellML equations
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL cmfe_CellMLEquations_Initialise(CellMLEquations,Err)
  CALL cmfe_Solver_CellMLEquationsGet(Solver,CellMLEquations,Err)
  !Add in the CellML environement
  CALL cmfe_CellMLEquations_CellMLAdd(CellMLEquations,CellML,CellMLIndex,Err)

  !Get the third solver  
  !Get the CellML equations
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,3,Solver,Err)
  CALL cmfe_CellMLEquations_Initialise(CellMLEquations,Err)
  CALL cmfe_Solver_CellMLEquationsGet(Solver,CellMLEquations,Err)
  !Add in the CellML environement
  CALL cmfe_CellMLEquations_CellMLAdd(CellMLEquations,CellML,CellMLIndex,Err)

  !Finish the creation of the problem solver CellML equations
  CALL cmfe_Problem_CellMLEquationsCreateFinish(Problem,Err)

  !Start the creation of the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateStart(Problem,Err)
  !Get the second solver  
  !Get the solver equations
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,2,Solver,Err)
  CALL cmfe_SolverEquations_Initialise(SolverEquations,Err)
  CALL cmfe_Solver_SolverEquationsGet(Solver,SolverEquations,Err)
  !Set the solver equations sparsity
  !CALL cmfe_SolverEquationsSparsityTypeSet(SolverEquations,cmfe_SolverEquationsSparseMatrices,Err)
  CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquations,CMFE_SOLVER_SPARSE_MATRICES,Err)  
  !Add in the equations set
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,EquationsSet,EquationsSetIndex,Err)
!!! Not found in earlier versions
!! Add equations for buffers
CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,FluoEquationsSet,EquationsSetIndex,Err)
CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,CaFluoEquationsSet,EquationsSetIndex,Err)
CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,mBufEquationsSet,EquationsSetIndex,Err)
CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,CamBufEquationsSet,EquationsSetIndex,Err)
  !Finish the creation of the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateFinish(Problem,Err)

!*******************************************************************************
! YOU CAN CHANGE THIS to change the dirchlet boundary conditions on the two ends of the line

  !Create the equations set boundary conditions
  WRITE(*,*) 'Set up boundary conditions'
  BCNODES = [1,NUMBER_GLOBAL_X_ELEMENTS+1]
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)
  DO node=1,2
    NODE_NUMBER = BCNODES(node)
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NODE_NUMBER,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      CONDITION = CMFE_BOUNDARY_CONDITION_FIXED
      VALUE=0.0_CMISSRP
! To fix variable change deludeln to different type (see below)
      CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,DependentField, &
       & CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,CMFE_NO_GLOBAL_DERIV, &
       & NODE_NUMBER,1,CONDITION,VALUE,Err)

! Set buffer boundary conditions - no flux at boundary
        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,FluoField, &
            & CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,CMFE_NO_GLOBAL_DERIV, &
            & NODE_NUMBER,1,CONDITION,VALUE,Err)
        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,CaFluoField, &
            & CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,CMFE_NO_GLOBAL_DERIV, &
            & NODE_NUMBER,1,CONDITION,VALUE,Err)

! Not nec
!        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,mBufField, &
!            & CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,CMFE_NO_GLOBAL_DERIV, &
!            & NODE_NUMBER,1,CONDITION,VALUE,Err)
!        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,CamBufField, &
!            & CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,CMFE_NO_GLOBAL_DERIV, &
!            & NODE_NUMBER,1,CONDITION,VALUE,Err)
!        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,ZThreeField, &
!            & CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_NO_GLOBAL_DERIV, &
!            & NODE_NUMBER,1,CONDITION,VALUEa,Err)


      !Need to set cellml model to zero at the nodes at which value of ca has been fixed.
      CALL cmfe_Field_ParameterSetUpdateNode(CellMLModelsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,&
       & 1,1,NODE_NUMBER,1,0_CMISSIntg,Err) 
    ENDIF
  ENDDO
! Prevent nuclear buffers from diffusing out of the nucleus nucpos
BCNODES = [14,NUMBER_GLOBAL_X_ELEMENTS+1]
DO node=1,2
    NODE_NUMBER = BCNODES(node)
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NODE_NUMBER,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
        CONDITION = CMFE_BOUNDARY_CONDITION_FIXED
        VALUE=0.0_CMISSRP

        ! Set buffer boundary conditions - no flux at boundary
        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,mBufField, &
        & CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,CMFE_NO_GLOBAL_DERIV, &
        & NODE_NUMBER,1,CONDITION,VALUE,Err)
        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,CamBufField, &
        & CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,CMFE_NO_GLOBAL_DERIV, &
        & NODE_NUMBER,1,CONDITION,VALUE,Err)
        ENDIF
    ENDDO

!!! Added Stuff
! Why did I add this?
!DO node=1,NUMBER_GLOBAL_X_ELEMENTS+1
!    NODE_NUMBER = BCNODES(node)
!    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NODE_NUMBER,1,NodeDomain,Err)
!    IF(NodeDomain==ComputationalNodeNumber) THEN
!        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,TRPNField, &
!            & CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_NO_GLOBAL_DERIV, &
!            & NODE_NUMBER,1,CONDITION,0.06_CMISSRP,Err)
!        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,CaSRField, &
!            & CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_NO_GLOBAL_DERIV, &
!            & NODE_NUMBER,1,CONDITION,0.7_CMISSRP,Err)
!        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,ZOneField, &
!            & CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_NO_GLOBAL_DERIV, &
!            & NODE_NUMBER,1,CONDITION,1.0_CMISSRP,Err)
!        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,ZTwoField, &
!            & CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_NO_GLOBAL_DERIV, &
!            & NODE_NUMBER,1,CONDITION,0.009_CMISSRP,Err)
!        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,ZThreeField, &
!            & CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_NO_GLOBAL_DERIV, &
!            & NODE_NUMBER,1,CONDITION,0.003_CMISSRP,Err)
!    ENDIF
!    ENDDO
!!! Added Stuff nucpos
!Change model in nucleus
!YOU CAN CHANGE THIS TO CHANGE THE POSITION OF THE NUCLEUS
  DO node=14,NUMBER_GLOBAL_X_ELEMENTS+1
    NODE_NUMBER = BCNODES(node)
    CALL cmfe_Field_ParameterSetUpdateNode(CompartmentTypeField,cmfe_FIELD_U_VARIABLE_TYPE, &
        & cmfe_FIELD_VALUES_SET_TYPE,1,1,node,1,2.0_CMISSDP,Err)
    CALL cmfe_Field_ParameterSetUpdateNode(mBufField,cmfe_FIELD_U_VARIABLE_TYPE, &
        & cmfe_FIELD_VALUES_SET_TYPE,1,1,node,1,0.065_CMISSDP,Err)
  ENDDO

!Add copy of Hinch model every 2um
! Make sure closest copy is
DO node=2,13,4
    NODE_NUMBER = BCNODES(node)
        WRITE(*,*) "The nodes with a Hinch model attached include",node
        CALL cmfe_Field_ParameterSetUpdateNode(CompartmentTypeField,cmfe_FIELD_U_VARIABLE_TYPE, &
            & cmfe_FIELD_VALUES_SET_TYPE,1,1,node,1,0.0_CMISSDP,Err)
ENDDO


!*******************************************************************************

  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)

  !Solve the problem
  CALL cmfe_Problem_Solve(Problem,Err)

  EXPORT_FIELD=.TRUE.
  IF(EXPORT_FIELD) THEN
    CALL cmfe_Fields_Initialise(Fields,Err)
    CALL cmfe_Fields_Create(Region,Fields,Err)
    CALL cmfe_Fields_NodesExport(Fields,"CellMLSplitReactionDiffusion1D","FORTRAN",Err)
    CALL cmfe_Fields_ElementsExport(Fields,"CellMLSplitReactionDiffusion1D","FORTRAN",Err)
    CALL cmfe_Fields_Finalise(Fields,Err)

  ENDIF
  

  WRITE(*,'(A)') "Program successfully completed."
  

  STOP
  
END PROGRAM CELLMLSPLITREACTIONDIFFUSION1DEXAMPLE
