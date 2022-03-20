%% To-do

% a) PROBLEM CLEANUP
%     OK!  - Delete Element_DiffReact
%     OK!  - Delete ElasticProblem?
%     OK!          - Fix NewElasticProblemMicro
%     OK!          - Fix micro FemTests/FemTestsSuite
%     OK!          - Fix micro TopOptTests/TopOptTestsSuite
%     OK!          - Fix micro HomogenizationTests
%     OK!  - Delete DiffReact_Problem?
%       - Clean NewDiffReactProblem

% b) BOUNDARY CONDITIONS
%       - Merge BoundaryConditions and BoundaryConditionsApplier
%       - Simplify and clean NewBoundaryConditions

% c) OTHER CLEANUP

% c.1) ASSEMBLER
%       - Simplify Assembler
%       - Use dofsInElem at assembleMatrix()
%       - Use dofsInElem at assembleBMatrix()
%       - Explore alternate ways to assemble B and C
%           (LHSintegrator_StiffnessElasticStoredB)
%       - Delete globalConnec from LHSintegrator/Assembler

% c.2) CANTILEVER BEAM
%       - CantileverBeam to CantileverBeamMeshCreator
%       - Fix CantileverBeamMeshCreator for 2D meshes
%       - Fix CantileverBeamMeshCreator for 3D meshes
%       - StressComputer increase performance through ("vectorize"/bsxfun/assmelby+product)
%       - StrainComputer increase performance through ("vectorize"/bsxfun/assmelby+product)



% z) COMPARISON
%       - First examples of: 
%       - Comparing product: pagemtimes, istrjstreLoop, pagefun 
%       - Comparing assembly: accumarray and sparse (Assembler)
%       - Comparing commutative of (product + assembly) vs
%           (assembly + product)

%% Comments
%       - We are three problems away from deleting them entirely (as well
%         as DOFs and Elements):  
%           - Hyperelastic_Problem is not used at all
%           - Thermal_Problem is not used at all
%           - Stokes_Problem is used in one test