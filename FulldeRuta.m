%% To-do
% a) CLEANUP
%           - Stokes_Problem still uses BoundaryConditionsApplier
%                   - Left for legacy purposes
%       OK! - Rename New stuff to just the name
%                   - ElasticProblem, DiffReactProblem
%                   - ElasticProblemMicro, DiffReactProblemMicro
%       OK!    - Delete old DOFs and Elements
%                   - Stokes, Hyperelastic?
%       OK! - Move dofsInElem to LHSintegrator

% b) CLOSING
%           - Adapt ThermalProblem via NewDiffReactProblem
%       OK! - Add vars2print for micro

% c) RESULT VISUALIZATIONS
%           - Set up GiD for ThermalProblem and visualize results
%           - Set up GiD for Micro problem and visualize results

% d) DOFSINELEM
%           - 

% e) PERFORMANCE
%           - See below


% - Stokes encara fa servir bcapplier
% - Tancar tèrmic i micro
% - Visualitzar resultats amb gid
% - Mirar Strain computer

%% Report
% We are NOT using the sparse method properly. It is way better to generate
% arrays of indices and values, and then call the sparse function to
% assemble the matrix.


%% Old to-do
% a) BOUNDARY CONDITIONS
%     OK!  - Merge BoundaryConditions and BoundaryConditionsApplier
%     OK!  - Simplify and clean NewBoundaryConditions

% b) PROBLEM CLEANUP
%     OK!  - Clean NewDiffReactProblem
%     OK!  - Clean NewElasticProblemMicro
%     HMM  - Clean Thermal_Problem
%               - Created test_thermal
%               - "Material not yet implemented"
%               - The physics don't really pop up in NewDiffReactProblem
%     HMM  - Clean Hyperelastic_Problem
%               - Created test_hyperelastic
%               - Problem was not working properly
%               - Uses discontinued methods? Don't really know
%                 Isotropic2dHyperElasticMaterial

% c) OTHER CLEANUP

% c.1) ASSEMBLER
%       - Simplify Assembler
%     OK!  - Use dofsInElem at assembleMatrix()
%          - Use dofsInElem at assembleBMatrix()
%               - See c.2 
%          - Explore alternate ways to assemble B and C
%           (LHSintegrator_StiffnessElasticStoredB)
%          - Delete globalConnec from LHSintegrator/Assembler

% c.2) CANTILEVER BEAM
%     OK!  - CantileverBeam to CantileverBeamMeshCreator
%     OK!  - Fix CantileverBeamMeshCreator for 2D meshes
%     OK!  - Fix CantileverBeamMeshCreator for 3D meshes
%       - StressComputer increase performance
%           OK!  - Through "vectorize"
%           OK!  - Through bsxfun
%             - Through assembly+product
%       - StrainComputer increase performance
%             - Through "vectorize"
%             - Through bsxfun
%             - Through assembly+product
