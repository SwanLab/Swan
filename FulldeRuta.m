%% To-do


% a) PROJECT CHARTER

% b) PROBLEM CLEANUP
%       - Delete ElementDiffReact
%       - Delete ElasticProblem?
%       - Delete DiffReact_Problem?
%       - Clean NewDiffReactProblem

% c) BOUNDARY CONDITIONS
%       - Simplify and clean BoundaryConditions
%       - Simplify and clean BoundaryConditionsApplier

% c) OTHER CLEANUP
%       - Simplify Assembler
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