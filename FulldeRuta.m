%% To-do

% a) RENAMING
%       - s.pdim 'FILTER'. 
%    OK!   - pdim to nDim, nunkn to dimField, nFields, dimAllFields

% b) ASSEMBLER
%    OK!   - Accumarray and sparse only in Assembler. 
%       - BMatrixComputer uses Assembler.
%       - LHSintegrator_StiffnessElasticStoredB uses Assembler.
%    OK!   - LHSintegrator uses Assembler.

% c) EXAMPLES
%       - NewFemExamples as a class
%       - Following cleancode techniques

% d) DIFFREACT_PROBLEM
%       - DiffReact_Problem to NewDiffReactProblem
%       - Delete Element_DiffReact (%)
%       - Use it in FilterPDE

% e) TOPOPT
%       - FEM to NewFem in TopOpt 

% f) COMPARISON
%       - First examples of: 
%       - Comparing product: pagemtimes, istrjstreLoop, pagefun 
%       - Comparing assembly: accumarray and sparse (Assembler)
%       - Comparing commutative of (product + assembly) vs
%           (assembly + product)
