%% To-do

% a) RENAMING
%    OK!   - s.pdim 'FILTER'. 
%    OK!   - pdim to nDim, nunkn to dimField, nFields, dimAllFields

% b) ASSEMBLER
%    OK!   - Accumarray and sparse only in Assembler. 
%    MEH   - BMatrixComputer uses Assembler.
%    MEH   - LHSintegrator_StiffnessElasticStoredB uses Assembler.
%    OK!   - LHSintegrator uses Assembler.
%    BTW   - ForcesComputer uses Assembler. (-ish)
%    ???   - StrainComputer, StressComputer

% c) EXAMPLES
%    OK!   - NewFemExamples as a class
%    OK!   - Following cleancode techniques
%    BTW   - Created CantileverBeam, PerformanceTests

% d) DIFFREACT_PROBLEM
%    OK!   - DiffReact_Problem to NewDiffReactProblem
%    OK!   - Delete Element_DiffReact (%)
%    OK!   - Use it in FilterPDE
%    MEH   - Cleanup on NewDiffReact
%    BTW   - NewDiffReactProblemMicro also done

% e) TOPOPT
%    OK!   - FEM to NewFem in TopOpt 

% f) COMPARISON
%       - First examples of: 
%       - Comparing product: pagemtimes, istrjstreLoop, pagefun 
%       - Comparing assembly: accumarray and sparse (Assembler)
%       - Comparing commutative of (product + assembly) vs
%           (assembly + product)
