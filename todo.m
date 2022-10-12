%% To-do

% UML Projectors
% plot fgaussfun
% geometry to mesh (jacobian in mesh)
% BMatComp(dNdx), dndx coming from interp
% create gradient of function and compSymmVector
% projectors with unfittedmesh at rhs without it noticing it

% Pending:
% - Projectors/Filters with Triangle&Quads (figures)

% Jose: 
% Done! - Delete element loops of projectors
% Done! - H1Projection
% - Clean code of projectors + tests
% Done! - Use filters in these tests and compare Filters vs L2Projectors vs H1Projectors]

% Ton:
% Done* - ACCEPT PULL REQUESTS FOR NULLSPACE
% Done* - projectors without discont mesh (use field as dispatching in lhs
%         integrator assembly)
% Done! - quadrilateral elements working
% Done! - assert fgauss

% u -> diff(u) (strain) -> plot -> project p1 p0 p1d
% Create a function in p1function named function fGaussFun  = comptueGrad(quad)
% project this function to p1 p0 and p1d. check results. ideally: a
% quadratic function u = x^2 -> du/dy = 0, du/dx = 2x; only evaluated at
% gauss points. project 

%% Changelog
% + plot() now plots all dimensions

%% Notes (Ton)
% FGaussFun SHOULD assert. The discrepancy in quadrature results from the
% quadrature used in the projector -- perhaps dispatch thru incoming
% function?

% - Quadratic needed for projecting to P1Discontinuous (???) a
%   Characteristic Function

%% Next steps

% - QUESTION: what to do eg. CutMeshProvisionalQuadrilateral. Projector or
%             evaluate?



%% Backlog

% - Let Paraview print strain as a P1DiscontFunct