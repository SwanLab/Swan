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