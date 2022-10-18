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
% Train of thought regarding Geometry/Mesh/Interpolation

%   - Basically, everything is linked up together. Geometry is used as a
%   flexible tool to perform various actions
%   - It is also very much dynamic, the same geometry can be re-computed
%   many times in a single problem with different parameters

%   - When a Mesh is created, a Geometry_Volumetric is also created.
%   However, nothing else is calculated.
%   - The calculations are performed when calling computeGeometry(q,int).
%   Thus, there is no one Jacobian/dNdx: 
%       - Quad: constant (normals) /linear (linear int) / quadratic (quad
%       fields)
%       - Int: linear (mesh) / quadratic (fields)


%       - mesh.computeInvJac(q,int)

%% Next steps

% - QUESTION: what to do eg. CutMeshProvisionalQuadrilateral. Projector or
%             evaluate?



%% Backlog

% - Let Paraview print strain as a P1DiscontFunct