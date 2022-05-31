%% Updates
% - Field
%   BoundaryConditions no longer there, but they can be used to
%   translate.

% - LHSintegrators
%   They now only take Meshes (for dV) and Fields as inputs*

% - DimensionVariables
%   Gone, alongside with DimensionScalar and DimensionVector. Now included
%   in Fields as FieldDimensions.

% - StokesProblem
%   Fully cleaned up. Elements, DOFs and BCAppliers are no more*.

% - StokesProblem
%   For some reason, Transient is consistently way faster than Steady.

% - Writing
%   Picking up speed, path ahead clear but should focus on that


%% Backlog
% - Move Geometry to Mesh
% - Merge Field to FeFunction
%       - Mesh.coord as a FeFunction as well
% - Elegant constructors with RequiredProperties and OptionalProperties