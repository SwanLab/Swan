%% To-do

% Jose: 
% - H1Projection
% - Delete element loops and clean code of projectos
% - Use filters in these tests and compare Filters vs L2Projectors vs H1Projectors]

% Ton:
%       - ACCEPT PULL REQUESTS FOR NULLSPACE
%       - projectors without discont mesh (use field as dispatching in lhs
%         integrator assembly)
% Done! - quadrilateral elements working
% Done! - assert fgauss

% u -> diff(u) (strain) -> plot -> project p1 p0 p1d
% Create a function in p1function named function fGaussFun  = comptueGrad(quad)
% project this function to p1 p0 and p1d. check results. ideally: a
% quadratic function u = x^2 -> du/dy = 0, du/dx = 2x; only evaluated at
% gauss points. project 

%% Changelog
% + Added fType to functions to check required quadrature orders and stuff.
%   Probably could have been done with isa() or class(), but this seems to
%   be cleaner

%% Notes (Ton)
% FGaussFun SHOULD assert. The discrepancy in quadrature results from the
% quadrature used in the projector -- perhaps dispatch thru incoming
% function?

% - Quadratic needed for projecting to P1Discontinuous (???) a
%   Characteristic Function

%% Next steps

% - Volumes?


%% Backlog

% - Let Paraview print strain as a P1DiscontFunct