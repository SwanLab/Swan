%% PROJECTORS | To-do

% Done! (T)  Rename FeFunction to P1Function
% Done! (T)  Rename PiecewiseConstantFunction to P0Function
% Done! (T)  Make them extend FeFunction (new)

% Done! (J)  Create plot method for P0Function
% Done! (T)  Create plot method for P1Function
% Done! (J)  Check: create P1function, plot and project to P0

% Done! (T)  Create ProjectorP0toP1 following the other example
% Done! (T)  Check: create P0function, plot and project to P1 (new)

%% Comments
% Projector P1 to P0 should return a P0Function...
%   ... but a P0Function requires a Mesh
%   ... requiring a Mesh implies defining the domain of the function
%   ... but then again, it kind of has already been implied by giving it
%   certain values at certain locations

% SOLUTION:
% Heavy lifiting required on P0Function. Almost everything can be moved to
%   the projector, structure should now be similar to P1Function
% Thus, meshes *only* needed in the projectors. Functions will only need
%   connec + element type