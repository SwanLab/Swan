%% PROJECTORS | To-do

% Done! (T)  Rename FeFunction to P1Function
% Done! (T)  Rename PiecewiseConstantFunction to P0Function
% Done! (T)  Make them extend FeFunction (new)

% Done! (J)  Create plot method for P0Function
% Done! (T)  Create plot method for P1Function
% Done! (J)  Check: create P1function, plot and project to P0

% Done! (T)  Create ProjectorP0toP1 following the other example
% Done! (T)  Check: create P0function, plot and project to P1 (new)

%% Comments / Questions
% - We can now define the displacement in the ElasticProblem as a
%   P1Function.
%       - Doing so would break the rest of the code / require heavy
%         rewriting of the topology optimization part
%       - fElem value in P1Function is defined as nnodeEl x ndim x nElem
%
% - We *could* define the strain as a P0Function, but there are several
%   issues with that:
%       - Strain is a nstre x nelem x ngaus matrix. In the examples that we
%         have used, only ngaus = 1 is considered.
%       - Thus, the fElem value in P0Function is always defined as a
%         ndim x 1 x nelem matrix.
%       - Should it be this way? Could we maintain the current fElem
%         structure and simply interpolate at other Gauss points?
%       - A part from that, there is another issue regarding dimensions.
%         Even if we sail ahead, in 2D cases the matrix would have ndim = 3
%         (XX, YY, XY). However, in 3D cases, the order would be
%         (XX, YY, ZZ, XY, YZ, XZ): the third index can mean different 
%         stuff -> needed for Paraview.
%       - We may need to specify R^n. Should we always consider the 3D case
%         for mesh-related stuff?
%       - Perhaps this could be a first step to move shape functions there,
%         solving the uneasiness we had with BMatrixComputer.
%
% - Also, plotting in 3D is scuffed.