%% PART I: LHS generalization
% - Allow for both continuous (P1)
% and discontinuous (P1 and P0) Galerkin
% - Done:
%       - P0Func (discontinuous)
%       - P1Func (continuous)
% - To do:
%       - P1Func (discontinuous) (should be similar to P0, M)
%       - FgaussDiscFunc (discontinuous*)


%% PART II: RHS generalization
% Almost there
% - Create evaluate functions

%% Cases
% Six first cases as tests:
%  - P0 to P1 continuous - DONE!
%  - P0 to P1 discontinuous - DONE!
%  - P1 contiuous to P0 - DONE!
%  - P1 discontiuous to P0
%  - FgaussDiscontinuous to P0
%  - FgaussDiscontinuous to P1 continuous
%  - FgaussDiscontinuous to P1 discontinuous (*)

% Plot should only be in P1 Discontinuous
% Temporarily:
%  - Projector to p1 continuous
%  - projector to p0
%  - projector to p1 discont
%       - case p0: fRepeated
%       - case p1: fRepeated
%       - case fGaussDiscont: ?

%%

% Delete computeValueInCenterElement in P1Function
% Think about computeDiscontinousField in P1Function (should return a P1
%   discontinuous function. More dofs, nodes between elements are not coupled)
%   - computeFnodesByelem is really a p1 discont function. perhaps rename
%     to transform or sth like that. move to discont.

% remove computeDiscontinousField 

% Create tests for these functions
%   - White circle in a black square (create using levelset, compare volumes)
%   - Quadrilateral elements
%   - Volumes

% P0Function: rename fElem to fValues
%   - createDiscontinuousP0 is kind of a Discontinuous P1 function
% P1Function: computeFnodesByElem now public, make it return fElem

% Create P1Function_Discontinuous // 

% Strain should be a FGaussDiscontFunction or sth like that. It's not a P0
%   function, only true if ngaus = 1. These gauss functions need the
%   quadrature: when the projector projects, it will evaluate it in these
%   gauss points
% FGaussDiscontFunction.plotDiscontinuous()

% Paraview: P1ContFunction or P1DiscontFunction?
% Matlab wants P1DiscontFunction

% Check: P1DiscontFunction's mass matrix is diagonal

% ProjectorP1toP0: mass matrix should use lhsintegrator

% Let Paraview print strain as a P1DiscontFunct

% ProjectorFromFgausToP1Discont  (similar to other projectors, mass matrix
% using LHSintegrator and discontinuousMesh. M should be invertible and
% have a block-like structure)


%% Comments / Questions
% - 