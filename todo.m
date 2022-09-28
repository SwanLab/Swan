%% Projector Generalization

% Project analytical functions to P0 P1 P1d (eg. sin)
% Project L2 functions to P0 P1 P1d (eg. circle w/ heaviside)
% (CharacteristicFunction)

% Test with quadrilaterals as well

% NOTE: assumption made for FGauss to P1Discont

%%

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