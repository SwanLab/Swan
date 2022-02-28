%% Todo
% {{done}} Create example 2D not using load
% {{done}}-ish Create example 3D not using load
% eliminate istre,jstre loop 
% investigate how to efficiently multiply B,C,B
% {{done}} Use BmatrixComputer in LHSintegrator_StifnessElastic
% With large example compare Sparse vs Accumarray
% Element_DiffReact K, M, Mr with LHSintegrator
%         - K done (had to fix LHSintegrator_Stiffness)
%         - M done
%         - Mr done
%         NOTE: this change affects PlottingTests. Solved, though!
% DONE - Eliminate computeLHS from Integrator_Simple
% * after solving the issues with PlottingTests
% DONE - Eliminate computeLHS from integratorComposite
% * after solving the issues with PlottingTests

% NEW: fix LHSintegrator_Stiffness

%% Endgame
% 7. Force use integrator for assembly
% 8. LHS integrator must be composed by integrator_simple for assemnbly