classdef TargetParameters < handle
    
    properties (Access = public)
        constr_tol
        optimality_tol
        Vfrac
        epsilon
        epsilon_velocity
        epsilon_perimeter
        epsilon_isotropy
        stressNormExponent
        iStep
    end
    
end