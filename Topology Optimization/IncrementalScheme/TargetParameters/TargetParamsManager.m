classdef TargetParamsManager < handle
    
    properties (GetAccess = public, SetAccess = private)
        targetParams
    end
    
    properties (Access = private)
        nSteps
        
        volumeFrac
        constraintTol
        optimalityTol
        epsilon
        epsilonPer
        epsilonIso
    end
    
    methods (Access = public)
        
        function obj = TargetParamsManager(cParams)
            obj.init(cParams);
            obj.generateSequences(cParams);
            obj.update(1);
        end
        
        function update(obj,iStep)
            obj.computeValues(iStep);
            obj.assignValues();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.nSteps = cParams.nSteps;
            obj.targetParams = TargetParameters();
        end
        
        function generateSequences(obj,cParams)
            obj.volumeFrac = LinearSequence(1/obj.nSteps,1,obj.nSteps,cParams.VfracInitial,cParams.VfracFinal);
            %obj.volumeFrac = LogarithmicSequence(0.85,obj.nSteps,cParams.VfracInitial,cParams.VfracFinal);
            obj.constraintTol = LinearSequence(0,1,obj.nSteps,cParams.constrInitial,cParams.constrFinal);
            obj.optimalityTol = LinearSequence(0,1,obj.nSteps,cParams.optimalityInitial,cParams.optimalityFinal);
            obj.epsilon = LinearSequence(0,1,obj.nSteps,cParams.epsilonInitial,cParams.epsilonFinal);
            obj.epsilonPer = LogarithmicSequence(-0.9,obj.nSteps,cParams.epsilonPerInitial,cParams.epsilonPerFinal);
            obj.epsilonIso = LinearSequence(0,1,obj.nSteps,cParams.epsilonIsoInitial,cParams.epsilonIsoFinal);
        end
        
        function computeValues(obj,iStep)
            obj.volumeFrac.update(iStep);
            obj.constraintTol.update(iStep);
            obj.optimalityTol.update(iStep);
            obj.epsilon.update(iStep);
            obj.epsilonPer.update(iStep);
%             epsi = zeros(obj.nSteps,1);
%             for i = 1:obj.nSteps
%             obj.epsilonPer.update(i);
%             epsi(i,1) = obj.epsilonPer.value;
%             %obj.volumeFrac.update(i);
%             %epsi(i,1) = obj.volumeFrac.value;            
%             end
%             plot(epsi,'+')
        end
        
        function assignValues(obj)
            obj.targetParams.constr_tol = obj.constraintTol.value;
            obj.targetParams.optimality_tol = obj.optimalityTol.value;
            obj.targetParams.Vfrac = obj.volumeFrac.value;
            obj.targetParams.epsilon = obj.epsilon.value;
            obj.targetParams.epsilon_velocity = obj.epsilon.value;
            obj.targetParams.epsilon_perimeter = obj.epsilonPer.value;
            obj.targetParams.epsilon_isotropy = obj.epsilonIso.value;
        end
        
    end
    
end