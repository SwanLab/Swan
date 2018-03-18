classdef ScalarProduct < handle
    %ScalarProduct Summary of this class goes here
    %   Detailed explanation goes here
    
    properties %(Access = private)
        epsilon
        Ksmooth
        Msmooth
    end
    
    methods 
        function obj = ScalarProduct(problemID,epsilon)
            obj.epsilon = epsilon;
            % Implement DiffReact_Problem instead of Physical with diff-react element
            physProb = Physical_Problem(problemID,'DIFF-REACT');
            physProb.mesh.scale='MACRO';
            physProb.preProcess;
            obj.Ksmooth = physProb.element.computeStiffnessMatrix;
            obj.Msmooth = physProb.element.computeMassMatrix(2);
        end
    end
    
    methods
        function sp = computeSP(obj,f,g)
            sp = f'*(((obj.epsilon)^2)*obj.Ksmooth + obj.Msmooth)*g;
        end
        
        %% !! USE APPROPIATE TERMINOLOGY -- SP WITH ONLY M IS CALLED...? !!
        function sp = computeSP_M(obj,f,g)
            sp = f'*(obj.Msmooth)*g;
        end
        
        function sp = computeSP_K(obj,f,g)
            sp = f'*(obj.Ksmooth)*g;
        end
    end
end

