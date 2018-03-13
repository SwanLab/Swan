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
            physProb = Physical_Problem(problemID,'DIFF-REACT');
%             physProb.mesh.ptype = 'DIFF-REACT';
            physProb.preProcess;
            %% !! Change how Ksmooth & Msmooth are computed !!
            [obj.Ksmooth, obj.Msmooth] = physProb.computeKM(2);
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

