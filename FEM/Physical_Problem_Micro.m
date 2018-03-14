classdef Physical_Problem_Micro < Physical_Problem
    %Physical_Problem Summary of this class goes here
    %   Detailed explanation goes here
    
    %% Public GetAccess properties definition =============================
    properties (GetAccess = public, SetAccess = private)        
    end
    
    %% Private properties definition ======================================
    properties (Access = private)
    end
    
    %% Public methods definition ==========================================
    methods (Access = public)
        function obj = Physical_Problem_Micro(problemID)
            obj@Physical_Problem(problemID);
            
            % Just to match Ferran's code%%%%%%%%%%%%%%%%%%%
            props.mu=0.375;
            props.kappa = 0.75;
            obj.material = obj.material.setProps(props);
%             obj.element.material = obj.element.material.setProps(props);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end       
        
        function [Chomog,tstrain,tstress] = computeChomog(obj)      
            
           % obj.variables = PhysicalVars_Elastic_2D_Micro(obj.dof.ndof);
            vstrain=diag(ones(obj.dim.nstre,1));
            Chomog =  zeros(obj.dim.nstre,obj.dim.nstre);
            tstrain = zeros(obj.dim.nstre,obj.geometry.ngaus,obj.dim.nstre,obj.mesh.nelem);
            tstress = zeros(obj.dim.nstre,obj.geometry.ngaus,obj.dim.nstre,obj.mesh.nelem);
            for istre=1:obj.dim.nstre
                obj.element.vstrain = vstrain(istre,:);
                obj.computeVariables();
                Chomog(:,istre) = obj.variables.stress_homog;
                tstrain(istre,:,:,:) = obj.variables.strain;
                tstress(istre,:,:,:) = obj.variables.stress;
            end
            obj.variables.Chomog = Chomog;
            obj.variables.tstrain = tstrain;
            obj.variables.tstress = tstress;
        end
    end
end

