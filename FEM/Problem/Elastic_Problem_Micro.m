classdef Elastic_Problem_Micro < Elastic_Problem
    %Elastic_Problem Summary of this class goes here
    %   Detailed explanation goes here
    
    %% Public GetAccess properties definition =============================
    properties (GetAccess = public, SetAccess = private)        
    end
    
    %% Private properties definition ======================================
    properties (Access = private)
    end
    
    %% Public methods definition ==========================================
    methods (Access = public)
        function obj = Elastic_Problem_Micro(problemID)
            obj@Elastic_Problem(problemID);
            obj.dof = DOF_ElasticMicro(problemID,obj.geometry,obj.mesh);
            
            % Just to match Ferran's code%%%%%%%%%%%%%%%%%%%
            props.mu=0.375;
            props.kappa = 0.75;
            obj.material = obj.material.setProps(props);
%             obj.element.material = obj.element.material.setProps(props);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end       
        
        function [Chomog,tstrain,tstress] = computeChomog(obj)      
            obj.element.quadrature.computeQuadrature('LINEAR');
            obj.element.interpolation_u.computeShapeDeriv(obj.element.quadrature.posgp)
            obj.element.geometry.computeGeometry(obj.element.quadrature,obj.element.interpolation_u);
           % obj.variables = PhysicalVars_Elastic_2D_Micro(obj.dof.ndof);
            vstrain = diag(ones(obj.element.nstre,1));
            Chomog =  zeros(obj.element.nstre,obj.element.nstre);
            tstrain = zeros(obj.element.nstre,obj.element.quadrature.ngaus,obj.element.nstre,obj.element.nelem);
            tstress = zeros(obj.element.nstre,obj.element.quadrature.ngaus,obj.element.nstre,obj.element.nelem);
            for istre=1:obj.element.nstre
                obj.element.vstrain = vstrain(istre,:);
                obj.computeVariables;
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

