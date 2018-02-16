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
        end
        
        function preProcess(obj)
            %props.mu=0.375;
            %props.kappa = 0.75;
            %obj.material = obj.material.setProps(props);            
            obj.bc = BC_Micro(obj.dim.nunkn,obj.problemID,obj.mesh.coord,obj.mesh.ptype,obj.dim.ndim);
            
            obj.dof = DOF.create(obj.geometry.nnode,obj.mesh.connec,obj.dim.nunkn,obj.mesh.npnod,obj.bc,obj.mesh.scale);
            obj.element = Element.create(obj.mesh,obj.geometry,obj.material,obj.bc,obj.dof,obj.dim);%Element_Elastic_Micro;
           % obj.variables = PhysicalVars_Elastic_2D_Micro(obj.dof.ndof);
            obj.solver = Solver.create(obj.mesh.ptype);
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

