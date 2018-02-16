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
        
        function computeVariables(obj,vstrain)
%             obj.element.computeLHS(obj.dim.nunkn,obj.dim.nstre,obj.mesh.nelem,obj.geometry,obj.material);
%             obj.element.computeRHS(obj.dim.nunkn,obj.dim.nstre,obj.mesh.nelem,obj.geometry.nnode,obj.material,obj.bc,obj.dof.idx,obj.geometry,vstrain);
%             
%             % Assembly
%             [obj.LHS,obj.RHS] = obj.Assemble(obj.element,obj.geometry.nnode,obj.dim.nunkn,obj.dof,obj.bc);
%                        
%             % Solver            
%             data.pnodes = obj.bc.pnodes;
%             data.nunkn = obj.dim.nunkn;
%             obj.solver.setSolverVariables(data);
%             sol = obj.solver.solve(obj.variables.d_u,obj.LHS,obj.RHS,obj.dof);
%             
%             obj.variables = obj.variables.computeVars(sol,obj.dim,obj.geometry,obj.mesh.nelem,obj.dof.idx,obj.element,obj.material,vstrain);
%           
            obj.element.vstrain = vstrain;
             
            tol   = 1e-6;
            x = zeros(length(obj.dof.vF),1);
            % Compute r & dr
            [r,dr] = obj.element.computeResidual(x);
            while dot(r,r) > tol
                inc_x = obj.solver.solve(dr,-r);
                x_new = x + inc_x;
                % Compute r & dr
                [r,dr] = obj.element.computeResidual(x_new);
                x = x_new;
            end
            
            obj.variables = obj.element.computeVars(x);
        end
        
        function [Chomog,tstrain,tstress] = computeChomog(obj)
           % obj.variables = PhysicalVars_Elastic_2D_Micro(obj.dof.ndof);
            vstrain=diag(ones(obj.dim.nstre,1));
            Chomog =  zeros(obj.dim.nstre,obj.dim.nstre);
            tstrain = zeros(obj.dim.nstre,obj.geometry.ngaus,obj.dim.nstre,obj.mesh.nelem);
            tstress = zeros(obj.dim.nstre,obj.geometry.ngaus,obj.dim.nstre,obj.mesh.nelem);
            for istre=1:obj.dim.nstre
                obj.computeVariables(vstrain(istre,:));
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

