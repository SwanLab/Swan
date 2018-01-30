classdef Physical_Problem_Micro < Physical_Problem
    %Physical_Problem Summary of this class goes here
    %   Detailed explanation goes here
    
    %% Public GetAccess properties definition =============================
    properties (GetAccess = public, SetAccess = private)
        StressHomog
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
            props.mu=0.375;
            props.kappa = 0.75;
            obj.material = obj.material.setProps(props);
            
            obj.bc = BC_Micro(obj.dim.nunkn,obj.problemID,obj.mesh.coord,obj.mesh.ptype,obj.dim.ndim);
            obj.dof = DOF(obj.geometry.nnode,obj.mesh.connec,obj.dim.nunkn,obj.mesh.npnod,obj.bc.fixnodes);
            obj.element = Element_Elastic_Micro;
            obj.variables = PhysicalVars_Elastic_2D_Micro(obj.dof.ndof);
            obj.solver = Solver.create(obj.mesh.scale);
        end
        
        function computeVariables(obj,vstrain)
            obj.element.computeLHS(obj.dim.nunkn,obj.dim.nstre,obj.mesh.nelem,obj.geometry,obj.material);
            obj.element.computeRHS(obj.dim.nunkn,obj.dim.nstre,obj.mesh.nelem,obj.geometry.nnode,obj.material,obj.bc,obj.dof.idx,obj.geometry,vstrain);
            
            % Assembly
            [obj.LHS,obj.RHS] = obj.Assemble(obj.element,obj.geometry.nnode,obj.dim.nunkn,obj.dof,obj.bc);
            
            %              comparison1(obj.LHS);
            %              comparison2(obj.RHS);
            
            % Solver            
            data.pnodes = obj.bc.pnodes;
            data.nunkn = obj.dim.nunkn;
            obj.solver.setSolverVariables(data);
            sol = obj.solver.solve(obj.variables.d_u,obj.LHS,obj.RHS,obj.dof);
            
            %             comparison3(sol);
            
            obj.variables = obj.variables.computeVars(sol,obj.dim,obj.geometry,obj.mesh.nelem,obj.dof.idx,obj.element,obj.material,vstrain);
            
            obj.StressHomog = obj.variables.stress_homog;
        end
        
        function obj = computeChomog(obj)
            obj.variables = PhysicalVars_Elastic_2D_Micro(obj.dof.ndof);
            vstrain=diag(ones(obj.dim.nstre,1));
            obj.variables.Chomog =  zeros(obj.dim.nstre,obj.dim.nstre);
            for istre=1:obj.dim.nstre
                obj.computeVariables(vstrain(istre,:));
                obj.variables.Chomog(:,istre) = obj.variables.stress_homog;
                obj.variables.tstrain(istre,:,:,:) = obj.variables.strain;
                obj.variables.tstress(istre,:,:,:) = obj.variables.stress;
            end
        end
    end
end

