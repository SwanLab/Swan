classdef Physical_Problem<handle
    %Physical_Problem Summary of this class goes here
    %   Detailed explanation goes here
    
    %% Public GetAccess properties definition =============================
    properties (GetAccess = public, SetAccess = private)
        variables
    end
    
    %% Restricted properties definition ===================================
    properties (GetAccess = ?Postprocess, SetAccess = private)
        mesh
        geometry
        dim
        RHS
        LHS
    end
    
    %% Private properties definition ======================================
    properties (Access = private)
        bc
        material
        element
        dof
        solver
        physicalVars
        problemID
    end
    
    %% Public methods definition ==========================================
    methods (Access = public)
        function obj=Physical_Problem(problemID)
            obj.problemID=problemID;
        end
        function preProcess(obj)
            obj.mesh = Mesh(obj.problemID);
            
            % Create Objects
            obj.dim = DIM(obj.mesh.ptype,obj.mesh.pdim);
            obj.geometry = Geometry(obj.mesh);
            
            switch obj.mesh.ptype
                case 'ELASTIC'
                    % !! IT HAS BEEN ASSUMED THAT THERE'S ONLY ISOTROPIC MATERIALS. 
                    %    THIS HAS TO BE CHANGED FOR THE OPT TOP PROBLEM !!
                    switch obj.mesh.pdim
                        case '2D'
                            obj.material = Material_Elastic_ISO_2D(obj.mesh.nelem);
                            obj.physicalVars = PhysicalVars_Elastic_2D;
                            obj.element = Element_Elastic();
                            obj.element.B = B2;
                        case '3D'
                            obj.material = Material_Elastic_ISO_3D(obj.mesh.nelem);
                            obj.physicalVars = PhysicalVars_Elastic_3D;
                            obj.element = Element_Elastic();
                            obj.element.B = B3;
                    end
                    
                case 'THERMAL'
                    error('Still not implemented.')
                otherwise
                    error('Invalid ptype.')
            end
            obj.bc = BC(obj.dim.nunkn,obj.problemID);
            obj.dof = DOF(obj.geometry.nnode,obj.mesh.connec,obj.dim.nunkn,obj.mesh.npnod,obj.bc.fixnodes);
            obj.solver = Solver_Analytical;
        end
        
        function computeVariables(obj)
            % Create Element_Elastic object
            obj.element.computeLHS(obj.dim.nunkn,obj.dim.nstre,obj.mesh.nelem,obj.geometry,obj.material);
            obj.element.computeRHS(obj.dim.nunkn,obj.mesh.nelem,obj.geometry.nnode,obj.bc,obj.dof.idx);
            
            % Assembly
            [obj.LHS,obj.RHS] = obj.Assemble(obj.element,obj.geometry.nnode,obj.dim.nunkn,obj.dof);
            
            % Solver
            sol = obj.solver.solve(obj.LHS,obj.RHS,obj.dof,obj.bc.fixnodes);
            obj.variables=obj.physicalVars.computeVars(sol,obj.dim,obj.geometry.nnode,obj.mesh.nelem,obj.geometry.ngaus,obj.dof.idx,obj.element,obj.material);
        end
        
        function postProcess(obj)
            iter = 1; % static
            postprocess = Postprocess();
            postprocess.ToGid(obj.problemID,obj,iter);
            postprocess.ToGidPost(obj.problemID,obj,iter);
        end
        
        function setMatProps(obj,props)
            obj.material.setProps(props);
        end
    end
    
    %% Private methods definition =========================================
    methods (Access = private, Static)
        function [LHS,RHS] = Assemble(element,nnode,nunkn,dof)
            
            % Compute LHS
            LHS = sparse(dof.ndof,dof.ndof);
            for i = 1:nnode*nunkn
                for j = 1:nnode*nunkn
                    a = squeeze(element.LHS(i,j,:));
                    LHS = LHS + sparse(dof.idx(i,:),dof.idx(j,:),a,dof.ndof,dof.ndof);
                end
            end
            
            % Compute RHS
            RHS = sparse(dof.ndof,1);
            for i = 1:length(dof.idx(:,1)) % nnode*nunkn
                b = squeeze(element.RHS(i,1,:));
                ind = dof.idx(i,:);
                RHS(ind) = b;
            end
        end
    end
end

