classdef Physical_Problem < FEM
    %Physical_Problem Summary of this class goes here
    %   Detailed explanation goes here
    
    %% Public GetAccess properties definition =============================
    properties (GetAccess = public, SetAccess = public)
        variables
        mesh
        dim
        dof
        bc
        problemID
        element
    end
    
    
    %% Restricted properties definition ===================================
    properties (GetAccess = {?Postprocess,?Physical_Problem_Micro}, SetAccess = protected)
        material
        
    end
    
    %% Public methods definition ==========================================
    methods (Access = public)
        function obj = Physical_Problem(problemID)
            obj.problemID = problemID;
            obj.mesh = Mesh(obj.problemID);
            obj.dim = DIM(obj.mesh.ptype,obj.mesh.pdim);
            obj.geometry = Geometry(obj.mesh);
            obj.material = Material.create(obj.mesh.ptype,obj.mesh.pdim,obj.mesh.nelem,obj.mesh.connec,obj.geometry.cartd,obj.geometry.nnode,obj.mesh.coord);
            obj.dof = DOF(problemID,obj.geometry.nnode,obj.mesh.connec,obj.dim.nunkn,obj.mesh.npnod,obj.mesh.scale);    
        end
        
        function preProcess(obj)
            obj.element = Element.create(obj.mesh,obj.geometry,obj.material,obj.bc,obj.dof,obj.dim);
            obj.solver = Solver.create();
        end
        
        function computeVariables(obj)
            tol   = 1e-6;
            x0 = zeros(length(obj.dof.free),1);
            % Compute r & dr
            [r,dr] = obj.element.computeResidual(x0);
            while dot(r,r) > tol
                inc_x = obj.solver.solve(dr,-r);
                x = x0 + inc_x;
                % Compute r & dr
                [r,dr] = obj.element.computeResidual(x);
                x0 = x;
            end
            
            obj.variables = obj.element.computeVars(x);
        end
        
        function print(obj)
            iter = 1; % static
            postprocess = Postprocess_PhysicalProblem();
            results.physicalVars = obj.variables;
            postprocess.print(obj,obj.problemID,iter,results);
        end
        
        function postProcess(obj)
            %    ToDo
            % Inspire in TopOpt
            
        end
        
        function setDof(obj,dof)
            obj.dof = dof;
        end
        
        function setMatProps(obj,props)
            obj.element.material = obj.material.setProps(props);
        end
        
        function [K, M] = computeKM(obj,job)
                        
            % !! Hyper-mega-ultra provisional !!
            
            dim_smooth.nnode = obj.geometry.nnode;
            dim_smooth.nunkn = 1;
            dim_smooth.nstre = 2;
            
            mesh_smooth = obj.mesh;
            mesh_smooth.ptype = 'DIFF-REACT';
            mesh_smooth.scale = 'MACRO';
            
            bc_smooth = obj.bc;
            bc_smooth.fixnodes = [];
            
            dof_smooth = DOF(obj.problemID,obj.geometry.nnode,obj.mesh.connec,...
                dim_smooth.nunkn,obj.mesh.npnod,mesh_smooth.scale);
            
            dof_smooth.neumann = [];
            dof_smooth.dirichlet = [];
            dof_smooth.neumann_values = [];
            dof_smooth.dirichlet_values = [];
            dof_smooth.periodic_free = [];
            dof_smooth.periodic_constrained = [];
            dof_smooth.constrained = [];
            dof_smooth.free = setdiff(1:dof_smooth.ndof,dof_smooth.constrained);
            
            element_smooth = Element.create(mesh_smooth,obj.geometry,obj.material,obj.bc,dof_smooth,dim_smooth);
            
            [K] = element_smooth.computeStiffnessMatrix;
            [M] = element_smooth.computeMassMatrix(job);
        end
    end
    
    %% Private methods definition =========================================
    methods (Access = protected, Static)
        %         function [LHS,RHS] = Assemble(element,nnode,nunkn,dof,bc)
        %             % Compute LHS
        %             LHS = sparse(dof.ndof,dof.ndof);
        %             for i = 1:nnode*nunkn
        %                 for j = 1:nnode*nunkn
        %                     a = squeeze(element.LHS(i,j,:));
        %                     LHS = LHS + sparse(dof.idx(i,:),dof.idx(j,:),a,dof.ndof,dof.ndof);
        %                 end
        %             end
        %             LHS = 1/2 * (LHS + LHS');
        %             % Compute RHS
        %             RHS = zeros(dof.ndof,1);
        %             for i = 1:nnode*nunkn
        %                 b = squeeze(element.Fext(i,1,:));
        %                 ind = dof.idx(i,:);
        %                 RHS = RHS + sparse(ind,1,b',dof.ndof,1);
        %             end
        %
        %             %Compute Global Puntual Forces
        %             if ~isempty(bc.iN)
        %                 FextPoint = zeros(dof.ndof,1);
        %                 FextPoint(bc.iN) = bc.neunodes(:,3);
        %                 RHS = RHS + FextPoint;
        %             end
        %         end
    end
end

