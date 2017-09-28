classdef Physical_Problem<handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        mesh
        bc
        dim
        variables
    end
    
    methods
        function obj = preProcess(obj,filename)
            obj.dim.nunkn = 3;
            obj.mesh = Mesh(filename);
            obj.bc = BC(obj.dim.nunkn);
        end
        
        function obj = computeVariables(obj)
            % Create Geometry and DOF objects
            geometry = Geometry(obj.mesh);
            dof = DOF(geometry.nnode,obj.mesh.connec,obj.dim.nunkn,obj.mesh.npnod,obj.bc.displacements);
            
            % Create Element_Elastic object
            element = Element_Elastic();
            element.computeLHS(obj.dim.nunkn,obj.mesh.nelem,geometry);
            element.computeRHS(obj.dim.nunkn,obj.mesh.nelem,geometry.nnode,obj.bc,dof.idx);
            
            % Assembly
            [LHS,RHS] = Assemble.Compute(element,geometry.nnode,obj.dim.nunkn,dof);
            
            % Solver
            d_u = zeros(dof.ndof,1);
            d_u = Solver.analytical(d_u,LHS,RHS',dof.vR,dof.vL,obj.bc.displacements);
            obj.variables.displacement = d_u;
            
            % preguntar
            obj.dim.nnode = geometry.nnode;
            obj.dim.ngaus = geometry.ngaus;
        end
        
        function obj = postProcess(obj,filename)
            iter = 1; % static
            postprocess = Postprocess();
            postprocess.ToGid(filename, obj,iter);
            postprocess.ToGidPost(filename,obj,iter);
        end
    end
    
end

