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
            obj.mesh = Mesh(filename);
            switch obj.mesh.ptype
                case 'ELASTIC'
                    switch obj.mesh.pdim
                        case '2D'
                            obj.dim.ndim=2;
                            obj.dim.nunkn=2;
                            obj.dim.nstre=3;      
                        case '3D'
                            obj.dim.ndim=3;
                            obj.dim.nunkn=3;
                            obj.dim.nstre=6;      
                    end
            end 
            obj.bc = BC(obj.dim.nunkn,filename);
        end
        
        function obj = computeVariables(obj)
            % Create Geometry and DOF objects
            geometry = Geometry(obj.mesh);
            dof = DOF(geometry.nnode,obj.mesh.connec,obj.dim.nunkn,obj.mesh.npnod,obj.bc.displacements);
            
            % Create Element_Elastic object
            element = Element_Elastic();
            element.computeLHS(obj.dim.nstre,obj.dim.nunkn,obj.mesh.nelem,geometry);
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

