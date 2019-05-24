classdef Elastic_Problem < FEM
    %Elastic_Problem Summary of this class goes here
    % Detailed explanation goes here
    
    %% Public GetAccess properties definition =============================
    properties (GetAccess = public, SetAccess = public)
    end
    
    
    %% Public methods definition ==========================================
    methods (Access = public)
        function obj = Elastic_Problem(fileName)
            obj.problemData.fileName = fileName;
            obj.mesh = Mesh_GiD(fileName); % Mesh defined twice, but almost free
            obj.createGeometry(obj.mesh);
            obj.dof = DOF_Elastic(fileName,obj.geometry,obj.mesh);
        end
        
        function preProcess(obj)
            cParams.ptype = obj.mesh.ptype;
            cParams.pdim  = obj.mesh.pdim;
            cParams.nelem = obj.geometry(1).interpolation.nelem;
            cParams.geometry = obj.geometry;
            cParams.mesh  = obj.mesh;            
            material = Material.create(cParams);
            obj.element = Element_Elastic.create(obj.mesh,obj.geometry,material,obj.dof);
            obj.solver = Solver.create;
        end
        
        function computeVariables(obj)
            Kred = obj.element.computeLHS;
            fext_red = obj.element.computeRHS;
            u = obj.solver.solve(Kred,fext_red);
            obj.variables = obj.element.computeVars(u);
        end
        
%         function print(obj)
%             postprocess = Postprocess_PhysicalProblem;
%             results.physicalVars = obj.variables;
%             postprocess.print(obj,obj.problemID,results);
%         end
        
        function postProcess(obj)
            % ToDo
            % Inspire in TopOpt
        end
        
        % !! THIS SHOULD BE DEFINED BY THE USER !!
        function createGeometry(obj,mesh)
            obj.geometry = Geometry(mesh,'LINEAR');
        end
    end
end

