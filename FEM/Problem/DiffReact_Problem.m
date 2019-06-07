classdef DiffReact_Problem < FEM
    
    properties
        material
    end
    
    properties (Access = protected)
        isRobinTermAdded
    end
    
    methods (Access = public)
        
        function obj = DiffReact_Problem(cParams)
            
            if isfield(cParams,'isRobinTermAdded')
                obj.isRobinTermAdded = cParams.isRobinTermAdded;
            else
                obj.isRobinTermAdded = false;
            end
            
            if ischar(cParams)
                obj.setupFromGiDFile(cParams);
            elseif isstruct(cParams)
                if isfield(cParams,'mesh')
                    obj.setupFromMesh(cParams);
                else
                    obj.setupFromGiDFile(cParams.fileName);
                end
            else
                error('Invalid input type');
            end
            obj.problemData.ptype = 'DIFF-REACT';
            obj.setScale();
        end
        
        function preProcess(obj)
            obj.createGeometry(obj.mesh);
            obj.setDOFs();
            obj.setElement();
            obj.solver = Solver.create();
        end
        
        function computeVariables(obj,x)
            if obj.isRobinTermAdded
                LHS = obj.element.computeLHS();
                x_reg = obj.solver.solve(LHS,x);
                obj.variables.x = x_reg;
            else
                bc = obj.element.getBcApplier();
                x_red  = bc.full_vector_2_reduced_vector(x);
                LHS = obj.element.computeLHS();
                x_reg = obj.solver.solve(LHS,x_red);
                obj.variables.x = bc.reduced_vector_2_full_vector(x_reg);
            end
        end
        
        
        function obj = setEpsilon(obj,epsilon)
            obj.element.setEpsilon(epsilon);
        end
        
        function createGeometry(obj,mesh)
            obj.geometry = Geometry(mesh,'LINEAR');
        end
        
    end
    
    methods (Access = protected)
        
        function setElement(obj)
            isRobinTermAdded = obj.isRobinTermAdded;
            obj.element = Element_DiffReact(obj.mesh,obj.geometry,obj.material,obj.dof,obj.problemData.scale,isRobinTermAdded);
        end
        
        function setDOFs(obj)
            obj.dof = DOF_DiffReact(obj.geometry);
        end
        
        function setScale(obj)
            obj.problemData.scale = 'MACRO';
        end
        
    end
    
    methods (Access = private)
        
        function setupFromGiDFile(obj,fileName)
            obj.inputReader.read(fileName);
            obj.createMesh();
            obj.problemData.fileName = fileName;
        end
        
        function setupFromMesh(obj,inputData)
            obj.mesh = inputData.mesh;
            if isfield(inputData,'fileName')
                obj.problemData.fileName = inputData.fileName;
            end
        end
        
    end
    
end
