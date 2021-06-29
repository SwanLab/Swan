classdef DiffReact_Problem < FEM
    
    properties
        material
    end
    
    properties (Access = protected)
        isRobinTermAdded
        bcApplierType
       interp  
       boundaryMesh
    end
    
    methods (Access = public)
        
        function obj = DiffReact_Problem(cParams)
            
            if isfield(cParams,'isRobinTermAdded')
                obj.isRobinTermAdded = cParams.isRobinTermAdded;
            else
                obj.isRobinTermAdded = false;
            end
            
            if isfield(cParams,'bcApplierType')
                obj.bcApplierType = cParams.bcApplierType;
            else
                obj.bcApplierType = '';
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
            obj.createInterpolation();
        end
        
        function preProcess(obj)
            if contains(class(obj.mesh),'Total')
                obj.createGeometry(obj.mesh.innerMeshOLD);
            else
                obj.createGeometry(obj.mesh);
            end
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
                x_red  = bc.fullToReducedVector(x);
                LHS = obj.element.computeLHS();
                x_reg = obj.solver.solve(LHS,x_red);
                obj.variables.x = bc.reducedToFullVector(x_reg);
            end
        end
        

        function obj = setEpsilon(obj,epsilon)
            obj.element.setEpsilon(epsilon);
        end
        
        function createGeometry(obj,mesh)
            s.mesh = mesh;
            obj.geometry = Geometry.create(s);
        end
        
    end
    
    methods (Access = protected)
        
        function setElement(obj)
            isRobinTermAdded = obj.isRobinTermAdded;
            bcType = obj.bcApplierType;
            obj.element = Element_DiffReact(obj.mesh,obj.geometry,...
                obj.material,obj.dof,obj.problemData.scale,isRobinTermAdded,bcType,obj.interp,obj.boundaryMesh);
        end
        
        function setDOFs(obj)
            obj.dof = DOF_DiffReact(obj.mesh,obj.interp);
        end
        
        function setScale(obj)
            obj.problemData.scale = 'MACRO';
        end
        
    end
    
    methods (Access = private)
        
        function createInterpolation(obj)
            if contains(class(obj.mesh),'Total')
                obj.interp{1} =Interpolation.create(obj.mesh.innerMeshOLD,'LINEAR');
            else
                obj.interp{1} =Interpolation.create(obj.mesh,'LINEAR');
            end
            
        end
        
        function setupFromGiDFile(obj,fileName)
            obj.inputReader.read(fileName);
            obj.createMesh();
            obj.createBoundaryMesh(fileName);
            obj.problemData.fileName = fileName;
        end
        
        function setupFromMesh(obj,s)
            obj.mesh = s.mesh;
            if isfield(s,'fileName')
                obj.problemData.fileName = s.fileName;
                obj.createBoundaryMesh(s.fileName);
            end
            
        end
        
        function createBoundaryMesh(obj,fileName)
            run(fileName);
            if exist('External_border_nodes','var') && ~isempty(External_border_nodes)
                s.borderNodes    = External_border_nodes;
                if exist('External_border_elements','var') 
                    s.borderElements = External_border_elements;
                else
                    s.borderElements = [];
                end
                s.backgroundMesh = obj.mesh;
                s.type = 'FromData';
                b = BoundaryMeshCreator.create(s);
                obj.boundaryMesh = b.create();
            else
                s.backgroundMesh = obj.mesh;
                s.dimension = 1:obj.mesh.ndim;
                s.type = 'FromReactangularBox';
                bC = BoundaryMeshCreator.create(s);
                obj.boundaryMesh = bC.create();
            end
        end
        
        
    end
    
end
