classdef LevelSetCreator < handle
    
    properties (Access = public)
        x
    end
    
    properties (Access = protected)
        mesh
        optimizerName
        hole_value
        nodeCoord
        levelSet
    end
    
    properties (Access = private)
        ini_design_value
        optimizer
        scalar_product
    end
    
    methods (Access = public)
        
        function xVal = compute_initial_design(obj)
            obj.computeInitialLevelSet();
            xVal = obj.x;
            % !! PROVISIONAL !!
            if strcmp(obj.optimizerName,'SLERP') %|| strcmp(optimizer,'HAMILTON-JACOBI')
                sqrt_norma = obj.scalar_product.computeSP(xVal,xVal);
                xVal = xVal/sqrt(sqrt_norma);
            end
                       
        end
        
    end
    
    methods (Access = public, Static)
        
        function obj = create(settings,mesh,epsilon)
            factory = LevelSetFactory();
            input.settings = settings;
            input.mesh     = mesh;
            input.epsilon  = epsilon;
            
            obj = factory.create(settings.initial_case,input);
            
        end
    end
    
    methods (Access = protected)
        
        function obj = compute(obj,input)
            obj.mesh = input.mesh;
            obj.optimizerName = input.settings.optimizer;
            obj.scalar_product = ScalarProduct(input.settings.filename,input.epsilon);
            obj.createNodalCoordinates(input.mesh);
            obj.computeInitialValue()
        end
    end
    
    methods (Access = private)
        
        function createNodalCoordinates(obj,mesh)
            coord.x = mesh.coord(:,1);
            coord.y = mesh.coord(:,2);
            if strcmp(mesh.pdim,'3D')
                coord.z = mesh.coord(:,3);
            end
            obj.nodeCoord = coord;
        end
        
        function computeInitialValue(obj)
            obj.setInitialValues();
            geometry = Geometry(obj.mesh,'LINEAR');
            obj.x = obj.ini_design_value*ones(geometry.interpolation.npnod,1);
        end
        
        function setInitialValues(obj)
            switch obj.optimizerName
                case {'SLERP', 'PROJECTED SLERP'}
                    obj.ini_design_value = -1.015243959022692;
                    obj.hole_value = 0.507621979511346;
                case 'HAMILTON-JACOBI'
                    obj.ini_design_value = -0.1;
                    obj.hole_value = 0.1;
                otherwise
                    obj.ini_design_value = 1;
                    obj.hole_value = 0;
            end
        end
        
    end
    
    methods (Abstract, Access = protected)
        x = computeInitialLevelSet(obj)
    end
    
end

