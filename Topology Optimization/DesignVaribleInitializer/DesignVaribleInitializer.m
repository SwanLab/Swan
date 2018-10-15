classdef DesignVaribleInitializer < handle
    properties
        x
        ini_design_value
        hole_value
        mesh
        optimizer
        scalar_product
    end
    
    methods (Static)
        function obj = create(settings,mesh,epsilon)
            switch settings.initial_case
                case 'circle'
                    obj = DesignVaribleInitializer_Circle(settings,mesh,epsilon);
                case 'sphere'
                    obj = DesignVaribleInitializer_Sphere(settings,mesh,epsilon);
                case 'cylinder'
                    obj = DesignVaribleInitializer_Cylinder(settings,mesh,epsilon);
                case 'horizontal'
                    obj = DesignVaribleInitializer_Horizontal(settings,mesh,epsilon);
                case 'square'
                    obj = DesignVaribleInitializer_Square(settings,mesh,epsilon);
                case 'feasible'
                    obj = DesignVaribleInitializer_Feasible(settings,mesh,epsilon);
                case 'rand'
                    obj = DesignVaribleInitializer_Random(settings,mesh,epsilon);
                case 'holes'
                    obj = DesignVaribleInitializer_Holes(settings,mesh,epsilon);
                case 'full'
                    obj = DesignVaribleInitializer_Full(settings,mesh,epsilon);
                case 'orientedFiber'
                    %To be changed! Not appropiate way of creating
                    %sublcasses
                    %obj = DesignVaribleInitializer_orientedFiber(settings,mesh,epsilon,direction);
                otherwise
                    error('Invalid initial value of design variable.');
            end
        end
    end
    
    methods
        function obj = DesignVaribleInitializer(settings,mesh,epsilon)
            geometry = Geometry(mesh,'LINEAR');
            obj.mesh = mesh;
            obj.optimizer = settings.optimizer;
            
            obj.setValues;
            obj.x = obj.ini_design_value*ones(geometry.interpolation.npnod,1);
            obj.scalar_product = ScalarProduct(settings.filename,epsilon);
        end
        
        function x = compute_initial_design(obj)
            x = obj.compute_initial_x;
            
            % !! PROVISIONAL !!
            if strcmp(obj.optimizer,'SLERP') %|| strcmp(optimizer,'HAMILTON-JACOBI')
                sqrt_norma = obj.scalar_product.computeSP(x,x);
                x = x/sqrt(sqrt_norma);
            end
        end
    end
    
    methods (Abstract)
        x = compute_initial_x(obj)
    end
    
    methods (Access = private)
        function setValues(obj)
            switch obj.optimizer
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
end

