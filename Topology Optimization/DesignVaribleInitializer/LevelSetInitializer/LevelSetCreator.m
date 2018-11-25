classdef LevelSetCreator < handle
    
    properties (Access = protected)
        x
    end
    
    properties (Access = protected)
        mesh
        optimizerName
        nodeCoord
        levelSet
        lsSize
        ndim
    end
    
    properties (Access = private)
        optimizer
        scalar_product
    end
    
    methods (Access = public)
        
        function compute_initial_design(obj)
            obj.computeLevelSet();
            %xVal = obj.x;
            xVal = obj.levelSet;
            % !! PROVISIONAL !!
            if strcmp(obj.optimizerName,'SLERP') %|| strcmp(optimizer,'HAMILTON-JACOBI')
                sqrt_norma = obj.scalar_product.computeSP(xVal,xVal);
                xVal = xVal/sqrt(sqrt_norma);
                obj.levelSet = xVal;
            end
            
        end
        
        function x = getValue(obj)
            x = obj.levelSet;
        end
        
    end
    
    methods (Access = public, Static)
        
        function obj = create(settings,mesh,epsilon)
            factory = LevelSetFactory();
            input.settings = settings;
            input.mesh     = mesh;
            input.epsilon  = epsilon;
            input.nHoles   = settings.N_holes;
            input.rHoles   = settings.R_holes;
            input.phaseHoles = settings.phase_holes;
            input.warningHoleBC = settings.warningHoleBC;
            input.ndim = mesh.ndim;
            input.coord = mesh.coord;
            input.m = settings.widthSquare;
            input.m1 = settings.widthH;
            input.m2 = settings.widthV;
            obj = factory.create(settings.initial_case,input);
        end
    end
    
    methods (Access = protected)
        
        function obj = compute(obj,input)
            obj.mesh = input.mesh;
            obj.lsSize = size(input.coord(:,1));
            obj.ndim   = input.ndim;
            obj.scalar_product = ScalarProduct(input.settings.filename,input.epsilon);
            obj.optimizerName = input.settings.optimizer;
            obj.createNodalCoordinates(input.coord);
        end
        
        function computeDesignVariable(obj)
            phi = obj.levelSet;
            switch obj.optimizerName
                case {'SLERP','HAMILTON-JACOBI','PROJECTED SLERP'}
                    obj.x = phi;
                otherwise
                    obj.x = 1 - heaviside(phi);
            end
        end
    end
    
    methods (Access = private)
        
        function createNodalCoordinates(obj,coord)
            obj.nodeCoord = coord;
        end
        
    end
    
    methods (Abstract, Access = protected)
        x = computeLevelSet(obj)
    end
    
end

