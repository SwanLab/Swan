classdef DensityCreatorByLevelSet < DensityCreator
    
    properties (Access = private)
        levelSet
        Mesh        
        Settings
        direction
        NumFibers
        filter
    end
    
    methods (Access = public)
      
      function obj = DensityCreatorByLevelSet(mesh,settings,dir,numFibers)
          obj.init(mesh,settings,dir,numFibers)
          obj.createLevelSet();
          obj.createFilter();
          obj.computeDensity()
      end
      
      function ls = getLevelSet(obj)
          ls = obj.levelSet;
      end

    end
    
    methods (Access = private)
        
        function init(obj,Mesh,Settings,dir,NumFibers)
           obj.Mesh = Mesh;     
           obj.Settings = Settings;
           obj.direction = dir;
           obj.NumFibers = NumFibers;
        end
        
        function createLevelSet(obj)
            epsilon = obj.Mesh.mean_cell_size;
            LS_initializer = DesignVaribleInitializer_orientedFiber(...
                obj.Settings,obj.Mesh,epsilon,...
                obj.direction,obj.NumFibers);
            LS_initializer.compute_initial_design();
            obj.levelSet = LS_initializer.x;
        end
        
        function createFilter(obj)
            dim = obj.Settings.ptype;
            fileName = obj.Settings.filename;
            obj.filter = Filter_P1_LevelSet_2D(fileName,dim);
            obj.filter.preProcess();
        end
        
        function computeDensity(obj)
            obj.density = obj.filter.getP0fromP1(obj.levelSet);            
        end
        
    end
    
end

