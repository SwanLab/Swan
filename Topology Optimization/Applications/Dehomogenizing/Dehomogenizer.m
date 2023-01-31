classdef Dehomogenizer < handle
    
    properties (Access = public)
       
    end
    
    properties (Access = private)
        dilation
        phi
        remesher
        levelSet
        uMesh
        nCell
        epsilon        
    end
    
    properties (Access = private)
        cellLevelSetParams
        backgroundMesh
        nCells
        theta
        mesh
    end
    
    methods (Access = public)
        
        function obj = Dehomogenizer(cParams)
            obj.init(cParams)            
        end

     function ls = compute(obj)           
            obj.computeDilation();
            obj.createMapping();
            obj.createLevelSet(); 
            nC = length(obj.nCells);
            ls = cell(nC,1);
            for iCell = 1:nC
                obj.nCell = obj.nCells(iCell);
                obj.createEpsilon();   
                obj.levelSet.computeLs(obj.epsilon);
                ls{iCell} = obj.levelSet.getValue();
             end            
     end

    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.backgroundMesh     = cParams.backgroundMesh;
            obj.nCells             = cParams.nCells;
            obj.theta              = cParams.theta;
            obj.cellLevelSetParams = cParams.cellLevelSetParams;            
            obj.mesh               = cParams.mesh;
            obj.remesher           = cParams.remesher;
        end
        
        function computeDilation(obj)
            s.theta = obj.theta;
            s.mesh  = obj.mesh;
            dC = DilationFieldComputer(s);
            d  = dC.compute();
            %dC.plot();
            obj.dilation = d;
        end      

        function createMapping(obj)
            s.mesh     = obj.mesh;
            s.theta    = obj.theta;
            s.dilation = obj.dilation;
            c = ConformalMappingComputer(s);
            phiV = c.compute();
            % c.plot();
            obj.phi = phiV;
        end
        
        function createEpsilon(obj)
            L = obj.mesh.computeCharacteristicLength();
            obj.epsilon = L/obj.nCell;
        end

        function createLevelSet(obj)
            s.coord  = obj.backgroundMesh.coord;            
            s.type   = 'periodicAndOriented';            
            s.backgroundMesh   = obj.backgroundMesh;
            s.mesh     = obj.mesh;
            s.remesher = obj.remesher;
            s.ndim     = 2;            
            s.phi      = obj.phi;            
            s.dilation = obj.dilation;
          %  s.epsilon = obj.epsilon;
            s.cellLevelSetParams = obj.cellLevelSetParams;
            lSet = LevelSetCreator.create(s);            
            obj.levelSet = lSet;   
        end             
        
    end
    
end