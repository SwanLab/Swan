classdef testVigdergauzMicroStructure < testShowingError & testTopOptComputation
    
    properties (Access = private)
        designVariable
        unfittedMesh
        volume
        numericalVolume
    end
    
    properties (Access = protected)
        testName = 'test_VigergauzMicroStructure';
        tol
    end
    
    methods (Access = public)
        
        function obj = testVigdergauzMicroStructure()
            obj.init();
            obj.createUnfittedMesh();
            obj.computeFractionVolume();
        end
        
    end
    
    methods (Access = protected)
        
        function selectComputedVar(obj)
            obj.designVariable = obj.topOpt.designVariable;
        end
        
        function computeError(obj)
            v = obj.volume;
            nv = obj.numericalVolume;
            obj.error = abs(v - nv)/v;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj)
            obj.tol = 2*1e-1;
            obj.volume = 0.6;
        end
        
        function createUnfittedMesh(obj)
            meshBackground = obj.topOpt.designVariable.mesh;
            interpolation = Interpolation.create(meshBackground,'LINEAR');
            s.unfittedType = 'INTERIOR';
            s.meshBackground = meshBackground;
            s.interpolationBackground = interpolation;
            s.includeBoxContour = false;
            cParams = SettingsMeshUnfitted(s);
            obj.unfittedMesh = UnfittedMesh(cParams);
            obj.unfittedMesh.compute(obj.designVariable.value);
        end
        
        function computeFractionVolume(obj)
            s.mesh = obj.unfittedMesh;
            s.type = 'COMPOSITE';
            s = obj.createInteriorParams(s,s.mesh);            
            integrator = Integrator.create(s);
            nodalF = ones(size(obj.designVariable.value));            
            xP = integrator.integrateAndSum(nodalF);                                    
            obj.numericalVolume = sum(xP);
        end
        
       function cParams = createInteriorParams(obj,cParams,mesh)
            cParamsInnerCut = obj.createInnerCutParams(mesh);
            cParams.compositeParams{1} = cParamsInnerCut;
            if mesh.innerMesh.nelem ~= 0
            cParamsInner = obj.createInnerParams(mesh);
            cParams.compositeParams{2} = cParamsInner;
            end
        end        
        
        function cParams = createInnerParams(obj,mesh)
            cParams.mesh = mesh.innerMesh;
            cParams.type = 'SIMPLE';
            cParams.globalConnec = mesh.globalConnec;
            cParams.npnod = mesh.innerMesh.npnod;
            cParams.backgroundMesh = obj.topOpt.designVariable.mesh;
            cParams.innerToBackground = mesh.backgroundFullCells;
        end        
        
        function cParams = createInnerCutParams(obj,mesh)
            cParams.mesh = mesh.innerCutMesh; 
            cParams.type = 'CutMesh';
            cParams.meshBackground = obj.topOpt.designVariable.mesh;
          
        end               
        
    end
    
    
end