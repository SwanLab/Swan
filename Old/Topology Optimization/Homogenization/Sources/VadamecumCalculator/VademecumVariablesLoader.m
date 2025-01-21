classdef VademecumVariablesLoader < handle
    
    properties (Access = public)
        Ctensor
        Ptensor
        density
    end
    
    properties (Access = private)
        fileName
        vadVariables
        interpolator
    end
    
    methods (Access = public)
        
        function obj = VademecumVariablesLoader(cParams)
            obj.init(cParams)
        end
        
        function load(obj)
            obj.loadVademecumVariables();
            obj.createInterpolator();
            obj.computeConstitutiveFromVademecum();
            obj.computePtensorFromVademecum();
            obj.computeDensityFromVademecum();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.fileName = cParams.fileName;
        end
        
        function loadVademecumVariables(obj)
            matFile   = [obj.fileName,'.mat'];
            file2load = fullfile('Vademecums',matFile);
            v = load(file2load);
            obj.vadVariables = v.d;
        end
        
        function createInterpolator(obj)
            sM.x = obj.vadVariables.domVariables.mxV;
            sM.y = obj.vadVariables.domVariables.myV;
            sI.mesh = StructuredMesh(sM);
            obj.interpolator = Interpolator(sI);
        end
        
        function computeConstitutiveFromVademecum(obj)
            s.vadVariables = obj.vadVariables;
            s.interpolator = obj.interpolator;
            ct = ConstitutiveTensorFromVademecum(s);
            obj.Ctensor = ct;
        end
        
        function computeDensityFromVademecum(obj)
            s.vadVariables = obj.vadVariables;
            s.interpolator = obj.interpolator;
            dt = DensityFromVademecum(s);
            obj.density = dt;
        end
        
        function computePtensorFromVademecum(obj)
            s.vadVariables = obj.vadVariables;
            s.interpolator = obj.interpolator;
            pt = AmplificatorTensorFromVademecum(s);
            obj.Ptensor = pt;
        end
        
    end
    
end