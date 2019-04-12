classdef ConstitutiveTensorFromVademecum < handle
    
    properties (Access = private)
        fileName
        loadVariables
        domain
        CtensorValues
        interpolator
    end
    
    methods (Access = public)
        
        function obj = ConstitutiveTensorFromVademecum(cParams)
            obj.init(cParams);
            obj.loadVademecumVariables();
            obj.obtainVademecumDomain();            
            obj.obtainConstitutiveTensorValues();
            obj.createInterpolator();
        end
        
        function C = computeCtensor(obj,x)
            obj.interpolator.setValues(x(:,1),x(:,2));            
            nstre = size(obj.CtensorValues,1);            
            C = zeros(nstre,nstre,size(x,1));
            for i = 1:nstre
                for j = 1:nstre
                    cv = squeeze(obj.CtensorValues(i,j,:,:));
                    c = obj.interpolate(cv);
                    C(i,j,:) = c;
                end
            end
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.fileName = cParams.fileName;            
        end
        
        function loadVademecumVariables(obj)
            matFile   = [obj.fileName,'.mat'];
            file2load = fullfile('Output',obj.fileName,matFile);
            v = load(file2load);
            obj.loadVariables = v.d;
        end
        
        function obtainVademecumDomain(obj)
            obj.domain.mxV = obj.loadVariables.domVariables.mxV;
            obj.domain.myV = obj.loadVariables.domVariables.myV;
        end
        
        function obtainConstitutiveTensorValues(obj)
            var = obj.loadVariables.variables;
            for imx = 1:length(obj.domain.mxV)
                for imy = 1:length(obj.domain.myV)
                    C(:,:,imx,imy) = var{imx,imy}.Ctensor;
                end
            end
            obj.CtensorValues = C;
        end
        
        function z = interpolate(obj,zG)            
            z = obj.interpolator.interpolate(zG);
        end
        
        function createInterpolator(obj)
            sM.x = obj.domain.mxV;
            sM.y = obj.domain.myV;            
            sI.mesh = StructuredMesh(sM);    
            obj.interpolator = Interpolator(sI);            
        end
        
    end
    
    
end
