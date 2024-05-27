classdef VademecumVariablesLoader < handle
    
    properties (Access = public)
        Ctensor
        Ptensor
        density
        structuredMesh        
    end
    
    properties (Access = private)
        fileName
        vademecumVariables
    end
    
    methods (Access = public)
        
        function obj = VademecumVariablesLoader(cParams)
            obj.init(cParams)
            obj.loadVademecumVariables();
            obj.createStructuredMesh();
            obj.computeConstitutive();        
            obj.computeDensityFromVademecum();            
            %obj.computePtensorFromVademecum();            
        end        
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.fileName = cParams.fileName;
        end
        
        function loadVademecumVariables(obj)
            fName = [obj.fileName,'WithAmplificators'];
            matFile   = [fName,'.mat'];
            file2load = fullfile('Vademecums',matFile);
            v = load(file2load);
            obj.vademecumVariables = v.d;
        end

        function createStructuredMesh(obj)
            s.x = obj.vademecumVariables.domVariables.mxV;
            s.y = obj.vademecumVariables.domVariables.myV;
            m = StructuredMesh(s);
            obj.structuredMesh = m;
        end        
        
        function computeConstitutive(obj)
            v     = obj.vademecumVariables;
            sMesh = obj.structuredMesh;
            for imx = 1:sMesh.nx
                for imy = 1:sMesh.ny
                    C(:,:,imx,imy) = v.variables{imx,imy}.('Ctensor');
                end
            end
            m = sMesh.mesh;
            for i = 1:size(C,1)
                for j = 1:size(C,2)
                    Cij = squeeze(C(i,j,:,:));
                    CijF = LagrangianFunction.create(m, 1, 'P1');
                    CijF.fValues  = Cij(:);
                    obj.Ctensor{i,j} = CijF;
                end
            end
        end

        function computeDensityFromVademecum(obj)
            v     = obj.vademecumVariables;
            sMesh = obj.structuredMesh;
            for imx = 1:sMesh.nx
                for imy = 1:sMesh.ny
                    rho(imx,imy) = v.variables{imx,imy}.('volume');
                end
            end
            m = sMesh.mesh;            
            rhoF = LagrangianFunction.create(m, 1, 'P1');
            rhoF.fValues = rho(:);
            obj.density  = rhoF;
        end

        function computePtensorFromVademecum(obj)
        %    s.vadVariables = obj.vademecumVariables;
        %    s.interpolator = obj.interpolator;
            pt = AmplificatorTensorFromVademecum(s);
            obj.Ptensor = pt;
        end             
        

        
    end
    
end