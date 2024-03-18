classdef VolumeFunctionalMicroParams < handle

    properties (Access = private)
        structuredMesh
        density  
        totalVolume        
        microParams        
    end

    properties (Access = private)
        fileName
        filter
        mesh
        volumeTarget
    end

    methods (Access = public)
        function obj = VolumeFunctionalMicroParams(cParams)
            obj.init(cParams);
            [mx,my,rho] = obj.loadVademecum();
            obj.createStructuredMesh(mx,my);         
            obj.createDensityFunction(rho)
            obj.createTotalVolume();
        end

        function [J,dJ] = computeFunctionAndGradient(obj,m)
            obj.microParams = m;
            J  = obj.computeFunction();
            dJm = obj.computeGradient();
            dJ = [];
            for ivar = 1:numel(dJm)
                dJ = [dJ;dJm{ivar}.fValues];                
            end              
        end
    end

    methods (Access = private)

        function init(obj,cParams)
           obj.mesh         = cParams.mesh;
           obj.fileName     = cParams.fileName;
           obj.filter       = cParams.filter;        
           obj.volumeTarget = cParams.volumeTarget;
        end

         function [mxV,myV,rho] = loadVademecum(obj)
            fName = [obj.fileName,'WithAmplificators'];
            matFile   = [fName,'.mat'];
            file2load = fullfile('Vademecums',matFile);
            v = load(file2load);
            var = v.d;
            mxV = var.domVariables.mxV;
            myV = var.domVariables.myV;
            rho = zeros(length(mxV),length(myV));
             for imx = 1:length(mxV)
                 for imy = 1:length(myV)
                     rho(imx,imy) = var.variables{imx,imy}.('volume');
                 end
             end
         end 

        function createStructuredMesh(obj,mxV,myV)
            s.x = mxV;
            s.y = myV;
            m = StructuredMesh(s);
            obj.structuredMesh = m;
        end       

        function createDensityFunction(obj,rhoV)
            m   = obj.structuredMesh.mesh;
            rho = LagrangianFunction.create(m, 1, 'P1');
            rho.fValues = rhoV(:);
            obj.density = rho;
        end

        function createTotalVolume(obj)
            m  = obj.mesh;            
            obj.totalVolume = m.computeVolume;
        end        

        function J = computeFunction(obj)
          s.operation = @(xV) obj.evaluateDensity(xV);
          rho         = DomainFunction(s);              
          m           = obj.mesh;
          volume      = Integrator.compute(rho,m,'QUADRATIC');
          tV          = obj.totalVolume;
          J           = volume/tV - obj.volumeTarget;
        end

        function dJ = computeGradient(obj)
          nVar = numel(obj.microParams.fun);
          dJ = cell(nVar,1);          
          for iVar = 1:nVar
            s.operation = @(xV) obj.evaluateVolumeDerivative(xV,iVar);
            s.ndimf = 1;
            dJ{iVar} =  DomainFunction(s);
          end
          dJ = obj.filter.compute(dJ,'QUADRATIC');
        end        

        function rho = evaluateDensity(obj,xV)
            [mL,cells] = obj.obtainLocalCoord(xV);
            nGaus = size(xV,2);
            ndimf = 1;
            nElem = obj.mesh.nelem;
            rho   = zeros(ndimf,nGaus,nElem);
            rhoV  = obj.density.sampleFunction(mL,cells);
            rho(1,:,:) = reshape(rhoV,nGaus,[]);
        end    

        function drho = evaluateDensityGradient(obj,xV,dir)
            [mL,cells] = obj.obtainLocalCoord(xV);
            nGaus = size(xV,2);
            nElem = obj.mesh.nelem;
            nDim  = obj.mesh.ndim;
            drhoV = obj.density.sampleGradient(mL,cells);
            drho  = squeezeParticular(drhoV,1);
            drho  = reshape(drho,nDim,nGaus,nElem);
            drho  = squeezeParticular(drho(dir,:,:),1);
        end        

        function [mL,cells] = obtainLocalCoord(obj,xV)
            mx = obj.microParams.fun{1};
            my = obj.microParams.fun{2};
            mxG = mx.evaluate(xV);
            myG = my.evaluate(xV);
            mG(:,1) = mxG(:);
            mG(:,2) = myG(:);
            sM      = obj.structuredMesh;
            [mL,cells] = sM.obtainLocalFromGlobalCoord(mG);
        end        

        function dJ = evaluateVolumeDerivative(obj,xV,iVar)
          drho = obj.evaluateDensityGradient(xV,iVar);
          tV   = obj.totalVolume;
          dJ   = drho/tV;
        end


    end

    methods (Static, Access = public)
        function title = getTitleToPlot()
            title = 'Volume';
        end
    end
end

