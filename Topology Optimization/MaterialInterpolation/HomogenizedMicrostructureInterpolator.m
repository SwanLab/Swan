classdef HomogenizedMicrostructureInterpolator < Material
    
    properties (Access = private)
        fileName
        structuredMesh
        Ctensor
        microParams
    end
    
    methods (Access = public)
        
        function obj = HomogenizedMicrostructureInterpolator(cParams)
            obj.init(cParams);
            obj.Ctensor        = cParams.Ctensor;
            obj.structuredMesh = cParams.structuredMesh; 
        end

        function C = obtainTensor(obj)
          s.operation = @(xV) obj.evaluate(xV);
          C = DomainFunction(s);            
        end

        function dC = obtainTensorDerivative(obj)
          nVar = numel(obj.microParams);
          dC   = cell(nVar,1);
            for iVar = 1:nVar
                s.operation = @(xV) obj.evaluateGradient(xV,iVar);
                s.ndimf = 1;
                dC{iVar} =  DomainFunction(s);
            end
        end

        function setDesignVariable(obj,x)
            obj.microParams = x;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
           obj.microParams = cParams.microParams;
           obj.Ctensor     = cParams.Ctensor;
           obj.structuredMesh = cParams.structuredMesh;
        end

        
        function C = evaluate(obj,xV)
            [mL,cells] = obj.obtainLocalCoord(xV);
            nGaus = size(xV,2);
            nElem = obj.microParams{1}.mesh.nelem;
            nStre = size(obj.Ctensor,1); 
          %  nDofs = size(mL,2);
            C  = zeros(nStre,nStre,nGaus,nElem);
            for i = 1:nStre
                for j = 1:nStre
                    Cv = obj.Ctensor{i,j}.sampleFunction(mL,cells);
                    Cij(1,1,:,:) = reshape(Cv,nGaus,[]);
                    C(i,j,:,:)   = Cij(1,1,:,:);
                end
            end
        end

        function dCt = evaluateGradient(obj,xV,dir)
            [mL,cells] = obj.obtainLocalCoord(xV);
            nGaus = size(xV,2);
            nElem = obj.microParams{1}.mesh.nelem;
            nDim  = obj.microParams{1}.mesh.ndim;
            nStre = size(obj.Ctensor,1);
            %  nDofs = size(mL,2);
            dC  = zeros(nDim,nStre,nStre,nGaus,nElem);
            for i = 1:nStre
                for j = 1:nStre
                    dCv = obj.Ctensor{i,j}.sampleGradient(mL,cells);
                    dCv = squeezeParticular(dCv,1);
                    dCij(:,1,1,:,:)  = reshape(dCv,nDim,nGaus,nElem);
                    dC(:,i,j,:,:) = dCij;
                end
            end     
            dCt = squeezeParticular(dC(dir,:,:,:,:),1);
        end


        function [mL,cells] = obtainLocalCoord(obj,xV)
            mx = obj.microParams{1};
            my = obj.microParams{2};
            mxG = mx.evaluate(xV);
            myG = my.evaluate(xV);
            mG(:,1) = mxG(:);
            mG(:,2) = myG(:);
            [mL,cells] = obj.structuredMesh.obtainLocalFromGlobalCoord(mG);
        end
        
    end
    
    methods (Access = private, Static)
        
        function xR = reshapeDesignVariable(x)
            xV = x.value;
            nV = x.nVariables;
            nx = length(xV)/nV;
            xR = cell(nV,1);
            for ivar = 1:nV
                i0 = nx*(ivar-1) + 1;
                iF = nx*ivar;
                xR{ivar} = xV(i0:iF);
            end
        end
        
    end
    
end
