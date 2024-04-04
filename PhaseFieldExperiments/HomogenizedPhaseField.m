classdef HomogenizedPhaseField < handle
    
    properties (Access = private)
        fileName
        structuredMesh
        Ctensor
        microParams
    end
    
    methods (Access = public)
        
        function obj = HomogenizedPhaseField(cParams)
            obj.init(cParams)
            [mxV, C] = obj.loadVademecum();
            obj.createStructuredMesh(mxV);
            obj.createCtensorFunction(C);
        end

        function C = obtainTensor(obj)
            s.operation = @(xV) obj.evaluate(xV);
            s.ndimf = 9;
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

        function ddC = obtainTensor2Derivative(obj)
        end

        function setDesignVariable(obj,x)
            obj.microParams = x;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
           obj.fileName    = cParams.fileName;
           obj.microParams = cParams.microParams;
        end

        function [mxV, C] = loadVademecum(obj)
            matFile = [obj.fileName,'.mat'];
            file2load = fullfile('Vademecums',matFile);
            v = load(file2load);
            for ilV = 1:length(mxV)
                C(:,:,)
            end
        end

        function createStructuredMesh(obj,mxV)
            s.x = mxV;
            m = StructuredMesh1D(s);
            obj.structuredMesh = m;
        end

        function C = evaluate(obj,xV)
            [mL,cells] = obj.obtainLocalCoord(xV);
            nGaus = size(xV,2);
            nElem = obj.microParams{1}.mesh.nelem;
            nStre = size(obj.Ctensor,1);
            C = zeros(nStre,nStre,nGaus,nElem);
            for i = 1:nStre
                for j = 1:nStre
                    Cv = obj.Ctensor{i,j}.sampleFunction(mL,cells);
                    Cij(1,1,:,:) = reshape(Cv,nGaus,[]);
                    C(i,j,:,:)   = Cij(1,1,:,:);
                end
            end
        end

        function [mL, cells] = obtainLocalCoord(obj,xV)
            mx = obj.microParams{1};
            mG = mx.evaluate(xV);
            [mL, cells] = obj.structuredMesh.obtainLocalFromGlobalCoord(mG);

        end
        
    end
    
end