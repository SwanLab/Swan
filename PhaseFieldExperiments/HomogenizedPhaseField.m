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
            [mxV, C] = obj.loadVademecum(); %%% Pending %%%
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
                s.ndimf = 9;
                dC{iVar} =  DomainFunction(s);
            end
        end

        function d2C = obtainTensorSecondDerivative(obj)
            nVar = numel(obj.microParams);
            d2C   = cell(nVar,1);
            for iVar = 1:nVar
                s.operation = @(xV) obj.evaluateHessian(xV,iVar);
                s.ndimf = 9;
                d2C{iVar} =  DomainFunction(s);
            end
        end

        function obj = setDesignVariable(obj,u,phi,t)
            obj.microParams{1} = phi;
        end             

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.fileName       = cParams.fileName;
            obj.microParams{1} = cParams.microParams;
        end

        function [mxV, C] = loadVademecum(obj)
            fName = [obj.fileName];
            matFile   = [fName,'.mat'];
            file2load = fullfile('VademecumDamage',matFile);
            v = load(file2load);
            mxV = v.alpha;
            C   = v.mat;
        end

        function createStructuredMesh(obj,mxV)
            s.x = mxV; %% [x1 x2 ...] %%%
            m = StructuredMesh1D(s);
            obj.structuredMesh = m;
        end

        function createCtensorFunction(obj,C)
            m = obj.structuredMesh.mesh;
            for i = 1:size(C,1)
                for j = 1:size(C,2)
                    Cij = squeeze(C(i,j,:));
                    CijF = LagrangianFunction.create(m, 1,'P1');
                    CijF.fValues = Cij;
                    Cij = CijF.project('P2');                                        
                    obj.Ctensor{i,j} = Cij;
                end
            end
        end

        function C = evaluate(obj,xV)
            [mL,cells] = obj.obtainLocalCoord(xV);
            nGaus = size(xV,2);
            nElem = obj.microParams{1}.mesh.nelem;
            nStre = size(obj.Ctensor,1);
            C = zeros(nStre,nStre,nGaus,nElem);
            for i = 1:nStre
                for j = 1:nStre
                    Cv = obj.Ctensor{i,j}.sampleFunction(mL',cells');
                    Cij(1,1,:,:) = reshape(Cv,nGaus,[]);
                    C(i,j,:,:)   = Cij(1,1,:,:);
                end
            end
        end

        function dCt = evaluateGradient(obj,xV,dir)
            [mL,cells] = obj.obtainLocalCoord(xV);
            nGaus = size(xV,2);
            nElem = obj.microParams{1}.mesh.nelem;
            nDim  = obj.microParams{1}.ndimf;
            nStre = size(obj.Ctensor,1);
            %  nDofs = size(mL,2);
            dC  = zeros(nDim,nStre,nStre,nGaus,nElem);
            for i = 1:nStre
                for j = 1:nStre
                    dCv = obj.Ctensor{i,j}.sampleGradient(mL',cells');
                    dCv = squeezeParticular(dCv,1);
                    dCij(:,1,1,:,:)  = reshape(dCv,nDim,nGaus,nElem);
                    dC(:,i,j,:,:) = dCij;
                end
            end
            dCt = squeezeParticular(dC(dir,:,:,:,:),1);
        end

        function d2Ct = evaluateHessian(obj,xV,dir)
            [mL,cells] = obj.obtainLocalCoord(xV);
            nGaus = size(xV,2);
            nElem = obj.microParams{1}.mesh.nelem;
            nDim  = obj.microParams{1}.ndimf;
            nStre = size(obj.Ctensor,1);
            d2C  = zeros(nDim,nStre,nStre,nGaus,nElem);
            for i = 1:nStre
                for j = 1:nStre
                    d2Cv = obj.Ctensor{i,j}.sampleHessian(mL',cells');
                    d2Cv = squeezeParticular(d2Cv,1);
                    d2Cij(:,1,1,:,:)  = reshape(d2Cv,nDim,nGaus,nElem);
                    d2C(:,i,j,:,:) = d2Cij;
                end
            end
            d2Ct = squeezeParticular(d2C(dir,:,:,:,:),1);
        end

        function [mL, cells] = obtainLocalCoord(obj,xV)
            mx = obj.microParams{1};
            mG = squeeze(mx.evaluate(xV));
            mG = reshape(mG,numel(mG),[]);
            [mL, cells] = obj.structuredMesh.obtainLocalFromGlobalCoord(mG);

        end

    end

end