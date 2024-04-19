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
            C{1} = DomainFunction(s);
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

            %v.phi(end+1) = 1;
            %v.mat(:,:,end+1) = zeros(3,3);

            mxV = v.phi;
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
                    %Cij = CijF.project('P2');                                        
                    obj.Ctensor{i,j} = CijF;
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
            nStre = size(obj.Ctensor,1);
            dCt = zeros(nStre,nStre,nGaus,nElem);
            for i = 1:nStre
                for j = 1:nStre
                    dCtensorDomF = Grad(obj.Ctensor{i,j});
                    dCtensor = dCtensorDomF.project('P1',obj.structuredMesh.mesh);

                    % for k=1:10
                    %     ss.filterType = 'LUMP';
                    %     ss.mesh       = obj.structuredMesh.mesh;
                    %     ss.trial      = LagrangianFunction.create(obj.structuredMesh.mesh,1,'P1');
                    %     filter        = Filter.create(ss);
                    %     dCtensor       = filter.compute(dCtensor,'QUADRATIC');
                    % end

                    dCv = dCtensor.sampleFunction(mL',cells');
                    dCij(1,1,:,:) = reshape(dCv,nGaus,[]);
                    dCt(i,j,:,:)   = dCij(1,1,:,:);
                end
            end
        end

        function d2Ct = evaluateHessian(obj,xV,dir)
            [mL,cells] = obj.obtainLocalCoord(xV);
            nGaus = size(xV,2);
            nElem = obj.microParams{1}.mesh.nelem;
            nStre = size(obj.Ctensor,1);
            d2Ct = zeros(nStre,nStre,nGaus,nElem);
            for i = 1:nStre
                for j = 1:nStre
                    dCtensorDomF = Grad(obj.Ctensor{i,j});
                    dCtensor = dCtensorDomF.project('P1',obj.structuredMesh.mesh);
                    d2CtensorDomF = Grad(dCtensor);
                    d2Ctensor = d2CtensorDomF.project('P1',obj.structuredMesh.mesh);

                    for k=1:10
                        ss.filterType = 'LUMP';
                        ss.mesh       = obj.structuredMesh.mesh;
                        ss.trial      = LagrangianFunction.create(obj.structuredMesh.mesh,1,'P1');
                        filter        = Filter.create(ss);
                        d2Ctensor       = filter.compute(d2Ctensor,'QUADRATIC');
                    end

                    d2Cv = d2Ctensor.sampleFunction(mL',cells');
                    d2Cij(1,1,:,:) = reshape(d2Cv,nGaus,[]);
                    d2Ct(i,j,:,:)   = d2Cij(1,1,:,:);
                end
            end
        end

        function [mL, cells] = obtainLocalCoord(obj,xV)
            mx = obj.microParams{1};
            mG = squeeze(mx.evaluate(xV));
            mG = reshape(mG,numel(mG),[]);
            [mL, cells] = obj.structuredMesh.obtainLocalFromGlobalCoord(mG);

        end

    end

end