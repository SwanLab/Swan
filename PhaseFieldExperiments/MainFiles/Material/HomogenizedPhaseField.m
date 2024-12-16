classdef HomogenizedPhaseField < handle

    properties (Access = private)
        fileName
        structuredMesh
        Ctensor
        microParams
        mesh

        degradation

        Gc
    end

    methods (Access = public)

        function obj = HomogenizedPhaseField(cParams)
            obj.init(cParams)
            [mxV, C] = obj.loadVademecum();
            obj.createStructuredMesh(mxV);
            obj.createCtensorFunction(C);
            obj.computeFunctionsAndDerivatives(mxV,C);
        end

        function C = obtainTensor(obj,phi)
            obj.setDesignVariable(phi)
            s.operation = @(xV) obj.evaluate(xV);
            s.ndimf = 6;
            s.mesh  = obj.mesh;
            C{1} = DomainFunction(s);
        end

        function dC = obtainTensorDerivative(obj,phi)
            obj.setDesignVariable(phi)
            nVar = numel(obj.microParams);
            dC   = cell(nVar,1);
            for iVar = 1:nVar
                s.operation = @(xV) obj.evaluateGradient(xV,iVar);
                s.ndimf = 6;
                s.mesh  = obj.mesh;
                dC{iVar} =  DomainFunction(s);
            end
        end

        function d2C = obtainTensorSecondDerivative(obj,phi)
            obj.setDesignVariable(phi)
            nVar = numel(obj.microParams);
            d2C   = cell(nVar,1);
            for iVar = 1:nVar
                s.operation = @(xV) obj.evaluateHessian(xV,iVar);
                s.ndimf = 6;
                s.mesh  = obj.mesh;
                d2C{iVar} =  DomainFunction(s);
            end
        end

        function setDesignVariable(obj,phi)
            obj.microParams{1} = phi;
        end             

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.fileName  = cParams.fileName;
            obj.mesh      = cParams.mesh;
        end

        function [mxV, C] = loadVademecum(obj)
            fName = [obj.fileName];
            matFile   = [fName,'.mat'];
            file2load = fullfile('VademecumDamage',matFile);
            v = load(file2load);
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
                    CijF.setFValues(Cij);
                    %Cij = CijF.project('P2');                                        
                    obj.Ctensor{i,j} = CijF;
                end
            end
        end

        function C = evaluate(obj,xV)
            %[mL,cells] = obj.obtainLocalCoord(xV);
            nGaus = size(xV,2);
            nElem = obj.microParams{1}.mesh.nelem;
            nStre = size(obj.Ctensor,1);
            C = zeros(nStre,nStre,nGaus,nElem);

            phiV = obj.microParams{1}.evaluate(xV);
            for i = 1:nStre
                for j = 1:nStre
                    % Cv = obj.Ctensor{i,j}.sampleFunction(mL',cells');
                    % Cij(1,1,:,:) = reshape(Cv,nGaus,[]);
                    % C(i,j,:,:)   = Cij(1,1,:,:);
                    
                    %C(i,j,:,:) = double(subs(obj.fun{i,j},phiV));
                    C(i,j,:,:) = obj.degradation.fun{i,j}(phiV);
                end
            end
        end

        function dCt = evaluateGradient(obj,xV,dir)         
            %[mL,cells] = obj.obtainLocalCoord(xV);
            nGaus = size(xV,2);
            nElem = obj.microParams{1}.mesh.nelem;
            nStre = size(obj.Ctensor,1);
            dCt = zeros(nStre,nStre,nGaus,nElem);

            phiV = obj.microParams{1}.evaluate(xV);
            for i = 1:nStre
                for j = 1:nStre
                    % dCtensorDomF = Grad(obj.Ctensor{i,j});
                    % dCtensor = dCtensorDomF.project('P1',obj.structuredMesh.mesh);
                    % 
                    % % for k=1:10
                    % %     ss.filterType = 'LUMP';
                    % %     ss.mesh       = obj.structuredMesh.mesh;
                    % %     ss.trial      = LagrangianFunction.create(obj.structuredMesh.mesh,1,'P1');
                    % %     filter        = Filter.create(ss);
                    % %     dCtensor       = filter.compute(dCtensor,'QUADRATIC');
                    % % end
                    % 
                    % dCv = dCtensor.sampleFunction(mL',cells');
                    % dCij(1,1,:,:) = reshape(dCv,nGaus,[]);
                    % dCt(i,j,:,:)   = dCij(1,1,:,:);

                    %dCt(i,j,:,:) = double(subs(obj.dfun{i,j},phiV)); WITH
                    %SYM
                    dCt(i,j,:,:) = obj.degradation.dfun{i,j}(phiV);
                    
                end
            end
        end

        function d2Ct = evaluateHessian(obj,xV,dir)
            %[mL,cells] = obj.obtainLocalCoord(xV);
            nGaus = size(xV,2);
            nElem = obj.microParams{1}.mesh.nelem;
            nStre = size(obj.Ctensor,1);
            d2Ct = zeros(nStre,nStre,nGaus,nElem);

            phiV = obj.microParams{1}.evaluate(xV);
            for i = 1:nStre
                for j = 1:nStre
                    % dCtensorDomF = Grad(obj.Ctensor{i,j});
                    % dCtensor = dCtensorDomF.project('P1',obj.structuredMesh.mesh);
                    % d2CtensorDomF = Grad(dCtensor);
                    % d2Ctensor = d2CtensorDomF.project('P1',obj.structuredMesh.mesh);
                    % 
                    % for k=1:10
                    %     ss.filterType = 'LUMP';
                    %     ss.mesh       = obj.structuredMesh.mesh;
                    %     ss.trial      = LagrangianFunction.create(obj.structuredMesh.mesh,1,'P1');
                    %     filter        = Filter.create(ss);
                    %     d2Ctensor       = filter.compute(d2Ctensor,'QUADRATIC');
                    % end
                    % 
                    % d2Cv = d2Ctensor.sampleFunction(mL',cells');
                    % d2Cij(1,1,:,:) = reshape(d2Cv,nGaus,[]);
                    % d2Ct(i,j,:,:)   = d2Cij(1,1,:,:);

                    %d2Ct(i,j,:,:) = double(subs(obj.ddfun{i,j},phiV)); 
                    d2Ct(i,j,:,:) = obj.degradation.ddfun{i,j}(phiV);
                end
            end
        end

        function [mL, cells] = obtainLocalCoord(obj,xV)
            mx = obj.microParams{1};
            mG = squeeze(mx.evaluate(xV));
            mG = reshape(mG,numel(mG),[]);
            [mL, cells] = obj.structuredMesh.obtainLocalFromGlobalCoord(mG);

        end

        function computeFunctionsAndDerivatives(obj,mxV,C)
            x = reshape(mxV,length(mxV),[]);
            y = C;

            fun   = cell(3,3);
            dfun  = cell(3,3);
            ddfun = cell(3,3);
            for i=1:3
                for j=1:3
                    % f = fit(x,squeeze(y(i,j,:)),'poly9');
                    % fun{i,j} = poly2sym(coeffvalues(f));
                    coeffs = polyfix(x,squeeze(y(i,j,:)),9,[0,1],[squeeze(y(i,j,1)),0]);
                    fun{i,j} = poly2sym(coeffs);
                    dfun{i,j} = diff(fun{i,j});
                    ddfun{i,j} = diff(dfun{i,j});
                    if all(coeffs)
                        obj.degradation.fun{i,j} = matlabFunction(fun{i,j});
                        obj.degradation.dfun{i,j} = matlabFunction(dfun{i,j});
                        obj.degradation.ddfun{i,j} = matlabFunction(ddfun{i,j});
                    else
                        obj.degradation.fun{i,j} = @(x) x-x;
                        obj.degradation.dfun{i,j} = @(x) x-x;
                        obj.degradation.ddfun{i,j} = @(x) x-x;
                    end

                end
            end
        end

    end

end