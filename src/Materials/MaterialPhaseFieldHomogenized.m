classdef MaterialPhaseFieldHomogenized < handle

    properties (Access = private)
        fileName
        structuredMesh
        Ctensor
        mesh
        phi
        degradation
        Gc
    end

    methods (Access = public)

        function obj = MaterialPhaseFieldHomogenized(cParams)
            obj.init(cParams)
            [mxV, C] = obj.loadVademecum();
            obj.computeFunctionsAndDerivatives(mxV,C);
        end

        function C = obtainTensor(obj,phi)
            obj.phi = phi;
            s.operation = @(xV) obj.evaluate(xV);
            s.ndimf = 6;
            s.mesh  = obj.mesh;
            C = DomainFunction(s);
        end

        function dC = obtainTensorDerivative(obj,phi)
            obj.phi = phi;
            s.operation = @(xV) obj.evaluateGradient(xV);
            s.ndimf = 6;
            s.mesh  = obj.mesh;
            dC =  DomainFunction(s);
        end

        function d2C = obtainTensorSecondDerivative(obj,phi)
            obj.phi = phi;
            s.operation = @(xV) obj.evaluateHessian(xV);
            s.ndimf = 6;
            s.mesh  = obj.mesh;
            d2C =  DomainFunction(s);
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
            file2load = fullfile('PFVademecum','Degradation',matFile);
            v = load(file2load);
            mxV = v.phi;
            C   = v.mat;
        end

        function C = evaluate(obj,xV)
            nStre = 3;
            nGaus = size(xV,2);
            nElem = obj.phi.mesh.nelem;
            C = zeros(nStre,nStre,nGaus,nElem);
            phiV = obj.phi.evaluate(xV);
            for i = 1:nStre
                for j = 1:nStre
                    C(i,j,:,:) = obj.degradation.fun{i,j}(phiV);
                end
            end
        end

        function dCt = evaluateGradient(obj,xV)
            nStre = 3;
            nGaus = size(xV,2);
            nElem = obj.phi.mesh.nelem;
            dCt = zeros(nStre,nStre,nGaus,nElem);
            phiV = obj.phi.evaluate(xV);
            for i = 1:nStre
                for j = 1:nStre
                    dCt(i,j,:,:) = obj.degradation.dfun{i,j}(phiV);
                end
            end
        end

        function d2Ct = evaluateHessian(obj,xV)
            nStre = 3;
            nGaus = size(xV,2);
            nElem = obj.phi.mesh.nelem;
            d2Ct = zeros(nStre,nStre,nGaus,nElem);
            phiV = obj.phi.evaluate(xV);
            for i = 1:nStre
                for j = 1:nStre
                    d2Ct(i,j,:,:) = obj.degradation.ddfun{i,j}(phiV);
                end
            end
        end

        function computeFunctionsAndDerivatives(obj,mxV,C)
            x = reshape(mxV,length(mxV),[]);
            y = C;
            nStre = size(C,1);

            fun   = cell(3,3);
            dfun  = cell(3,3);
            ddfun = cell(3,3);
            for i=1:nStre
                for j=1:nStre
                    degreePoly = 9;
                    fixedPointX = [0,1];
                    fixedPointY = [squeeze(y(i,j,1)),0];
                    coeffs = polyfix(x,squeeze(y(i,j,:)),degreePoly,fixedPointX,fixedPointY);
                    fun{i,j} = poly2sym(coeffs);
                    dfun{i,j} = diff(fun{i,j});
                    ddfun{i,j} = diff(dfun{i,j});
                    if all(coeffs)
                        obj.degradation.fun{i,j}   = matlabFunction(fun{i,j});
                        obj.degradation.dfun{i,j}  = matlabFunction(dfun{i,j});
                        obj.degradation.ddfun{i,j} = matlabFunction(ddfun{i,j});
                    else
                        obj.degradation.fun{i,j}   = @(x) x-x;
                        obj.degradation.dfun{i,j}  = @(x) x-x;
                        obj.degradation.ddfun{i,j} = @(x) x-x;
                    end

                end
            end
        end

    end

end