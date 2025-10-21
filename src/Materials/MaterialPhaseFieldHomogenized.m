classdef MaterialPhaseFieldHomogenized < handle

    properties (Access = private)
        fileName
        degradation
        mesh
        young
    end

    methods (Access = public)

        function obj = MaterialPhaseFieldHomogenized(cParams)
            obj.init(cParams)
            obj.loadVademecum();
        end

        function C = obtainTensor(obj,phi)
            fun = obj.degradation.fun;
            s.operation = @(xV) obj.evaluate(phi,fun,xV);
            s.ndimf = 6;
            s.mesh  = obj.mesh;
            C = DomainFunction(s);
        end

        function dC = obtainTensorDerivative(obj,phi)
            fun = obj.degradation.dfun;
            s.operation = @(xV) obj.evaluate(phi,fun,xV);
            s.ndimf = 6;
            s.mesh  = obj.mesh;
            dC =  DomainFunction(s);
        end

        function d2C = obtainTensorSecondDerivative(obj,phi)
            fun = obj.degradation.ddfun;
            s.operation = @(xV) obj.evaluate(phi,fun,xV);
            s.ndimf = 6;
            s.mesh  = obj.mesh;
            d2C =  DomainFunction(s);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.fileName  = cParams.fileName;
            obj.mesh      = cParams.mesh;
            obj.young     = cParams.young;
        end

        function loadVademecum(obj)
            fName = [obj.fileName];
            matFile   = [fName,'.mat'];
            file2load = fullfile('TOVademecum','Interpolation',matFile);
            v = load(file2load);
            if isfield(v,'Interpolation')
                E = obj.young;
                nStre = size(v.Interpolation.fun,1);
                for i=1:nStre
                    for j=1:nStre
                        for k=1:nStre
                            for l=1:nStre
                                obj.degradation.fun{i,j,k,l} = @(x) E.*v.Interpolation.fun{i,j,k,l}(x);
                                obj.degradation.dfun{i,j,k,l} = @(x) E.*v.Interpolation.dfun{i,j,k,l}(x);
                                obj.degradation.ddfun{i,j,k,l} = @(x) E.*v.Interpolation.ddfun{i,j,k,l}(x);
                            end
                        end
                    end
                end
            else
                phi = v.phi; mat = v.mat;
                DHF = DamageHomogenizationFitter();
                [f,df,ddf] = DHF.computePolynomialFitting(9,phi,mat);
                obj.degradation.fun = f;
                obj.degradation.dfun = df;
                obj.degradation.ddfun = ddf;
            end
        end

        function C = evaluate(~,phi,fun,xV)
            nStre = size(fun,1);
            nGaus = size(xV,2);
            nElem = phi.fun.mesh.nelem;
            C = zeros(2,2,2,2,nGaus,nElem);
            phiV = phi.evaluate(xV);
            for i = 1:nStre
                for j = 1:nStre
                    for k=1:nStre
                        for l=1:nStre
                            C(i,j,k,l,:,:) = fun{i,j,k,l}(phiV);
                        end
                    end
                end
            end
        end

    end

end