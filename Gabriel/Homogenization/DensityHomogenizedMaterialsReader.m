classdef DensityHomogenizedMaterialsReader < handle

    properties (Access = private)
        fileName
        homogenization
        mesh
        young
    end

    methods (Access = public)

        function obj = DensityHomogenizedMaterialsReader(cParams)
            obj.init(cParams)
            obj.loadVademecum();
        end

        function C = obtainTensor(obj,rho)
            fun = obj.homogenization.fun;
            s.operation = @(xV) obj.evaluate(rho,fun,xV);
            s.ndimf = 6;
            s.mesh = obj.mesh;
            C = DomainFunction(s);
        end

        function dC = obtainTensorDerivative(obj,rho)
            fun = obj.homogenization.dfun;
            s.operation = @(xV) obj.evaluate(rho,fun,xV);
            s.ndimf = 6;
            s.mesh = obj.mesh;
            dC = DomainFunction(s);
        end

        function d2C = obtainTensorSecondDerivative(obj,rho)
            fun = obj.homogenization.ddfun;
            s.operation = @(xV) obj.evaluate(rho,fun,xV);
            s.ndimf = 6;
            s.mesh = obj.mesh;
            d2C = DomainFunction(s);
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            obj.fileName = cParams.fileName;
            obj.mesh = cParams.mesh;
            obj.young = cParams.young;
        end

        function loadVademecum(obj)
            fName = [obj.fileName];
            matFile = [fName,'.mat'];
            file2load = fullfile('HMVademecum','Homogenization',matFile);

            v = load(file2load);

            if isfield(v,'homogenization')

                E = obj.young;
                nStre = size(v.homogenization.fun,1);

                for i=1:nStre
                    for j=1:nStre
                        for k=1:nStre
                            for l=1:nStre

                                obj.homogenization.fun{i,j,k,l} = ...
                                    @(x) E.*v.homogenization.fun{i,j,k,l}(x);

                                obj.homogenization.dfun{i,j,k,l} = ...
                                    @(x) E.*v.homogenization.dfun{i,j,k,l}(x);

                                obj.homogenization.ddfun{i,j,k,l} = ...
                                    @(x) E.*v.homogenization.ddfun{i,j,k,l}(x);

                            end
                        end
                    end
                end

            else

                rho = v.rho;
                mat = v.mat;

                [f,df,ddf] = DensityHomogenizationFitter.computePolynomial(8,rho,mat);

                obj.homogenization.fun = f;
                obj.homogenization.dfun = df;
                obj.homogenization.ddfun = ddf;

            end

        end

        function C = evaluate(~,rho,fun,xV)

            nStre = size(fun,1);
            nGaus = size(xV,2);
            nElem = rho.mesh.nelem;

            C = zeros(2,2,2,2,nGaus,nElem);

            rhoV = rho.evaluate(xV);

            for i=1:nStre
                for j=1:nStre
                    for k=1:nStre
                        for l=1:nStre

                            C(i,j,k,l,:,:) = fun{i,j,k,l}(rhoV);

                        end
                    end
                end
            end

        end

    end

end