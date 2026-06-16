classdef HomogenizedMicroDensityFixedInterpolator < handle

    properties (Access = private)
        bField
        rhoField
        degradation
        mesh
        young
        fileName
    end

    methods (Access = public)

        function obj = HomogenizedMicroDensityFixedInterpolator(cParams)
            obj.init(cParams);
            obj.loadVademecum();
        end

        function C = obtainTensor(obj)
            fun = obj.degradation.fun;
            s.operation = @(xV) obj.evaluate(fun, xV, 'value');
            s.ndimf = 6;
            s.mesh  = obj.mesh;
            C = DomainFunction(s);
        end

        function dC = obtainTensorDerivative(obj)
            fun = obj.degradation.dfun;
            s.ndimf = 6;
            s.mesh  = obj.mesh;
            s.operation  = @(xV) obj.evaluate(fun, xV, 'db');
            dC{1} = DomainFunction(s);
            s.operation  = @(xV) obj.evaluate(fun, xV, 'drho');
            dC{2} = DomainFunction(s);
        end

        function d2C = obtainTensorSecondDerivative(obj)
            fun = obj.degradation.ddfun;
            s.ndimf = 6;
            s.mesh  = obj.mesh;
            s.operation  = @(xV) obj.evaluate(fun, xV, 'db2');
            d2C{1,1} = DomainFunction(s);
            s.operation  = @(xV) obj.evaluate(fun, xV, 'drho2');
            d2C{2,2} = DomainFunction(s);
            s.operation  = @(xV) obj.evaluate(fun, xV, 'dbrho');
            d2C{1,2} = DomainFunction(s);
            d2C{2,1} = d2C{1,2};
        end

        function setDesignVariable(obj, x)
            obj.bField   = x{1};
            obj.rhoField = x{2};
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.mesh     = cParams.mesh;
            obj.young    = cParams.young;
            obj.fileName = cParams.fileName;
        end

        function loadVademecum(obj)
            matFile   = [obj.fileName, '.mat'];
            file2load = fullfile('TOVademecum', 'Interpolation', matFile);
            v = load(file2load);
            E = obj.young;
            nStre = size(v.Interpolation.fun, 1);
            
            for i = 1:nStre
                for j = 1:nStre
                    for k = 1:nStre
                        for l = 1:nStre
                            
                            obj.degradation.fun{i,j,k,l} = @(b,rho) E .* v.Interpolation.fun{i,j,k,l}(b,rho);

                            obj.degradation.dfun{i,j,k,l} = @(b,rho) obj.scaleGrad(v.Interpolation.dfun{i,j,k,l}(b,rho), E);

                            obj.degradation.ddfun{i,j,k,l} =@(b,rho) obj.scaleHess(v.Interpolation.ddfun{i,j,k,l}(b,rho), E);
                        end
                    end
                end
            end
        end

        function C = evaluate(obj, fun, xV, component)
            nStre = size(fun, 1);
            nGaus = size(xV, 2);
            nElem = obj.mesh.nelem;
            C = zeros(2, 2, 2, 2, nGaus, nElem);

            bV   = obj.bField.evaluate(xV);
            rhoV = obj.rhoField.evaluate(xV);

            bFlat   = reshape(bV,   nGaus*nElem, 1);
            rhoFlat = reshape(rhoV, nGaus*nElem, 1);

            for i = 1:nStre
                for j = 1:nStre
                    for k = 1:nStre
                        for l = 1:nStre
                            if isempty(fun{i,j,k,l})
                                continue
                            end
                            raw = fun{i,j,k,l}(bFlat, rhoFlat);
                            switch component
                                case 'value'
                                    val = raw;
                                case 'db'
                                    val = raw{1};
                                case 'drho'
                                    val = raw{2};
                                case 'db2'
                                    val = raw{1,1};
                                case 'drho2'
                                    val = raw{2,2};
                                case 'dbrho'
                                    val = raw{1,2};
                            end
                            C(i,j,k,l,:,:) = reshape(val, nGaus, nElem);
                        end
                    end
                end
            end
        end

    end

    methods (Access = private, Static)

        function dJ = scaleGrad(dJ_in, E)
            
            if iscell(dJ_in)
                dJ = cell(1, 2);
                dJ{1} = E .* dJ_in{1};
                dJ{2} = E .* dJ_in{2};
            else
                dJ = dJ_in;
            end
        end

        function ddJ = scaleHess(ddJ_in, E)
            
            if iscell(ddJ_in)
                ddJ = cell(2, 2);
                ddJ{1,1} = E .* ddJ_in{1,1};
                ddJ{2,2} = E .* ddJ_in{2,2};
                ddJ{1,2} = E .* ddJ_in{1,2};
                ddJ{2,1} = E .* ddJ_in{2,1};
            else
                ddJ = ddJ_in;
            end
        end

    end

end