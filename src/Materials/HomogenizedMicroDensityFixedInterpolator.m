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
            CH = obj.evaluateHomogenized();
            r  = Expand(obj.rhoField{1}, 4);
            C  = (r.^3 + 1e-3*(1-r).^3) .* CH;
        end

        function dC = obtainTensorDerivative(obj)
            CH     = obj.evaluateHomogenized();
            dCH_db = obj.evaluateHomogenizedDerivative();
            r      = Expand(obj.rhoField{1}, 4);

            
            penalty_b = r.^3 + 1e-3*(1-r).^3;
            dC{1} = dCH_db .* penalty_b;

           
            dC{2} = CH .* (3*r.^2 - 3*1e-3*(1-r).^2);
        end

        function setDesignVariable(obj, x)
            % x = {{b_filtrado}, {rho_filtrado}}
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
            E     = obj.young;
            nStre = size(v.Interpolation.fun, 1);
            for i = 1:nStre
                for j = 1:nStre
                    for k = 1:nStre
                        for l = 1:nStre
                            obj.degradation.fun{i,j,k,l} = ...
                                @(x) E .* v.Interpolation.fun{i,j,k,l}(x);
                            obj.degradation.dfun{i,j,k,l} = ...
                                @(x) E .* v.Interpolation.dfun{i,j,k,l}(x);
                        end
                    end
                end
            end
        end

        function CH = evaluateHomogenized(obj)
            fun = obj.degradation.fun;
            s.operation = @(xV) obj.evaluateCH(fun, xV);
            s.ndimf = 6;
            s.mesh  = obj.mesh;
            CH = DomainFunction(s);
        end

        function dCH = evaluateHomogenizedDerivative(obj)
            fun = obj.degradation.dfun;
            s.operation = @(xV) obj.evaluateCH(fun, xV);
            s.ndimf = 6;
            s.mesh  = obj.mesh;
            dCH = DomainFunction(s);
        end

        function C = evaluateCH(obj, fun, xV)
            nStre = size(fun, 1);
            nGaus = size(xV, 2);
            nElem = obj.mesh.nelem;
            C     = zeros(2, 2, 2, 2, nGaus, nElem);
            phiV  = obj.bField{1}.evaluate(xV);
            for i = 1:nStre
                for j = 1:nStre
                    for k = 1:nStre
                        for l = 1:nStre
                            C(i,j,k,l,:,:) = fun{i,j,k,l}(phiV);
                        end
                    end
                end
            end
        end

    end
end