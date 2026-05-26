classdef HomogenizedMicroDensityInterpolator < handle
    % Interpolador único que combina C^h(b) com SIMP (ρ)
    
    properties (Access = private)
        bField          % campo b (cisalhamento)
        rhoField        % campo ρ (densidade)
        fileName
        degradation     % funções C^h(b), dC^h/db
        mesh
        young
    end

    methods (Access = public)

        function obj = HomogenizedMicroDensityInterpolator(cParams)
            obj.init(cParams)
            obj.loadVademecum();
        end

        function C = obtainTensor(obj)
            % 1. Obter C^h(b)
            CH = obj.evaluateHomogenized();
            
            % 2. Obter ρ
            r = Expand(obj.rhoField{1}, 4);
            
            % 3. Combinar com SIMP
            rho3 = r.^3;
            one_minus_rho3 = (1 - r).^3;
            C = rho3 .* CH + 1e-3 * one_minus_rho3 .* CH;
        end

        function dC = obtainTensorDerivative(obj)
            % Derivada em relação a b
            dCH_db = obj.evaluateHomogenizedDerivative();
            r = Expand(obj.rhoField{1}, 4);
            
            dC_db = dCH_db .* (r.^3 + 1e-3 * (1 - r).^3);
            
            % Derivada em relação a ρ
            CH = obj.evaluateHomogenized();
            d_penalty = 3 * r.^2;
            d_void_penalty = -3 * 1e-3 * (1 - r).^2;
            dC_drho = CH .* (d_penalty + d_void_penalty);
            
            % Retornar como cell array (b, ρ)
            dC{1} = dC_db;
            dC{2} = dC_drho;
        end

        function setDesignVariable(obj, x)
            % x{1} = b, x{2} = ρ
            obj.bField = x{1};
            obj.rhoField = x{2};
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.fileName = cParams.fileName;
            obj.mesh = cParams.mesh;
            obj.young = cParams.young;
        end

        function loadVademecum(obj)
            fName = [obj.fileName];
            matFile = [fName, '.mat'];
            file2load = fullfile('TOVademecum', 'Interpolation', matFile);
            v = load(file2load);
            
            E = obj.young;
            nStre = size(v.Interpolation.fun, 1);
            
            for i = 1:nStre
                for j = 1:nStre
                    for k = 1:nStre
                        for l = 1:nStre
                            obj.degradation.fun{i,j,k,l} = @(x) E .* v.Interpolation.fun{i,j,k,l}(x);
                            obj.degradation.dfun{i,j,k,l} = @(x) E .* v.Interpolation.dfun{i,j,k,l}(x);
                        end
                    end
                end
            end
        end

        function CH = evaluateHomogenized(obj)
            fun = obj.degradation.fun;
            s.operation = @(xV) obj.evaluateCH(fun, xV);
            s.ndimf = 6;
            s.mesh = obj.mesh;
            CH = DomainFunction(s);
        end

        function dCH = evaluateHomogenizedDerivative(obj)
            fun = obj.degradation.dfun;
            s.operation = @(xV) obj.evaluateCH(fun, xV);
            s.ndimf = 6;
            s.mesh = obj.mesh;
            dCH = DomainFunction(s);
        end

        function C = evaluateCH(obj, fun, xV)
            nStre = size(fun, 1);
            nGaus = size(xV, 2);
            nElem = obj.mesh.nelem;
            C = zeros(2, 2, 2, 2, nGaus, nElem);
            phiV = obj.bField{1}.evaluate(xV);  % ← b, não ρ!
            for i = 1:nStre
                for j = 1:nStre
                    for k = 1:nStre
                        for l = 1:nStre
                            C(i, j, k, l, :, :) = fun{i, j, k, l}(phiV);
                        end
                    end
                end
            end
        end

    end

end