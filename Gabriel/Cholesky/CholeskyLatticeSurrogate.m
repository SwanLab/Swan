classdef CholeskyLatticeSurrogate < handle

    properties (Access = private)
        model
        mandelToComponents
        n
        diagonal
        lower
        densityFactor
    end

    methods (Access = public)

        function obj = CholeskyLatticeSurrogate(model,file,s)
            d     = load(file,'columns');
            A     = mandelMap(d.columns(s.data.yFeatures));
            obj.model = model;
            obj.mandelToComponents = A\eye(size(A,1));
            obj.n = sqrt(size(A,1));
            [obj.diagonal,obj.lower] = choleskyLayout(obj.n);
            obj.densityFactor = s.output.densityFactor;
        end

        function Y = computeOutputValues(obj,X)
            L = obj.buildFactor(obj.model.computeOutputValues(X));
            Y = obj.toComponents(obj.factor(X).*pagemtimes(L,'none',L,'transpose'));
        end

        function [Y,dY] = computeOutputValuesAndGradient(obj,X)
            [Z,dZ] = obj.model.computeOutputValuesAndGradient(X);
            L  = obj.buildFactor(Z);
            P  = pagemtimes(L,'none',L,'transpose');
            f  = obj.factor(X);
            df = obj.factorDerivative(X);
            Y  = obj.toComponents(f.*P);
            dY = zeros(size(dZ,1),size(Y,2),size(Y,1));
            for i = 1:size(dZ,1)
                dL = obj.buildFactorDerivative(L,reshape(dZ(i,:,:),size(dZ,2),[])');
                dP = pagemtimes(dL,'none',L,'transpose');
                dM = f.*(dP + pagetranspose(dP)) + df(:,:,:,i).*P;
                dY(i,:,:) = reshape(obj.toComponents(dM)',1,size(Y,2),[]);
            end
        end

        function dY = computeGradient(obj,X)
            [~,dY] = obj.computeOutputValuesAndGradient(X);
        end

        function E = computeError(obj,X,Y)
            E = norm(obj.computeOutputValues(X)-Y,'fro')/max(norm(Y,'fro'),eps);
        end

    end

    methods (Access = private)

        function f = factor(obj,X)
            f = ones(1,1,size(X,1));
            if obj.densityFactor
                f = reshape(X(:,1),1,1,[]);
            end
        end

        function df = factorDerivative(obj,X)
            df = zeros(1,1,size(X,1),size(X,2));
            if obj.densityFactor
                df(:,:,:,1) = 1;
            end
        end

        function L = buildFactor(obj,Z)
            nD = numel(obj.diagonal);
            L  = zeros(obj.n*obj.n,size(Z,1));
            L(obj.diagonal,:) = exp(Z(:,1:nD))';
            L(obj.lower,:)    = Z(:,nD+1:end)';
            L  = reshape(L,obj.n,obj.n,[]);
        end

        function dL = buildFactorDerivative(obj,L,dZ)
            nD = numel(obj.diagonal);
            Lf = reshape(L,obj.n*obj.n,[]);
            dL = zeros(size(Lf));
            dL(obj.diagonal,:) = Lf(obj.diagonal,:).*dZ(:,1:nD)';
            dL(obj.lower,:)    = dZ(:,nD+1:end)';
            dL = reshape(dL,obj.n,obj.n,[]);
        end

        function Y = toComponents(obj,M)
            Y = (obj.mandelToComponents*reshape(M,obj.n*obj.n,[]))';
        end

    end

end