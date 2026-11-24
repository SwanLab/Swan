classdef SymmetricLatticeSurrogate < handle

    properties (Access = private)
        model
        elements
        columns
    end

    methods (Access = public)

        function obj = SymmetricLatticeSurrogate(model,file,s)
            d              = load(file,'columns','settings');
            components     = SymmetricLatticeSurrogate.parseComponents(d.columns(s.data.yFeatures));
            % transformation = CellTransformation(d.settings.transformation);
            if isfield(s,'transformation')
                transformation = s.transformation;
            else
                transformation = CellTransformation(d.settings.transformation);
            end
            generators     = cellfun(@(name) transformation.getSymmetry(name),d.settings.dataset.symmetries,'UniformOutput',false);
            obj.model      = model;
            obj.elements   = SymmetricLatticeSurrogate.buildGroup(generators,components);
            obj.columns    = [];
            if isfield(s,'symmetryColumns')
                obj.columns = s.symmetryColumns;
            end
        end

        function Y = computeOutputValues(obj,X)
            Y = 0;
            for e = obj.elements
                Y = Y + obj.model.computeOutputValues(obj.mapInput(e,X))*e.inverseOutputMap';
            end
            Y = Y/numel(obj.elements);
        end

        function [Y,dY] = computeOutputValuesAndGradient(obj,X)
            Y  = 0;
            dY = 0;
            for e = obj.elements
                [Ye,dYe] = obj.model.computeOutputValuesAndGradient(obj.mapInput(e,X));
                dX = pagemtimes(obj.inputJacobian(e,X),'transpose',dYe,'none');
                Y  = Y  + Ye*e.inverseOutputMap';
                dY = dY + pagemtimes(dX,e.inverseOutputMap');
            end
            Y  = Y/numel(obj.elements);
            dY = dY/numel(obj.elements);
        end

        function dY = computeGradient(obj,X)
            [~,dY] = obj.computeOutputValuesAndGradient(X);
        end

        function E = computeError(obj,X,Y)
            E = norm(obj.computeOutputValues(X)-Y,'fro')/max(norm(Y,'fro'),eps);
        end

    end

    methods (Access = private)

        function c = parameterColumns(obj,X)
            % c = 2:size(X,2);
            if isempty(obj.columns)
                c = 2:size(X,2);
            else
                c = obj.columns;
            end
        end

        function Xe = mapInput(obj,e,X)
            c       = obj.parameterColumns(X);
            Xe      = X;
            Xe(:,c) = e.mapParams(X(:,c));
        end

        function J = inputJacobian(obj,e,X)
            c       = obj.parameterColumns(X);
            nInputs = size(X,2);
            J = repmat(eye(nInputs),1,1,size(X,1));
            J(c,c,:) = e.mapJacobian(X(:,c));
        end

    end

    methods (Access = private, Static)

        function elements = buildGroup(generators,components)
            nOutputs = size(components,1);
            elements = struct('mapParams',@(p) p, ...
                              'mapJacobian',@(p) repmat(eye(size(p,2)),1,1,size(p,1)), ...
                              'outputMap',eye(nOutputs), ...
                              'inverseOutputMap',eye(nOutputs));
            for iG = 1:numel(generators)
                h = generators{iG};
                M = SymmetricLatticeSurrogate.vectorMap(h,components);
                composed = arrayfun(@(g) SymmetricLatticeSurrogate.compose(h,M,g),elements,'UniformOutput',false);
                elements = [elements,composed{:}]; %#ok<AGROW>
            end
        end

        function e = compose(h,M,g)
            e.mapParams        = @(p) h.mapParams(g.mapParams(p));
            e.mapJacobian      = @(p) pagemtimes(h.mapParamsJacobian(g.mapParams(p)),g.mapJacobian(p));
            e.outputMap        = M*g.outputMap;
            e.inverseOutputMap = e.outputMap\eye(size(M,1));
        end

        function M = vectorMap(h,components)
            nOutputs = size(components,1);
            M = zeros(nOutputs);
            for k = 1:nOutputs
                v      = zeros(nOutputs,1);
                v(k)   = 1;
                C      = h.mapTensor(SymmetricLatticeSurrogate.vectorToTensor(v,components));
                M(:,k) = SymmetricLatticeSurrogate.tensorToVector(C,components);
            end
        end

        function C = vectorToTensor(v,components)
            n = max(components(:));
            C = zeros(n,n,n,n);
            for k = 1:size(components,1)
                c = components(k,:);
                for idx = [c; c([2 1 3 4]); c([1 2 4 3]); c([2 1 4 3]); c([3 4 1 2]); c([4 3 1 2]); c([3 4 2 1]); c([4 3 2 1])]'
                    C(idx(1),idx(2),idx(3),idx(4)) = v(k);
                end
            end
        end

        function v = tensorToVector(C,components)
            v = arrayfun(@(k) C(components(k,1),components(k,2),components(k,3),components(k,4)),(1:size(components,1))');
        end

        function components = parseComponents(names)
            names      = cellstr(names);
            components = cell2mat(cellfun(@(name) name(2:end) - '0',names(:),'UniformOutput',false));
        end

    end

end