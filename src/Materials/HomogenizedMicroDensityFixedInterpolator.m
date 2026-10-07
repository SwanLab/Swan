classdef HomogenizedMicroDensityFixedInterpolator < handle

    properties (Access = private)
        designFields
        degradation
        mesh
        young
        fileName
        parameterNames
        nVariables
    end

    methods (Access = public)

        function obj = HomogenizedMicroDensityFixedInterpolator(cParams)
            obj.init(cParams);
            obj.loadVademecum();
        end

        function C = obtainTensor(obj)
            obj.validateDesignVariables();
            fun = obj.degradation.fun;
            s.operation = @(xV) obj.evaluate(fun, xV, []);
            s.ndimf = 6;
            s.mesh = obj.mesh;
            C = DomainFunction(s);
        end

        function dC = obtainTensorDerivative(obj)
            obj.validateDesignVariables();
            fun = obj.degradation.dfun;
            dC = cell(1, obj.nVariables);

            for q = 1:obj.nVariables
                s.ndimf = 6;
                s.mesh = obj.mesh;
                s.operation = obj.createEvaluationOperation(fun, q);
                dC{q} = DomainFunction(s);
            end
        end

        function d2C = obtainTensorSecondDerivative(obj)
            obj.validateDesignVariables();

            if ~obj.hasSecondDerivatives()
                error( ...
                    'HomogenizedMicroDensityFixedInterpolator:SecondDerivativeUnavailable', ...
                    'Second derivatives are not available in the loaded vademecum.');
            end

            fun = obj.degradation.ddfun;
            d2C = cell(obj.nVariables, obj.nVariables);

            for q = 1:obj.nVariables
                for r = 1:obj.nVariables
                    s.ndimf = 6;
                    s.mesh = obj.mesh;
                    s.operation = obj.createEvaluationOperation(fun, [q r]);
                    d2C{q,r} = DomainFunction(s);
                end
            end
        end

        function setDesignVariable(obj, x)

            if ~iscell(x) || numel(x) ~= obj.nVariables
                error( ...
                    'HomogenizedMicroDensityFixedInterpolator:InvalidDesignVariables', ...
                    'Expected %d design variables.', obj.nVariables);
            end

            obj.designFields = x;
        end

        function names = getParameterNames(obj)
            names = obj.parameterNames;
        end

    end

    methods (Access = private)

        function init(obj, cParams)

            obj.mesh = cParams.mesh;
            obj.young = cParams.young;
            obj.fileName = cParams.fileName;
            obj.designFields = {};

            if isfield(cParams, 'parameterNames')
                obj.parameterNames = ...
                    HomogenizedMicroDensityFixedInterpolator.normalizeNames( ...
                    cParams.parameterNames);
            elseif isfield(cParams, 'nVariables')
                obj.parameterNames = arrayfun( ...
                    @(k) sprintf('p%d', k), ...
                    1:cParams.nVariables, ...
                    'UniformOutput', false);
            else
                obj.parameterNames = {};
            end

            obj.nVariables = numel(obj.parameterNames);

        end

        function loadVademecum(obj)

            matFile = [obj.fileName, '.mat'];
            file2load = fullfile('TOVademecum', 'Interpolation', matFile);
            v = load(file2load);

            if ~isfield(v, 'Interpolation')
                error( ...
                    'HomogenizedMicroDensityFixedInterpolator:InvalidVademecum', ...
                    'The file %s does not contain Interpolation.', file2load);
            end

            interpolation = v.Interpolation;

            obj.resolveParameterMetadata(interpolation);

            if ~isfield(interpolation, 'fun') || ...
               ~isfield(interpolation, 'dfun')
                error( ...
                    'HomogenizedMicroDensityFixedInterpolator:InvalidVademecum', ...
                    'Interpolation must contain fun and dfun.');
            end

            E = obj.young;
            nStre = size(interpolation.fun, 1);

            obj.degradation.fun = cell(size(interpolation.fun));
            obj.degradation.dfun = cell(size(interpolation.fun));
            obj.degradation.ddfun = cell(size(interpolation.fun));

            hasDdfun = isfield(interpolation, 'ddfun') && ...
                ~isempty(interpolation.ddfun);

            for i = 1:nStre
                for j = 1:nStre
                    for k = 1:nStre
                        for l = 1:nStre

                            f = interpolation.fun{i,j,k,l};

                            if ~isempty(f)
                                obj.degradation.fun{i,j,k,l} = ...
                                    @(varargin) E .* f(varargin{:});
                            end

                            df = interpolation.dfun{i,j,k,l};

                            if ~isempty(df)
                                obj.degradation.dfun{i,j,k,l} = ...
                                    @(varargin) ...
                                    HomogenizedMicroDensityFixedInterpolator.scaleCell( ...
                                    df(varargin{:}), E);
                            end

                            if hasDdfun

                                ddf = interpolation.ddfun{i,j,k,l};

                                if ~isempty(ddf)
                                    obj.degradation.ddfun{i,j,k,l} = ...
                                        @(varargin) ...
                                        HomogenizedMicroDensityFixedInterpolator.scaleCell( ...
                                        ddf(varargin{:}), E);
                                end

                            end

                        end
                    end
                end
            end

        end

        function resolveParameterMetadata(obj, interpolation)

            namesFromFile = {};

            if isfield(interpolation, 'parameterNames')
                namesFromFile = ...
                    HomogenizedMicroDensityFixedInterpolator.normalizeNames( ...
                    interpolation.parameterNames);
            elseif isfield(interpolation, 'cfg') && ...
                   isfield(interpolation.cfg, 'parameterNames')
                namesFromFile = ...
                    HomogenizedMicroDensityFixedInterpolator.normalizeNames( ...
                    interpolation.cfg.parameterNames);
            end

            if isempty(obj.parameterNames)

                if isempty(namesFromFile)
                    error( ...
                        'HomogenizedMicroDensityFixedInterpolator:MissingParameterMetadata', ...
                        ['The vademecum must contain parameterNames, or ', ...
                         'cParams.parameterNames must be provided.']);
                end

                obj.parameterNames = namesFromFile;

            elseif ~isempty(namesFromFile) && ...
                   ~isequal(obj.parameterNames, namesFromFile)

                error( ...
                    'HomogenizedMicroDensityFixedInterpolator:ParameterOrderMismatch', ...
                    'The requested parameter order does not match the vademecum.');

            end

            obj.nVariables = numel(obj.parameterNames);

            if obj.nVariables < 1
                error( ...
                    'HomogenizedMicroDensityFixedInterpolator:InvalidParameterMetadata', ...
                    'At least one design variable is required.');
            end

        end

        function operation = createEvaluationOperation(obj, fun, derivativeIndex)
            operation = @(xV) obj.evaluate(fun, xV, derivativeIndex);
        end

        function C = evaluate(obj, fun, xV, derivativeIndex)

            nStre = size(fun, 1);
            nGaus = size(xV, 2);
            nElem = obj.mesh.nelem;

            C = zeros(2, 2, 2, 2, nGaus, nElem);

            args = cell(1, obj.nVariables);

            for q = 1:obj.nVariables
                values = obj.designFields{q}.evaluate(xV);
                args{q} = reshape(values, nGaus*nElem, 1);
            end

            for i = 1:nStre
                for j = 1:nStre
                    for k = 1:nStre
                        for l = 1:nStre

                            f = fun{i,j,k,l};

                            if isempty(f)
                                continue
                            end

                            raw = f(args{:});

                            if isempty(derivativeIndex)
                                val = raw;
                            elseif numel(derivativeIndex) == 1
                                val = raw{derivativeIndex};
                            else
                                val = raw{derivativeIndex(1), derivativeIndex(2)};
                            end

                            C(i,j,k,l,:,:) = reshape(val, nGaus, nElem);

                        end
                    end
                end
            end

        end

        function validateDesignVariables(obj)

            if ~iscell(obj.designFields) || ...
               numel(obj.designFields) ~= obj.nVariables
                error( ...
                    'HomogenizedMicroDensityFixedInterpolator:DesignVariablesNotSet', ...
                    'Call setDesignVariable with %d fields before evaluation.', ...
                    obj.nVariables);
            end

        end

        function tf = hasSecondDerivatives(obj)

            tf = false;
            fun = obj.degradation.ddfun;

            for i = 1:numel(fun)
                if ~isempty(fun{i})
                    tf = true;
                    return
                end
            end

        end

    end

    methods (Access = private, Static)


        function names = normalizeNames(inputNames)

            if isstring(inputNames)
                names = cellstr(inputNames);
            elseif ischar(inputNames)
                names = {inputNames};
            elseif iscell(inputNames)
                names = inputNames;
            else
                error( ...
                    'HomogenizedMicroDensityFixedInterpolator:InvalidParameterNames', ...
                    'parameterNames must be char, string or cell array.');
            end

            for i = 1:numel(names)
                names{i} = char(names{i});
            end

        end

        function out = scaleCell(in, E)

            if ~iscell(in)
                out = E .* in;
                return
            end

            out = cell(size(in));

            for i = 1:numel(in)
                if ~isempty(in{i})
                    out{i} = E .* in{i};
                end
            end

        end

    end

end
