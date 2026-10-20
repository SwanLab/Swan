classdef LatticeHomogenizedMaterialsReader < handle

    properties (Access = private)
        fileName
        problemNN
        mesh
        young
        parameterNames
        nFeatures
    end

    methods (Access = public)

        function obj = LatticeHomogenizedMaterialsReader(cParams)
            obj.init(cParams);
            obj.loadVademecum();
        end

        function C = obtainTensor(obj,designVariables)
            designVariables = obj.formatDesignVariables(designVariables);
            obj.validateDesignVariables(designVariables);

            s.operation = @(xV) obj.evaluateTensor(designVariables,xV);
            s.ndimf     = 6;
            s.mesh      = obj.mesh;

            C = DomainFunction(s);
        end

        function dC = obtainTensorDerivative(obj,designVariables,iFeature)
            designVariables = obj.formatDesignVariables(designVariables);
            obj.validateDesignVariables(designVariables);
            obj.checkFeatureIndex(iFeature);

            s.operation = @(xV) obj.evaluateTensorDerivative( ...
                designVariables,iFeature,xV);
            s.ndimf = 6;
            s.mesh  = obj.mesh;

            dC = DomainFunction(s);
        end

        function dC = obtainTensorDerivatives(obj,designVariables)
            designVariables = obj.formatDesignVariables(designVariables);
            obj.validateDesignVariables(designVariables);

            dC = cell(obj.nFeatures,1);

            for iFeature = 1:obj.nFeatures
                s.operation = @(xV) obj.evaluateTensorDerivative( ...
                    designVariables,iFeature,xV);
                s.ndimf = 6;
                s.mesh  = obj.mesh;

                dC{iFeature} = DomainFunction(s);
            end
        end

        function names = getParameterNames(obj)
            names = obj.parameterNames;
        end

    end

    methods (Access = private)

        function init(obj,cParams)
            requiredFields = {'fileName','mesh','young','parameterNames'};

            for iField = 1:numel(requiredFields)
                fieldName = requiredFields{iField};

                if ~isfield(cParams,fieldName)
                    error('LatticeHomogenizedMaterialsReader:MissingParameter', ...
                        'Missing parameter: %s.',fieldName);
                end
            end

            obj.fileName      = cParams.fileName;
            obj.mesh          = cParams.mesh;
            obj.young         = cParams.young;
            obj.parameterNames = obj.normalizeNames(cParams.parameterNames);
            obj.nFeatures     = numel(obj.parameterNames);

            if ~isscalar(obj.young) || ...
                    ~isfinite(obj.young) || ...
                    obj.young <= 0
                error('LatticeHomogenizedMaterialsReader:InvalidYoung', ...
                    'young must be a positive scalar.');
            end

            if obj.nFeatures < 1
                error('LatticeHomogenizedMaterialsReader:InvalidParameters', ...
                    'At least one parameter must be provided.');
            end
        end

        function loadVademecum(obj)
            [folder,name,extension] = fileparts(obj.fileName);

            if isempty(extension)
                file = [name,'.mat'];
            else
                file = [name,extension];
            end

            if isempty(folder)
                file2load = fullfile( ...
                    'Gabriel','HMVademecum', ...
                    'Homogenization', ...
                    file);
            else
                file2load = obj.fileName;

                if isempty(extension)
                    file2load = [file2load,'.mat'];
                end
            end

            if ~isfile(file2load)
                error('LatticeHomogenizedMaterialsReader:FileNotFound', ...
                    'Lattice NN vademecum file not found: %s.',file2load);
            end

            v = load(file2load);

            if ~isfield(v,'problemNN')
                error('LatticeHomogenizedMaterialsReader:InvalidVademecum', ...
                    'The lattice NN vademecum must contain problemNN.');
            end

            obj.problemNN = v.problemNN;

            if isfield(v,'parameterNames')
                namesFromFile = obj.normalizeNames(v.parameterNames);

                if ~isequal(namesFromFile,obj.parameterNames)
                    error('LatticeHomogenizedMaterialsReader:ParameterOrderMismatch', ...
                        'The requested parameter order does not match the vademecum.');
                end
            end
        end

        function designVariables = formatDesignVariables(~,designVariables)
            if isempty(designVariables)
                error('LatticeHomogenizedMaterialsReader:EmptyDesignVariables', ...
                    'At least one design variable must be supplied.');
            end

            if ~iscell(designVariables)
                try
                    designVariables = num2cell(designVariables);
                catch
                    error('LatticeHomogenizedMaterialsReader:InvalidDesignVariables', ...
                        'designVariables must be a cell array or an object array.');
                end
            end

            designVariables = designVariables(:);
        end

        function validateDesignVariables(obj,designVariables)
            if numel(designVariables) ~= obj.nFeatures
                error('LatticeHomogenizedMaterialsReader:InvalidNumberOfDesignVariables', ...
                    'Expected %d design variables and received %d.', ...
                    obj.nFeatures,numel(designVariables));
            end

            for iFeature = 1:obj.nFeatures
                if ~ismethod(designVariables{iFeature},'evaluate')
                    error('LatticeHomogenizedMaterialsReader:InvalidDesignVariable', ...
                        'Design variable %d cannot be evaluated.',iFeature);
                end
            end
        end

        function checkFeatureIndex(obj,iFeature)
            if ~isscalar(iFeature) || ...
                    ~isfinite(iFeature) || ...
                    iFeature < 1 || ...
                    iFeature > obj.nFeatures || ...
                    iFeature ~= floor(iFeature)

                error('LatticeHomogenizedMaterialsReader:InvalidFeatureIndex', ...
                    'iFeature must be an integer between 1 and %d.', ...
                    obj.nFeatures);
            end
        end

        function X = evaluateDesignVariables(obj,designVariables,xV)
            nGaus   = size(xV,2);
            nElem   = obj.mesh.nelem;
            nPoints = nGaus*nElem;

            X = zeros(nPoints,obj.nFeatures);

            for iFeature = 1:obj.nFeatures
                value = designVariables{iFeature}.evaluate(xV);

                if numel(value) ~= nPoints
                    error('LatticeHomogenizedMaterialsReader:UnexpectedFieldSize', ...
                        'Unexpected number of values for design variable %d.', ...
                        iFeature);
                end

                X(:,iFeature) = reshape(value,[],1);
            end

            if any(~isfinite(X(:)))
                error('LatticeHomogenizedMaterialsReader:NonFiniteInput', ...
                    'Non-finite design variable values detected.');
            end
        end

        function C = evaluateTensor(obj,designVariables,xV)
            nGaus = size(xV,2);
            nElem = obj.mesh.nelem;

            X = obj.evaluateDesignVariables(designVariables,xV);

            Y = obj.problemNN.computeOutputValues(X);
            Y = obj.young.*Y;

            obj.checkOutput(Y,size(X,1));

            C = obj.buildTensor(Y,nGaus,nElem);
        end

        function dC = evaluateTensorDerivative( ...
                obj,designVariables,iFeature,xV)

            nGaus = size(xV,2);
            nElem = obj.mesh.nelem;

            X = obj.evaluateDesignVariables(designVariables,xV);

            dY = obj.problemNN.computeGradient(X);

            nPoints = size(X,1);

            if size(dY,1) ~= obj.nFeatures
                error('LatticeHomogenizedMaterialsReader:UnexpectedJacobianFeatures', ...
                    'Unexpected number of NN Jacobian input features.');
            end

            if size(dY,2) ~= 6 || size(dY,3) ~= nPoints
                error('LatticeHomogenizedMaterialsReader:UnexpectedJacobianSize', ...
                    'Unexpected NN Jacobian dimensions.');
            end

            dYFeature = reshape( ...
                dY(iFeature,:,:), ...
                6,nPoints)';

            dYFeature = obj.young.*dYFeature;

            obj.checkOutput(dYFeature,nPoints);

            dC = obj.buildTensor(dYFeature,nGaus,nElem);
        end

        function checkOutput(~,Y,nPoints)
            if size(Y,1) ~= nPoints || size(Y,2) ~= 6
                error('LatticeHomogenizedMaterialsReader:UnexpectedOutputSize', ...
                    'Unexpected lattice NN output dimensions.');
            end

            if any(~isfinite(Y(:)))
                error('LatticeHomogenizedMaterialsReader:NonFiniteOutput', ...
                    'Non-finite lattice NN output detected.');
            end
        end

        function C = buildTensor(~,Y,nGaus,nElem)
            nPoints = nGaus*nElem;

            if size(Y,1) ~= nPoints || size(Y,2) ~= 6
                error('LatticeHomogenizedMaterialsReader:InvalidTensorData', ...
                    'Tensor data must have size nPoints x 6.');
            end

            C = zeros(2,2,2,2,nGaus,nElem);

            C1111 = reshape(Y(:,1),1,1,1,1,nGaus,nElem);
            C2222 = reshape(Y(:,2),1,1,1,1,nGaus,nElem);
            C1122 = reshape(Y(:,3),1,1,1,1,nGaus,nElem);
            C1112 = reshape(Y(:,4),1,1,1,1,nGaus,nElem);
            C2212 = reshape(Y(:,5),1,1,1,1,nGaus,nElem);
            C1212 = reshape(Y(:,6),1,1,1,1,nGaus,nElem);

            C(1,1,1,1,:,:) = C1111;
            C(2,2,2,2,:,:) = C2222;

            C(1,1,2,2,:,:) = C1122;
            C(2,2,1,1,:,:) = C1122;

            C(1,1,1,2,:,:) = C1112;
            C(1,1,2,1,:,:) = C1112;
            C(1,2,1,1,:,:) = C1112;
            C(2,1,1,1,:,:) = C1112;

            C(2,2,1,2,:,:) = C2212;
            C(2,2,2,1,:,:) = C2212;
            C(1,2,2,2,:,:) = C2212;
            C(2,1,2,2,:,:) = C2212;

            C(1,2,1,2,:,:) = C1212;
            C(1,2,2,1,:,:) = C1212;
            C(2,1,1,2,:,:) = C1212;
            C(2,1,2,1,:,:) = C1212;
        end

    end

    methods (Static, Access = private)

        function names = normalizeNames(inputNames)
            if isstring(inputNames)
                names = cellstr(inputNames);
            elseif ischar(inputNames)
                names = {inputNames};
            elseif iscell(inputNames)
                names = inputNames;
            else
                error('LatticeHomogenizedMaterialsReader:InvalidParameterNames', ...
                    'parameterNames must be char, string or cell array.');
            end

            names = names(:)';

            for i = 1:numel(names)
                names{i} = char(names{i});
            end
        end

    end

end