classdef HomogenizedMicrostructureBOnlyInterpolator < handle
    

    properties (Access = private)
        bField
        fileName
        degradation
        mesh
        young
    end

    methods (Access = public)

        function obj = HomogenizedMicrostructureBOnlyInterpolator(cParams)

            obj.init(cParams);
            obj.loadVademecum();

        end

        function C = obtainTensor(obj)

            fun = obj.degradation.fun;

            s.operation = @(xV) obj.evaluate(fun,xV);
            s.ndimf     = 6;
            s.mesh      = obj.mesh;

            C = DomainFunction(s);

        end

        function dC = obtainTensorDerivative(obj)

            fun = obj.degradation.dfun;

            s.operation = @(xV) obj.evaluate(fun,xV);
            s.ndimf     = 6;
            s.mesh      = obj.mesh;

           
            dC{1} = DomainFunction(s);

        end

        function d2C = obtainTensorSecondDerivative(obj)

            fun = obj.degradation.ddfun;

            s.operation = @(xV) obj.evaluate(fun,xV);
            s.ndimf     = 6;
            s.mesh      = obj.mesh;

            
            d2C{1,1} = DomainFunction(s);

        end

        function setDesignVariable(obj,x)

            

            if iscell(x)

                if isempty(x)
                    error('The design-variable cell is empty.');
                end

                obj.bField = x{1};

            else

                obj.bField = x;

            end

        end

    end

    methods (Access = private)

        function init(obj,cParams)

            obj.fileName = cParams.fileName;
            obj.mesh     = cParams.mesh;
            obj.young    = cParams.young;
            obj.bField   = [];

            
            if isfield(cParams,'density')

                x = cParams.density;

                if iscell(x)

                    if ~isempty(x)
                        obj.bField = x{1};
                    end

                elseif isa(x,'LagrangianFunction')

                    obj.bField = x;

                end

            end

        end

        function loadVademecum(obj)

            matFile = [obj.fileName,'.mat'];

            file2load = fullfile( ...
                'TOVademecum', ...
                'Interpolation', ...
                matFile);

            if ~exist(file2load,'file')
                error('Homogenization file not found:\n%s',file2load);
            end

            v = load(file2load);

            if ~isfield(v,'Interpolation')
                error([ ...
                    'The file %s does not contain the variable ', ...
                    '"Interpolation".'],file2load);
            end

            requiredFields = {'fun','dfun','ddfun'};

            for iField = 1:numel(requiredFields)

                fieldName = requiredFields{iField};

                if ~isfield(v.Interpolation,fieldName)
                    error( ...
                        'Interpolation.%s is missing in %s.', ...
                        fieldName,file2load);
                end

            end

            E = obj.young;

            interpolation = v.Interpolation;

            nStre = size(interpolation.fun,1);

            obj.degradation.fun   = cell(size(interpolation.fun));
            obj.degradation.dfun  = cell(size(interpolation.dfun));
            obj.degradation.ddfun = cell(size(interpolation.ddfun));

            for i = 1:nStre
                for j = 1:nStre
                    for k = 1:nStre
                        for l = 1:nStre

                            if ~isempty(interpolation.fun{i,j,k,l})

                                fLocal = interpolation.fun{i,j,k,l};

                                obj.degradation.fun{i,j,k,l} = ...
                                    @(b) E.*fLocal(b);

                            end

                            if ~isempty(interpolation.dfun{i,j,k,l})

                                dfLocal = interpolation.dfun{i,j,k,l};

                                obj.degradation.dfun{i,j,k,l} = ...
                                    @(b) E.*dfLocal(b);

                            end

                            if ~isempty(interpolation.ddfun{i,j,k,l})

                                ddfLocal = interpolation.ddfun{i,j,k,l};

                                obj.degradation.ddfun{i,j,k,l} = ...
                                    @(b) E.*ddfLocal(b);

                            end

                        end
                    end
                end
            end

            fprintf('\n');
            fprintf('B-only homogenization loaded successfully.\n');
            fprintf('File: %s\n',file2load);
            fprintf('Young scaling: %.6e\n',E);

        end

        function C = evaluate(obj,fun,xV)

            if isempty(obj.bField)
                error([ ...
                    'The b design field has not been assigned. ', ...
                    'Call setDesignVariable before evaluating the tensor.']);
            end

            nStre = size(fun,1);
            nGaus = size(xV,2);
            nElem = obj.mesh.nelem;

            C = zeros( ...
                2,2,2,2, ...
                nGaus,nElem);

            bV = obj.bField.evaluate(xV);

            bFlat = reshape( ...
                real(bV), ...
                nGaus*nElem,1);

            for i = 1:nStre
                for j = 1:nStre
                    for k = 1:nStre
                        for l = 1:nStre

                            if isempty(fun{i,j,k,l})
                                continue
                            end

                            value = fun{i,j,k,l}(bFlat);

                            value = reshape( ...
                                real(value), ...
                                nGaus,nElem);

                            C(i,j,k,l,:,:) = value;

                        end
                    end
                end
            end

        end

    end

end