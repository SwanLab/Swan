classdef TrainingDataset < handle

    properties (Access = private)
        cases
        parameterNames
        parameterValues
        directory
        saveMAT
        saveCSV
    end

    methods (Access = public)

        function obj = TrainingDataset(s)
            obj.init(s);
        end

        function samples = compute(obj)
            nSamples = numel(obj.cases);
            samples = cell(nSamples, 1);

            for i = 1:nSamples
                fprintf('Computing sample %d/%d\n', i, nSamples);

                training = TrainingPreconditioning(obj.cases{i});
                data = training.compute();

                data.parameterNames  = obj.parameterNames;
                data.parameterValues = obj.parameterValues(i,:);

                samples{i} = data;
            end
        end

        function export(obj, samples)
            if ~(obj.saveMAT || obj.saveCSV)
                return
            end

            if ~isfolder(obj.directory)
                mkdir(obj.directory);
            end

            if obj.saveMAT
                obj.exportMAT(samples);
            end

            if obj.saveCSV
                obj.exportCSV(samples);
            end
        end

    end

    methods (Access = private)

        function init(obj, s)
            obj.cases           = s.cases;
            obj.parameterNames  = s.parameterNames;
            obj.parameterValues = s.parameterValues;
            obj.directory       = s.directory;
            obj.saveMAT         = s.saveMAT;
            obj.saveCSV         = s.saveCSV;

            assert(~isempty(obj.cases), ...
                'At least one sample is required.');

            assert(size(obj.parameterValues, 1) == numel(obj.cases), ...
                'There must be one parameter row per sample.');

            assert(size(obj.parameterValues, 2) == ...
                numel(obj.parameterNames), ...
                'There must be one name per parameter column.');
        end

        function exportMAT(obj, samples)
            for i = 1:numel(samples)
                data = samples{i};
                fileName = 'sample';
                for j = 1:numel(data.parameterNames)
                    name  = data.parameterNames{j};
                    value = num2str(data.parameterValues(j), '%.8g');
                    value = strrep(value, '.', '_');
                    fileName = [fileName, '_', name, '_', value];
                end
                filePath = fullfile(obj.directory, [fileName, '.mat']);
                save(filePath, '-struct', 'data');
            end
        end

        function exportCSV(obj, samples)
            nSamples = numel(samples);
            first = samples{1};
            mesh = first.mesh;

            ndim   = mesh.ndim;
            nModes = size(first.T, 2);
            nK     = nModes*(nModes + 1)/2; % Upper triangular elem including diag

            tBlocks = cell(nSamples, 1);
            kRows   = zeros(nSamples, ...
                size(obj.parameterValues, 2) + nK);

            for i = 1:nSamples
                data = samples{i};
                mesh = data.mesh;
                T    = data.T;
                K    = data.Kcoarse;
                p    = data.parameterValues;
                nNodes = size(mesh.coord, 1);

                assert(mesh.ndim == ndim && size(T, 2) == nModes, ...
                    'CSV samples must have the same dimension and number of modes.');

                assert(size(T, 1) == nNodes*ndim, ...
                    'T must contain ndim displacement components per node.');

                assert(isequal(size(K), [nModes, nModes]), ...
                    'Kcoarse dimensions must match the number of modes.');

                % T: one row per node, components grouped by mode
                tValues = zeros(nNodes, ndim*nModes);

                for k = 1:nModes
                    columns = (k-1)*ndim + (1:ndim);
                    tValues(:, columns) = reshape(T(:,k), ndim, []).';
                end

                tBlocks{i} = [repmat(p, nNodes, 1), mesh.coord, tValues];

                % K: upper triangle, traversed row by row
                kValues = zeros(1, nK);
                position = 1;

                for row = 1:nModes
                    count = nModes - row + 1;
                    kValues(position:position + count - 1) = ...
                        K(row, row:nModes);
                    position = position + count;
                end

                kRows(i,:) = [p, kValues];
            end

            tData = vertcat(tBlocks{:});

            writematrix(tData, ...
                fullfile(obj.directory, 'DataT.csv'));

            writematrix(kRows, ...
                fullfile(obj.directory, 'DataK.csv'));
        end

    end

end