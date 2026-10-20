classdef MultiField < DesignVariable

    properties (GetAccess = public, SetAccess = private)
        fields
    end

    properties (Access = private)
        plotting
        plotter
    end

    methods (Access = public)

        function obj = MultiField(cParams)
            obj.validateInput(cParams);
            obj.fields     = cParams.fields;
            obj.nVariables = numel(obj.fields);
            cParams.fun    = obj.createConcatenatedFunction(obj.fields);
            obj.init(cParams);
            obj.createPlotter(cParams);
        end

        function fun = obtainDomainFunction(obj)
            fun = obj.fields;
        end

        function update(obj,value)
            value  = value(:);
            nTotal = obj.computeNumberDofs();

            if numel(value) ~= nTotal
                error('MultiField:InvalidNumberOfValues', ...
                    'Expected %d values and received %d.', ...
                    nTotal,numel(value));
            end

            offset = 0;

            for i = 1:obj.nVariables
                field   = obj.fields{i};
                nValues = numel(field.fValues);
                indices = offset + (1:nValues);
                values  = reshape(value(indices),size(field.fValues));

                field.setFValues(values);

                offset = offset + nValues;
            end

            obj.fun.setFValues(value);
        end

        function fixedDofs = getFixedDofs(obj)
            fixedDofs = obj.isFixed;
        end

        function n = obtainNumberVariables(obj)
            n = obj.nVariables;
        end

        function plot(obj)
            if obj.plotting
                obj.plotter.plot();
            end
        end

    end

    methods (Access = private)

        function createPlotter(obj,cParams)
            if isfield(cParams,'plotting')
                obj.plotting = cParams.plotting;
            else
                obj.plotting = false;
            end

            if obj.plotting
                obj.plotter = Plotter.create(obj);
            end
        end

        function n = computeNumberDofs(obj)
            n = 0;

            for i = 1:obj.nVariables
                n = n + numel(obj.fields{i}.fValues);
            end
        end

    end

    methods (Access = private, Static)

        function validateInput(cParams)
            if ~isfield(cParams,'fields')
                error('MultiField:MissingFields', ...
                    'cParams.fields must be provided.');
            end

            if ~iscell(cParams.fields)
                error('MultiField:InvalidFields', ...
                    'cParams.fields must be a cell array.');
            end

            if isempty(cParams.fields)
                error('MultiField:EmptyFields', ...
                    'At least one field must be provided.');
            end

            for i = 1:numel(cParams.fields)
                if ~isprop(cParams.fields{i},'fValues')
                    error('MultiField:InvalidField', ...
                        'Each field must contain fValues.');
                end
            end
        end

        function fun = createConcatenatedFunction(fields)
            values = [];

            for i = 1:numel(fields)
                values = [values; fields{i}.fValues(:)];
            end

            fun = fields{1}.copy();
            fun.setFValues(values);
        end

    end

end