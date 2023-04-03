classdef DesignFieldTester < handle
    properties (Access = private)
        iterations
        tolerateError
        designFields
    end

    methods (Access = public)
        function obj = DesignFieldTester()
            obj.iterations = 3;
            obj.tolerateError = 1e-10;
            obj.createFields();
        end
        function testFilter(obj)
            obj.designFields.filter();
            results = obj.designFields.filteredField;
            %Validate Results
            s.results = results;
            B = FieldFilterTester(3);
            B.loadResults(s);
            B.validate;
        end
        function testProjector(obj)
            obj.testFilter();
            obj.designFields.project();
            results = obj.designFields.projectedField;
            %Validate Results
            s.results = results;
            B = FieldProjectorTester(3);
            B.loadResults(s);
            B.validate;
        end
        function testProjectedFieldDerivator(obj)
            obj.testProjector();
            obj.designFields.deriveProjectedField();
            results = obj.designFields.derivedProjectedField;
            %Validate Results
            s.results = results;
            B = ProjectedFieldDerivatorTester(3);
            B.loadResults(s);
            B.validate;
        end
        function testFilteredFieldDerivator(obj)
%            obj.designFields.deriveFilteredField();
            results = obj.designFields.derivedFilteredField;
            %Validate Results
            s.results = results;
            B = FilteredFieldDerivatorTester(3);
            B.loadResults(s);
            B.validate;
        end
    end
    methods (Access = private)
        function createFields(obj)
            %Load Initial Paramaters

            file = fullfile("DensityBasedProjection",'Test','Data','projectorParameters.mat');            
            a = load(file);
            s.projectorParameters = a.projectorParameters;

            file = fullfile("DensityBasedProjection",'Test','Data','mesh.mat');            
            a = load(file);
            s.mesh = a.mesh;

            file = fullfile("DensityBasedProjection",'Test','Data','field.mat');            
            a = load(file);
            s.field = a.field;

            file = fullfile("DensityBasedProjection",'Test','Data','filterParameters.mat');            
            a = load(file);
            s.filterParameters = a.filterParameters;

            %Create the DesignFieldObject
            obj.designFields = DesignField(s);
        end 

    end
end

