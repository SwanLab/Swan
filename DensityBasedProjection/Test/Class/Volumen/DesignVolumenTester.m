classdef DesignVolumenTester < handle
    properties (Access = private)
        iterations
        tolerateError
        designVolumen
        data
    end

    methods (Access = public)
        function obj = DesignVolumenTester()
            obj.iterations = 3;
            obj.tolerateError = 1e-10;
            obj.createVolumen();
        end
        function testVolumen(obj)
            obj.designVolumen.computeVolumen();
            results = obj.designVolumen.volumen;
            %Validate Results
            s.results = results;
            B = VolumenComputerTester(3);
            B.loadResults(s);
            B.validate;
        end
        function testVolumenDerivator(obj)
            obj.designVolumen.deriveVolumen();
            results = obj.designVolumen.derivedVolumen;
            %Validate Results
            s.results = results;
            B = VolumenDerivatorTester(3);
            B.loadResults(s);
            B.validate;
        end
        function testVolumenFraction(obj)
            obj.designVolumen.computeVolumenFraction(obj.data.D,obj.data.I);
            results = obj.designVolumen.volumenFraction;
            %Validate Results
            s.results = results;
            B = VolumenFractionTester(3);
            B.loadResults(s);
            B.validate;
        end
    end
    methods (Access = private)
        function createVolumen(obj)
            %Load Initial Paramaters

            file = fullfile("DensityBasedProjection",'Test','Data','filterParameters.mat');            
            a = load(file);
            s.filterParameters = a.filterParameters;

            file = fullfile("DensityBasedProjection",'Test','Data','mesh.mat');            
            a = load(file);
            s.mesh = a.mesh;

            file = fullfile("DensityBasedProjection",'Test','Data','designFields.mat');            
            a = load(file);
            s.designField= a.designFields;
            s.designField.deriveProjectedField
            file = fullfile("DensityBasedProjection",'Test','Data','volfracD.mat');            
            a = load(file);
            s.volumenFraction = a.volfracD;
            file = fullfile("DensityBasedProjection",'Test','Data','D.mat');            
            a = load(file);
            obj.data.D = a.D;
            file = fullfile("DensityBasedProjection",'Test','Data','I.mat');            
            a = load(file);
            obj.data.I = a.I;

            %Create the DesignFieldObject
            obj.designVolumen = DesignVolumen(s);
            obj.designVolumen.computeVolumenFraction(obj.data.D,obj.data.I);            
        end 

    end
end

