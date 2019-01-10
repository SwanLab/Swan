classdef DensityGaussResultsPrinter < ResultsPrinter ...
        

    properties (Access = protected)
        simulationStr = 'DensityGauss';
    end

    properties (Access = private)
        fieldName = 'RegularizedDensity';
        headPrinter = GaussHeadPrinter;
    end

    methods (Access = public)

        function obj = DensityGaussResultsPrinter(d)
            obj.init(d); 
            obj.printHeader();
        end

    end

    methods (Access = protected)

        function printHeader(obj)
            d.fileID = obj.fileID;
            d.gaussDescriptor = obj.gaussDescriptor;
            d.etype = obj.etype;
            d.ngaus = obj.ngaus;
            d.ndim  = obj.ndim;
            d.posgp = obj.posgp;
            obj.headPrinter.print(d);
        end

        function printResults(obj)
            dens = obj.fields;
            dS = obj.createScalarGaussDataBase(dens, obj.fieldName,'OnGaussPoints');
            ScalarGaussPrinter(dS);
        end

    end

end