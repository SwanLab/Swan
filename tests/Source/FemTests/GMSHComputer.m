classdef GMSHComputer < handle

    properties (Access = public)
        computation
        variables
    end

    properties (Access = private)
        testName
    end

    methods (Access = public)
        function obj = GMSHComputer(cParams)
            obj.testName = cParams.testName;
        end

        function compute(obj)
            filePath = fullfile('tests','Source','ReadingFilesTests','ReadingFiles','testReadingGmsh.msh');
            reader = GmsReader(filePath);
            reader.read();
            readData = reader.getDataBase();
            readData.isElemInThisSet = readData.isElemInThisSet{1};
            obj.variables = readData;
            obj.computation = obj;
        end
    end
end