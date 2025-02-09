classdef FieldRepresenterPrinter < handle

    properties (Access = protected)
        fileID
        fieldName
        iter
        fieldType
        fieldPosition
        simulationStr
        fieldValues
        nComp
    end

    methods (Access = public, Abstract)
        printResultsLineHeader(obj)
        printFieldLines(obj)
    end


end