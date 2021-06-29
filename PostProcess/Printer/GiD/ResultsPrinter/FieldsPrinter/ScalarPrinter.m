classdef ScalarPrinter < FieldPrinter
    
    properties (Access = protected)
        fieldType = 'Scalar';
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = protected)
        
        function print(obj)
            obj.printResultsLineHeader()
            obj.printValuesLine();
            obj.printFieldLines();
            obj.printEndValuesLine();
        end       
        
     function printFieldLines(obj)
            obj.fieldRepresenter.printFieldLines()
        end
        
        function printResultsLineHeader(obj)
           obj.fieldRepresenter.printResultsLineHeader()
        end         
        
    end
        
    
end