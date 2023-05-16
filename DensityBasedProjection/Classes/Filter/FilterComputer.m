classdef FilterComputer < handle 
    properties (Access = public)
        filteredField
    end
    properties (Access = private)
        filterParameters
        
    end
    methods (Access = public)
        function obj = FilterComputer(cParams)
            obj.inputData(cParams);
        end

        function filteredField = compute(obj,field)
            filteredField = (obj.filterParameters.H*field(:))./obj.filterParameters.Hs;
        end
    end
    methods (Access = private)
        function inputData(obj,cParams)
            obj.filterParameters =cParams.filterParameters;
        end
        function filter(obj,field)
            
        end 
    end
end