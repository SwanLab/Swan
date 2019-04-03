classdef SettingsVademecumCellVariablesCalculator < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsVademecumCellVariablesCalculator'
    end
    
    properties (Access = public)
        fileName
        freeFemFileName
        mxMin
        mxMax
        myMin
        myMax
        nMx
        nMy
        outPutPath
        freeFemSettings
        print
    end
    
    methods (Access = public)
        
        function obj = SettingsVademecumCellVariablesCalculator(varargin)
            switch nargin
                case 1
                    obj.loadParams(varargin{1});
            end
        end
        
    end
end