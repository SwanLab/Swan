classdef SettingsVademecumCellVariablesCalculator < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsVademecumCellVariablesCalculator.json'
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
        smoothingExponentSettings
    end
    
    methods (Access = public)
        
        function obj = SettingsVademecumCellVariablesCalculator(varargin)
            obj.freeFemSettings = SettingsFreeFemMeshGenerator();            
            switch nargin
                case 1
                    obj.loadParams(varargin{1});
            end
        end
        
    end
end