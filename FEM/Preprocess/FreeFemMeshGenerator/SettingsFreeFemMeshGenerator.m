classdef SettingsFreeFemMeshGenerator < AbstractSettings
    
    properties (Access = protected)
        defaultParamsName = 'paramsFreeFemMeshGenerator.json'
    end
    
    properties (Access = public)
        mxV 
        myV 
        hMax
        qNorm
        intElements
        borderElements
        fileName
        printingDir
        freeFemFileName
    end
    
    methods (Access = public)
        
        function obj = SettingsFreeFemMeshGenerator(varargin)
            switch nargin
                case 1
                    obj.loadParams(varargin{1});
            end
        end
        
    end
    
end