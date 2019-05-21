classdef AbstractSettings < handle
    
    properties (GetAccess = public, SetAccess = private)
        loadedFile
    end
    
    properties (Access = protected, Abstract)
        defaultParamsName
    end
    
    properties (GetAccess = protected, SetAccess = private)
        cParams
    end
    
    methods (Access = protected)
        
        function obj = AbstractSettings()
            obj.loadParams(obj.defaultParamsName);
        end
        
        function loadParams(obj,p)
            if ~isempty(p)
                if isstruct(p)
                    obj.cParams = p;
                    obj.loadedFile = 'Settings loaded from struct.';
                else
                    switch obj.getFileType(p)
                        case {'','.m'}
                            error('Only JSON files accepted!');
                            obj.loadParamsFromMatlabScript(p);
                        case '.json'
                            obj.loadParamsFromJSON(p);
                        otherwise
                            error('Invalid extension');
                    end
                end
                obj.assignParams();
            end
        end
        
    end
    
    methods (Access = private)
        
        function loadParamsFromJSON(obj,paramsFileName)
            obj.loadedFile = paramsFileName;
            obj.cParams = jsondecode(fileread(paramsFileName));
        end
        
        function loadParamsFromMatlabScript(obj,paramsFileName)
            obj.loadedFile = paramsFileName;
            run(paramsFileName);
            obj.cParams = who;
            obj.clearConstructionParams();
            s = struct;
            for i = 1:length(obj.cParams)
                param = obj.cParams{i};
                if isprop(obj,param)
                    s.(param) = eval(param);
                else
                    obj.warnOfInvalidConstructionParams(param);
                end
            end
            obj.cParams = s;
        end
        
        function assignParams(obj)
            fields = fieldnames(obj.cParams);
            for i = 1:length(fields)
                param = fields{i};
                if isprop(obj,param)
                    obj.(param) = obj.cParams.(param);
                else
                    obj.warnOfInvalidConstruction(param);
                end
            end
        end
        
        function clearConstructionParams(obj)
            obj.removeVar('paramsFileName');
            obj.removeVar('obj');
        end
        
        function warnOfInvalidConstructionParams(obj,param)
            warning([param ' is not a property of ' class(obj)]);
        end
        
        function removeVar(obj,var)
            pos = strcmp(obj.cParams,var);
            obj.cParams(pos) = [];
        end
        
    end
    
    methods (Access = private, Static)
        
        function ext = getFileType(fileName)
            [~,~,ext] = fileparts(fileName);
        end
        
    end
    
end