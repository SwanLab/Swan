classdef AbstractSettings < handle
    
    properties (GetAccess = public, SetAccess = private)
        loadedFile
        cParams
    end
    
    properties (Access = protected, Abstract)
        defaultParamsName
    end
    
    properties (GetAccess = protected, SetAccess = private)
        customParams
    end
    
    methods (Access = protected)
        
        function obj = AbstractSettings()
            obj.loadParams(obj.defaultParamsName);
        end
        
        function loadParams(obj,p)
            if ~isempty(p)
                if isstruct(p)
                    obj.customParams = p;
                else
                    switch obj.getFileType(p)
                        case {'','.m'}
                            error('THIS IS NO JSON FILE!');
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
            obj.customParams = jsondecode(fileread(paramsFileName));
        end
        
        function loadParamsFromMatlabScript(obj,paramsFileName)
            obj.loadedFile = paramsFileName;
            run(paramsFileName);
            obj.customParams = who;
            obj.clearCustomParams();
            s = struct;
            for i = 1:length(obj.customParams)
                param = obj.customParams{i};
                if isprop(obj,param)
                    s.(param) = eval(param);
                else
                    obj.warnOfInvalidCustomParams(param);
                end
            end
            obj.customParams = s;
        end
        
        function assignParams(obj)
            fields = fieldnames(obj.customParams);
            for i = 1:length(fields)
                param = fields{i};
                if isprop(obj,param)
                    obj.(param) = obj.customParams.(param);
                else
                    obj.warnOfInvalidCustomParams(param);
                end
            end
        end
        
        function clearCustomParams(obj)
            obj.removeVar('paramsFileName');
            obj.removeVar('obj');
        end
        
        function warnOfInvalidCustomParams(obj,param)
            warning([param ' is not a property of ' class(obj)]);
        end
        
        function removeVar(obj,var)
            pos = strcmp(obj.customParams,var);
            obj.customParams(pos) = [];
        end
        
    end
    
    methods (Access = private, Static)
        
        function ext = getFileType(fileName)
            [~,~,ext] = fileparts(fileName);
        end
        
    end
    
end