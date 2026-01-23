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
    
    methods (Access = public)
        
        function loadParams(obj,p)
            if ~isempty(p)
                if isstruct(p)
                    obj.loadFromStruct(p);
                elseif isobject(p)
                    obj.loadFromObject(p);
                else
                    obj.loadFromFile(p);
                end
            end
        end
        
        function s = getParams(obj)
            s = struct;
            props = properties(obj);
            for i = 1:length(props)
                prop = props{i};
                s.(prop) = obj.(prop);
            end
        end
        
    end
    
    methods (Access = protected)
        
        function obj = AbstractSettings()
            obj.loadParams(obj.defaultParamsName);
        end
        
    end
    
    methods (Access = private)
        
        function loadFromStruct(obj,s)
            obj.cParams = s;
            obj.loadedFile = 'Settings loaded from struct.';
            obj.assignParams();
        end
        
        function loadFromObject(obj,obj2)
            s = obj2.getParams();
            obj.loadFromStruct(s);
        end
        
        function loadFromFile(obj,f)
            switch obj.getFileType(f)
               case {'','.m'}
                   error('Only JSON files accepted!');
                   obj.loadParamsFromMatlabScript(f);
               case '.json'
                   obj.loadParamsFromJSON(f);
               otherwise
                   error('Invalid extension');
            end
            %f = [f,'.json'];
            %obj.loadParamsFromJSON(f);
            obj.loadedFile = f;
            obj.assignParams();
        end
        
        function loadParamsFromJSON(obj,fileName)
            obj.cParams = jsondecode(fileread(fileName));
        end
        
        function loadParamsFromMatlabScript(obj,fileName)
            run(fileName);
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
                    obj.warnOfInvalidConstructionParams(param);
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