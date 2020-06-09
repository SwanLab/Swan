classdef ClassWritter < handle
    
    properties (Access = private)
        lines
        string2Print
    end
    
    properties (Access = private)
        fileName
    end
    
    methods (Access = public)
        
        function obj = ClassWritter(cParams)
            obj.init(cParams)
            obj.computeLines();
            obj.computeString();
            obj.createNewFile();
            obj.applySmartIndent();
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.fileName = cParams;
        end
        
        function computeLines(obj)
            fName = obj.fileName;
            l = cell(0);
            l{end+1} = ['classdef ',fName,' < handle\n'];
            l{end+1} = '\n';
            l{end+1} = 'properties (Access = public)\n';
            l{end+1} = '\n';
            l{end+1} = 'end';
            l{end+1} = '\n\n';
            l{end+1} = 'properties (Access = private)\n';
            l{end+1} = '\n';
            l{end+1} = 'end';
            l{end+1} = '\n\n';
            l{end+1} = 'properties (Access = private)\n';
            l{end+1} = '\n';
            l{end+1} = 'end';
            l{end+1} = '\n\n';
            l{end+1} = 'methods (Access = public)\n';
            l{end+1} = '\n';
            l{end+1} = ['function obj = ',fName,'(cParams)\n'];
            l{end+1} = 'obj.init(cParams)\n';
            l{end+1} = '\n';
            l{end+1} = 'end';
            l{end+1} = '\n\n';
            l{end+1} = 'end';
            l{end+1} = '\n\n';
            l{end+1} = 'methods (Access = private)\n';
            l{end+1} = '\n';
            l{end+1} = 'function init(obj,cParams)\n';
            l{end+1} = '\n';
            l{end+1} = 'end';
            l{end+1} = '\n\n';
            l{end+1} = 'end';
            l{end+1} = '\n\n';
            l{end+1} = 'end';
            obj.lines = l;
        end
        
        function computeString(obj)
            str2Print = strcat(obj.lines{:});
            str = sprintf(str2Print);
            obj.string2Print = str;
        end
        
        function createNewFile(obj)
            editorService = com.mathworks.mlservices.MLEditorServices;
            editorApp = editorService.getEditorApplication();
            editorApp.newEditor(obj.string2Print);
        end
        
    end
    
    methods (Access = private, Static)
        
        function applySmartIndent()
            h = matlab.desktop.editor.getActive;
            h.smartIndentContents            
        end
        
    end
    
end