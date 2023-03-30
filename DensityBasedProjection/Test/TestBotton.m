classdef TestBotton < handle
    properties (Access = private)
        name
        size
        botton
    end

    methods (Access = public)
        function obj = TestBotton(name,size)
            obj.name = name;
            obj.size = size;
            obj.createButton();
        end
    end
    methods (Access = private)
        function createButton(obj)
            obj.botton = uicontrol('Style','pushbutton','String',obj.name,'Position',obj.size,'Callback',@obj.Callback);
        end
    end
    methods (Static)
        function Callback(hObject,eventdata)
            designFieldTester = DesignFieldTester;
            designFieldTester.TestFilter;
        end

    end
end