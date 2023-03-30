close all
clc
clear all
f = figure('Position',[200 200 200 150]);

% Create the buttons
btnGlobalTest = uicontrol('Style','pushbutton','String','Run Global Test','Position',[50 50 150 30],'BackgroundColor', 'r','Callback',@btn1Callback);
btnFilterTest = uicontrol('Style','pushbutton','String','Run Filter Test','Position',[50 90 150 30],'BackgroundColor', [0.8 0.8 0.8],'Callback',@btn2Callback);
btnProjectorTest = uicontrol('Style','pushbutton','String','Run Projector Test','Position',[50 130 150 30],'BackgroundColor', [0.8 0.8 0.8],'Callback',@btn3Callback);
btnFilteredFieldDerivatorTest = uicontrol('Style','pushbutton','String','Run Projected Field Derivator Test','Position',[50 170 180 30],'BackgroundColor', [0.8 0.8 0.8],'Callback',@btn4Callback);
btnProjectorFieldDerivatorTest =  uicontrol('Style','pushbutton','String','Run Filtered Field Derivator Test','Position',[50 210 180 30],'BackgroundColor', [0.8 0.8 0.8],'Callback',@btn5Callback);

function btn1Callback(hObject,eventdata)
  GlobalTest(3)
  disp('-----------')
end

function btn2Callback(hObject,eventdata)
designFieldTester = DesignFieldTester;
designFieldTester.testFilter;
  disp('-----------')
end

function btn3Callback(hObject,eventdata)
designFieldTester = DesignFieldTester;
designFieldTester.testProjector;
  disp('-----------')


end
function btn4Callback(hObject,eventdata)
designFieldTester = DesignFieldTester;
designFieldTester.testFilteredFieldDerivator;
  disp('-----------')
end

function btn5Callback(hObject,eventdata)
designFieldTester = DesignFieldTester;
designFieldTester.testProjectedFieldDerivator;
  disp('-----------')
end
