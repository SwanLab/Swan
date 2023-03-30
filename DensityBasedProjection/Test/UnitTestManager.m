close all
clc
clear all
f = figure('Position',[200 200 250 450]);

% Create the buttons
btnGlobalTest = uicontrol('Style','pushbutton','String','Run Global Test','Position',[50 50 150 30],'BackgroundColor', 'r','Callback',@btn1Callback);
btnFilterTest = uicontrol('Style','pushbutton','String','Run Filter Test','Position',[50 90 150 30],'BackgroundColor', [0.8 0.8 0.8],'Callback',@btn2Callback);
btnProjectorTest = uicontrol('Style','pushbutton','String','Run Projector Test','Position',[50 130 150 30],'BackgroundColor', [0.8 0.8 0.8],'Callback',@btn3Callback);
btnFilteredFieldDerivatorTest = uicontrol('Style','pushbutton','String','Run Projected Field Derivator Test','Position',[50 170 180 30],'BackgroundColor', [0.8 0.8 0.8],'Callback',@btn4Callback);
btnProjectorFieldDerivatorTest =  uicontrol('Style','pushbutton','String','Run Filtered Field Derivator Test','Position',[50 210 180 30],'BackgroundColor', [0.8 0.8 0.8],'Callback',@btn5Callback);

btnFEMComputerTest =  uicontrol('Style','pushbutton','String','Run FEM Test','Position',[50 250 180 30],'BackgroundColor', [0.3569 0.8118 0.9569],'Callback',@btn6Callback);
btnDisplacementComputerTest =  uicontrol('Style','pushbutton','String','Run Displacement Test','Position',[50 290 180 30],'BackgroundColor', [0.3569 0.8118 0.9569],'Callback',@btn7Callback);
btnForceComputerTest =  uicontrol('Style','pushbutton','String','Run Force Test','Position',[50 330 180 30],'BackgroundColor',[0.3569 0.8118 0.9569],'Callback',@btn8Callback);
btnPenalizerTest =  uicontrol('Style','pushbutton','String','Run Penalize Test','Position',[50 370 180 30],'BackgroundColor', [0.3569 0.8118 0.9569],'Callback',@btn9Callback);
btnStifnessMatrixComputerTest =  uicontrol('Style','pushbutton','String','Run StifnessMatrix Test','Position',[50 410 180 30],'BackgroundColor', [0.3569 0.8118 0.9569],'Callback',@btn10Callback);

function btn1Callback(hObject,eventdata)
  GlobalTest(3)
  disp('-----------')
end

function btn2Callback(hObject,eventdata)
B = DesignFieldTester;
B.testFilter;
  disp('-----------')
end

function btn3Callback(hObject,eventdata)
B = DesignFieldTester;
B.testProjector;
  disp('-----------')


end
function btn4Callback(hObject,eventdata)
B = DesignFieldTester;
B.testFilteredFieldDerivator;
  disp('-----------')
end

function btn5Callback(hObject,eventdata)
B = DesignFieldTester;
B.testProjectedFieldDerivator;
  disp('-----------')
end
function btn6Callback(hObject,eventdata)
B = FEMcomputerTester(3);
B.compute;
B.validate;
  disp('-----------')
end

function btn7Callback(hObject,eventdata)
B = DisplacementComputerTester(3);
B.compute;
B.validate;
  disp('-----------')
end
function btn8Callback(hObject,eventdata)
B = ForceComputerTester(3);
B.compute;
B.validate;
  disp('-----------')
end
function btn9Callback(hObject,eventdata)
B = PenalizerTester(3);
B.compute;
B.validate;
  disp('-----------')
end
function btn10Callback(hObject,eventdata)
B = StifnessMatrixComputerTester(3);
B.compute;
B.validate;
  disp('-----------')
end



