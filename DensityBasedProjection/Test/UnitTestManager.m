close all
clc
clear all
f = figure('Position',[200 60 280 650]);

% Create the buttons
btnGlobalTest = uicontrol('Style','pushbutton','String','Run Global Test','Position',[50 10 150 30],'BackgroundColor', [1.0000 0.1020 0.1020],'Callback',@btn1Callback);
btnOptimizerTest =  uicontrol('Style','pushbutton','String','Run Optimizer Test','Position',[50 50 150 30],'BackgroundColor', [1.0000 0.8000 0],'Callback',@btnOptCallback);
btnFilterTest = uicontrol('Style','pushbutton','String','Run Filter Test','Position',[50 90 150 30],'BackgroundColor', [0.8 0.8 0.8],'Callback',@btn2Callback);
btnProjectorTest = uicontrol('Style','pushbutton','String','Run Projector Test','Position',[50 130 150 30],'BackgroundColor', [0.8 0.8 0.8],'Callback',@btn3Callback);
btnFilteredFieldDerivatorTest = uicontrol('Style','pushbutton','String','Run Projected Field Derivator Test','Position',[50 170 180 30],'BackgroundColor', [0.8 0.8 0.8],'Callback',@btn4Callback);
btnProjectorFieldDerivatorTest =  uicontrol('Style','pushbutton','String','Run Filtered Field Derivator Test','Position',[50 210 180 30],'BackgroundColor', [0.8 0.8 0.8],'Callback',@btn5Callback);

btnFEMComputerTest =  uicontrol('Style','pushbutton','String','Run FEM Test','Position',[50 250 180 30],'BackgroundColor', [0 0.5843 0.7137],'Callback',@btn6Callback);
btnDisplacementComputerTest =  uicontrol('Style','pushbutton','String','Run Displacement Test','Position',[50 290 180 30],'BackgroundColor', [0.3569 0.8118 0.9569],'Callback',@btn7Callback);
btnForceComputerTest =  uicontrol('Style','pushbutton','String','Run Force Test','Position',[50 330 180 30],'BackgroundColor',[0.3569 0.8118 0.9569],'Callback',@btn8Callback);
btnPenalizerTest =  uicontrol('Style','pushbutton','String','Run Penalize Test','Position',[50 370 180 30],'BackgroundColor', [0.3569 0.8118 0.9569],'Callback',@btn9Callback);
btnStifnessMatrixComputerTest =  uicontrol('Style','pushbutton','String','Run StifnessMatrix Test','Position',[50 410 180 30],'BackgroundColor', [0.3569 0.8118 0.9569],'Callback',@btn10Callback);


btnCostFieldDerivator =  uicontrol('Style','pushbutton','String','Run Cost Field Derivator Test','Position',[50 450 180 30],'BackgroundColor', [0.56, 0.93, 0.56],'Callback',@btn11Callback);

btnVolumen =  uicontrol('Style','pushbutton','String','Run Volumen Test','Position',[50 490 180 30],'BackgroundColor', [0.8039 0.5647 0.5137],'Callback',@btn12Callback);
btnVolumenFraction =  uicontrol('Style','pushbutton','String','Run Fraction Volumen Test','Position',[50 530 180 30],'BackgroundColor', [0.8039 0.5647 0.5137],'Callback',@btn13Callback);
btnDerivedVolumen =  uicontrol('Style','pushbutton','String','Run Derived Volumen Test','Position',[50 570 180 30],'BackgroundColor', [0.8039 0.5647 0.5137],'Callback',@btn14Callback);



function btn1Callback(hObject,eventdata)
  clc
  GlobalTest(3);
  disp('-----------')
end

function btn2Callback(hObject,eventdata)
clc
B = DesignFieldTester;
B.testFilter;
  disp('-----------')
end

function btn3Callback(hObject,eventdata)
clc
B = DesignFieldTester;
B.testProjector;
  disp('-----------')


end
function btn4Callback(hObject,eventdata)
clc
B = DesignFieldTester;
B.testFilteredFieldDerivator;
  disp('-----------')
end

function btn5Callback(hObject,eventdata)
clc
B = DesignFieldTester;
B.testProjectedFieldDerivator;
  disp('-----------')
end
function btn6Callback(hObject,eventdata)
clc
B = CostComputerTester;
B.testFEM;
  disp('-----------')
end

function btn7Callback(hObject,eventdata)
clc
B = DisplacementComputerTester(3);
B.compute;
B.validate;
  disp('-----------')
end
function btn8Callback(hObject,eventdata)
clc
B = ForceComputerTester(3);
B.compute;
B.validate;
  disp('-----------')
end
function btn9Callback(hObject,eventdata)
clc
B = PenalizerTester(3);
B.compute;
B.validate;
  disp('-----------')

end
function btn10Callback(hObject,eventdata)
clc
B = StifnessMatrixComputerTester(3);
B.compute;
B.validate;
  disp('-----------')
end
function btn11Callback(hObject,eventdata)
B = CostComputerTester;
B.testCostDerivator;
  disp('-----------')
end
function btn12Callback(hObject,eventdata)
clc
B = DesignVolumenTester;
B.testVolumen;
  disp('-----------')
end
function btn13Callback(hObject,eventdata)
clc
B = DesignVolumenTester;
B.testVolumenFraction;
  disp('-----------')
end
function btn14Callback(hObject,eventdata)
clc
B = DesignVolumenTester;
B.testVolumenDerivator;
  disp('-----------')
end
function btnOptCallback (hObject,eventdata)
clc
B = OptimizerTester;
B.validate;
  disp('-----------')
end 

