clc
clear all
close all
beta = 1e-5;
%matValues = zeros(49,5);
contador = 1;
for radio = [1] %0.14 0.28 0.43 0.57 0.72 0.86 1.01 1.73 2.02 2.88
        for eta = [0.5] %0.15 0.2 0.25 0.3 0.35 0.4 0.45
            %sound(sin(2*pi*800*linspace(0, 0.1, 0.1 * 44100)),44100)
            close all
            iter = 150;
            disp(['Eta:',num2str(eta*100),'  Radio:',num2str(radio),'  Iter:',num2str(iter)]);
            scaleParameters.eta = eta;
            scaleParameters.beta = beta;
            scaleParameters.rad = radio;
            save('scaleParameters.mat', 'scaleParameters')
            s.testName = 'test_arturo';
            t = TopOptComputer(s);
            t.compute();
            %sound(sin(2*pi*600*linspace(0, 0.5, 0.5 * 44100)),44100)
            nameFile = ['Eta',num2str(eta*100),'R',num2str(radio*100),'h','Beta',num2str(beta*100),'Iter',num2str(iter),'Morph']
            rutaGraf = ['C:\Users\artur\Documents\GitHub\SWAM\Swan\LenghtScaleAnalysis\',nameFile,'.png'];
            rutaArch = ['C:\Users\artur\Documents\GitHub\SWAM\Swan\LenghtScaleAnalysis\',nameFile,'.m'];
            complianceValue = (t.computation.constraint.shapeFunctions{3}.value+t.computation.cost.shapeFunctions{1}.value);
            volumenValue = (1+t.computation.constraint.shapeFunctions{4}.value)*t.computation.constraint.shapeFunctions{4}.vTarget;
            title(['Compliance (adim): ',num2str(complianceValue),'     Volume (adim): ',num2str(volumenValue)])
            saveas(2, rutaGraf);
     %       matValues(contador,:) = [eta,radio,beta,complianceValue,volumenValue];
     %       contador = contador+1;
        end 
end 
%sound(sin(2*pi*300*linspace(0, 2.5, 2.5 * 44100)),44100)