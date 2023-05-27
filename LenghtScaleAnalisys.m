clc
clear all
close all
radio = 'R1h';
for beta = [0.1 0.25 0.5 0.75 1 1.25 1.5 1.75 2]
        for eta = [0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45]
            %sound(sin(2*pi*800*linspace(0, 0.1, 0.1 * 44100)),44100)
            close all
            iter = 150;
            disp(['Eta:',num2str(eta*100),'  Beta:',num2str(beta),'  Iter:',num2str(iter)]);
            scaleParameters.eta = eta;
            scaleParameters.beta = beta;
            save('scaleParameters.mat', 'scaleParameters')
            s.testName = 'test_arturo';
            t = TopOptComputer(s);
            t.compute();
            %sound(sin(2*pi*600*linspace(0, 0.5, 0.5 * 44100)),44100)
            nameFile = ['Eta',num2str(eta*100),'Beta',num2str(beta*100),radio,'Iter',num2str(iter),'Morph']
            rutaGraf = ['C:\Users\artur\Documents\GitHub\SWAM\Swan\LenghtScaleAnalysis\',nameFile,'.png'];
            rutaArch = ['C:\Users\artur\Documents\GitHub\SWAM\Swan\LenghtScaleAnalysis\',nameFile,'.m'];
            complianceValue = (t.computation.constraint.shapeFunctions{3}.value+t.computation.cost.shapeFunctions{1}.value);
            volumenValue = (1+t.computation.constraint.shapeFunctions{4}.value)*t.computation.constraint.shapeFunctions{4}.vTarget;
            title(['Compliance (adim): ',num2str(complianceValue),'     Volumen (adim): ',num2str(volumenValue)])
            saveas(2, rutaGraf);
        end 
end 
%sound(sin(2*pi*300*linspace(0, 2.5, 2.5 * 44100)),44100)