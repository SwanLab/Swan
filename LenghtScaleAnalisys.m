clc
clear all
close all
beta = 1.5;
iter = 150;
radioVec = [0.28 0.57 0.86 1.15 1.44 1.73 2.02]; %0.28 0.57 0.86 1.15 1.44 1.73 2.02
etaVec = [0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45]; %0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45       0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
matValues = zeros(length(etaVec)*length(radioVec),7);
contador = 1;

for radio = radioVec
        for eta = etaVec
            %sound(sin(2*pi*800*linspace(0, 0.1, 0.1 * 44100)),44100)
            close all
            
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
            rutaArch = ['C:\Users\artur\Documents\GitHub\SWAM\Swan\LenghtScaleAnalysis\',nameFile,'.mat'];
            complianceValue = (t.computation.constraint.shapeFunctions{2}.value+t.computation.cost.shapeFunctions{1}.value);
            complianceEroded = (t.computation.constraint.shapeFunctions{1}.value+t.computation.cost.shapeFunctions{1}.value);
            complianceDilated = (t.computation.constraint.shapeFunctions{3}.value+t.computation.cost.shapeFunctions{1}.value);            
            volumenValue = (1+t.computation.constraint.shapeFunctions{4}.value)*t.computation.constraint.shapeFunctions{4}.vTarget;
            %complianceValue = t.computation.cost.shapeFunctions{1}.value;
            %volumenValue = t.computation.constraint.shapeFunctions{1}.value;
            title(['Compliance (adim): ',num2str(complianceValue),'     Volume (adim): ',num2str(volumenValue)])
            saveas(2, rutaGraf);
            matValues(contador,:) = [eta,radio,beta,complianceValue,complianceEroded,complianceDilated,volumenValue];
            contador = contador+1;

%             p.mesh    = t.computation.designVariable.mesh;
%             p.fValues = t.computation.designVariable.value;
%             s.testName = 'test_arturo';
%             s.filename = s.testName;
%             Result = P1Function(p);
%             Result.print(s);
%             saveStruc.topOpt = t;
%             %saveStruc.GID = Result;
%             %save(rutaArch,'-struct','saveStruc');
        end 
end 
save('C:\Users\artur\Documents\GitHub\SWAM\Swan\LenghtScaleAnalysis\Cantiliver3BridgeTest2DataMorph.mat','matValues')
%sound(sin(2*pi*300*linspace(0, 2.5, 2.5 * 44100)),44100)