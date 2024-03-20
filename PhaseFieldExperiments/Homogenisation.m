type0 = "Square";

[density,matHomog] = runCases(type0);
 C = createFunctions(density,matHomog);
 plot(C);

 function [density,matHomog] = runCases(type0)
     if type0 == "Circle"
         max = 0.5;
     elseif type0 == "Square"
         max = 0.94;
     end
     squareLength = linspace(0,max,20);
     matHomog = zeros(3,3,length(squareLength));
     density = zeros(length(squareLength),1);
    
     for i=1:length(squareLength)
         l = squareLength(i);
         if l<0.015
             type = "Full";
             l = 0;
         else
             type = type0;
         end
         mat = Tutorial02p2FEMElasticityMicro(l,type);
    
         figure(1)
         cla reset
         mat.mesh.plot
    
    
         matHomog(:,:,i) = mat.stateProblem.Chomog;
         if type == "Circle"
             holeArea = pi*l^2;
         elseif type == "Square"
             holeArea = l^2;
         elseif type == "Full"
             holeArea = 0;
         end
         density(i,1) = mat.mesh.computeVolume;
         %density(i,1) = 1-2*l;
     end
 end

function C = createFunctions(density,matHomog)
C = {};
for i=1:3
    for j=1:3
        sM.coord = density;
        sM.connec = 1:length(density)-1 + [1,0];
        %sM.kFace = 0;
        s.mesh = Mesh.create(sM);
        s.fValues = squeeze(matHomog(i,j,:));
        s.order = 'P1';
        C{i,j} = LagrangianFunction(s);
    end
end
end

function plot(C)
figure(2)
for i=1:3
    for j=1:3
        k = 3*(i-1)+j;
        subplot(3,3,k)
        C{i,j}.plot;
        set ( gca, 'XDir', 'reverse' ) 
        title(['C',num2str(i),num2str(j)]);
    end
end
end











