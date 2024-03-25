close all

type0 = "Square";

[density,matHomog] = runCases(type0);
 Chomog = createFunctions(density,matHomog);
 Ciso = computeIsotropicDamage();
 plot(Ciso,Chomog);

 function [density,matHomog] = runCases(type0)
     if type0 == "Circle"
         max = 0.5;
     elseif type0 == "Square"
         max = 0.94;
     elseif type0 == "Full"
         max = 1;
     end
     squareLength = linspace(0,max,30);
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
             area = mat.mesh.computeVolume;
             holeArea = 1 - area;
             diameter = 2*l;
             perimeter = 2*pi*l/4;
         elseif type == "Square"
             area = mat.mesh.computeVolume;
             holeArea = 1 - area;
             diameter = l;
             perimeter = 4*l/4;
         elseif type == "Full"
             area = mat.mesh.computeVolume;
             holeArea = 1 - area;
             diameter = 0;
             perimeter = 0;
         end
         density(i,1) = perimeter;
     end
 end

 function  Ciso = computeIsotropicDamage()
    C = zeros(3,3);
    E = 1;
    v = 0.3;
    constant = E/(1-v^2);

    C(1,1) = constant; 
    C(1,2) = constant*v;
    C(2,1) = constant*v;
    C(2,2) = constant;
    C(3,3) = constant*(1-v)/2;

    x = linspace(0,1,30);
    for i=1:3
        for j=1:3
            sM.coord = x';
            sM.connec = [1:length(x)-1]' + [0,1];
            s.mesh = Mesh.create(sM);
            s.fValues = [(1-x).^2*C(i,j)]';
            s.order = 'P1';
            Ciso{i,j} = LagrangianFunction(s);
        end
    end
 end

function Chomog = createFunctions(density,matHomog)
for i=1:3
    for j=1:3
        sM.coord = density;
        sM.connec = [1:length(density)-1]' + [0,1];
        s.mesh = Mesh.create(sM);
        s.fValues = squeeze(matHomog(i,j,:));
        s.order = 'P1';
        Chomog{i,j} = LagrangianFunction(s).project('P2');
    end
end
end

function plot(Ciso,Chomog)
figure(2)
for i=1:3
    for j=1:3
        k = 3*(i-1)+j;
        subplot(3,3,k)
        Chomog{i,j}.plot;
        %set ( gca, 'XDir', 'reverse' ) %% ONLY FOR AREA DENSITY
        hold on
        Ciso{i,j}.plot
        hold off
        title(['C',num2str(i),num2str(j)]);
    end
end
end











