PFH = PhaseFieldHomogenisator();
% [mat,phi] = PFH.computeHomogMaterial("Circle","Perimeter",100);
% save('CircleMicroDamagePerimeter','mat','phi')
% [mat,phi] = PFH.computeHomogMaterial("Circle","Area",100);
% save('CircleMicroDamageArea','mat','phi')
% 
[mat,phi]  = PFH.computeHomogMaterial("Square","Perimeter",100);
save('SquareMicroDamagePerimeter','mat','phi')
% [mat,phi]  = PFH.computeHomogMaterial("Square","Area",100);
% save('SquareMicroDamageArea','mat','phi')

% [mat,phi] = PFH.computeIsotropicMaterial("AT1",100);
% save('IsoMicroDamage','mat','phi')
% [~,AT2] = PFH.computeIsotropicMaterial("AT2",30);

%%% PLOT %%%
close all
figure(2)
tiledlayout(3,3)
for i=1:3
    for j=1:3
        k = 3*(i-1)+j;
        nexttile
        hold on

        plot(squeeze(P_Circle,C_P_Circle(i,j,:)),'.','Color','#0072BD');
        plot(squeeze(A_Circle,C_A_Circle(i,j,:)),'.','Color','#D95319');

        plot(squeeze(P_Square,C_P_Square(i,j,:)),'--','Color','#0072BD');
        plot(squeeze(A_Square,C_A_Square(i,j,:)),'--','Color','#D95319');

        plot(squeeze(AT1,C_Iso(i,j,:)),'-','Color','#7E2F8E');
        plot(squeeze(AT2,C_Iso(i,j,:)),'-','Color','#77AC30');

        xlabel(['C',num2str(i),num2str(j)]);
        ylabel("$\alpha$",'Interpreter','latex');
    end
end
% leg = legend('Circle (Perimeter)', ...
%              'Circle (Area)', ...
%              'Circle (Length)', ...
%              'Analytical (AT1)', ...
%              'Analytical (AT2)','Orientation', 'Vertical');

% leg = legend('Square (Perimeter)', ...
%              'Square (Area)', ...
%              'Square (Length)', ...
%              'Analytical (AT1)', ...
%              'Analytical (AT2)','Orientation', 'Vertical');

leg = legend('Circle (Perimeter)', ...
             'Circle (Area)', ...
             'Circle (Length)', ...
             'Square (Perimeter)', ...
             'Square (Area)', ...
             'Analytical (AT1)', ...
             'Analytical (AT2)','Orientation', 'Vertical');
leg.Layout.Tile = 'east';

function plotFun(C,alpha,op1,op2,op3)
    switch C.order
        case "P1"
            x = alpha;
            y = C.fValues;
            p = plot(x,y);
            p.LineStyle = op1;
            p.Color = op2;
            p.Marker = op3;
        case "P2"
            [x,sortIdx] = sort(C.getCoord());
            y = C.fValues(sortIdx);
            p = plot(x,y);
            p.LineStyle = op1;
            p.Color = op2;
            p.Marker = op3;
    end
end










