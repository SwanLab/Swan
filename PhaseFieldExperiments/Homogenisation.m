% [C_P_Circle,P_Circle] = PFH.computeHomogMaterial("Circle","Perimeter",30);
% [C_A_Circle,A_Circle] = PFH.computeHomogMaterial("Circle","Area",30);
% [C_L_Circle,L_Circle] = PFH.computeHomogMaterial("Circle","Diameter",30);
% 
% [C_P_Square,P_Square] = PFH.computeHomogMaterial("Square","Perimeter",30);
% [C_A_Square,A_Square] = PFH.computeHomogMaterial("Square","Area",30);
% [C_L_Square,L_Square] = PFH.computeHomogMaterial("Square","Diameter",30);
% 
% [C_Iso,AT1] = PFH.computeIsotropicMaterial("AT1",30);
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

        plotFun(C_P_Circle{i,j},P_Circle,"none","#0072BD",".");
        plotFun(C_A_Circle{i,j},A_Circle,"none","#D95319",".");
        plotFun(C_L_Circle{i,j},L_Circle,"none","#EDB120",".");

        plotFun(C_P_Square{i,j},P_Square,"--","#0072BD","none");
        plotFun(C_A_Square{i,j},A_Square,"--","#D95319","none");
        %plotFun(C_L_Square{i,j},L_Square,"--","#EDB120","none");

        plotFun(C_Iso{i,j},AT1,"-","#7E2F8E","none")
        plotFun(C_Iso{i,j},AT2,"-","#77AC30","none")

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
             'Square (Perimeter/Length)', ...
             'Square (Area)', ...
             'Analytical (AT1)', ...
             'Analytical (AT2)','Orientation', 'Vertical');
leg.Layout.Tile = 'east';

function plotFun(C,alpha,op1,op2,op3)
    switch C.order
        case "P1"
            x = alpha;
            y = C.fValues;
            p = plot(y,x);
            p.LineStyle = op1;
            p.Color = op2;
            p.Marker = op3;
        case "P2"
            [x,sortIdx] = sort(C.getCoord());
            y = C.fValues(sortIdx);
            p = plot(y,x);
            p.LineStyle = op1;
            p.Color = op2;
            p.Marker = op3;
    end
end










