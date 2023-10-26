clear all;
close all;

path = 'C:\Users\berna\Desktop\TFG\CASES\Circular Column\';
cases = [24;23;25;26;28;27]; %[1;2;5;13;15;16;19;23];
nCases = length(cases);

for iCase=1:nCases
    n = num2str(cases(iCase));
    name = append('C',n,'\C',n,'_monitoring2.fig');
    fig = openfig(append(path,name));
    axObjs = fig.Children;
    nPlots = size(axObjs,1);
    dataObjs = axObjs(nPlots).Children;
    data{iCase,1} = dataObjs.YData;
    switch nPlots
        case 7
            data{iCase,2} = 'MMA';
        case 8
            data{iCase,2} = 'fmincon';
    end
end

figure

% for iCase=1:nCases
%     plot(data{iCase},'DisplayName',append('C',num2str(cases(iCase))))
%     hold on
% end
for iCase=1:nCases
    switch data{iCase,2}
        case 'fmincon'
            plot(data{iCase},'b','DisplayName',append('C',num2str(cases(iCase))))
            hold on
        case 'MMA'
            plot(data{iCase},'r','DisplayName',append('C',num2str(cases(iCase))))
            hold on
    end
end
grid on
grid minor
legend('Location','southeast')
xlabel('nIterations')
