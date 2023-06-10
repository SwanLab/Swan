
close all
compliance = zeros(7,9);
complianceEroded = zeros(7,9);
complianceDilated = zeros(7,9);
differenceDilate = zeros(7,9);
differenceErode = zeros(7,9);

contador = 1;
for i  = 1:9:63
    compliance(contador,:) = matValues(i:(i+8),4);
    complianceDilated(contador,:) = matValues(i:(i+8),5);
    complianceEroded(contador,:) = matValues(i:(i+8),6);
    contador = contador+1;
end 

eta = (0.5-matValues(1:9,1))';
colores_rgb = [
    1 0 0;       % Rojo
    1 0.5 0;     % Naranja
    1 1 0;       % Amarillo
    0 1 0;       % Verde claro
    0 1 1;       % Cian
    0 0.5 1;     % Azul claro
    0 0 1;       % Azul
    0.5 0 0.5;   % Morado
    1 0.5 0.5    % Rosa
];

for i = 1:size(compliance, 1)
    differenceDilate(i, :) = (complianceDilated(i, :)-compliance(i, :))./(compliance(i, :));
    differenceErode(i, :)= (compliance(i, :)- complianceEroded(i,:))./(compliance(i, :));
end 
figure(1);
hold on;
for i = 1:size(compliance, 1)
%   plot(eta,compliance(i, :),'-','color', colores_rgb(i,:), 'LineWidth', 2);
%    plot(eta,complianceEroded(i, :),'--','color', colores_rgb(i,:), 'LineWidth', 2);
%   plot(eta,complianceDilated(i, :),':','color', colores_rgb(i,:), 'LineWidth', 2);
   plot(eta,differenceErode(i, :)*100,':','color', colores_rgb(i,:), 'LineWidth', 2);
   plot(eta,differenceDilate(i, :)*100,'-','color', colores_rgb(i,:), 'LineWidth', 2);
end
xlabel('Eta');
ylabel('Increment[%]');
grid on;
%legend('R1h','R2h','R3h','R4h','R5h','R6h','R7h');

% linearRegresionDilate = polyfit([1;2;3;4;5;6;7], differenceDilate(:,1), 1);
% linearRegresionEroded = polyfit([1;2;3;4;5;6;7], differenceErode(:,1), 1);
% 
% linearRegresionDilatePlot  = ([1;2;3;4;5;6;7]*linearRegresionDilate(1))+linearRegresionDilate(2);
% linearRegresionErodePlot  = ([1;2;3;4;5;6;7]*linearRegresionEroded(1))+linearRegresionEroded(2);
% 
% figure(2)
% hold on;
% plot([1;2;3;4;5;6;7],differenceDilate(:,1),'-','color', [1 0 0], 'LineWidth', 2);
% plot([1;2;3;4;5;6;7],differenceErode(:,1),'-','color', [0.5 0 0.5], 'LineWidth', 2);
% plot([1;2;3;4;5;6;7],linearRegresionDilatePlot(:,1),':','color', [1 0 0], 'LineWidth', 2);
% plot([1;2;3;4;5;6;7],linearRegresionErodePlot(:,1),':','color', [0.5 0 0.5], 'LineWidth', 2);

% close all
% volumen = zeros(7,9);
% contador = 1;
% for i  = 1:9:63
%     volumen(contador,:) = matValues(i:(i+8),7);
%     contador = contador+1;
% end 
% 
% eta = (0.5-matValues(1:9,1))';
% colores_rgb = [
%     1 0 0;       % Rojo
%     1 0.5 0;     % Naranja
%     1 1 0;       % Amarillo
%     0 1 0;       % Verde claro
%     0 1 1;       % Cian
%     0 0.5 1;     % Azul claro
%     0 0 1;       % Azul
%     0.5 0 0.5;   % Morado
%     1 0.5 0.5    % Rosa
% ];
% 
% figure;
% hold on;
% for i = 1:size(volumen, 1)
%     plot(eta,volumen(i, :),'-','color', colores_rgb(i,:), 'LineWidth', 2);
% end
% xlabel('Eta');
% ylabel('Volume');
% grid on;
% legend('R1h','R2h','R3h','R4h','R5h','R6h','R7h');

%matValues(contador,:) = [eta,radio,beta,complianceValue,complianceEroded,complianceDilated,volumenValue];
