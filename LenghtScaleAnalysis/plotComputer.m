
compiance = zeros(7,7);
contador = 1;
for i  = 1:7:49
    compiance(contador,:) = matValues(i:(i+6),4);
    contador = contador+1;
end 
eta = matValues(1:7,1)';
colores_rgb = [
    1 0 0;       % Rojo
    1 0.5 0;     % Naranja
    1 1 0;       % Amarillo
    0 1 0;       % Verde claro
    0 1 1;       % Cian
    0 0.5 1;     % Azul claro
    0 0 1        % Azul
];


figure;
hold on;
for i = 1:size(compiance, 2)
    plot(eta,compiance(i, :), 'color', colores_rgb(i,:), 'LineWidth', 2);
end
xlabel('Compiance');
ylabel('Eta');
title('Compiance value respect different eta values for different radius of influence');
grid on;
legend('R1h','R2h','R3h','R4h','R5h','R6h','R7h');
%[0.28 0.57 0.86 1.1547 1.4434 1.7321 2.0207]
% Mostrar la figura