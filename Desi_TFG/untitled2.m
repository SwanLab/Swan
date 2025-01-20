% Especifica la ruta del archivo Excel y el nombre de salida CSV
excelFile = 'C:\Users\desir\Desktop\Machine_Learning_Desi\datos.xlsx';  % Cambia la ruta al archivo Excel
csvFile = 'C:\Users\desir\Desktop\Machine_Learning_Desi\datos5.csv';  % Ruta del archivo CSV de salida


% Lee el archivo Excel como una matriz (sin encabezados)
data = readmatrix(excelFile);  % O usa xlsread(excelFile) si prefieres esta opción

% Si la primera columna debe ser excluida (por ejemplo, es un identificador), elimínala
data = data(:, 2:end);  % Excluir la primera columna

% Escribe los datos en un archivo CSV
writematrix(data, csvFile);

disp(['El archivo Excel se ha convertido a CSV y guardado en: ', csvFile]);