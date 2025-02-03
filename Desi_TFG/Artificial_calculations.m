%% Introduce variables
Waterpl_Area = input('Waterpl_Area(m^2): ');
Watted_Area = input('Watted_Area(m^2): ');
V_Total = input('Ship velocity (m/s): ');
Holtrop = input('Fung power (Kw): ');

%% Calculations
Caa = 0.001 * (Waterpl_Area / Watted_Area);
A_Cf = (105 * (((150 * 10^-6) / 9.74)^(1/3)) - 0.64) * (10^-3);
Rn = (V_Total * 9.74) / (1.05 * 10^-6);
Cfs = 0.075 / ((log10(Rn) - 2)^2);
k = -0.095 + 25.6 * (0.629 / (((9.74 / 3.627)^2) * sqrt(3.627 / 0.74)));
Cts = (1 + k) * Cfs + A_Cf + Caa;
Rts = 0.5 * 1028 * (V_Total^2) * Watted_Area * Cts;
Prop = Rts * V_Total;

% Prop (kw)
Prop = Prop / 1000;
Prop = Prop + Holtrop;

fprintf('Prop (kw): %.4f\n', Prop);

%% Determine rpm
img = imread('Graph2.png');
gray_img = rgb2gray(img);
binary_img = gray_img < 100;
y_input = Prop;
y_pixel = size(img, 1) - round((y_input / 50) * size(img, 1));

for x_pixel = 1:size(img, 2) 
    if binary_img(y_pixel, x_pixel) == 1 
        x_value = 1200 + (x_pixel / size(img, 2)) * (3000 - 1200); 
        disp(['Revolutions (rpm): ', num2str(x_value)]);
        break;
    end
end

%% Determine consumption
img2 = imread('Graph1.png');
gray_img2 = rgb2gray(img2);
binary_img2 = gray_img2 < 100;
x_input = x_value;
x_pixel = round((x_input - 1200) / (3000 - 1200) * size(img2, 2));

for y_pixel2 = size(img2, 1):-1:1
    if binary_img2(y_pixel2, x_pixel) == 1 
        y_value2 = ((size(img2, 1) - y_pixel2) / size(img2, 1)) * 18;
        disp(['Fuel consumption (l/h): ', num2str(y_value2)]);
        break;
    end
end