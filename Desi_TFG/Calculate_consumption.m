%% Load chart
img = imread('Graph1.png');

% Convert to grayscale
gray_img = rgb2gray(img);

% Set threshold
binary_img = gray_img < 100;

%% Ask for x (propulsion) and calculate y (consumption)
x_input = input('Introduce the propulsion value: ');

% Normalize x
x_pixel = round((x_input - 1200) / (3000 - 1200) * size(img, 2));

% Find Y
for y_pixel = size(img, 1):-1:1
    if binary_img(y_pixel, x_pixel) == 1  
        y_value = ((size(img, 1) - y_pixel) / size(img, 1)) * 18;
        disp(['Prop (w) = ', num2str(x_input)]);
        disp(['Consumption (l/h): ', num2str(y_value)]);
        break;
    end
end

