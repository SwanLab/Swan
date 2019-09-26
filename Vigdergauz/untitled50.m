xmin = 0;
xmax = pi/2*0.9;
x = linspace(0,xmax,100);
y = tan(x);
plot(x,y);

theta = 0.01;

mxmy_MxMax_Max = 1/theta; 
mxmy_MxMax_Min = pi/4/theta;

mxmy_MyMax_Max = 4/pi*theta;
mxmy_MyMax_Min = theta;

hold on

plot([xmin xmax],[mxmy_MxMax_Max mxmy_MxMax_Max])
plot([xmin xmax],[mxmy_MxMax_Min mxmy_MxMax_Min])
plot([xmin xmax],[mxmy_MyMax_Max mxmy_MyMax_Max])
plot([xmin xmax],[mxmy_MyMax_Min mxmy_MyMax_Min])