function y = sinc(x)
if x==0
	y = 1;
else
	y = sin(pi*x)/(pi*x);
end
