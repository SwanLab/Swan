function plothull(U, domain)
%PLOTHULL Plot weighting functions
%	plothull(U)
%	plothull(U, domain)
%
%	U      - colums are weighting functions
%	domain - x axis (default is 1:size(U,1))

%If multiple weighting functions, plot all
if iscell(U)
	P=length(U);
	for i=1:P
		%If nargin==1, plot only shape of weighting functions
		if nargin==1
			x=1:size(U{i},1);
		%Else plot exactly
		else
			x=linspace(domain(i,1),domain(i,2),size(U{i},1));
		end
		% TODO: optional setting?
		figure
%		h = figure;
%		set(h, 'WindowStyle', 'docked');
		plot(x,U{i},'LineWidth',2);
		axis tight
	end
	return
end

%Plot only given weighting function
if nargin==1
	x=1:size(U,1);
else
	x=linspace(domain(1,1),domain(1,2),size(U,1));
end
figure
plot(x,U,'LineWidth',2);
axis tight
