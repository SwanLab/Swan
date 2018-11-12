function V = compute_initial_radius(paint)
b = 0.5;
V = 1-integral(@(theta) compute_integrant(theta,b),0,2*pi,'RelTol',0,'AbsTol',1e-15);
painting(paint)


end

function painting(paint)
if paint
N = 100;
x = linspace(0,1,N);
y = linspace(0,1,N);

for ix = 1:N
    for jx = 1:N
        r = sqrt((x(ix)-x0)^2 + (y(jx)-y0)^2);
        theta = atan((y(ix)-y0)/(x(ix)-x0));
        f(ix,jx) = phi_function(r,theta,b);
    end
end

thet = [0:0.01:2*pi];
radius = abs(compute_radius(thet,b));

figure(1)
plot(thet*180/pi,radius)
xlabel('$\theta$','Interpreter','LaTex')
ylabel('r','Interpreter','Tex')
title(['$R(\theta)$'],'Interpreter','LaTex')
figure(2)
surf(x,y,f')
xlabel('x','Interpreter','Tex')
ylabel('y','Interpreter','Tex')
zlabel('$\Phi$','Interpreter','LaTex')

figure(3)
contour(x,y,f,[0 0])
grid on
axis equal
end


end




function f = compute_integrant(theta,b)
r = compute_radius(theta,b);
f = r'.^2/2;
end



function radius = compute_radius(theta,b)
options = optimoptions('fsolve','Display','Off','TolX',1e-15,'TolFun',1e-15);
radius = zeros(length(theta),1);
for itheta = 1:length(theta)
    radius(itheta,1) = (fzero(@(r) phi_function(r,theta(itheta),b),0,options));
end
end


function f = phi_function(r,theta,b)

f = ((cos(pi*r*cos(theta))).^2).*((cos(pi*r*sin(theta))).^2)-(b);

end
