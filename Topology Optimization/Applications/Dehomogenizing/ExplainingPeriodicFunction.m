function ExplainingPeriodicFunction
x   = linspace(-10,10,1000);
N   = 5;
eps = (max(x)-min(x))/N;
m   = linspace(0.01,0.99,20);
types = {'FLOOR','COS'};
for it = 1:length(types)
    for im = 1:length(m)
        mOr(im,it) = computeRealInclusionSize(x,N,eps,m(im),types{it});
    end
end
figure(100)
plot(m,m,'-+');
hold on
plot(m,mOr,'-+')
axis([0 1 0 1])
end


function mOr = computeRealInclusionSize(x,N,eps,m,type)


y = (x-min(x)-eps)/eps;%-eps/4;
fp = peridicFunction(y,type);
r  = (fp-0.5)/(m*0.5) - 1;

xi = zeros(size(r));
xi(r>0) = 1;


[(r>0)',x']

index = find((r(1:end-1).*r(2:end))'<0);
index = [1;index;length(x)];

for i = 1:length(index)-1

    xa = x(index(i));
    xb = x(index(i+1));
    mm(i) = xb - xa;
    isMat(i) = xi(index(i+1))>0;
end


si  = (reshape(mm(2:end),2,[])');
isM = (reshape(isMat(2:end),2,[])');
[N,size(si,1)]

mOr = mode(si(~isM))/eps;
figure()
plotTopology(x,index)
figure(300)
plot(x,xi,'+-')
end

function fp = peridicFunction(y,type)

switch type
    case 'COS'
        T = pi;
        fp = abs(cos(T*(y))).^2;
    case 'FLOOR'
        fp = y - floor(y);
end
end

function plotTopology(x,index)

for i = 2:length(index)-1
 xmin = x(index(i));
 xmax = x(index(i+1));
 m = createMesh(xmin,xmax);
 if mod(i,2) == 0 
    m.plotSolidColor('black')
 else
    m.plotSolidColor('white') 
 end
end
close all
end

function m = createMesh(xmin,xmax)
N = 2;
ymin = 0;
ymax = 0.05;
xv = linspace(xmin,xmax,N);
yv = linspace(ymin,ymax,N);
            [X,Y] = meshgrid(xv,yv);
            s.coord(:,1) = X(:);
            s.coord(:,2) = Y(:);
              [F,V] = mesh2tri(X,Y,zeros(size(X)),'f');
              s.coord  = V(:,1:2);
              s.connec = F;

%            s.connec = delaunay(s.coord);
            m = Mesh.create(s);


end