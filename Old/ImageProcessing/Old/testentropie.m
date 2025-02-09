%% testentropie


u0 = imresize(double(imread('data/einstein.png')),.1)/255; %% 300x400
%u0 = double(u0);
imagesc(u0);

[M N K] = size(u0);
MN = M*N;
I=reshape([1:M*N],M,N);
east=[I(:,2:end), I(:,end)];
north=[I(2:end,:); I(end,:)];
D1 = sparse(I,east,1,MN,MN) -speye(MN,MN);
D2 = sparse(I,north,1,MN,MN) -speye(MN,MN);
%D = [D1 ; D2];


DK1=D1;
DK2=D2;
for i=2:K
    DK1 = blkdiag(DK1,D1);
    DK2 = blkdiag(DK2,D2);
end
DK = [DK1; DK2];

phi=zeros(M,N,K); %% useless
eta=zeros(MN*2*K,1);

delta=.1;
gamma = 1;

tau = delta/(4*gamma^2);
for i=1:1000
    phi = reshape(DK'*eta,M,N,K);
    phi = exp(-(50*u0+gamma*phi)/delta);
    phi = bsxfun(@rdivide,phi,sum(phi,3));
    eta = eta + tau*gamma*DK*phi(:);
    normeta= max(1,hypot(eta(1:MN*K),eta(MN*K+1,end)));
    eta = eta./[normeta;normeta];
    if mod(i,20)==0
        imagesc(phi); drawnow();
    end   
end

