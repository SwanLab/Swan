function stX=stdata(X)
%把样本数据集标准化、规范化 具有均值0标准差1
stX=zeros(size(X));%初始化
[m,n]=size(X);
for i=1:n
    mX=mean(X(:,i));%样本每个特征（每一列）的均值
    sX=std(X(:,i));%样本每个特征（每一列）的方差
    for j=1:m
        stX(j,i)=(X(j,i)-mX)/sX;%标准化、规范化
    end
end
end