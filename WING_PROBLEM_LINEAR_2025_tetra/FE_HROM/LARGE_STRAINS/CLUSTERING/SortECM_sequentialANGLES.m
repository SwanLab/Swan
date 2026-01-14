function NEW_ORDER_clusters  =  SortECM_sequentialANGLES(BasisFint_cluster)
%  % Criterion based on the principla angles formed by the column
%   % spaces of BasisFint_cluster
% JAHO, 18-MArch-2022
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/CLUSTERING_HROM/08_HOMOG_2D_1T/Test_ECMsequence.mlx
if nargin == 0
    load('tmp.mat')
end


%NEW_ORDER_clusters = zeros(1,length(BasisFint_cluster))  ;
n = length(BasisFint_cluster) ;

T = zeros(n) ;


for i = 1:n
    disp(['icluster =',num2str(i)])
    for j = i+1:n
        
        [UU,SS,VV] = SVDT(BasisFint_cluster{i}'*BasisFint_cluster{j}) ;
        % T_{ij}  \defeq  \dfrac{\sum_{k=1}^{p} \S^{ij}_k}{p}
        T(i,j) = sum(SS)/length(SS) ;
        
    end
end

T = T + T' ;

t = sum(T,1)/n;
figure(784)
hold on
xlabel('Cluster')
ylabel('t: Degree of alignment')
plot(t)

[tPRIME,NEW_ORDER_clusters  ]= sort(t,'descend') ; 
