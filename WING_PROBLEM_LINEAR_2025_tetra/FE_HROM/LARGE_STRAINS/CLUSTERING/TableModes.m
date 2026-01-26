function NameColums = TableModes(TableInfoModesCluster)

if nargin == 0
    load('tmp.mat')
    
end
format short 
%format bank
disp('Writing tex files...')
nclusters = length(TableInfoModesCluster.DISP) ; 
NameRows = cell(1,nclusters) ;
for i=1:nclusters
    NameRows{i} = ['k=',num2str(i)] ;
end
varTypes = cell(size(NameRows)) ;
varTypes(:) = {'int64'} ;
varNames = fieldnames(TableInfoModesCluster)' ;

INFO  = zeros(nclusters,length(varNames)) ; 
for j  =1:length(varNames)
    INFO(:,j) = TableInfoModesCluster.(varNames{j}) ; 
end

T2 = array2table(INFO,'VariableNames',varNames,'RowNames',NameRows);
table2latex(T2,'ModesAndPoints') ;
disp(T2)

