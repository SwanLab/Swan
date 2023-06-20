function [GlobalC] = AssemblyGlobalLagrangeCondensed(ACinputs)
    subC = ACinputs.subC;

    boundaries = size(subC,1);

    GlobalC = -subC{1,1};
    for i=2:boundaries
        columnC = [  subC{i-1,2};
                    -subC{i,1}  ];
        GlobalC = [               GlobalC                      sparse(size(GlobalC,1),size(columnC,2));
                     sparse(size(columnC,1),size(GlobalC,2))                   columnC                ];
    end
    columnC = sparse(size(subC{end,2},1), size(GlobalC,2));
    GlobalC = [ GlobalC ;
                columnC ];
end