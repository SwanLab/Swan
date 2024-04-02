function saveDisplacements(uFun,name)
    title = strcat('./AA_Results/displacements_',name);
    uFun.print(title);
end