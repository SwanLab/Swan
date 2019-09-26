function FunctionalProgrammingExample

employees{1} = {'Oriol',1000};
employees{2} = {'Marc',4000};
employees{3} = {'Alex',3000};

increase = 500;

happierEmployees = increaseSalaries(employees,increase);

happierEmployees{1}
happierEmployees{2}
happierEmployees{3}

end

function employee = increaseSalary(employee,increase)
employee{2} = employee{2} + increase;
end

function employees = increaseSalaries(employees,increase)
nEmployees = numel(employees);
for ie = 1:nEmployees
    employees{ie} = increaseSalary(employees{ie},increase);
end
end
