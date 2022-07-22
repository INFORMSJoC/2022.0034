function RunInstances(folder_name)

file_names = dir(folder_name + "*.mat");

objective_table = cell(length(file_names), 5);
time_table = cell(length(file_names), 5);

for instanceIndex = 1 : length(file_names)
    split_instance = strsplit(file_names(instanceIndex).name, '.');
    instance_tag = split_instance{1};
    
    fprintf("Processing Instance %s\n", instance_tag)
    
    load(strcat(folder_name, file_names(instanceIndex).name));
    
    %% Call Algorithm
    [final_solutions, ~, computation_times, groups] = CoMax(objective, feasible_set);
    
    %% Process Solution Info   
    final_values = zeros(1, size(final_solutions, 2));
    for i = 1 : size(final_solutions, 2)
        final_values(i) = objective.f(final_solutions(:, i));
    end
    
    [best_value, best_index] = max(final_values);
    best_solution = final_solutions(:, best_index);
    
    fprintf('------------------------------------------------\n')
    fprintf('Best Value found: %.3f \n', best_value);
    fprintf('Total time spent: %.3f \n', sum(computation_times))
    fprintf('------------------------------------------------\n')
    
    [objective_table(instanceIndex, :), time_table(instanceIndex, :)] = FillSummaryTable(instance_tag, final_values, groups, computation_times);
end

ObjectiveValues = cell2table(objective_table, 'VariableNames', {'Instance', 'Box', 'Analytic Center', 'Random', 'MVIE'})
ComputationTimes = cell2table(time_table, 'VariableNames', {'Instance', 'Box', 'Analytic Center', 'Random', 'MVIE'})

end

function [ objective_entry, time_entry ] = FillSummaryTable(instance_tag, values, groups, computation_times)
    
objective_entry = cell(1, 5);
time_entry = cell(1, 5);

objective_entry{1} = instance_tag;
time_entry{1} = instance_tag;

best_values = zeros(1, length(groups) - 1);

for i = 2 : length(groups)
    objective_entry{i} = sprintf('%.3f', max(values(groups(i-1)+1 : groups(i))));
    time_entry{i} = sprintf('%.3f', computation_times(i-1));
end

for i = length(groups) + 1 : 5
    objective_entry{i} = '-';
    time_entry{i} = '-';
end

end
