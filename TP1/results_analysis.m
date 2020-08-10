% pega todos os arquivos de log
filePattern = fullfile(pwd, 'execution_log_*.csv');
files = dir(filePattern);
num_files = size(files,1);
results_cell = cell(num_files,5);

figure()
for i = 1:num_files
    data = readtable(files(i).name);
    best_solutions = data.Var1;
    num_generations = data.Var2; 
    results_cell{i,1} = mean(best_solutions);
    results_cell{i,2} = mean(num_generations);
    
    % pega os parâmetros da execução pelo titulo 
    exec_config = files(i).name; 
    exec_config = strrep(exec_config,'execution_log_generations-1000','');
    exec_config = strrep(exec_config,'_iterations-200.csv','');
    exec_config = strrep(exec_config,'_popsize-','pop: ');
    exec_config = strrep(exec_config,'_parent-perc',', par. %: ');
    exec_config = strrep(exec_config,'_mutationperc-',', mut. %: ');
    exec_config = strrep(exec_config,'_crossover_perc',', cross. %: ');
    results_cell{i,3} = num2str(exec_config);
    results_cell{i,4} = best_solutions;
    results_cell{i,5} = num_generations;
    
    
    
    % plotting
    scatter(mean(num_generations), mean(best_solutions), 'filled'); % media de gerações e melhores soluções
    hold on
end

% plotting information
xlabel('Média do número de gerações');
ylabel('Media do melhor fitness');
legend(results_cell{:,3})
title('Todas as disposições testadas')


sorted_results = sortrows(results_cell, [1 2]);
num_results = 10;
best_results = sorted_results(1:num_results,:);

display(best_results);

% plot do top 10 de soluções (pelo menor numero de gerações)
% figure()
% for i = 1:num_results
%     scatter(best_results{i,1}, best_results{i,2},'filled')
%     hold on
% end
% xlabel('Média do número de gerações');
% ylabel('Media do melhor fitness');
% legend(best_results{:,3})
% title('Todas as disposições testadas')

