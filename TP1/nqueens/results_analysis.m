% pega todos os arquivos de log
filePattern = fullfile(pwd, 'execution_log_*.csv');
files = dir(filePattern);
num_files = size(files,1);
results_cell = cell(num_files,5);
id1 = 1;
id2 = id1;
id3 = id1;
id4 = id1;
plot_by_populations = Inf*ones(num_files/4, 2, 4);
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
    exec_config = strrep(exec_config,'_crossover_perc-',', cross. %: ');
    results_cell{i,3} = exec_config;
    results_cell{i,4} = best_solutions;
    results_cell{i,5} = num_generations;
    
    if (contains(exec_config, 'pop: 20,'))
        plot_by_populations(id1, :, 1) = [mean(num_generations), mean(best_solutions)];
        id1 = id1 + 1;
    elseif (contains(exec_config, 'pop: 50,'))
        plot_by_populations(id2, :, 2) = [mean(num_generations), mean(best_solutions)];
        id2 = id2 + 1;
    elseif(contains(exec_config, 'pop: 100,'))
        plot_by_populations(id3, :, 3) = [mean(num_generations), mean(best_solutions)];
        id3 = id3 + 1;
    elseif(contains(exec_config, 'pop: 200,'))
        plot_by_populations(id4, :, 4) = [mean(num_generations), mean(best_solutions)];
        id4 = id4 + 1;
    end
    
    % plotting
    scatter(mean(num_generations), mean(best_solutions),'filled'); % media de gerações e melhores soluções
    hold on
end

% plotting information
xlabel('Média do número de gerações');
ylabel('Media do melhor fitness');
legend(results_cell{:,3})
title('Todas as disposições testadas')

% plot por população
figure()
for i = 1:size(plot_by_populations,3)
     % plotting
    scatter(plot_by_populations(:, 1, i), plot_by_populations(:, 2, i),'filled'); % media de gerações e melhores soluções
    hold on
end
xlabel('Média do número de gerações');
ylabel('Media do melhor fitness');
legend('População: 20', 'População: 50', 'População: 100', 'População: 200')
title('Disposições testadas por população')

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

