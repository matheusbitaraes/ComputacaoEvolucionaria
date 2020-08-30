% definições gerais
num_exec = 200; %numero de execuções para cada disposição de parametros
dim = 8; % dimensão do tabuleiro (dim x dim)
should_plot = false;
n_generations = 1000; % numero de gerações

% descomentar para multiplas execuções
% Ns = [20, 50, 100, 200]; % tamanhos da população a serem testados
% parent_percs = [0.3, 0.2]; % porcentagens de individuos aleatorios selecionados para selecao de pais
% mutation_percs = [0.1, 0.05]; % porcentagens de chance de mutação
% crossover_percs = [0.5, 0.7, 0.9];% probabilidade de crossover


Ns = [200]; % tamanhos da população a serem testados
parent_percs = [0.3]; % porcentagens de individuos aleatorios selecionados para selecao de pais
mutation_percs = [0.05]; % porcentagens de chance de mutação
crossover_percs = [0.7];% probabilidade de crossover

% iteracoes para as possiveis solucoes
for idn = 1: length(Ns)
    for idc = 1: length(crossover_percs)
        for idp = 1: length(parent_percs)
            for idm = 1: length(mutation_percs)
                N = Ns(idn);
                mutation_perc = mutation_percs(idm);
                crossover_perc = crossover_percs(idc);
                parent_perc = parent_percs(idp);

                execution_log = Inf * ones(num_exec, 2);
                for n = 1:num_exec
                    % criação da populacão inicial
                    population = zeros(N, dim + 1); % vetor N x dim, contendo todos os N individuos + o numero de colisoes na ultima coluna
                    for i = 1:N
                        individual = randperm(dim); % indivíduo da população
                        num_colisions = fitness_nq(individual); % numero de colisões
                        population(i,:) = [individual, num_colisions];
                    end

                    % criação de variáveis para geração de gráficos no final
                    best_fitness = zeros(1, n_generations);
                    average_fitness = zeros(1, n_generations);
                    average_parent_fitness = zeros(1, n_generations);
                    average_children_fitness = zeros(1, n_generations);

                    for i = 1:n_generations

                        sorted_population = sortrows(population, dim + 1);

                        % selecão dos pais e cruzamento
                        children = zeros (N, dim + 1); % alocando memoria para os filhos
                        for j = 1:2:N-1

                            % seleção dos pais: seleciona os 2 melhores pais
                            % dentre uma selecão aleatoria de individuos
                            mixed_population = population(randperm(N),:); % mistura a população aleatoriamente
                            random_individuals = sortrows(mixed_population(1:round(parent_perc * N), :), dim + 1);
                            selected_parents = random_individuals(1:2, :);

                             % realiza o cruzamento de acordo com a probabilidade estipulada
                            if(rand() > crossover_perc)
                                children(j:j+1, 1:dim) =...
                                    CutAndCrossfill_Crossover(selected_parents(:, 1:dim));
                            else
                                children(j:j+1, 1:dim) = selected_parents(:, 1:dim);
                            end

                            % mutação: pega um dos dois filhos e troca os valores de duas posições aleatoriamente
                            if (rand() < mutation_perc)
                                child_index = j + round(rand()); % escolha aleatoria entre um dos filhos
                                indices = randperm(dim,2); % pega os indices a serem trocados

                                % realiza a troca dos indices
                                mutated_child = children(child_index,:);
                                aux = mutated_child(indices(1));
                                mutated_child(indices(1)) = mutated_child(indices(2));
                                mutated_child(indices(2)) = aux;
                                children(child_index,:) = mutated_child;
                            end
                        end

                        % calculo de fitness para os filhos
                        for k = 1:size(children, 1)
                            children(k, dim + 1) =...
                                fitness_nq(children(k, 1:dim)); % calculando os novos fitness
                        end

                        % integra filhos à população e ordena por fitness
                        new_population = sortrows([population; children], dim + 1);

                        % armazena melhor solução
                        best_solution = new_population(1,:);

                        % armazena valores de fitness para gráficos
                        best_fitness(i) = best_solution(dim + 1);
                        average_fitness(i) = mean(sorted_population(:, dim + 1));
                        average_parent_fitness(i) = mean(selected_parents(:, dim + 1));
                        average_children_fitness(i) = mean(children(:, dim + 1));

                        % criterio de parada: se tivermos um fitness minimo (igual à zero),
                        % temos uma solução ótima e podemos parar
                        if (best_solution(dim + 1) == 0)
                            break;
                        end

                        % sobreviventes: ordena pelo fitness e seleciona N melhores individuos
                        population = new_population(1:N,:);
                    end

                    execution_log(n,1) = best_solution(dim + 1); % melhor fitnes
                    execution_log(n,2) = i; % numero de geracoes

                    %plot de uma unica solução
                    if(should_plot)
                        figure()
                        hold on
                        title('Fitness Analysis Over Generations')
                        xlabel('Generations')
                        ylabel('Fitness')
                        lw = 1.5;
                        plot(best_fitness(1:i), 'LineWidth', lw)
                        plot(average_fitness(1:i), 'LineWidth', lw)
                        plot(average_parent_fitness(1:i), 'LineWidth', lw)
                        plot(average_children_fitness(1:i), 'LineWidth', lw)
                        legend('Best Fitness','Average Population Fitness',...
                            'Average Parent Fitness','Average Children Fitness')
                    end
                end

                %             informação sobre a execucao
                filename = sprintf('execution_log_generations-%d_popsize-%d_parent-perc%.2f_mutationperc-%.2f_crossover_perc-%.2f_iterations-%d.csv',...
                    n_generations, N, parent_perc, mutation_perc, crossover_perc, num_exec);
                csvwrite(filename, execution_log)

            end
        end    
    end
end
