% definições gerais
dim = 8; % dimensão do tabuleiro (dim x dim)
N = 20; % tamanho da população
parent_perc = 0.3; % porcentagem de pais selecionados para cruzamento
n_generations = 100; % numero de gerações
mutation_perc = 0.1; % porcentagem de chance de mutação

% criação da populacão inicial
population = zeros(N, dim + 1); % vetor N x dim, contendo todos os N individuos + o numero de colisoes na ultima coluna
for i = 1:N
    individual = randperm(dim); % indivíduo da população
    num_colisions = fitness_nq(individual); % numero de colisões
    population(i,:) = [individual, num_colisions];
end

% criação de variáveis para geração de gráficos no final
best_fitness = zeros(1,n_generations);
average_fitness = zeros(1,n_generations);
average_children_fitness = zeros(1,n_generations);

for i = 1:n_generations
    
    % selecão dos pais (os 40% com melhor fitness irão para o cruzamento)
    sorted_population = sortrows(population, dim + 1);
    num_parents = round(N*parent_perc);
    selected_parents = sorted_population(1:num_parents,:);
    
    % cruzamento
    children = zeros (num_parents, dim + 1); % alocando memoria para os filhos
    for j = 1:2:num_parents-1
        children(j:j+1,1:dim) = CutAndCrossfill_Crossover(selected_parents(j:j+1,1:dim));
        
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
        children(k, dim + 1) = fitness_nq(children(k, 1:dim)); % calculando os novos fitness
    end

    % integra filhos à população e ordena por fitness
    new_population = sortrows([population; children], dim + 1);
    
    % armazena melhor solução
    best_solution = new_population(1,:);
    
    % armazena valores de fitness para gráficos
    best_fitness(i) = best_solution(dim + 1);
    average_fitness(i) = mean(new_population(:, dim + 1));
    average_children_fitness(i) = mean(children(:, dim + 1));
    
    % criterio de parada: se tivermos um fitness minimo (igual à zero),
    % temos uma solução ótima e podemos parar
    if (best_solution(dim + 1) == 0)
        break; 
    end
    
    % sobreviventes: ordena pelo fitness e seleciona N melhores individuos
    population = new_population(1:N,:);
end
plot(best_fitness(1:i))
hold on
plot(average_fitness(1:i))
plot(average_children_fitness(1:i))
legend('Best Fitness','Average Population Fitness','Average Children Fitness')

best_solution
i

% plotar melhor solução