function [x, f] = matheus_bitaraes_rastrigin(nvar, ncal)
%funcao da segunda tarefa do TP2
%   [x, f] = matheus_bitaraes_rastrigin(nvar, ncal)
%   x: vetor das variáveis de decisão do melhor indivíduo
%   f: falor da função em x
%   nvar: numero de variáveis de decisão
%   ncal: numero total de chamadas da funcao de calculo de fitness do
%   problema

%% definição de parâmetros e hiper-parâmetros
pop_size = 100;
num_parents = pop_size; % numero de pais a serem selecionados igual ao numero da populaçao
x_inf = -5.25; % limite inferior de x
x_sup = 5.25; % limite superior de x
li = 10; % numero de bits para representar os valores
k = round(pop_size * 0.4); % quantidade de individuos no torneio
pc = 0.8; % probabilidade de cruzamento por individuo
pm = 0.2; % probabilidade de mutação de cada indivíduo
pb = 0.5; % probabilidade de bitflip de cada binario de um individuo a ser mutado

%% Implementação
fit_col = li + 2; % coluna da matriz população onde serão armazenados os valores de fitness
realn_col = li + 1; % coluna da matriz população onde serão armazenados os numero reais

% proteção caso a função seja chamada sem argumentos
if (nargin < 2)
   nvar = 10;
   ncal = 10000;
end

% definição da função Rastrigin (é uma função com vários minimos locais)
% a função recebe um x [numero de variaveis x 1] e retorna a avaliação da
% funcão
rastr_func = @(x)(10*size(x,1) + sum(x.^2 - 10*cos(x.*2*pi())));

% função de decodificação binario para decimal(referencia: material
% disponibilizado no moodle)
bin_decode = @(x_bin) (...
    x_inf + (x_sup - x_inf) .* ...
    sum(x_bin.*repmat(pow2(li-1:-1:0), [size(x_bin,1), 1, size(x_bin,3)]), 2) ./...
    (pow2(li)-1));

% inicializa a população com valores binários aleatórios. A população é de
% dimensão [ numero de variaveis x numero de bits binários + 2 x tamanho da populacao ] 
pop = round(rand(nvar, li, pop_size));


% coloca população decodificada para numeros reais na coluna li + 1
pop(:, realn_col, :) = bin_decode(pop(:,1:li,:));

% cálculo inicial do fitness
num_eval = 0; 
for i = 1:pop_size
    % atualiza o numero de chamadas à função rastrigin
    num_eval = num_eval + 1;

    % coloca o fitness calculado na coluna li + 2 (fit_col) e
    % atribui valor negativo, para que o problema se torne de
    % maximização
    pop(:, fit_col, i) = -rastr_func(pop(:,realn_col,i)); 
end
    
% inicio do loop geracional. O modelo de população será o modelo
% geracional, ou seja, toda populacao é substituída pelos filhos gerados
while num_eval <= ncal
    
    % seleção dos pais, cruzamento e mutação
    children = zeros(nvar, li + 2, num_parents);
    for i = 1:2:num_parents
        
        % seleção dos pai. Torneio ou roleta, com igual probabilidade
        if (rand() > 0.5)
            selected_parents = tournament(pop, fit_col, 2, k);
        else
            selected_parents = roulette(pop, fit_col, 2);
        end
        
        % cruzamento: 1 ponto de corte por variável
        if (pc < rand())
            for idv = 1:nvar
                cut_point = round(li*rand()); % ponto de crossover da variável
                children(idv,1:li,i) = [selected_parents(idv,1:cut_point,1)...
                    selected_parents(idv,cut_point + 1:li,2)];
                children(idv,1:li,i+1) = [selected_parents(idv,1:cut_point,2)...
                    selected_parents(idv,cut_point + 1:li,1)];
            end
        else
            children(:,:,i:i+1) = selected_parents;
        end
        
            
        % mutação bit flip
        for idm = i:i+1
            if (pm < rand())
                children(:,:,idm) = bitflip(children(:,:,idm), pb);
            end
        end
    end
    
    % coloca população decodificada para numeros reais na coluna li + 1
    children(:, realn_col, :) = bin_decode(children(:,1:li,:));

    % calculo do fitness dos filhos
    for i = 1:num_parents
        
        % atualiza o numero de chamadas à função rastrigin
        num_eval = num_eval + 1;

        % coloca o fitness calculado na coluna li + 2 (fit_col) e
        % atribui valor negativo, para que o problema se torne de
        % maximização
        children(:, fit_col, i) = -rastr_func(children(:,realn_col,i)); 
    end
    
    % seleção: ordena pelo fitness e reduz a população ao tamanho original
    over_pop = cat(3, pop, children);
    fits(:,1) = over_pop(1, fit_col, :); % pega os fitness dos participantes
    fits(:,2) = 1:size(over_pop,3); % guarda indices antes da ordenação
    sorted_fits = sortrows(fits, 1, 'descend');
    pop = over_pop(:, :, sorted_fits(1:pop_size,2));
end

% fim do loop geracional
x = pop(:, realn_col, 1);
f = -pop(1, fit_col, 1);
end

function parents = tournament(pop, fit_col, num_parents, k)

% torneio
pop_size = size(pop, 3);

% ordena aleatoriamente os individuos e pega os k primeiros
participants = pop(:, :, randperm(pop_size, k));
fits(:,1) = participants(1, fit_col, :); % pega os fitness dos participantes
fits(:,2) = 1:k; % guarda indices antes da ordenação
sorted_fits = sortrows(fits, 1, 'descend');
parents = participants(:, :, sorted_fits(1:num_parents,2));

end

function parents = roulette(pop, fit_col, num_parents)

parents = pop(:, :, 1:num_parents) * 0;

% roleta
fits(:,1) = pop(1, fit_col, :);
positive_fits = abs(sum(fits)) + fits;
roleta = cumsum(positive_fits ./ sum(positive_fits)); 
for id = 1:num_parents
    resultado_roleta = rand();
    idx = 1;
    while roleta(idx) < resultado_roleta
        idx = idx + 1;
    end
    parents(:, :, id) = pop(:, :, idx);
end

end

function mutated = bitflip(children, pb)
mutated = children * 0;

% itera nas variaveis
nvar = size(children, 1);
li = size(children, 2) - 2;
for i = 1:nvar
    var_mutation = (rand(1, li) <= pb);
    mutated(i, 1:li) = abs(children(i, 1:li) - var_mutation);
end
end