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
