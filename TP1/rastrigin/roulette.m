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
