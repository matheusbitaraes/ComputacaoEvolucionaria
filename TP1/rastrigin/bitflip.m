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