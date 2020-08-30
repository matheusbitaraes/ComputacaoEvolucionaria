% código para executar multiplas vezes a função rastrigin e plotar a
% distribuição

num_exec = 100; % numero de execuções
executions = zeros(num_exec, 1);

for i = 1:num_exec
    [~, optimal] = matheus_bitaraes_rastrigin(10,10000);
    executions(i) = optimal;
end
interval = 0:100;
mu = mean(executions);
sigma = std(executions);

figure()
% histogram(executions, 10)
hold on
title('Executions distribution')
xlabel('f*');
ylabel('occurency');
norm = fitdist(executions,'normal');
y = pdf(norm, interval);
plot(interval, y)
plot([mu mu],[0 y(round(mu))], 'r')
text(mu + 10, y(round(mu)), sprintf('average = %f', mu))