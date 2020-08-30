function plot_rastrigin(x,y)
% definição da função Rastrigin (é uma função com vários minimos locais)
% a função recebe um x [numero de variaveis x 1] e retorna a avaliação da
% funcão
rastr_func = @(x)(10*size(x,1) + sum(x.^2 - 10*cos(x.*2*pi())));

if (nargin < 1)
    x = -6:0.1:6;
    y = -6:0.1:6;
    z = zeros(1,length(x));
    for i = 1:length(x)
        for j = 1:length(x)
            z(i,j) = rastr_func([x(i);y(j)]);
        end
    end
    
    figure(1)
    contour(x,y,z)
    hold on
    
elseif (nargin == 2)
    figure(1)
    plot(x,y,'-*','LineWidth', 2)
end

    title('Rastrigin Fuction (2 variables)')
    xlabel('X1');
    ylabel('X2');
end