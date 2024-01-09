function [X_hat] = optim_bss(Graph_cell, x)
    % Number of graphs
    k = numel(Graph_cell);
    % Dimensionality of the signal
    n = numel(x);
    
    % CVX optimization
    cvx_begin
        variable x_k(n, k)
        
        % Objective function (loss)
        loss = 0;
        for i = 1:k
            L = Graph_cell{i}.L;
            loss = loss + x_k(:, i)' * L * x_k(:, i);
        end
        minimize(loss);
        
        % Constraints
        subject to
            x == sum(x_k, 2);
            sum(x_k, 1) == zeros(1, k);
    cvx_end

    X_hat = x_k;
end
