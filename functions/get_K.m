function K = get_K(A, B, epsilon, alpha, delta, c)
    [n, n] = size(A);
    [n, m] = size(B); 
    % Variables de d�cision
    K = sdpvar(m, n, 'full'); % La matrice K � optimiser
    w = sdpvar(n, n, 'full');
    z = sdpvar(m, n, 'full');

    Constraints = [];
    
    for i = 1:n
        for j = 1:n
             Constraints = [Constraints, w(i, j) >= A(i, j)  - sum(B(i, :) * K(:, j))]; % leq ou geq ????
             Constraints = [Constraints, w(i, j) >= -A(i, j) + sum(B(i, :) * K(:, j))];
             Constraints = [Constraints, sum(w(i, :)) <= 1 - alpha / epsilon];  
        end      
    end
    
%     for i = 1:n
%         Constraints = [Constraints, sum(w(i, :)) <= 1 - alpha / epsilon];
%     end
     
    for l = 1:m
        for j = 1:n
            Constraints = [Constraints, z(l, j) >= K(l, j)];
            Constraints = [Constraints, z(l, j) >= -K(l, j)];
            Constraints = [Constraints, sum(z(l, :)) <= delta * c / epsilon];
        end
    end
    

    Objective = 0; % Par exemple, minimiser t
    
    % Options pour le solveur
    options = sdpsettings('solver', 'linprog', 'verbose', 0);
    
    % R�soudre le probl�me de programmation lin�aire
    diagnostics = optimize(Constraints, Objective, options);
    
    if diagnostics.problem == 0
        K = value(K); % Retourner la solution
    else
        K = []; % Pas de solution trouv�e
        disp('Le probl�me n�a pas de solution.');
    end
end

