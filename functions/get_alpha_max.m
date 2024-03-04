function [alpha_max, epsilon_opt, K] = get_alpha_max(c, A, B, delta)
    n = size(A, 1);
    m = size(B, 2);
    k = sdpvar(m, n, 'full');
    w = sdpvar(n, n, 'full');
    z = sdpvar(m, n, 'full');
    y = sdpvar(1);
    epsilon = sdpvar(1);

    constraints = [];
    for i = 1:n
        for j = 1:n
            constraints = [constraints, w(i, j) >= A(i, j) * epsilon - sum(B(i, :) * k(:, j) * c)];
            constraints = [constraints, w(i, j) >= -A(i, j) * epsilon + sum(B(i, :) * k(:, j) * c)];
            constraints = [constraints, -epsilon + sum(w(i, :)) <= y];
        end
    end

    for l = 1:m
        for j = 1:n
            constraints = [constraints, z(l, j) >= k(l, j)];
            constraints = [constraints, z(l, j) >= -k(l, j)];
            constraints = [constraints, sum(z(l, :)) <= delta];
        end
    end

    objective = y;

    ops = sdpsettings('verbose', 1);
    diagnostics = optimize(constraints, objective, ops);

    if diagnostics.problem == 0
        Kprime = value(k);
        epsilon_opt = value(epsilon);
        K = Kprime / epsilon_opt * c; % Kprime = K * epsilon / c
        alpha_max = epsilon_opt * (1 - norm(A - B * K, 'inf'));
    else
        disp('Le problème n’a pas été résolu');
    end
end