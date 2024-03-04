function K = best_k(A, B, delta, epsilon)

n = size(A, 1);
m = size(B, 2);
k = sdpvar(m, n, 'full');
w = sdpvar(n, n, 'full');
z = sdpvar(m, n, 'full');
y = sdpvar(1);

constraints = [];
for i = 1:n
    for j = 1:n
        constraints = [constraints, w(i, j) >= A(i, j) - sum(B(i, :) * k(:, j))];
        constraints = [constraints, w(i, j) >= -A(i, j) + sum(B(i, :) * k(:, j))];
        constraints = [constraints, sum(w(i, :)) <= y];
    end
end

for l = 1:m
    for j = 1:n
        constraints = [constraints, z(l, j) >= k(l, j)];
        constraints = [constraints, z(l, j) >= -k(l, j)];
        constraints = [constraints, sum(z(l, :)) <= delta / epsilon];
    end
end

objective = y;

ops = sdpsettings('verbose', 1);
diagnostics = optimize(constraints, objective, ops);

if diagnostics.problem == 0
    K = value(k);
%     disp('Optimal value of k:');
%     disp(k_value);
else
    disp('Le problème n’a pas été résolu');
end




