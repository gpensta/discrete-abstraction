alpha = .035; epsilon = 0.04; max_u = 6 ; % epsilon size of a target set, alpha : The distance between a point in R^n and the center of a target.
grid_step = 2 * alpha; % Any point in R^2 is within a distance less than alpha from a point in the grid, points of the grid are targets.
xlimi = ceil(.8/grid_step) * grid_step * [-1, 1]; ylimi = ceil(.8/grid_step) * grid_step * [-1, 1]; % To center grid on (0, 0)
[x, y] = meshgrid(xlimi(1):grid_step:xlimi(2), ylimi(1):grid_step:ylimi(2));
c = 0.8; % what is allowed for tracking, what is left is for ref control.
g = -9.81; l = 1.0; m = 1; beta = 1; dt = .1;
f = @(x, u) [x(1) + dt * x(2); x(2) + dt * (-g/l*sin(x(1)) - beta/(m*l^2)*x(2) + u)]; % pendulum state equations

U = zonotope([0; 0], max_u * [1 0; 0 1]); % total control allowed

x0 = [1 * grid_step; 0]; % start point, should be a point of the grid.
target = [0.; 0.]; % should also be a pt of the grid

Z = zonotope(x0, epsilon * [0 1; 1 0]);

[a, b] = linearize(f, [0; 0], [0]);
A = eye(2);
B = [];
for k  = 1:2
    A = A * a;
    B = [a^(k-1) * b B];
end

K = best_k(A, B, c * max_u, epsilon); % Find the value of K that minimizes ||A - BK||, subject to the constraint ||K|| < c * max_u / epsilon.
track = -K*Z + K*x0; % Inputs used for tracking.
ref_zonotope = minkDiff(U, track); % Inputs left used for ref controls
 
Acl = A - B*K; % Closed loop 
disp("x0")
disp(x0);
disp("target")
disp(target)
beta = (((1 - norm(Acl, 'inf')) * epsilon - alpha)) / norm(B, 'inf');
disp("beta > 0");
disp(beta);
disp("norm Acl < 1");
disp(norm(Acl, 'inf'));
% disp("||Acl|| epsilon + ||B||beta < epsilon - alpha");
% disp(norm(Acl, 'inf') * epsilon + norm(B, 'inf') * beta);

% Find if there exists an ur \in ref_zonotope solution to 
%|| A x0 + B ur - target || < alpha.

max_u1 = interval(ref_zonotope).sup(1);
max_u2 = interval(ref_zonotope).sup(2);
u = sdpvar(2, 1);

Constraints = [u(1) <= max_u1, u(1) >= -max_u1, 
               u(2) <= max_u2, u(2) >= -max_u2, 
               norm(A * x0 + B * u - target, 'inf') <= alpha];
           
diagnostics = optimize(Constraints, []);

Bw = k_step_backward(zonotope([0; 0], epsilon * [1 0; 0 1]), U, A, B, 20);

if diagnostics.problem == 0 % Means there is a solution.
    disp('Can reach target in one step.');
    optim_u = value(u);
    disp(optim_u);
    disp("norm(A * x  + B * optim_u - target, 'inf') <= alpha");
    disp(norm(A * x0 + B * optim_u - target, 'inf'));
    F = (A - B*K)* zonotope(x0, epsilon * [0 1; 1 0]) + B * K * x0 + B * zonotope(optim_u, beta * [0 1; 1 0]); % F(T_0, ur)
    Fw = (A - B*K)*Z + B*K*x0  + B * ref_zonotope;  % F(T_0, U_r)
    validXj = A * x0 + B * ref_zonotope + zonotope([0; 0], alpha * [0 1; 1 0]);
    figure; hold on;
    axis equal; 
    xlim(xlimi);
    ylim(ylimi);
    grid on;
    scatter(x(:), y(:), 'k', 'x', 'SizeData', 5);
    center = F.center;
    scatter(center(1), center(2), 'k', 'x', 'SizeData', 40, 'DisplayName', '');
    plot(F, [1,2], 'b', 'linewidth', .5, 'FaceColor', 'none', 'DisplayName', 'F(x0 + ball(epsilon), ur + ball(beta))');
    plot(zonotope(target, epsilon * [0 1; 1 0]), [1,2], 'r', 'linewidth', .5,'FaceColor', 'none', 'DisplayName', 'target + ball(epsilon)');
    %plot(zonotope(target, alpha * [0 1; 1 0]), [1,2], 'r--', 'linewidth', .5,'FaceColor', 'none', 'DisplayName', 'target + ball(alpha)');
    plot(Fw, [1,2], 'b--', 'linewidth', .5,'FaceColor', 'none', 'DisplayName', 'F(x0 + Ball(epsilon), Ur)');
%     plot(Fw.box, [1,2], 'k--', 'linewidth', .5,'FaceColor', 'none', 'DisplayName', '');
    plot(zonotope(x0, epsilon * [0 1; 1 0]), [1,2], 'b', 'linewidth', .5,'FaceColor', 'none', 'DisplayName', '');
    plot(validXj, [1,2], 'g', 'linewidth', .5,'FaceColor', 'none', 'DisplayName', 'F(x0, Ur) + ball(alpha)');
    plot(Bw, [1,2], 'k', 'linewidth', .5,'FaceColor', 'none');

else
    disp("Can't reach target in one step.")
    Fw = (A - B*K)*Z + B*K*x0  + B * ref_zonotope;  % F(T_0, U_r)
    validXj = A * x0 + B * ref_zonotope + zonotope([0; 0], alpha * [0 1; 1 0]);
    figure; hold on;
    axis equal; 
    xlim(xlimi);
    ylim(ylimi);
    grid on;
    scatter(x(:), y(:), 'k', 'x', 'SizeData', 5);
    scatter(x0(1), x0(2), 'k', 'x', 'SizeData', 40, 'DisplayName', '');

    plot(zonotope(target, epsilon * [0 1; 1 0]), [1,2], 'r', 'linewidth', .5,'FaceColor', 'none');
%     plot(zonotope(target, alpha * [0 1; 1 0]), [1,2], 'g--', 'linewidth', .5,'FaceColor', 'none');
%     plot(Fw, [1,2], 'b--', 'linewidth', .5,'FaceColor', 'none');
%     plot(Fw.box, [1,2], 'k--', 'linewidth', .5,'FaceColor', 'none');
    plot(zonotope(x0, epsilon * [0 1; 1 0]), [1,2], 'b', 'linewidth', .5,'FaceColor', 'none');
    plot(validXj, [1,2], 'g', 'linewidth', .5,'FaceColor', 'none');
    plot(Bw, [1,2], 'k', 'linewidth', .5,'FaceColor', 'none');
end


function res = k_step_backward(x, U, A, B, k)
    for i = 1:k
        x = inv(A) *  (x + (-B * U));
    end
    res = x;
end







