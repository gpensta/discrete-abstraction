% Epsilon max optim with a fixed alpha 
% Diff lin points
% Pb de temps de calcul, du au calcul de la perturbation
clear all
alpha = .1;  % 
grid_step = alpha; % chaque point est au plus à une distance alpha d'un point de la grille.
max_x = .8; 
xlimi = ceil(max_x/grid_step) * grid_step * [-1, 1]; ylimi = ceil(max_x/grid_step) * grid_step * [-1, 1]; % To center grid on (0, 0)
delta_max = 10; 
U_tot = zonotope([0; 0], delta_max * [1 0; 0 1]); % total control allowed
[x, y] = meshgrid(xlimi(1):grid_step:xlimi(2), ylimi(1):grid_step:ylimi(2));
beta = 0;

num_steps = 2; % pour avoir rk(A) = n.


%%


% xi = [0; 0]; % Vanderpol 
% u_star = 0;  
% [A_d, B_d] = get_A_B_vdp(xi, u_star, dt);
% A = eye(2);
% B = [];
% for k  = 1:num_steps
%     A = A * A_d;
%     B = [A_d^(k-1) * B_d B];
% end

% pendulum
[A, B] = get_A_B_pendulum([0; 0]);


g1 = size(x, 1); % taille de la grille -> g1 x g2
g2 = size(x, 2); % 
length_g = g1 * g2;

p = size(B, 2);
n = size(A, 1);

epsilon_max = sdpvar(1, 1); % 

t = sdpvar(n, n, length_g, 'full'); % var intermediaire 
v = sdpvar(p, n, length_g, 'full'); % -v < k < v
epsilon = sdpvar(length_g, 1, 'full');
K = sdpvar(p, n, length_g, 'full');

I = [];

% Quelles sont les transitions g -> h possibles ? 
for g = 1:length_g
    [i, j] = ind2coord(g, g2); % coord dans la grille pour le pt g. 
    xg = [x(i, j); y(i, j)];
    [A, B] = get_A_B_pendulum(xg);
    F = A * xg + B *zonotope([0; 0], .6 * delta_max *  [1 0; 0 1]);  %   + zonotope([0; 0], alpha * [1 0; 0 1]); 
    for h = 1:length_g
        [i, j] = ind2coord(h, g2);
        xh = [x(i, j); y(i, j)];
        if contains(F, xh)
            I = [I [g; h]];  % transitions g -> h
        end        
   end
end

blck_transitions = [];
[omegag, etag, ~, ~] = linearization_parameters([xlimi(2); ylimi(2)], alpha, delta_max); % omega : coeffdir, etag : ordonnée origine. 

res = optimizeAndCheck(I, blck_transitions, A, B, delta_max, alpha, beta, x, y, K, t, v, epsilon, epsilon_max, g2, p, n, omegag, etag);
% disp("transitions valides trouvées");
% disp(res);

G = zeros(length_g);
for i = 1:size(res, 2)
    G(res(1, i), res(2, i)) = 1;
end 

filename = sprintf('/p2p/pendulum4');

scc_raw = tarjanSCC(G); % algorithme pour extraire les composantes fortement connectées. 
scc = [];

for i = 1:length(scc_raw)
    if length(scc_raw{i}) > 1
        scc = [scc; {scc_raw{i}}];
    end
end

sccGraphviz(G, grid_step, xlimi(1), ylimi(1), scc, filename);


function res = optimizeAndCheck(I, blck_transitions, A, B, delta_max, alpha, beta, x, y, K, t, v, epsilon, epsilon_max, g2, p, n, omegag, etag)
    I_prev = I;
    I(:, blck_transitions) = [];
    num_transitions = size(I, 2);
    fprintf("\n");
    disp("Nombre de transition");
    disp(num_transitions);
    
    delta_gh = sdpvar(num_transitions, 1, 'full'); % chaque transition a son delta. 
    epsilon_gh = sdpvar(num_transitions, 1, 'full');
    u_gh = sdpvar(p, num_transitions, 'full');

    constraints = [];

    for i = 1:num_transitions % parcours des transition
    progress(i); % disp progression
    g = I(1, i);
    h = I(2, i);
    epsilong = epsilon(g, 1); 
    epsilonh = epsilon(h, 1);
    [ii, jj] = ind2coord( g, g2);
    xg = [x(ii, jj); y(ii, jj)];
    [A, B] = get_A_B_pendulum(xg);
%     [omegag, etag, A, B] = linearization_parameters(xg, alpha, delta_max); % omega : coeffdir, etag : ordonnée origine. 
    [ii, jj] = ind2coord( h, g2);
    xh = [x(ii, jj); y(ii, jj)];
    
    kg = K(:, :, g);
    ugh = u_gh(:, i);
    deltagh = delta_gh(i, 1);
    epsilongh = epsilon_gh(i, 1);

    tg = t(:, :, g);
    vg = v(:, :, g);
    
    constraints = [constraints, .5 * alpha <= epsilong];
    constraints = [constraints, epsilong <= epsilon_max];
    constraints = [constraints, .5 * alpha <= epsilonh];
    constraints = [constraints, epsilonh <= epsilon_max];
%     constraints = [constraints, 0 <= epsilongh];

    for r = 1:n % row 
        for c = 1:n % col  || A xg + B u_gh  + (A eps_g - B K) Ball_1 + B Ball_beta + Ball_1 (eta + omega * eps_g) - xh || < eps_h
            constraints = [constraints, tg(r, c) >= 0];
            constraints = [constraints, -tg(r, c) <= A(r, c) * epsilong - B(r, :) * kg(:, c) <= tg(r, c)];
        end
        constraints = [constraints, (A(r, :) * xg  + B(r, :) * ugh - xh(r, 1) + 1 * sum(tg(r, :)) + (etag + omegag * epsilong) <= epsilonh):num2str(i)]; % labelise les contraintes liées aux transitions g->h
        constraints = [constraints, (A(r, :) * xg  + B(r, :) * ugh - xh(r, 1) - 1 * sum(tg(r, :)) - (etag + omegag * epsilong) >= -epsilonh):num2str(i)]; 
    end
    constraints = [constraints, 0 <= deltagh];
    for r = 1:p % row 
        for c = 1:n % col
           constraints = [constraints, 0 <= vg( r, c)];
           constraints = [constraints, -vg(r, c) <= kg(r, c) <= vg(r, c)];
        end
       constraints = [constraints,  (sum(vg(r, :)) + ugh(r, 1) <= deltagh - beta):num2str(i)];
       constraints = [constraints,  (-sum(vg(r, :)) + ugh(r, 1) >= -deltagh + beta):num2str(i)];    
    end
    constraints = [constraints, deltagh <= delta_max];
end
fprintf("\n");
disp("size(constraints)");
disp(size(constraints));
% Fonction objectif

Objective = epsilon_max; % minimiser epsilon max pour avoir une abstraction fine. 

% Résolution du problème
options = sdpsettings('solver', 'linprog', 'verbose', 0); 
sol = optimize(constraints, Objective, options);

% Vérification des résultats
if sol.problem == 0
    disp('Solution trouvée :');
    disp('epsilon_max = ');
    disp(value(epsilon_max));
    res = I;
else
    if (size(I_prev, 2) == size(I, 2)) && (size(I, 2) < 30)
        res = "no solution";
        disp(res);
        return;
    end
    disp('Echec de la résolution du problème');
    disp(yalmiperror(sol.problem));
    new_blck_transitions = []; 
    for i = 1:size(I, 2)
        progress(i);
        violation = check(constraints(num2str(i))); % constraints("1") correspond à la transition g->h = 1;
        if  min(violation) < -1e-5
            new_blck_transitions = [new_blck_transitions, [i; min(violation)]];
        end
    end
    disp("length(new_blck_transitions)");
    disp(length(new_blck_transitions));
    [~, idx] = sort(new_blck_transitions(2, :));
    new_blck_transitions = new_blck_transitions(:, idx);
    new_blck_transitions = new_blck_transitions(1, 1:ceil(size(new_blck_transitions,2)*.5)); % on enleve .3 de la taille des transitions bloquantes.
    
    res = optimizeAndCheck(I, new_blck_transitions, A, B, delta_max, alpha, beta, x, y, K, t, v, epsilon, epsilon_max, g2, p, n, omegag, etag);
end
end

function [i, j] = ind2coord(g, n)
    % Convert unique index g to (i, j) coordinates
    % n is the number of columns in the grid
    i = floor((g - 1) / n) + 1;
    j = g - (i - 1) * n;
end

function [A, B] = get_A_B_pendulum(x_star)
    g = -9.81; l = 1.0; m = 1;  dt = .1;
    f = @(x, u) [x(1) + dt * x(2); x(2) + dt * (-g/l*sin(x(1)) - 1/(m*l^2)*x(2) + u)]; % pendulum state equations
    [a, b] = linearize(f, x_star, [0]);
            A = eye(2);
            B = [];
            for k  = 1:2
                A = A * a;
                B = [a^(k-1) * b B];
            end  
end


function g = c2i(i, j, n)
    % Convert (i, j) coordinates to unique index g
    % n is the number of columns in the grid
    g = (i - 1) * n + j;
end

function [coeff_dir_w, offset_w, A, B] = linearization_parameters(xg, alpha, delta)
    g = -9.81;
    l = 1;
    beta = 1; 
    m = 1; 
    dt = .1;
    num_steps = 2;

    f = @(x, u) [x(1) + dt * x(2); 
                 x(2) + dt * (-g/l*sin(x(1)) - beta/(m*l^2)*x(2) + u)];

    fn = @(x, u_list) apply_fn_iteratively(f, x, u_list);

    % Define range of epsilon
    epsilons = alpha:alpha:3*alpha;
    n = size(xg, 1); 
    W_sups = zeros(size(epsilons));

    for i = 1:length(epsilons)
        epsilon = epsilons(i);
        S = zonotope(xg, epsilon * eye(n));
        U = zonotope(zeros(num_steps, 1), delta * eye(num_steps));
        
        [A, B, W] = get_multistep_w(fn, S, U);
        W_sups(i) = max(max(abs(W.sup), abs(W.inf)));
    end
    
    p = polyfit(epsilons, W_sups, 1);
    affine_fit = polyval(p, epsilons);

    % Ensure the affine function is always above sup W by adding a margin
    margin = max(W_sups - affine_fit);
    p(2) = p(2) + margin;  % Adjust the intercept

    % Extract the coefficients
    coeff_dir_w = p(1);  % Slope
    offset_w = p(2);     % Intercept
end






