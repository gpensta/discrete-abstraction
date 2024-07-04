% Epsilon max optim with a fixed alpha 
% Diff lin points
% Pb de temps de calcul, du au calcul de la perturbation, essayer avec le
% calcul en un point. 
clear all
alpha = .1;  % 
grid_step = alpha; % chaque point est au plus � une distance alpha d'un point de la grille.
max_x = .7  ; 
xlimi = ceil(max_x/grid_step) * grid_step * [-1, 1]; ylimi = ceil(max_x/grid_step) * grid_step * [-1, 1]; % To center grid on (0, 0)
delta_max = 1.8; 
[x, y] = meshgrid(xlimi(1):grid_step:xlimi(2), ylimi(1):grid_step:ylimi(2));
beta = 0;

% num_steps = 2; % pour avoir rk(A) = n.


%%
[A, B] = get_A_B_vdp([0.; 0.]); 


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
    [A, B] = get_A_B_vdp(xg);
    F = A * xg + B *zonotope(zeros(size(B, 2), 1), 1 * delta_max * eye(size(B, 2))); 
    for h = 1:length_g
        [i, j] = ind2coord(h, g2);
        xh = [x(i, j); y(i, j)];
        if contains(F, xh)
            I = [I [g; h]];  % transitions g -> h
        end        
   end
end

blck_transitions = [];
omegag = 0; 
etag = 0.;

res = optimizeAndCheck(I, blck_transitions, A, B, delta_max, alpha, beta, x, y, K, t, v, epsilon, epsilon_max, g2, p, n, omegag, etag);
% disp("transitions valides trouv�es");
% disp(res);

G = zeros(length_g);
for i = 1:size(res, 2)
    G(res(1, i), res(2, i)) = 1;
end 

filename = sprintf('/p2p/vdp3');

scc_raw = tarjanSCC(G); % algorithme pour extraire les composantes fortement connect�es. 
scc = [];

for i = 1:length(scc_raw)
    if length(scc_raw{i}) > 3
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
    [A, B] = get_A_B_vdp(xg);

    [ii, jj] = ind2coord( h, g2);
    xh = [x(ii, jj); y(ii, jj)];
    
    kg = K(:, :, g);
    ugh = u_gh(:, i);
    deltagh = delta_gh(i, 1);
    epsilongh = epsilon_gh(i, 1);

    tg = t(:, :, g);
    vg = v(:, :, g);
    
    constraints = [constraints, alpha * .5 <= epsilong];
    constraints = [constraints, epsilong <= epsilon_max];
%     constraints = [constraints, 0 <= epsilongh];

    for r = 1:n % row 
        for c = 1:n % col  || A xg + B u_gh  + (A eps_g - B K) Ball_1 + B Ball_beta + Ball_1 (eta + omega * eps_g) - xh || < eps_h
            constraints = [constraints, tg(r, c) >= 0];
            constraints = [constraints, -tg(r, c) <= A(r, c) * epsilong - B(r, :) * kg(:, c) <= tg(r, c)];
        end
        constraints = [constraints, (A(r, :) * xg  + B(r, :) * ugh - xh(r, 1) + 1 * sum(tg(r, :)) + (etag + omegag * epsilong) <= epsilonh):num2str(i)]; % labelise les contraintes li�es aux transitions g->h
        constraints = [constraints, (A(r, :) * xg  + B(r, :) * ugh - xh(r, 1) - 1 * sum(tg(r, :)) - (etag + omegag * epsilong) >= -epsilonh):num2str(i)]; 
    end
    constraints = [constraints, 0 <= deltagh];
    for r = 1:p % row 
        for c = 1:n % col
           constraints = [constraints, 0 <= vg(r, c)];
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

% R�solution du probl�me
options = sdpsettings('solver', 'linprog', 'verbose', 0); 
sol = optimize(constraints, Objective, options);

% V�rification des r�sultats
if sol.problem == 0
    disp('Solution trouv�e :');
    disp('epsilon_max = ');
    disp(value(epsilon_max));
    res = I;
else
    if (size(I_prev, 2) == size(I, 2)) && (size(I, 2) < 30)
        res = "no solution";
        disp(res);
        return;
    end
    disp('Echec de la r�solution du probl�me');
    disp(yalmiperror(sol.problem));
    new_blck_transitions = []; 
    for i = 1:size(I, 2)
        progress(i);
        violation = check(constraints(num2str(i))); % constraints("1") correspond � la transition g->h = 1;
        if  min(violation) < -1e-5
            new_blck_transitions = [new_blck_transitions, [i; min(violation)]];
        end
    end
    disp("length(new_blck_transitions)");
    disp(length(new_blck_transitions));
    [~, idx] = sort(new_blck_transitions(2, :));
    new_blck_transitions = new_blck_transitions(:, idx);
    new_blck_transitions = new_blck_transitions(1, 1:ceil(size(new_blck_transitions,2)*.25)); % on enleve .3 de la taille des transitions bloquantes.
    
    res = optimizeAndCheck(I, new_blck_transitions, A, B, delta_max, alpha, beta, x, y, K, t, v, epsilon, epsilon_max, g2, p, n, omegag, etag);
end
end

function [i, j] = ind2coord(g, n)
    % Convert unique index g to (i, j) coordinates
    % n is the number of columns in the grid
    i = floor((g - 1) / n) + 1;
    j = g - (i - 1) * n;
end



function g = c2i(i, j, n)
    % Convert (i, j) coordinates to unique index g
    % n is the number of columns in the grid
    g = (i - 1) * n + j;
end

function [A_s, B_s] = get_A_B_vdp(x_star)
    num_steps = 5;
    deltaT = 0.04;
    A = [0, 2; -0.8 - 5*x_star(1)*x_star(2), 2 - 2.5*x_star(1)^2];
%   A = [0, 2; -0.8 - 20*x_star(1)*x_star(2), 2 - 10*x_star(1)^2];
    B = [0; 1];
    % Utilisation de la fonction c2d pour obtenir les matrices discr�tes A_d et B_d
    sys = ss(A, B, eye(size(A)), 0); % Cr�ation d'un syst�me d'�tat-space continu ss(A, B, C, D)
    sys_d = c2d(sys, deltaT, 'zoh'); % Discr�tisation du syst�me
    % Extraction des matrices A_d et B_d du syst�me discr�tis�
    A_d = sys_d.A;
    B_d = sys_d.B;
    A_s = eye(2);
    B_s = [];
        for k  = 1:num_steps
        A_s = A_s * A_d; 
        B_s = [A_d^(k-1) * B_d B_s];
    end
end




