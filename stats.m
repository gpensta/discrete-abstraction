% Supposons que G est votre matrice d'adjacence

% Nombre de transitions
num_transitions = nnz(G);

outdegree = sum(G, 2);

min_outdegree = min(outdegree);
max_outdegree = max(outdegree);
avg_outdegree = mean(outdegree);


fprintf('Nombre de transitions : %d\n', num_transitions);
fprintf('Degrés de sortie minimum : %d\n', min_outdegree);
fprintf('Degrés de sortie maximum : %d\n', max_outdegree);
fprintf('Degrés de sortie moyen : %.2f\n', avg_outdegree);


eps = value(epsilon);
k_val = value(K);

% myw = omegag * alpha + etag;

max_ABK_inf_norm = -inf;
min_ABK_inf_norm = inf;
nb_states = size(G, 2);
ABK = ones(nb_states, 1); 

for g = 1:nb_states
    [i, j] = ind2coord(g, size(x, 2)); % coord dans la grille pour le pt g. 
    xg = [x(i, j); y(i, j)];
%     [A, B] = get_A_B_pendulum(xg);
    [A, B] = get_A_B_vdp(xg);
    eps_g = eps(g);
    kg_prime = k_val(:, :, g);
    kg = kg_prime / eps_g;
    ABK_inf_norm = norm(A - B * kg, inf);
    ABK(g) = norm(A - B * kg, inf);
    max_ABK_inf_norm = max(max_ABK_inf_norm, ABK_inf_norm);
    min_ABK_inf_norm = min(min_ABK_inf_norm, ABK_inf_norm);
end

disp("min(eps)");   
disp(min(eps));
disp("max(eps)");
disp(max(eps));
% disp("nb_states")
% disp(nb_states);
disp("max_ABK_inf_norm");
disp(max_ABK_inf_norm);
disp("min_ABK_inf_norm");
disp(min_ABK_inf_norm);
avgABK = mean(ABK(~isnan(ABK)));
disp("avgABK");
disp(avgABK);

% for i = 1:835
%    epss =  ABK(I(1, i)) * alpha + norm(B, inf) * beta + myw;
%    if epss > .05
%       disp(i); 
%       break;
%    end
% end


g = I(1, i);

[i, j] = ind2coord(g, size(x, 2)); % coord dans la grille pour le pt g. 
xg = [x(i, j); y(i, j)];

% disp("xg");
% disp(xg);



function [A_s, B_s] = get_A_B_vdp(x_star)
    num_steps = 5;
    deltaT = 0.04;
    A = [0, 2; -0.8 - 5*x_star(1)*x_star(2), 2 - 2.5*x_star(1)^2];
%   A = [0, 2; -0.8 - 20*x_star(1)*x_star(2), 2 - 10*x_star(1)^2];
    B = [0; 1];
    % Utilisation de la fonction c2d pour obtenir les matrices discrï¿½tes A_d et B_d
    sys = ss(A, B, eye(size(A)), 0); % Crï¿½ation d'un systï¿½me d'ï¿½tat-space continu ss(A, B, C, D)
    sys_d = c2d(sys, deltaT, 'zoh'); % Discrï¿½tisation du systï¿½me
    % Extraction des matrices A_d et B_d du systï¿½me discrï¿½tisï¿½
    A_d = sys_d.A;
    B_d = sys_d.B;
    A_s = eye(2);
    B_s = [];
    for k  = 1:num_steps
        A_s = A_s * A_d; 
        B_s = [A_d^(k-1) * B_d B_s];
    end
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

