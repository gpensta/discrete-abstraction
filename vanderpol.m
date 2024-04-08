clear 
close all
addpath('./CORA_2022')
addpath('./functions')

% Vanderpol oscillator 
% Dimensions of state domain : (x1, x2) \in [-1, 1]^2, inputs \in [-1, 1]. 

%% *************************** Dynamics ***********************************

% f_u =  @(t,x,u)(-[ -2*x(2,:) ; 0.8*x(1,:) + 10*x(1,:).^2.*x(2,:) - 2*x(2,:) + u] );
% n = 2;
% m = 1; % number of control inputs

%% ************************** Discretization ******************************

% deltaT = 0.01;
% %Runge-Kutta 4
% k1 = @(t,x,u) (  f_u(t,x,u) );
% k2 = @(t,x,u) ( f_u(t,x + k1(t,x,u)*deltaT/2,u) );
% k3 = @(t,x,u) ( f_u(t,x + k2(t,x,u)*deltaT/2,u) );
% k4 = @(t,x,u) ( f_u(t,x + k1(t,x,u)*deltaT,u) );
% f_ud = @(t,x,u) ( x + (deltaT/6) * ( k1(t,x,u) + 2*k2(t,x,u) + 2*k3(t,x,u) + k4(t,x,u)  )   );

%% ************************** Linearization *******************************

% Définition des points d'équilibre x_star et u_star
% x_star = [.7; .7]; % Comme dans le papier
% u_star = 0;  
% 
% % Matrices A et B linéarisées évaluées au point d'équilibre
% A = [0, 2; -0.8 - 20*x_star(1)*x_star(2), 2 - 10*x_star(1)^2];
% B = [0; 1];
% 
% % Fonction dynamique linéarisée continue
% f_c_lin = @(t,x,u) A*(x-x_star) + B*(u-u_star);

%% ************************** Discretization ******************************

% Définissez les matrices A, B, ainsi que le pas de temps deltaT
    
% deltaT = 0.025; % Pas de temps pour la discrétisation
% ss_sup = 1.1;

% % Utilisation de la fonction c2d pour obtenir les matrices discrètes A_d et B_d
% sys = ss(A, B, eye(size(A)), 0); % Création d'un système d'état-space continu ss(A, B, C, D)
% sys_d = c2d(sys, deltaT, 'zoh'); % Discrétisation du système
% 
% % Extraction des matrices A_d et B_d du système discrétisé
% A_d = sys_d.A;
% B_d = sys_d.B;


%% ************************** Alpha Max ***********************************

delta = 10;
c = 0.9;         
num_steps = 5;
deltaT = 0.04; % Pas de temps pour la discrétisation
ss_sup = 1.5;

% A_s = eye(2);
% B_s = [];
% for k  = 1:num_steps
%     A_s = A_s * A_d;
%     B_s = [A_d^(k-1) * B_d B_s];
% end
% 
% [alpha, epsilon, K] = get_alpha_max(c, A_s, B_s, delta);
% 
% disp("alpha_max");
% disp(alpha);
% disp("epsilon_opt");
% disp(epsilon);
% 
% ss = .5;
% ss_step = 2 * ss / 10;
% 
% [x, y] = meshgrid(-ss:ss_step:ss, -ss:ss_step:ss);

%% Get Alpha Min

[x, y] = meshgrid(-ss_sup:2 * ss_sup / 8 :ss_sup, -ss_sup:2 * ss_sup / 8 :ss_sup);
alphas = zeros(size(x));
epsilons = zeros(size(x));
normKmax = 0;
Kmax = ones(1);
for i = 1:size(x,1)
    for j = 1:size(x,2)
        x_star = [x(i,j); y(i,j)]; % Comme dans le papier
%         if (i == 9) && (j == 1)
%             disp("f");
%         end
        u_star = 0;  
        [A_d, B_d] = get_A_B_vdp(x_star, u_star, deltaT);
        A_s = eye(2);
        B_s = [];
        for k  = 1:num_steps
            A_s = A_s * A_d;
            B_s = [A_d^(k-1) * B_d B_s];
        end
%         disp(A_s);
%         disp(B_s);
        [alpha_max, epsilon_opt, K] = get_alpha_max(c, A_s, B_s, delta);
        alphas(i, j) = alpha_max;
        epsilons(i, j) = epsilon_opt;
        if norm(K, 'inf') > Kmax
            Kmax = K;
        else
        end
    end
end
% close all; 
% figure;
% surf(x, y, Kmaxs);
% xlabel('x_1');
% ylabel('x_2');
% zlabel('K_{max}')

alpha_star = min(alphas(:));
epsilon_star = min(epsilons(:)); 

% % Trouver le minimum et son indice dans la matrice
% [minVal, linearInd] = min(alphas(:));
% % Convertir l'indice linéaire en indices de souscription
% [row, col] = ind2sub(size(alphas), linearInd);


disp("alpha min");
disp(alpha_star);
disp("epsilon min");
disp(epsilon_star);
% disp("row, col, alpha");
% disp(row);
% disp(col);

[minVal, linearInd] = min(epsilons(:));

% Disturbances

% [x, y] = meshgrid(-ss_sup:2 * ss_sup / 6 :ss_sup, -ss_sup:2 * ss_sup / 6 :ss_sup);
% w_rad = zeros(size(x));
% for i = 1:size(x,1)
%     for j = 1:size(x,2)
%         x_star = [x(i,j); y(i,j)]; % Comme dans le papier
%         W = get_W_vdp(x_star, epsilon_star, delta, deltaT, num_steps);
%         w_rad(i, j) = .5 * (W.sup(2) -  W.inf(2)); % pq 2 ? 
% %         disp(W.center);
%     end
% end

% figure;
% surf(x, y, ws);
% xlabel('x_1');
% ylabel('x_2');
% zlabel('W')


% % Trouver le minimum et son indice dans la matrice
% [max_w_rad, linearInd] = max(w_rad(:));
% % Convertir l'indice linéaire en indices de souscription
% % [row, col] = ind2sub(size(alphas), linearInd);
% 
% % disp("max_w_rad");
% % disp(max_w_rad);
% % disp("max_w_rad / alpha_star");
% % disp(max_w_rad / alpha_star);
% 
alpha_star = 0.5 * alpha_star;



% Vérification sanitaire

for i = 1:size(x,1)
    for j = 1:size(x,2)
        x_star = [x(i,j); y(i,j)]; % Comme dans le papier
        u_star = 0;
        [A_d, B_d] = get_A_B_vdp(x_star, u_star, deltaT);
        A_s = eye(2);
        B_s = [];
        for k  = 1:num_steps
            A_s = A_s * A_d;
            B_s = [A_d^(k-1) * B_d B_s];
        end
        K = get_K(A_s, B_s, epsilon_star, alpha_star, delta, c);
        if (alpha_star - (1 - norm(A_s - B_s * K, 'inf')) * epsilon_star) > 1e-6
%             disp("alpha_min > (1 - norm(A_s - B_s * Kmax, 'inf')) * epsilon_min : not verified")
            disp("Problem");
            disp((1 - norm(A_s - B_s * K, 'inf')) * epsilon_star);
            break
        end
    end
end

disp("Sanity check passed");

%% ********************    Construction du graphviz  **********************************
% 
grid_step = 2 * alpha_star; % alpha = + grande dis entre pt et grid_pt, au pire diag / 2 = alpha , alpha = sqrt(2) / 2 grid_step...
xlimi = ceil(ss_sup/grid_step) * grid_step * [-1, 1]; ylimi = ceil(ss_sup/grid_step) * grid_step * [-1, 1]; % Pour centrer la grille sur (0, 0).
[x, y] = meshgrid(xlimi(1):grid_step:xlimi(2), ylimi(1):grid_step:ylimi(2));
U = zonotope(zeros(size(B_s, 2), 1), delta * eye(size(B_s, 2))); % Control set centré. 
G = zeros(size(x,1) * size(x, 2)); % state = j + i * size(x, 2), encodage de la matrice d'adjacence.


disp("size(G)");
disp(size(G));


if (size(G, 1) < 700)

    % disp("interval(ref_zonotope).sup")
    % disp(interval(ref_zonotope).sup)

    for i = 1:size(x, 1)
        for j = 1:size(x, 2)
            x0 = [x(i, j); y(i, j)];
            [A_d, B_d] = get_A_B_vdp(x0, 0, deltaT);
            A_s = eye(2);
            B_s = [];
            for k  = 1:num_steps
                A_s = A_s * A_d; 
                B_s = [A_d^(k-1) * B_d B_s];
            end

            K = get_K(A_s, B_s, epsilon_star, alpha_star, delta, c);
            Z = zonotope(zeros(size(A_s, 1), 1), epsilon_star * [0 1; 1 0]); 
            track = -K*Z;
            ref_zonotope = minkDiff(U, track);
            
%             W = get_W_vdp(x0, epsilon_star, delta, deltaT, num_steps);
%             wc = W.center;

            Fw = A_s * x0 + 0.025 * B_s * ref_zonotope + zonotope([0; 0], alpha_star * [0 1; 1 0]);  

            % Une fois l'équation A* Ball(eps) + B Ball(beta) \subset B(eps - alpha)
            % ||A|| eps + ||B|| beta + w <= eps - alpha ~> 
            %(1 - |A|) eps - alpha -  w > 0
            % Pour x_0 -> x_j iff |A x_0 + B u_r  - x_j|< alpha

            vert = Fw.box.vertices;
            ind_vert = zeros(2, 4);
            % Encapsulage du zononope dans une boite, qui nécessairement inclue
            % les x_j candidats.

            for k = 1:4
                vert = Fw.box.vertices;
                [closest_x, closest_y, row, col] = project_closest(x, y, grid_step, vert(1, k), vert(2, k)); % proj sur point de la grille
                ind_vert(1, k) = row; 
                ind_vert(2, k) = col;
            end

            % Trouver la plus grande grille possible

            indi_min = min(ind_vert(1, :));
            indi_max = max(ind_vert(1, :));
            indj_min = min(ind_vert(2, :));
            indj_max = max(ind_vert(2, :));

            for ii = indi_min:indi_max
                for jj = indj_min:indj_max
                    if contains(Fw, [x(ii,jj); y(ii, jj)]) % Condition pour l'existence d'une transition.
                        target = [x(ii, jj); y(ii, jj)];
                        G(j + size(x, 2) * (i-1), jj + size(x, 2) * (ii-1)) = 1; % Matrice d'adjacence
                    end  
                end 
            end
        end
    end         


    filename = sprintf('/vanderpol/alpha%.2fepsilon%.2fmaxu%dc%.2fss%.2fdeltaT%.2f', alpha_star, epsilon_star, delta, c, ss_sup, deltaT);

    scc_raw = tarjanSCC(G); % algorithme pour extraire les composantes fortement connectées. 
    scc = [];

    for i = 1:length(scc_raw)
        if length(scc_raw{i}) > 1
            scc = [scc; {scc_raw{i}}];
        end
    end

    sccGraphviz(G, grid_step, xlimi(1), ylimi(1), scc, filename); % script pour tracer le graph.
    
else
    disp("Graph is too big (> 900 x 900)")
end

%% Compute disturbances

% Meshgrid sur le state-space.

% [x, y] = meshgrid(-ss_sup:2 * ss_sup / 6 :ss_sup, -ss_sup:2 * ss_sup / 6 :ss_sup);
% w_diam = zeros(size(x));
% for i = 1:size(x,1)
%     for j = 1:size(x,2)
%         x_star = [x(i,j); y(i,j)]; % Comme dans le papier
%         W = get_W_vdp(x_star, epsilon_star, delta, deltaT, num_steps);
%         w_diam(i, j) = .5 * (W.sup(2) -  W.inf(2));
%     end
% end

% figure;
% surf(x, y, ws);
% xlabel('x_1');
% ylabel('x_2');
% zlabel('W')


% % Trouver le minimum et son indice dans la matrice
% [maxWDiam, linearInd] = max(w_diam(:));
% % Convertir l'indice linéaire en indices de souscription
% % [row, col] = ind2sub(size(alphas), linearInd);
% 
% disp("maxWDiam");
% disp(maxWDiam);


% For ss = .05, eps = 0.0594 => W_max = 0.0045 (7 %)

%% Functions


function W = get_W_vdp(x_star, epsilon, delta, deltaT, num_steps)
    n = 2; % Dimensions of the state space
    S = zonotope(x_star, epsilon * eye(n));
    U = zonotope(zeros(num_steps, 1), delta * eye(num_steps));
    % linearization_area = A_s * zonotope(zeros(n, 1), epsilon_opt * eye(n)) + B_s * zonotope(zeros(num_steps, 1), delta * eye(num_steps));
    f = @(x, u) [x(1) + deltaT * 2 * x(2); 
                 x(2) + deltaT * (-.8 * x(1) + 2 * x(2) - 10 * x(1) * x(1) * x(2) + u)];
    fn = @(x, u_list) apply_fn_iteratively(f, x, ones(num_steps, 1));
    [A, B, W] = get_multistep_w(fn, S, U);
end

function [A_d, B_d] = get_A_B_vdp(x_star, u_star, deltaT)
    A = [0, 2; -0.8 - 5*x_star(1)*x_star(2), 2 - 2.5*x_star(1)^2];
%      A = [0, 2; -0.8 - 20*x_star(1)*x_star(2), 2 - 10*x_star(1)^2];
    B = [0; 1];
    % Utilisation de la fonction c2d pour obtenir les matrices discrètes A_d et B_d
    sys = ss(A, B, eye(size(A)), 0); % Création d'un système d'état-space continu ss(A, B, C, D)
    sys_d = c2d(sys, deltaT, 'zoh'); % Discrétisation du système
    % Extraction des matrices A_d et B_d du système discrétisé
    A_d = sys_d.A;
    B_d = sys_d.B;
end

function K = get_K(A, B, epsilon, alpha, delta, c)
    [n, n] = size(A);
    [n, m] = size(B); 
    % Variables de décision
    K = sdpvar(m, n, 'full'); % La matrice K à optimiser
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
    
    % Résoudre le problème de programmation linéaire
    diagnostics = optimize(Constraints, Objective, options);
    
    if diagnostics.problem == 0
        K = value(K); % Retourner la solution
    else
        K = []; % Pas de solution trouvée
        disp('Le problème n’a pas de solution.');
    end
end