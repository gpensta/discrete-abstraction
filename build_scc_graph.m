clear;
% profile on
% alpha = 0.099; epsilon = 4 * alpha; delta = 20;

% OBJECTIF : Construire une reach set abstraction : F(T_i, U_i) \subset
% T_j, identifier sous graphs fortement connectés avec algorithme Tarjan,
% la supervision doit garantir de rester dans le SCC. 

% a = [1. 0.1;
%      .9810 0.9];
% b = [0.;
%      .1]; % Dynamique du pendule linéarisée pas de temps de .1.

% A = eye(2);
% B = [];
% for k  = 1:num_steps
%     A = A * a;
%     B = [a^(k-1) * b B];
% end

delta = 6; % U = [-delta, delta]
c = .6; % what is allowed for tracking, what is left is for ref control.
num_steps = 2;

[A_d, B_d] = get_A_B_pendulum([0.; 0.], .1); % deltaT = .01
A = eye(2);
B = [];
for k  = 1:num_steps
    A = A * A_d;
    B = [A_d^(k-1) * B_d B];
end


[alpha, epsilon, K] = get_alpha_max(c, A, B, delta); 
% + grd alpha, avec epsilon optimal, epsilon limite la taille de K, donc
% moins de ref control.
% alpha = 0.05;
epsilon = 1 * alpha; 
disp("alpha");
disp(alpha);
disp("epsilon");
disp(epsilon); 
grid_step = 2 * alpha; % alpha = + grande dis entre pt et grid_pt, au pire diag / 2 = alpha , alpha = sqrt(2) / 2 grid_step...
ss = .5; % size state-space
xlimi = ceil(ss/grid_step) * grid_step * [-1, 1]; ylimi = ceil(ss/grid_step) * grid_step * [-1, 1]; % Pour centrer la grille sur (0, 0).
[x, y] = meshgrid(xlimi(1):grid_step:xlimi(2), ylimi(1):grid_step:ylimi(2));
U = zonotope(zeros(size(B, 2), 1), delta * eye(size(B, 2))); % Control set centré. 
G = zeros(size(x,1) * size(x, 2)); % state = j + i * size(x, 2), encodage de la matrice d'adjacence.
disp("size(G)");
disp(size(G));

Z = zonotope(zeros(size(A, 1), 1), epsilon * [0 1; 1 0]); 
track = -K*Z;
ref_zonotope = minkDiff(U, track);
disp("interval(ref_zonotope).sup")
disp(interval(ref_zonotope).sup)

for i = 1:size(x, 1)
    for j = 1:size(x, 2)
        x0 = [x(i, j); y(i, j)];
%         [A_d, B_d] = get_A_B_pendulum(x0); % deltaT = .01
%         A = eye(2);
%         B = [];
%         for k  = 1:num_steps
%             A = A * A_d;
%             B = [A_d^(k-1) * B_d B];
%         end
%         [alpha, epsilon, K] = get_alpha_max(c, A, B, delta);
        Fw = A * x0 + B * ref_zonotope + zonotope([0; 0], alpha * [0 1; 1 0]);  
        
        % Une fois l'équation A* Ball(eps) + B Ball(beta) \subset B(eps - alpha)
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



filename = sprintf('/test3/alpha%.2fepsilon%.2fmaxu%dc%.2fss%.2f', alpha, epsilon, delta, c, ss);

scc_raw = tarjanSCC(G); % algorithme pour extraire les composantes fortement connectées. 
scc = [];

for i = 1:length(scc_raw)
    if length(scc_raw{i}) > 1
        scc = [scc; {scc_raw{i}}];
    end
end

sccGraphviz(G, grid_step, xlimi(1), ylimi(1), scc, filename); % script pour tracer le graph.
% profile off
% profile viewer
