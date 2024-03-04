clear;
% alpha = 0.099; epsilon = 4 * alpha; max_u = 20;

% OBJECTIF : Construire une reach set abstraction : F(T_i, U_i) \subset
% T_j, identifier sous graphs fortement connectés avec algorithme Tarjan,
% la supervision doit garantir de rester dans le SCC. 

alphas = [.1]; % [.1 .15 .2]; % max dist entre un point de la grille et n'importe quel point. 
max_us = [6]; % [10 20 30]; 
epsilons = [1]; %[2 3 4];
cs = [.8];  %[.7 .8 .9]; % Le controle est divisé en tracking (c * max_u) et en reference ((1 - c) * max_u)

for alpha_idx = 1:length(alphas)
    for max_u_idx = 1:length(max_us)
        for eps_idx = 1:length(epsilons)
            for c_idx = 1:length(cs)
                alpha = alphas(alpha_idx);
                epsilon = epsilons(eps_idx) * alpha;
                max_u = max_us(max_u_idx);
                grid_step = 2 * alpha; % alpha = + grande dis entre pt et grid_pt, au pire diag / 2 = alpha , alpha = sqrt(2) / 2 grid_step...
                xlimi = ceil(.8/grid_step) * grid_step * [-1, 1]; ylimi = ceil(.8/grid_step) * grid_step * [-1, 1]; % Pour centrer la grille sur (0, 0).
                [x, y] = meshgrid(xlimi(1):grid_step:xlimi(2), ylimi(1):grid_step:ylimi(2));
                c = cs(c_idx); % what is allowed for tracking, what is left is for ref control.
                g = -9.81; l = 1.0; m = 1; beta = 1; dt = .1;
                f = @(x, u) [x(1) + dt * x(2); x(2) + dt * (-g/l*sin(x(1)) - beta/(m*l^2)*x(2) + u)];
                U = zonotope([0; 0], max_u * [1 0; 0 1]);
                G = zeros(size(x,1) * size(x, 2)); % state = j + i * size(x, 2).
                for i = 1:size(x, 1)
                    for j = 1:size(x, 2)
                        x0 = [x(i, j); y(i, j)];
                        Z = zonotope(x0, epsilon * [0 1; 1 0]);
                        [a, b] = linearize(f, x0, [0]);
                        A = eye(2);
                        B = [];
                        for k  = 1:2
                            A = A * a;
                            B = [a^(k-1) * b B];
                        end
                        K = best_k(A, B, c * max_u, epsilon); % delta = u, epsilon
                        beta = (((1 - norm(A - B * K, 'inf')) * epsilon - alpha)) / norm(B, 'inf');
                        if beta > 0
                            track = -K*Z + K*x0;
                            ref_zonotope = minkDiff(U, track);
                            Fw = A * x0 + B * ref_zonotope + zonotope([0; 0], alpha * [0 1; 1 0]);  
                            vert = Fw.box.vertices;
                            ind_vert = zeros(2, 4);
                            for k = 1:4
                                vert = Fw.box.vertices;
                                [closest_x, closest_y, row, col] = project_closest(x, y, grid_step, vert(1, k), vert(2, k));
                                ind_vert(1, k) = row; 
                                ind_vert(2, k) = col;
                            end
                            indi_min = min(ind_vert(1, :));
                            indi_max = max(ind_vert(1, :));
                            indj_min = min(ind_vert(2, :));
                            indj_max = max(ind_vert(2, :));

                            for ii = indi_min:indi_max
                                for jj = indj_min:indj_max
                                    if contains(Fw, [x(ii,jj); y(ii, jj)])   
                                        target = [x(ii, jj); y(ii, jj)];
                                            G(j + size(x, 2) * (i-1), jj + size(x, 2) * (ii-1)) = 1;
                                    end  
                                end
                            end
                        end
                    end         
                end
               
                
                filename = sprintf('/test3/alpha%.2fepsilon%.2fmaxu%dc%.2f', alpha, epsilon, max_u, c);
                
                scc_raw = tarjanSCC(G);
                scc = [];
                
                for i = 1:length(scc_raw)
                    if length(scc_raw{i}) > 1
                        scc = [scc; {scc_raw{i}}];
                    end
                end
                sccGraphviz(G, grid_step, xlimi(1), ylimi(1), scc, filename);
            end
        end
    end
end





