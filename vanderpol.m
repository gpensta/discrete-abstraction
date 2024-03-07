clear 
close all
addpath('./CORA_2022')
addpath('./functions')


%% *************************** Dynamics ***********************************

f_u =  @(t,x,u)(-[ -2*x(2,:) ; 0.8*x(1,:) + 10*x(1,:).^2.*x(2,:) - 2*x(2,:) + u] );
n = 2;
m = 1; % number of control inputs


%% ************************** Discretization ******************************

deltaT = 0.01;
%Runge-Kutta 4
k1 = @(t,x,u) (  f_u(t,x,u) );
k2 = @(t,x,u) ( f_u(t,x + k1(t,x,u)*deltaT/2,u) );
k3 = @(t,x,u) ( f_u(t,x + k2(t,x,u)*deltaT/2,u) );
k4 = @(t,x,u) ( f_u(t,x + k1(t,x,u)*deltaT,u) );
f_ud = @(t,x,u) ( x + (deltaT/6) * ( k1(t,x,u) + 2*k2(t,x,u) + 2*k3(t,x,u) + k4(t,x,u)  )   );

%% ************************** Linearization *******************************

% D�finition des points d'�quilibre x_star et u_star
x_star = [.7; .7]; % Comme dan sle papier
u_star = 0;  

% Matrices A et B lin�aris�es �valu�es au point d'�quilibre
A = [0, 2; -0.8 - 20*x_star(1)*x_star(2), 2 - 10*x_star(1)^2];
B = [0; 1];

% Fonction dynamique lin�aris�e continue
f_c_lin = @(t,x,u) A*(x-x_star) + B*(u-u_star);

%% ************************** Discretization ******************************

% D�finissez les matrices A, B, ainsi que le pas de temps deltaT
    
deltaT = 0.01; % Pas de temps pour la discr�tisation

% Utilisation de la fonction c2d pour obtenir les matrices discr�tes A_d et B_d
sys = ss(A, B, eye(size(A)), 0); % Cr�ation d'un syst�me d'�tat-space continu ss(A, B, C, D)
sys_d = c2d(sys, deltaT, 'zoh'); % Discr�tisation du syst�me

% Extraction des matrices A_d et B_d du syst�me discr�tis�
A_d = sys_d.A;
B_d = sys_d.B;


%% ************************** Alpha Max ***********************************

delta = 1;
c = .9;         

num_steps = 2;

A_s = eye(2);
B_s = [];
for k  = 1:num_steps
    A_s = A_s * A_d;
    B_s = [A_d^(k-1) * B_d B_s];
end

[alpha_max, epsilon_opt, K] = get_alpha_max(c, A_s, B_s, delta);


disp("alpha_max");
disp(alpha_max);
disp("epsilon_opt");
disp(epsilon_opt);

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
        disp('Le probl�me n�a pas �t� r�solu');
    end
end


