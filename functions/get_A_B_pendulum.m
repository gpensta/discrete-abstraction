function [A_d, B_d] = get_A_B_pendulum(x_star, deltaT)
    g = -9.81; l = 1.; beta = 1.; m = 1.; 
    A = [0 1;  -(g/l) * cos(x_star(1)) (beta)/(m*l*l)];
    B = [0; 1];
    sys = ss(A, B, eye(size(A)), 0); % Création d'un système d'état-space continu ss(A, B, C, D)
    sys_d = c2d(sys, deltaT, 'zoh'); % Discrétisation du système
    A_d = sys_d.A;
    B_d = sys_d.B;
end
