function [Phi, Gamma] = getlinmodel_exp(x)
    dt = .05;
    A = [0 1; cos(x(1)) * 9.81 -1];
    B = [0 0; 0 1];
    sys = ss(A, B, [], []); % Créez un système d'état-espace continu
    sys_d = c2d(sys, dt); % Convertissez le système en un système discret avec un pas de temps dt
    Phi = sys_d.a; % Récupérez la matrice de transition d'état discrète
    Gamma = sys_d.b; % Récupérez la matrice d'entrée discrète
end
