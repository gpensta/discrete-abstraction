function [Phi, Gamma] = getlinmodel_exp(x)
    dt = .05;
    A = [0 1; cos(x(1)) * 9.81 -1];
    B = [0 0; 0 1];
    sys = ss(A, B, [], []); % Cr�ez un syst�me d'�tat-espace continu
    sys_d = c2d(sys, dt); % Convertissez le syst�me en un syst�me discret avec un pas de temps dt
    Phi = sys_d.a; % R�cup�rez la matrice de transition d'�tat discr�te
    Gamma = sys_d.b; % R�cup�rez la matrice d'entr�e discr�te
end
