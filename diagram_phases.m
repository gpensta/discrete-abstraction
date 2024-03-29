% Définition des matrices A et B

close all
clear;

[A, B] = get_A_B_pendulum([0.; 0.], .01);

% Définition de l'état initial x0
x0 = [-.5; 0.5]; % Exemple, à remplacer par vos valeurs

% Définition de plusieurs états initiaux x0
x0s = ((1 - (-1)) * rand(5, 2) - 1* ones(5, 2))'; % Chaque colonne est un état initial différent
x0s = [x0s [0; 0]];
% Nombre de pas de temps à simuler
N = 25;

% Création d'une grille de points dans l'espace d'état [-1, 1] x [-1, 1]
[x, y] = meshgrid(-1:0.1:1, -1:0.1:1);

% Calcul des vecteurs de phase pour chaque point de la grille
u = A(1,1)*x + A(1,2)*y;
v = A(2,1)*x + A(2,2)*y;

% Normalisation des vecteurs pour une meilleure visualisation
norm = sqrt(u.^2 + v.^2);
u = u ./ norm;
v = v ./ norm;

% Dessin du diagramme des phases
figure;
quiver(x, y, u, v);
hold on;

% Couleurs pour les différentes trajectoires
colors = lines(size(x0s, 2));

% Simulation et dessin de chaque trajectoire
for j = 1:size(x0s, 2)
    % Initialisation de la trajectoire pour l'état initial actuel
    trajectory = zeros(2, N+1);
    trajectory(:,1) = x0s(:,j);
    
    % Simulation de la trajectoire
    for k = 1:N
        trajectory(:,k+1) = A * trajectory(:,k); % + B*u si vous avez une entrée u
    end
    
    % Dessin de la trajectoire
    plot(trajectory(1,:), trajectory(2,:), '-', 'Color', colors(j,:), 'LineWidth', 2);
    plot(trajectory(1,1), trajectory(2,1), 'o', 'MarkerFaceColor', colors(j,:), 'MarkerEdgeColor', 'k'); % Point de départ
end

xlabel('x_1');
ylabel('x_2');
title('Pendules et trajectoires');
axis([-1 1 -1 1]);
axis square;
grid on;

