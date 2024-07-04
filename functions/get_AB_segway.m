function [A, B] = get_AB_segway(x_star, u_star)

    % Define the symbolic variables for states and inputs
    syms x1 x2 tau real
    
    % Paramètres du Segway
    m = 100; 
    M = 2; % Masse de la roue 
    r = .2; % Rayon des roues
    Iw = .6; % Moment d'inertie des roues
    Ib = 35; % Moment d'inertie du corps
    l = .2; % Longueur du pendule (du centre de masse à l'axe des roues)
    g = 9.81; % Accélération due à la gravité sur Terre

    % Calcul des paramètres dynamiques du modèle
    m11 = (m + M)*r^2 + Iw;
    m12 = m*l*r;
    m22 = m*l^2 + Ib;
    G_b = m*g*l;


    % Define the functions bar_M and M1 as given
    bar_M = m11*m22 - (m12*cos(x1))^2;
    M1 = m11 + m12*cos(x1);

    % Define the system dynamics f
    f = [x2; 
         (m12 -m12^2 * x2^2 * cos(x1) * sin(x1) + (m11 * G_b * sin(x1))) / bar_M ...
         - (M1 * tau) / bar_M]; % Comme dans le papier


    % Define the symbolic variables as vectors
    x_sym = [x1; x2];
    u_sym = tau;

    % Compute the Jacobian matrices symbolically
    A_sym = jacobian(f,  x_sym);
    B_sym = jacobian(f, u_sym);

    % Substitute the linearization point into the Jacobian matrices
    A = double(subs(A_sym, [x_sym; u_sym], [x_star; 0]));
    B = double(subs(B_sym, [x_sym; u_sym], [x_star; 0]));
end