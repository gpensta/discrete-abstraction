function [A, B] = get_lin_model(x)
    m = 1.;  g = -9.81; l = 1.; dt = .1; beta = 1;
    A = [1,    dt;
        - cos(x(1)) * dt * g / l ,  1 - dt * beta / (m * l * l)]; 
    B = [0 0;
         0 dt];
end