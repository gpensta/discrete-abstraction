function x_prime = nl_model(x, u)
    g = -9.81; l = 1.; m = 1; beta = 1; dt = .01;
    x1 = x(1); x2 = x(2);
    x_prime = [x1  + x2 * dt; x2 + dt * (-g / l * sin(x1) - beta / (m * l^2) * x2 + u)];
end

