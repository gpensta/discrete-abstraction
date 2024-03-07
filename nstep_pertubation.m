g = -9.81;
l = 1;
beta = 1; 
m = 1; 
dt = .1;

f = @(x, u) [x(1) + dt * x(2); 
             x(2) + dt * (-g/l*sin(x(1)) - beta/(m*l^2)*x(2) + u)];

fn = @(x, u_list) apply_fn_iteratively(f, x, u_list);


x_prime = fn([0.; 0.], [6, 6]);

% a = [1. 0.1;
%      .9810 0.9];
% b = [0.;
%      .1]; 
%  
% A = eye(2);
% B = [];
% for k  = 1:num_steps
%     A = A * a;
%     B = [a^(k-1) * b B];
% end

% disp(A);
% disp(B); 

epsilon = .0837;
delta = 6;

n = 2;
m = num_steps;



S = zonotope(zeros(n, 1), epsilon * eye(n));
U = zonotope(zeros(m, 1), delta * eye(m));

% profile on

[A, B, W] = get_multistep_w(fn, S, U);

disp(W.inf);
disp(W.sup);


% f = @(x) [x(1) + dt * x(2); x(2) + dt * (-g/l*sin(x(1)) - x(2) + x(4))]; 
% S = zonotope(zeros(n, 1), epsilon * eye(n));
% U = zonotope(zeros(2, 1), delta * eye(2));

% profile off
% profile viewer



