% -------------------------------------------------------------------------
% Solution for the advection equation
% Upwind in space + Implicit Euler in time
% -------------------------------------------------------------------------

a = 1;            
x_min = -1;         
x_max = 1;         
n_cells = 200;   
x = linspace(x_min, x_max, n_cells);

dx = (x_max-x_min)/(n_cells-1); 
C = 0.8; % Courant
dt = C*dx/a;

% Initial condition

a_elipse = 0.5;
z = -0.7;
delta = 0.005;
alpha = 10;
beta = log(2)/(36*delta^2);

G = @(x, beta, z) exp(-beta * (x - z).^2);

F = @(x, alpha, a_elipse) sqrt(max(1 - alpha^2 * (x - a_elipse).^2, 0));

phi_0 = zeros(1, n_cells);

for i = 1:n_cells
    if x(i) >= -0.8 && x(i) <= -0.6
        phi_0(i) = (1/6) * (G(x(i), beta, z - delta) + G(x(i), beta, z + delta) + 4 * G(x(i), beta, z));

    elseif x(i) >= -0.4 && x(i) <= -0.2
        phi_0(i) = 1;

    elseif x(i) >= 0 && x(i) <= 0.2
        phi_0(i) = 1 - abs(10 * (x(i) - 0.1));

    elseif x(i) >= 0.4 && x(i) <= 0.6
        phi_0(i) = (1/6) * (F(x(i), alpha, a_elipse - delta) + F(x(i), alpha, a_elipse + delta) + 4 * F(x(i), alpha, a_elipse));
   
    else
        phi_0(i) = 0;
    end
end

figure;
plot(x, phi_0, 'k-');
xlabel('x');
ylabel('\phi');
title('Initial Condition');


phi = phi_0;   
phi_new = phi;

figure;
plot(x, phi_0,'DisplayName','Initial Condition'); 
hold on;
eachplot = plot(x, phi,'DisplayName','Solution');
legend show;
xlabel('x'); 
ylabel('\phi');

lambda = a*dt/dx;         % CFL-like parameter for implicit scheme

% Create the coefficient matrix for the implicit scheme
A = (1 + lambda) * eye(n_cells) - lambda * diag(ones(n_cells-1,1), -1);
A(1,end) = -lambda;


% Upwind in space and implicit euler in time + periodic BC
for i = 0:dt:2
    
    phi_new = (A \ phi')';
    phi = phi_new;

    if mod(i,2) < dt
        set(eachplot, 'YData', phi,'DisplayName', ['dt = ', num2str(i)]);
        drawnow;
    end

end


