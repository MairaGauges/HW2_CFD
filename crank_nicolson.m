% -------------------------------------------------------------------------
% Solution for the advection equation
% Crank-Nicolson scheme
% -------------------------------------------------------------------------

method = 'crank-nicolson';

a = 1;            
x_min = -1;         
x_max = 1;         
n_cells = 200;   
x = linspace(x_min, x_max, n_cells);

dx = (x_max-x_min)/(n_cells-1); 
C = 0.8; % Courant
dt = C*dx/a;

phi = initialcondition(n_cells,x_min,x_max); 

phi_new = phi;

alpha = C/4;

A = eye(n_cells) + alpha * diag(ones(n_cells-1,1), 1) - alpha * diag(ones(n_cells-1,1), -1);
A(1, end) = -alpha;      
A(end, 1) = alpha; 

figure;
plot(x, phi,'DisplayName','Initial Condition','LineWidth', 0.8); 
hold on;
legend show;
xlabel('x'); 
ylabel('\phi');

for n = 0:dt:2

    phi_exp = phi + alpha * [phi(end); phi(1:end-1)] - alpha * [phi(2:end); phi(1)]; 

    phi_new = A \ phi_exp;

    phi = phi_new;

end

plot(x, phi,'DisplayName','Solution','LineStyle', '--','LineWidth', 1);

[diffusive_error,dispersive_error] = error_calculation(x_max,x_min,n_cells,C,a,method);
