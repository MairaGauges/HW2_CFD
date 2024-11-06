% -------------------------------------------------------------------------
% Solution for the advection equation
% Upwind in space + Explicit Euler in time
% -------------------------------------------------------------------------

method = 'explicit';

a = 1;            
x_min = -1;         
x_max = 1;         
n_cells = 200;   
x = linspace(x_min, x_max, n_cells);

dx = (x_max-x_min)/(n_cells-1); 
C = 1; % Courant
dt = C*dx/a;

phi = initialcondition(n_cells,x_min,x_max); 
phi_new = phi;

figure;
plot(x, phi,'DisplayName','Initial Condition','LineWidth', 0.8); 
hold on;
legend show;
xlabel('x'); 
ylabel('\phi');

for i = 0:dt:2

    phi_new = phi - a*dt/dx * (phi-[phi(end),phi(1:end-1)]);
    phi = phi_new;

end

plot(x, phi,'DisplayName','Solution','LineStyle', '--','LineWidth', 1);

[diffusive_error,dispersive_error] = error_calculation(x_max,x_min,n_cells,C,a,method);
