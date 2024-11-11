% -------------------------------------------------------------------------
% Solution for the advection equation
% Lax-Wendroff scheme
% -------------------------------------------------------------------------

method = 'lax-wendroff';

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

figure;
plot(x, phi,'DisplayName','Initial Condition','LineWidth', 0.8); 
hold on;
legend show;
xlabel('x'); 
ylabel('\phi');

for n = 0:dt:2

    for i = 2:n_cells-1

        phi_new(i) = phi(i) - (C/2) * (phi(i+1) - phi(i-1)) + ...
               (C^2/2) * (phi(i+1) - 2 * phi(i) + phi(i-1));

    end

    phi_new(1) = phi(1) - (C/2) * (phi(2) - phi(n_cells)) + ...
               (C^2/2) * (phi(2) - 2 * phi(1) + phi(n_cells));
    
    
    phi_new(n_cells) = phi(n_cells) - (C/2) * (phi(1) - phi(n_cells-1)) + ...
               (C^2/2) * (phi(1) - 2 * phi(n_cells) + phi(n_cells-1));

    %phi_new(n_cells-1);
    %phi_new(n_cells) = phi_new(2);

    phi = phi_new;


end

plot(x, phi,'DisplayName','Solution','LineStyle', '--','LineWidth', 1);

%[diffusive_error,dispersive_error] = error_calculation(x_max,x_min,n_cells,C,a,method);

