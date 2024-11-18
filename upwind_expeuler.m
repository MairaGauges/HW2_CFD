% -------------------------------------------------------------------------
% Solution for the advection equation
% Upwind in space + Explicit Euler in time
% -------------------------------------------------------------------------

method = 'explicit'; % label for the stability plot and error calculation

a = 1; % velocity can be changed to negative here           
x_min = -1;         
x_max = 1;         
n_cells = 200;   
x = linspace(x_min, x_max, n_cells);

dx = (x_max-x_min)/(n_cells-1); 
C_numbers = [0.1,0.5,0.8,1]; % Courant numbers, needs to be kept as a vector

phi_0 = initialcondition(n_cells,x_min,x_max); 
phi = phi_0;

figure;
plot(x, phi,'DisplayName','Initial Condition','LineWidth', 0.8); 
hold on;
legend show;
xlabel('x'); 
ylabel('\phi');
title('Solution for the advection equation using Upwind and Explicit Euler');

for C = C_numbers

    dt = abs(C*dx/a);
    phi = phi_0;
    
    for i = 0:dt:2

        if a > 0
            phi_new = phi - C * (phi-[phi(end);phi(1:end-1)]);

        else % Upwind needs to shift
            phi_new = phi - C * (phi-[phi(2:end); phi(1)]);

        end

        phi = phi_new;
    
    end

    plot(x, phi,'DisplayName',sprintf('C = %.1f', C),'LineStyle', '--','LineWidth', 1);

end

[diffusive_error,dispersive_error] = error_calculation(n_cells,C_numbers,method);
