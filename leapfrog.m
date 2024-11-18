% -------------------------------------------------------------------------
% Solution for the advection equation
% Central differences in space + Leap-frog in time
% -------------------------------------------------------------------------

method = 'leap-frog'; % label for the stability plot and error calculation

a = 1;            
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
title('Solution for the advection equation using Central Differences and Leap-frog');

for C = C_numbers

    dt = C*dx/a;
    phi_old = phi; % level n-1

    % Using Crank-Nicolson to compute the first level n

    alpha = C/4;

    A = eye(n_cells) + alpha * diag(ones(n_cells-1,1), 1) - alpha * diag(ones(n_cells-1,1), -1);
    A(1, end) = -alpha;      
    A(end, 1) = alpha; 

    phi_exp = phi + alpha * [phi(end); phi(1:end-1)] - alpha * [phi(2:end); phi(1)]; 
    
    phi_new = A \ phi_exp; % 1st time step

    for i = dt:dt:2 %starting one time step after
    
        phi_new = phi_old + C * ([phi(end);phi(1:end-1)]-[phi(2:end);phi(1)]);
        phi_old = phi;
        phi = phi_new;
    
    end
    
    plot(x, phi,'DisplayName',sprintf('C = %.1f', C),'LineStyle', '--','LineWidth', 1);

end

[diffusive_error,dispersive_error] = error_calculation(n_cells,C_numbers,method);
