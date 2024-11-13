% -------------------------------------------------------------------------
% Solution for the advection equation
% Upwind in space + Implicit Euler in time
% -------------------------------------------------------------------------

method = 'implicit';

a = 1;            
x_min = -1;         
x_max = 1;         
n_cells = 200;   
x = linspace(x_min, x_max, n_cells);

dx = (x_max-x_min)/(n_cells-1); 
C_numbers = [0.1, 0.5, 0.8, 1]; % Courant

% Initial condition
phi_0 = initialcondition(n_cells,x_min,x_max);   
phi = phi_0;

figure;
plot(x, phi,'DisplayName','Initial Condition','LineWidth', 1); 
hold on;
% eachplot = plot(x, phi,'DisplayName','Solution');
legend show;
xlabel('x'); 
ylabel('\phi');
title('Solution for the advection equation using Upwind and Implicit Euler');

for C = C_numbers

    dt = C*dx/a;
    phi = phi_0;

    A = (1 + C) * eye(n_cells) - C * diag(ones(n_cells-1,1), -1);
    A(1,end) = -C;
    
    % Upwind in space and implicit euler in time + periodic BC
    for i = 0:dt:2
        
        phi_new = A \ phi;
        phi = phi_new;
    
        % if mod(i,2) < dt
        %     set(eachplot, 'YData', phi,'DisplayName', ['dt = ', num2str(i)]);
        %     drawnow;
        % end
    
    end

    plot(x, phi,'DisplayName',sprintf('C = %.1f', C),'LineStyle', '--','LineWidth', 1);

end

[diffusive_error,dispersive_error] = error_calculation(n_cells,C_numbers,method);
