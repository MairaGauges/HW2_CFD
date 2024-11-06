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
C = 0.8; % Courant
dt = C*dx/a;

% Initial condition
phi = initialcondition(n_cells,x_min,x_max);   
phi_new = phi;

figure;
plot(x, phi,'DisplayName','Initial Condition','LineWidth', 1); 
hold on;
% eachplot = plot(x, phi,'DisplayName','Solution');
legend show;
xlabel('x'); 
ylabel('\phi');

lambda = a*dt/dx;

A = (1 + lambda) * eye(n_cells) - lambda * diag(ones(n_cells-1,1), -1);
A(1,end) = -lambda;

% Upwind in space and implicit euler in time + periodic BC
for i = 0:dt:2
    
    phi_new = (A \ phi')';
    phi = phi_new;

    % if mod(i,2) < dt
    %     set(eachplot, 'YData', phi,'DisplayName', ['dt = ', num2str(i)]);
    %     drawnow;
    % end

end

plot(x, phi,'DisplayName','Solution','LineStyle', '--','LineWidth', 1);

