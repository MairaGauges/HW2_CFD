% -------------------------------------------------------------------------
% Dispersive and dissipative error calculation
% -------------------------------------------------------------------------


function [diffusive_error,dispersive_error] = error_calculation(x_max,x_min,n_cells,C,a,method)

dx = (x_max-x_min)/(n_cells-1);
dt = C*dx/a;

% k_values = linspace(-(pi/dx-pi*(n_cells*dx)), pi/dx-pi*(n_cells*dx), n_cells);
k_values = linspace(-pi/dx,pi/dx,n_cells);

diffusive_error = zeros(size(k_values));
dispersive_error = zeros(size(k_values));
G_values = zeros(size(k_values));

figure;
hold on;

for idx = 1:length(k_values)

    k = k_values(idx); 

    G_exact = exp(-1i*a*k*dt);
     
   switch method
            case 'explicit'  % Explicit Euler with upwind scheme
                G_numerical = 1 - C + C * (cos(k*dx) - 1i*sin(k*dx));
    
            case 'implicit' % Implicit Euler with upwind scheme
                G_numerical = 1 / (1 + C * (1 - cos(k*dx) - 1i * sin(k*dx)));
            
    
            otherwise
                error_calculation('ups');
    end
    
    diffusive_error(idx) = abs(G_numerical)/abs(G_exact);

    dispersive_error(idx) = angle(G_numerical)/angle(G_exact);
    % dispersive_error(idx) = angle(G_numerical)/(-C*k*dx);

    G_values(idx) = G_numerical;
end


plot(real(G_values), imag(G_values));
hold on;
theta = linspace(0, 2*pi, n_cells); 
x_circle = cos(theta); 
y_circle = sin(theta); 
plot(x_circle, y_circle, 'k--', 'LineWidth', 1.5);

title('Stability diagram');
xlabel('Re(G)');
ylabel('Im(G)');
legend('Amplification factor', 'Stability Region');
axis equal;
grid on;

% For plotting
plot_k_values = k_values > 0 & k_values < pi / dx;

figure;
subplot(2,1,1);
plot((k_values(plot_k_values)*dx), diffusive_error(plot_k_values));
xlabel('Wavenumber');
ylabel('Diffusive Error');
title('Diffusive Error');

subplot(2,1,2);
plot((k_values(plot_k_values)*dx), dispersive_error(plot_k_values));
xlabel('Wavenumber');
ylabel('Dispersive Error');
title('Dispersive Error');

end
