% -------------------------------------------------------------------------
% Dispersive and dissipative error calculation
% -------------------------------------------------------------------------


function [diffusive_error,dispersive_error] = error_calculation(x_max,x_min,n_cells,C,a,method)

dx = (x_max-x_min)/(n_cells-1);
dt = C*dx/a;

k_values = linspace(0, pi/dx, n_cells);
% k_values = linspace(0, pi);

diffusive_error = zeros(size(k_values));
dispersive_error = zeros(size(k_values));

for idx = 1:length(k_values)
    k = k_values(idx); 

    % G_exact = exp(-1i*a*k*dt);
    
   switch method
            case 'explicit'  % Explicit Euler with upwind scheme
                G_numerical = 1 - C + C * (cos(k*dx) - 1i*sin(k*dx));
    
            case 'implicit' % Implicit Euler with upwind scheme
                G_numerical = 1 / (1 + C * (1 - cos(k*dx) - 1i * sin(k*dx)));
    
            otherwise
                error_calculation('ups');
    end
    
    % Calculate diffusive error as |G_numerical| - |G_exact|, ideally |G_exact| = 1
    % diffusive_error(idx) = abs(abs(G_numerical) - 1);
    diffusive_error(idx) = abs(G_numerical);
    
    % Calculate dispersive error as phase difference between G_numerical and G_exact
    % dispersive_error(idx) = abs(angle(G_numerical)-angle(G_exact));
    dispersive_error(idx) = angle(G_numerical)/(-C*k*dx);
end

figure;
subplot(2,1,1);
plot(k_values, diffusive_error);
xlabel('Wavenumber');
ylabel('Diffusive Error');
title('Diffusive Error');

subplot(2,1,2);
plot(k_values, dispersive_error);
xlabel('Wavenumber');
ylabel('Dispersive Error');
title('Dispersive Error');

end
