% -------------------------------------------------------------------------
% Dispersive and dissipative error calculation
% -------------------------------------------------------------------------

function [diffusive_error,dispersive_error] = error_calculation(n_cells,C_numbers,method)

theta_values = linspace(0,pi,n_cells);

figure;
subplot(2,1,1);
hold on;
xlabel('Wavenumber');
ylabel('Diffusive Error');
title('Diffusive Error');
grid on;

subplot(2,1,2);
hold on;
xlabel('Wavenumber');
ylabel('Dispersive Error');
title('Dispersive Error');
grid on;

G_all_values = cell(length(C_numbers), 1);

for C_idx = 1:length(C_numbers)

    C = C_numbers(C_idx);

    diffusive_error = zeros(size(theta_values));
    dispersive_error = zeros(size(theta_values));
    G_values = zeros(size(theta_values));


    for idx = 1:length(theta_values)
    
        theta = theta_values(idx);
    
         
       switch method
    
            case 'explicit'  % Explicit Euler with upwind scheme
                G_numerical = 1 - C + C * (cos(theta) - 1i*sin(theta));
            
            case 'implicit' % Implicit Euler with upwind scheme
                G_numerical = 1 / (1 + C*(1 - cos(theta) + 1i*sin(theta)));
            
            case 'lax-wendroff'
                G_numerical = 1 - 1i*C*sin(theta) + (C^2)*(cos(theta)-1);
            
           case 'crank-nicolson' % Crank-Nicolson with central differences 
                G_numerical = (1 - 1i*(C/2)*sin(theta)) / (1 + 1i*(C/2)*sin(theta));
            
           case 'leap-frog' % Leap-frog with central differences
                G_numerical = -C*1i*sin(theta)+sqrt(-C^2*(sin(theta)^2)+1);
                
            otherwise
                fprintf('Method not available.');
                break
        end
    
        diffusive_error(idx) = abs(G_numerical);
        dispersive_error(idx) = angle(G_numerical)/(-C*theta);
           
        G_values(idx) = G_numerical;
    
    end

    subplot(2,1,1);
    plot(theta_values, diffusive_error, 'DisplayName', sprintf('C = %.2f', C), 'LineWidth', 0.8);
    
    subplot(2,1,2);
    plot(theta_values, dispersive_error, 'DisplayName', sprintf('C = %.2f', C), 'LineWidth', 0.8);

    % Stability diagram plot
    theta_values_stability = linspace(0, 2*pi, 2*n_cells);
    G_values_stability = zeros(size(theta_values_stability));

    for idx = 1:length(theta_values_stability)
        theta = theta_values_stability(idx);
        
        switch method
            case 'explicit'
                G_values_stability(idx) = 1 - C + C * (cos(theta) - 1i*sin(theta));

            case 'implicit'
                G_values_stability(idx) = 1 / (1 + C*(1 - cos(theta) + 1i*sin(theta)));

            case 'lax-wendroff'
                G_values_stability(idx) = 1 - 1i*C*sin(theta) + (C^2)*(cos(theta)-1);

            case 'crank-nicolson'
                G_values_stability(idx) = (1 - 1i*(C/2)*sin(theta)) / (1 + 1i*(C/2)*sin(theta));

            case 'leap-frog'
                G_values_stability(idx) = -C*1i*sin(theta) + sqrt(1 - C^2*(sin(theta)^2));
  
        end
    end

    G_all_values{C_idx} = G_values_stability;

end

subplot(2,1,1);
legend show;
xt = xticks;
yt = yticks;
xticklabels(arrayfun(@(x) sprintf('%.2f', x), xt, 'UniformOutput', false));
yticklabels(arrayfun(@(y) sprintf('%.2f', y), yt, 'UniformOutput', false));

subplot(2,1,2);
legend show;
xt = xticks;
yt = yticks;
xticklabels(arrayfun(@(x) sprintf('%.2f', x), xt, 'UniformOutput', false));
yticklabels(arrayfun(@(y) sprintf('%.2f', y), yt, 'UniformOutput', false));

figure;
hold on;
alpha = linspace(0, 2*pi, n_cells); 
x_circle = cos(alpha); 
y_circle = sin(alpha); 
plot(x_circle, y_circle, 'k--', 'LineWidth', 1.5);

title('Stability diagram');
xlabel('Re(G)');
ylabel('Im(G)');
legend('Stability Region');
axis equal;
grid on;

for C_idx = 1:length(C_numbers)
    C = C_numbers(C_idx);
    G_values = G_all_values{C_idx};
    plot(real(G_values), imag(G_values),'DisplayName', sprintf('C = %.2f', C), 'LineWidth', 0.8); 
end

legend show;

end
