% -------------------------------------------------------------------------
% Dispersive and dissipative error calculation
% -------------------------------------------------------------------------


function [diffusive_error,dispersive_error] = error_calculation(n_cells,C,method)

% dx = (x_max-x_min)/(n_cells-1);
% dt = C*dx/a;

% k_values = linspace(-(pi/dx-pi*(n_cells*dx)), pi/dx-pi*(n_cells*dx), n_cells);
% k_values = linspace(-pi/dx,pi/dx,n_cells);
theta_values = linspace(0,pi,n_cells);

diffusive_error = zeros(size(theta_values));
dispersive_error = zeros(size(theta_values));
G_values = zeros(size(theta_values));
G_exact_values = zeros(size(theta_values));

for idx = 1:length(theta_values)

    % k = k_values(idx); 
    theta = theta_values(idx);

    % G_exact = exp(-1i*a*k*dt);
    G_exact = exp(-1i*C*theta);
     
   switch method
            case 'explicit'  % Explicit Euler with upwind scheme
                G_numerical = 1 - C + C * (cos(theta) - 1i*sin(theta));
    
            case 'implicit' % Implicit Euler with upwind scheme
                G_numerical = 1 / (1 + C * (1 - cos(theta) - 1i * sin(theta)));

            case 'lax-wendroff'
                G_numerical = 1 - 1i*C*sin(theta) - (C^2)*(1-cos(theta));

            case 'crank-nicolson'
                G_numerical = (1 - 1i*(C/2)*sin(theta))/(1 + 1i*(C/2)*sin(theta));
    
            otherwise
                fprintf('ups');
    end
    
    % diffusive_error(idx) = abs(G_numerical)/abs(G_exact);
    diffusive_error(idx) = abs(G_numerical);

    % phase_numerical = angle(G_numerical);
    % phase_numerical = mod(phase_numerical,2*pi);

    % phase_exact = angle(G_exact);
    % phase_exact = mod(phase_exact,2*pi);

    % if phase_exact == 0
    %     dispersive_error(idx) = phase_numerical-phase_exact;
    % else
        dispersive_error(idx) = angle(G_numerical)/(-C*theta);
    % end
       
    G_values(idx) = G_numerical;
    G_exact_values(idx) = G_exact;

end

% figure;
% plot(real(G_exact_values), imag(G_exact_values));
% axis equal;

figure;
plot(real(G_values), imag(G_values), 'b-','LineWidth', 0.8);
hold on;
plot(real(G_values), -imag(G_values),'b-','LineWidth', 0.8);
hold on;
alpha = linspace(0, 2*pi, n_cells); 
x_circle = cos(alpha); 
y_circle = sin(alpha); 
plot(x_circle, y_circle, 'k--', 'LineWidth', 1.5);

title('Stability diagram');
xlabel('Re(G)');
ylabel('Im(G)');
legend('Amplification factor','', 'Stability Region');
axis equal;
grid on;


figure;

subplot(2,1,1);
plot(theta_values(2:end-1), diffusive_error(2:end-1));
xlabel('Wavenumber');
ylabel('Diffusive Error');
title('Diffusive Error');

xt = xticks;
yt = yticks;
xticklabels(arrayfun(@(x) sprintf('%.2f', x), xt, 'UniformOutput', false));
yticklabels(arrayfun(@(y) sprintf('%.2f', y), yt, 'UniformOutput', false));

subplot(2,1,2);
plot(theta_values(2:end-1), dispersive_error(2:end-1));
xlabel('Wavenumber');
ylabel('Dispersive Error');
title('Dispersive Error');

xt = xticks;
yt = yticks;
xticklabels(arrayfun(@(x) sprintf('%.2f', x), xt, 'UniformOutput', false));
yticklabels(arrayfun(@(y) sprintf('%.2f', y), yt, 'UniformOutput', false));

end
