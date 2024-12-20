% -------------------------------------------------------------------------
% Initial condition
% -------------------------------------------------------------------------

function [phi_0] = initialcondition(n_cells,x_min,x_max)

x = linspace(x_min, x_max, n_cells);

a_elipse = 0.5;
z = -0.7;
delta = 0.005;
alpha = 10;
beta = log(2)/(36*delta^2);

G = @(x, beta, z) exp(-beta * (x - z).^2);

F = @(x, alpha, a_elipse) sqrt(max(1 - alpha^2 * (x - a_elipse).^2, 0));

phi_0 = zeros(n_cells,1);

for i = 1:n_cells
    if x(i) >= -0.8 && x(i) <= -0.6
        phi_0(i) = (1/6) * (G(x(i), beta, z - delta) + G(x(i), beta, z + delta) + 4 * G(x(i), beta, z));

    elseif x(i) >= -0.4 && x(i) <= -0.2
        phi_0(i) = 1;

    elseif x(i) >= 0 && x(i) <= 0.2
        phi_0(i) = 1 - abs(10 * (x(i) - 0.1));

    elseif x(i) >= 0.4 && x(i) <= 0.6
        phi_0(i) = (1/6) * (F(x(i), alpha, a_elipse - delta) + F(x(i), alpha, a_elipse + delta) + 4 * F(x(i), alpha, a_elipse));

    else
        phi_0(i) = 0;
    end
end

% uncomment if running only this function alone
% figure;
% plot(x, phi_0);
% xlabel('x');
% ylabel('\phi');
% title('Initial Condition');

end