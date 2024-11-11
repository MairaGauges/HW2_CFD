
a = 1;            
x_min = -1;         
x_max = 1;         
n_cells = 200;   
x = linspace(x_min, x_max, n_cells);

dx = (x_max-x_min)/(n_cells-1); 
C = 0.8; % Courant
dt = C*dx/a;

phi = initialcondition(n_cells,x_min,x_max); 

% wave 1
%phi(x>-0.5) = 0;

%wave 2
% phi(x>0) = 0;
% phi(x<-0.5) = 0;

%wave 3
% phi(x>0.25) = 0;
% phi(x<-0) = 0;

%wave 2
% phi(x<0.25) = 0;


%test constant initial condition
% phi = 2*ones(1,200);

phi_new = phi;

figure;
plot(x, phi,'DisplayName','Initial Condition','LineWidth', 0.8); 
hold on;
legend show;
xlabel('x'); 
ylabel('\phi');

for i = 0:dt:2
    [flux_lim,r] = flux_limiter_array(phi);
    
    %first order upwind
    F_low = phi-[phi(end),phi(1:end-1)];
    
    %second order upwind as higher order scheme
    %F_high = 0.5*(3*phi -4*[phi(end),phi(1:end-1)] + [phi(end-1:end), phi(1:end-2)]);

    %F according to lax-wendruff 
    F_high = 0.5*(1+C)*([phi(2:end), phi(1)]) + 0.5*(C-1)*([phi(end),phi(1:end-1)]) - C*phi;


    F = F_low + flux_lim.*(F_high-F_low);
    phi_new = phi - a*dt/dx * F;
    phi = phi_new;




% testing / plotting  how errors grow
%     if mod(count,100) == 0
%         plot(x, phi,'DisplayName','Solution','LineStyle', '--','LineWidth', 1);
%     end


end

[flux_lim,r] = flux_limiter_array(phi);


plot(x, phi,'DisplayName','Solution','LineStyle', '--','LineWidth', 1);
%plot(x, r,'DisplayName','r','LineStyle', '-','LineWidth', 1);
%plot(x, flux_lim,'DisplayName','flux limiter','LineStyle', '-','LineWidth', 1);

function [f_lim,r]= flux_limiter_array(phi)

r1 = (phi - [phi(end),phi(1:end-1)]);
r2=([phi(2:end), phi(1)] - phi); 

r = r1./r2;
r(r1 == r2) = 1;

% test minmod 
% f_lim = max(0,min(1,r));
% f_lim(isinf(r)) = 0;

%test van Leer
f_lim = 2*r ./ (r+1); 
f_lim(isinf(r)) = 0;
f_lim(r<0) = 0;

%Superbee
%f_lim = max(0, max(min(2*r,1),min(r,1)));
%f_lim(isinf(r)) = 0;

%design own flux limiter
% f_lim = 0;
% f_lim = max(0,min(r,1));

% f_lim = zeros(1,200);
% f_lim(r >= 0.1 & r <= 1.9) = 1;
f_lim = 1;

end 

