% -------------------------------------------------------------------------
% VERY SIMPLE FIRST APPROACH
% Upwind in space + Euler in time + Initial condition square wave
% -------------------------------------------------------------------------

a = 1;            
x_min = -1;         
x_max = 1;         
n_cells = 200;   
x = linspace(x_min, x_max, n_cells);

dx = (x_max-x_min)/(n_cells-1);  
dt = 0.01;

% Courant
C = a*dt/dx; % like this or should I define dt based on C that I know is stable?

% Initial condition
phi_0 = (x >= -0.6) & (x <= -0.4);
phi = phi_0;                           

% is 1 revolution 1 time step???
% Upwind in space and explicit euler in time + periodic BC
phi_new = phi - a*dt/dx * (phi-[phi(end),phi(1:end-1)]);

figure;
plot(x, phi_0); 
hold on;
plot(x, phi_new);
legend('Initial condition', 'After one revolution');
xlabel('x'); 
ylabel('\phi');





