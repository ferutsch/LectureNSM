% Finite Difference Method

% Felix Rutsch, February 2025

% cantilever beam

clearvars;

%% parameters
N = 10;             % number of points
L = 1;              % [m] length of cantilever
q = 0;              % [N/m] line load
F = 1;              % [N] load at tip of cantilever
E = 1;              % [N/m2] Youngs modulus

delta_x = L/(N-1);
x_axis = 0:delta_x:L;

I = 1*ones(1,N);   % [m4] moment of inertia
%I = 1.5-x_axis*1; % linear moment of inertia

%% solve

[Q, M, w] = FDM_solve_forward(N,L,q,F,E,I);


%% analytical solution
% for constant EI, F = 1, q(x) = 0
I_const = 1;
w_analytical = @(x) (F*L*x.^2/2 - F*x.^3/6)/(E*I_const);

%% plotting

subplot(3,1,1)
plot(x_axis,Q)
title('Querkraft')
ylabel('Q(x)')
xlabel('x')
grid on

subplot(3,1,2)
plot(x_axis,M)
title('Biegemoment')
ylabel('M(x)')
xlabel('x')
grid on

subplot(3,1,3)
plot(x_axis,w)
hold on
plot(x_axis,w_analytical(x_axis))
set(gca, 'YDir','reverse')
title('Durchbiegung')
ylabel('w(x)')
xlabel('x')
grid on



