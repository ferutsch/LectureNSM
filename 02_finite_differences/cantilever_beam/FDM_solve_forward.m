function [Q, M, w] = FDM_solve_forward(N,L,q,F,E,I)

delta_x = L/(N-1);
x_axis = 0:delta_x:L;

% coeficient pattern for forward difference method
forward_coefs = zeros(N-1,N);
for i=1:N-1
    forward_coefs(i,i) = -1;
    forward_coefs(i,i+1) = 1;
end

% coefficient pattern for second derivative
second_derivative_coefs = zeros(N-2,N);
for i = 1:N-2
    second_derivative_coefs(i,i) = 1;
    second_derivative_coefs(i,i+1) = -2;
    second_derivative_coefs(i,i+2) = 1;
end

%% calculate shear force

% vector
Q_vector = zeros(N,1);
Q_vector(1:N-1) = -q*delta_x;
Q_vector(N) = F; % BC: Q(x=L) = F

% coefficient matrix
Q_coefs = zeros(N,N);
Q_coefs(1:N-1,1:N) = forward_coefs;
Q_coefs(N,N) = 1; % BC: Q(x=L)

% solve
Q = Q_coefs\Q_vector;


%% calculate bending moment

% vector
M_vector = zeros(N,1);
M_vector(1:N-1) = Q(i)*delta_x;
M_vector(N) = 0; % BC: M(x=L) = 0

% coefficient matrix
M_coefs = zeros(N,N);
M_coefs(1:N-1,1:N) = forward_coefs;
M_coefs(N,N) = 1; % BC: M(x=L)

% solve
M = M_coefs\M_vector;

%% calculate beam deflection
% vector
w_vector = zeros(N,1);
for i=3:N 
    w_vector(i) = -M(i-1)/(E*I(i-1));
end
w_vector(1) = 0; % BC: w(x=0) = 0
w_vector(2) = 0; % BC: w'(x=0) = 0

% coefficient matrix
w_coefs = zeros(N,N);
w_coefs(3:N,1:N) = second_derivative_coefs/(delta_x^2);
w_coefs(1,1) = 1; % BC w(x=0)
w_coefs(2,1:3) = [-3, 4, -1]/(2*delta_x); % BC: w'(x=0) -> 2nd order forward difference

% solve
w = w_coefs\w_vector;

end

