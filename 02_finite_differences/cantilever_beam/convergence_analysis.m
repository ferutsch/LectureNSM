% convergence study

clearvars
close all

L = 1;              % [m] length of cantilever
q = 0;              % [N/m] line load
F = 1;              % [N] load at tip of cantilever
E = 1;              % [N/m2] Youngs modulus
I_const = 1;        % [m4] moment of inertia

stepnumber = [11, 51, 101, 501, 1001, 5001, 10001]';
delta_x = L./(stepnumber-1);

% analytical solution
% for constant EI, F = 1, q(x) = 0
w_function = @(x) (F*L*x.^2/2 - F*x.^3/6)/(E*I_const);
w_analyt = w_function(1);

error1 = zeros(size(stepnumber,1),1);
error2 = zeros(size(stepnumber,1),1);

for i=1:size(stepnumber,1)

    N = stepnumber(i);

    I = I_const*ones(N,1);  

    [~, ~, w1] = FDM_solve_forward(N,L,q,F,E,I);
    error1(i) = abs(w1(end)-w_analyt);

    [~, w2] = FDM_solve_2nd_order(N,L,q,F,E,I);
    error2(i) = abs(w2(end)-w_analyt);

end

loglog(delta_x,error1,'-o','DisplayName','forward difference')
hold on
loglog(delta_x,error2,'-o','DisplayName','2nd order central difference')
grid on
ylabel('error |w*-w| at x=L')
xlabel('stepsize \Delta x')
legend('Location','southeast')

%p1 = ( log(error_trpz(2)/error_trpz(4)) )/( log(stepsize(2)/stepsize(4)) )