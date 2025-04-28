%% planar mathematical pendulum
clearvars
close all

% Parameters
g = 9.81;               % gravity (m/s^2)
L = 1.0;                % length of pendulum (m)
delta_t = 0.05;         % time step
T = 10;                 % total time (s)
N = round(T/delta_t);   % number of time steps

% Preallocate
phi = zeros(1, N+1);
v = zeros(1, N+1);
t = 0:delta_t:T;

% Initial conditions
phi(1) = pi/4;  % 45 degrees
v(1) = 0;

method = 'explicit Euler';
linearized = false;

switch method
    case 'explicit Euler'
        for n=1:N
            phi(n+1) = phi(n) + v(n)*delta_t;

            if linearized
            v(n+1) = v(n) - g/L*phi(n)*delta_t;
            else
            v(n+1) = v(n) - g/L*sin(phi(n))*delta_t; 
            end
        end

    case 'implicit Euler'
    for n = 1:N
        % Initial guess for Newton iteration
        phi_guess = phi(n);
        v_guess = v(n);

        % TODO linearisierte variante 

        % Newton iteration
        for iter = 1:10
            % Residuals
            F1 = phi_guess - phi(n) - delta_t * v_guess;
            F2 = v_guess - v(n) + delta_t * (g/L) * sin(phi_guess);

            % Jacobian
            J = [1,        -delta_t;
                 delta_t*(g/L)*cos(phi_guess), 1];

            % Update
            delta = -J \ [F1; F2];
            phi_guess = phi_guess + delta(1);
            v_guess = v_guess + delta(2);

            if norm(delta) < 1e-10
                break;
            end
        end

        % Store updated values
        phi(n+1) = phi_guess;
        v(n+1) = v_guess;
    end 

    case 'symplectic Euler'
    for n=1:N

        phi(n+1) = phi(n) + v(n)*delta_t;

        if linearized
            v(n+1) = v(n) - g/L*phi(n+1)*delta_t;
        else
            v(n+1) = v(n) - g/L*sin(phi(n+1))*delta_t; 
        end
    end

    otherwise
        error('unknown method name')

end

% Plot
visualize.angle(phi, v, t)

% energy
visualize.energy(phi, v, t, L, g)

% animate
visualize.animate(phi, t, L, delta_t)

% compare with analytical solution
%visualize.analytical(phi,t,g,L)


%#ok<*UNRCH> % suppress warning unreachable code lines
