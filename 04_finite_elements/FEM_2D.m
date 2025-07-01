%% Code example 4: Finite Element Method

% Lecture 'Numerical Simulation Methods'
% Felix Rutsch, July 2025

% solve stationary diffusion equation with FEM

% y
% ^
% |
% o ---> x

clearvars
close all

%--- PARAMETERS (SI UNITS) ---
Lx = 1.0;              % [m]
Ly = 1.0;              % [m]
nx = 10;                % elements in x direction
ny = 10;                % elements in y direction

k   = 1;               % [W/(m K)], thermal conductivity 
%alpha = k/(rho*cp);   % [m^2/s], thermal diffusivity

% adapt heat source

%Q0 = 1000;
%src_size = 0.2; % [m] size of the heat source
%ncell = src_size/(2*hx);
%cx = ceil(Nx/2);
%cy = ceil(Ny/2);

%f = zeros(Nx,Ny);
%f( cx-ncell+1:cx+ncell, cy-ncell+1:cy+ncell ) = Q0;

%--- FLATTEN SOURCE ---
%Fv = reshape(f, Nx*Ny, 1);
%[nodes, edof] = mesh(Lx, Ly, nx, ny); 


% --- GENERATE FE MESH ---
[nodes, edof] = mesh_rect_Q4(Lx, Ly, nx, ny); % edof = "element degrees of freedom"
% move a single node by changing nodes vector on that position?

%--- ASSEMBLY ---
nnodes = (nx+1)*(ny+1);         % number of nodes
nelem = nx*ny;                  % number of elements

K = sparse(nnodes, nnodes);     % initialize stiffness matrix
F = zeros(nnodes, 1);           % initialize load vector

% assembly loop
for e = 1:nelem

    if e==45
        Fe = [1 1 1 1] * 12;
    else
        Fe = [0 0 0 0];
    end

    % element stiffness matrix
    [ke, fe] = elem_heat_Q4(Fe,k);

    % assembly
    K(edof(e,:), edof(e,:)) = K(edof(e,:), edof(e,:)) + ke;
    F(edof(e,:)) = F(edof(e,:)) + fe;
end

%--- BOUNDARY CONDITIONS ---
boundaryNodes = [];
for n=1:nnodes
    if nodes(1,n) == -Lx/2 || nodes(1,n) == Lx/2 || ...
       nodes(2,n) == -Ly/2 || nodes(2,n) == Ly/2          
    boundaryNodes = [boundaryNodes n];
    end
end
freeNodes = setdiff(1:nnodes, boundaryNodes);


%--- SOLVE ---
u = zeros(nnodes,1);
u(freeNodes) = K(freeNodes, freeNodes) \ F(freeNodes);


%--- VISUALIZE ---
u_matrix = reshape(u,[ny+1,nx+1]);
x = -0.5:0.1:0.5;
y = -0.5:0.1:0.5;
imagesc(x,y,u_matrix);

%visualize.heat_source(x,y,hx,hy,f)
%visualize.temperature_interp(x,y,hx,hy,phi)
%visualize.temperature(x,y,phi)


function [Ke, be] = elem_heat_Q4(Fe,kappa)

Xe = [0 0 1 1;
      0 1 1 0]; % node coordinates

qp = [-1  1 -1  1;
      -1 -1  1  1] / sqrt(3); % gauss points

qw = [ 1  1  1  1];

% Initalization
Ne = 4;                 % Size of Ke, be
Ke = zeros(Ne,Ne);      % Element stiffness matrix
be = zeros(Ne,1);       % Element force vector
qn = size(qp,2);        % Number of quadrature points
    
% Quadrature loop
for k = 1:qn

    % Quadrature point
    xi = qp(1,k);
    eta = qp(2,k);
    
    % Evaluation of parametric gradient shape functions Ni
    dNdXi = 0.25 * [ -1+eta, -1+xi;
                      1-eta, -1-xi;
                      1+eta,  1+xi;
                     -1-eta,  1-xi ];

    % Jacobian (parametric gradient of coordinate transformation)
    J = Xe * dNdXi;
    Jinv = inv(J);
    detJ = abs(det(J));

    % Physical gradient of shape functions Ni 
    dNdX = dNdXi * Jinv;     

    % Integrand evaluation for element stiffness matrix
    Kek = dNdX * kappa * dNdX' * detJ;
    Ke = Ke + qw(k) * Kek;

    
    % element load vector

    % Evaluation of shape functions Ni
    Neval = 0.25 * [ (1-xi)*(1-eta);
                     (1+xi)*(1-eta);
                     (1+xi)*(1+eta);
                     (1-xi)*(1+eta) ];

    % Integrand evaluation for element force vector
    bek = Neval * (Fe * Neval) * detJ;
    be = be + qw(k) * bek;
end

end
