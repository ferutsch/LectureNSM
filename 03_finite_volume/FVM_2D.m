clearvars

%--- PARAMETERS (SI UNITS) ---

% y
% ^
% |
% o ---> x

Lx = 1.0;              % [m]
Ly = 1.0;              % [m]
Nx = 100; Ny = 100;      % grid points
hx = Lx/Nx;            % [m]
hy = Ly/Ny;            % [m]

rho = 1.0e3;           % [kg/m^3], mass density, e.g. 1.0e3 for water
cp  = 4.18e3;          % [J/(kg K)], heat capacity, e.g. 4.18e3
k   = 0;            % [W/(m K)], thermal conductivity, e.g. 0.6 
%alpha = k/(rho*cp);    % [m^2/s], thermal diffusivity

v1 = 0.00001; v2 = 0.00001;   % [m/s]

method = 'UDS'; % for the convective term: UDS or CDS

%--- SPATIAL GRID ---
x = hx/2:hx:Lx-hx/2;  % cell‐center x
y = hy/2:hy:Lx-hy/2;  % cell‐center y

Q0 = 1000;
src_size = 0.2; % [m] size of the heat source
ncell = src_size/(2*hx);
cx = ceil(Nx/2);
cy = ceil(Ny/2);

f = zeros(Nx,Ny);
f( cx-ncell+1:cx+ncell, cy-ncell+1:cy+ncell ) = Q0;

%--- FLATTEN SOURCE ---
Fv = reshape(f, Nx*Ny, 1);

%--- ASSEMBLE A and b ---
N = Nx*Ny;
A = sparse(N,N);
b = zeros(N,1);
idx = @(i,j) (j-1)*Nx + i;

switch method
    case 'CDS' % central differencing scheme
    % CDS for both the advective and diffusive terms   
    aW = k/hx^2 + rho*cp * v1 / (2*hx);
    aE = k/hx^2 - rho*cp * v1 / (2*hx);
    aS = k/hy^2 + rho*cp * v2 / (2*hy);
    aN = k/hy^2 - rho*cp * v2 / (2*hy);
    aP = aW + aE + aS + aN;

    case 'UDS' % upwind differencing scheme
    % CDS for the diffusive term
    % UDS for the advective term (maximum operator)
    aW = k/hx^2 + max(0,rho*cp * v1 / hx);
    aE = k/hx^2 + max(0,-rho*cp * v1 / hx);
    aS = k/hy^2 + max(0,rho*cp * v2 / hy);
    aN = k/hy^2 + max(0,-rho*cp * v2 / hy);
    aP = aW + aE + aS + aN;
end

for j = 2:Ny-1
  for i = 2:Nx-1
    p = idx(i,j);
    
    A(p, idx(i-1,j  )) = -aW;
    A(p, idx(i+1,j  )) = -aE;
    A(p, idx(i  ,j-1)) = -aS;
    A(p, idx(i  ,j+1)) = -aN;
    A(p, idx(i,  j  )) = aP;
    b(p) = Fv(p);
  end
end


%--- DIRICHLET BC: phi=300 K on all boundaries ---
T0 = 300;
for i = 1:Nx
  for j = [1 Ny]
    p = idx(i,j);
    A(p,p) = 1; b(p) = T0;
  end
end
for j = 2:Ny-1
  for i = [1 Nx]
    p = idx(i,j);
    A(p,p) = 1; b(p) = T0;
  end
end

%--- SOLVE & RESHAPE ---
phi_vec = A \ b;
phi = reshape(phi_vec, [Nx, Ny]); 

%--- VISUALIZE ---
visualize.heat_source(x,y,hx,hy,f)
visualize.temperature_interp(x,y,hx,hy,phi)
visualize.temperature(x,y,phi)



