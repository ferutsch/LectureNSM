% -------------------------------------------------------------------------
% fem2d: Educational MATLAB code for finite element discretizations in 2D
% Cyber-Physical Simulation Group, TU Darmstadt
% Creator: Prof. Dr. Oliver Weeger
% -------------------------------------------------------------------------
%
% Mesh generation for a rectangle of size Lx*Ly with ex*ey elements
% Bilinear, 4-node quadrilateral Q4 mesh
%

function [X, elem_nodes] = mesh_rect_Q4(Lx,Ly,ex,ey)

nx = ex+1;
ny = ey+1;
X = zeros(2,nx*ny);
k = 0;
for j=0:ey
    xi2 = -0.5 + j/ey;
    for i=0:ex
        xi1 = -0.5 + i/ex;
        k = k+1;
        X(:,k) = [Lx*xi1; Ly*xi2];
    end
end

elem_nodes = zeros(ex*ey,4);
k = 0;
for j = 0:(ey-1)
    for i = 1:ex
        i0 = j*nx + i;
        k = k+1;
        elem_nodes(k,:) = [i0, i0+1, i0+1+nx, i0+nx];
    end
end

end