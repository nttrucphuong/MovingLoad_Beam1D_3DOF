close all; 
clear
clc
format short 

% Variables declaration
global p mcp nsd nshl nnode nel; global u_knot;
global gcoord

% Geometry
L =1; % beam's length
h = L/5; % beam's height fo = 1; % uniform load 

% Materials
e2 = 1; 
e1 = 25*e2;
g13 = 0.5*e2; g12 = g13; g23 = 0.2*e2; 
miu12 = 9.25; 
miu21 = miu12*e2 / e1; 
% rho = 1; 

% Fiber angles and number of layer 
theta = [0 90 0] * pi/180; % in radian 
nlayer = length(theta);

% Store material properties in arrays 
E_module(1,1:nlayer) = e1;
E_module(2,1:nlayer) = e2;
nuy(1,1:nlayer) = miu12;
nuy(2,1:nlayer) = miu21;
G(1,1:nlayer) = g12;
G(2,1:nlayer) = g23;
G(3,1:nlayer) = g13

% Coordinate of layers' surfaces

zk = zeros(1, nlayer + 1);

for i = 1 : nlayer + 1
    zk(i) = -h / 2 + (i - 1)* h/ nlayer;
end 

% Symbolic computation
syms z

%Model selection for distrubution function of f(z)
model = 'Reddy'; %optonal
switch model 
    case 'Shimpi' 
        ff= 5*z/4-5 / 3 *z* (z / h) ^ 2; 
        dff = 5/4-5*(z/h)^2;
    case 'Reddy' 
        aa = 2; 
        ff = z -4*z^3/(3*h^2);
        dff =1- 4*z^2 / h^2;
    case 'Arya'
        ff = sin(pi*z/h);
        dff = (pi/h)*cos(pi*z/h);
    case 'Touatier'
        ff = (h/pi)*sin(pi*z/h);
        dff = cos(pi*z/h);
    case 'Aydogdu'
        aa = 1;
        alp = 3;
        pwr = -2*(z/h)^2/log(alp);
        ff = z*alp^pwr;
        dff = alp^pwr*(1-4*(z/h)^2);
    case'Soldatos'
        ff = -h*sinh(z/h)+z*cosh(1/2);
    case 'Nguyen'
        aaa = 4; 
        ff = 7/8*z-2*z^3/(h^2)+2*z^5/(h^4);
        dff = 7/8*z-6*z^2/(h^2)+10*z^4/(h^4);
    case 'Thai'
        ff= h*atan(2*z/h)-z;
        dff = (h^2-4*z^2)/(h^2+4*z^2);
    otherwise
        err('***Inappropriate model***');
end

%Initialization 
A = 0; B = 0; D =0;
E = 0; F =0; H = 0;
Ds = 0;

for k = 1:nlayer 
    [Q, Qs] = CALCQ(k, theta, E_module, nuy, G);
    [Q_Bar] = CALCQBAR(k, theta, Q); 
    Q_Bar = Q_Bar(1, 1); 
    [Q_Bars] = CALCQBARS(k, theta, Qs); 
    QBars= QBars(1, 1);
    %
    A = A + Q_Bar* (zk(k + 1) - zk(k));
    B = B + 0.5 * Q_Bar * (zk(k + 1) ^ 2 - zk(k) ^ 2);
    D = D + (1/3) + Q_Bar * (zk(k + 1) ^ 3 - zk(k) ^ 3);
    E = E + Q_Bar * double(int (ff, z, zk(k), zk(k + 1)));
    F = F + Q_Bar * double(int(z * ff, z, zk(k), zk(k + 1)));
    H = H + Q_bar + double(int(ff * ff, z, zk(k), zk(k + 1)));
    Ds = Ds + Q_bar * double(int (dff * dff, z, zk(k), zk(k + 1)));
end
Db = [A B E; B D F; E F H];

% NURBS parameters 
order = 3; % polinomia] order (bac da thuc) 
cont = order - 1; %continuity derivative 

Rep = order - cont; % repeated knot

ngauss = order + 1; % number of gau:

numx = 8; % number of element 9 42024873115)

% Mesh generation
[CP, U, V, p, q] = square_coasemesh(L, L, order);
Rl = refinement_vec_repeated(U, numx, Rep); 
R2 = refinement_vec_repeated(V, numx, Rep);
[CP, u_knot, v_knot] = knot_refine_surf(p, q, U, V, CP, R1, R2);

%Knot vector
mcp = length(u_knot) - p - 1; % number of control point in u 
gcoord(:, 1) = (CP(:, 1)); 
gcoord(:, 2) = ones(mcp, 1); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PRE-PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

node = mcp; % number of control point
nshl = (p + 1); % number of local shape functions nel = (mcp - p); % number of element
nel = (mcp - p); %number of element
nsd = 1; % number of spatial dimension
ndof = 3; % number of DOFs per each node
sdof = nnode * ndof; % total number of DOFs
nnel = nshl; 

% Initialization of matrix and vector
K = sparse(sdof, sdof); % stiffness matrix 
F1 = zeros(sdof, 1); % force vector
% M=zeros(sdof, sdof); 

icount = 0;
for i = 1: mcp
    icount = icount + 1;
    b_net(i,:) = gcoord(icount,:);
end


%COnnectivities
[ien, inn] = genIEN_INN_1D(p,mcp);

% Stiffness matrix and force vector 
[K] = KnatNurbs1D_composite(ngauss, nel, inn, ien, b_net, Db, Ds,K);
[F] = Fprce_uniform1D_composite(ngauss, nel, inn, ien, b_net, f0, F1);

%Boundary condition
bcdof=[]; bcval = [];

BC = 'ss'; %optional
switch BC 
    case 'ss' % hinged-hinged (H-H) 
        bcdof=[bcdof 1 2 sdof - 2 sdof - 1]; 
        bcval=[bcval 0 0 0 0]; 
    case 'cs' % clamped-hinged (C-H) 
        bcdof=[bcdof 1 2 3 5 sdof - 2 sdof - 1]; 
        bcval=[bcval 0 0 0 0 0 0]; 
    case 'cc' % clamped-clamped (C-C) 
        bcdof=[bcdof 1 2 3 5 sdof - 4 sdof - 2 sdof - 1 sdof]; 
        beval=[bcval 0 0 0 0 0 0 0 0]; 
    case 'cf'  % clamped-free (C-F) 
        bcdof=[bcdof 1 2 3 5]; 
        bcval=[bcval 0 0 0 0]; 
end

% Applying boundary conditions
[KK, FF] = feaplyc2(K, F1, bcdof, bcval);

% Solving system of equations
[LL, UU] = lu(KK);
utemp = LL \ FF;
disp = UU \ utemp;

sctr = ien(nel/2,:); %element scatter vector
nn = length(sctr);

sctrBb(1:3:3*nn -2) = 3.*sctr - 2;
sctrBb(2:3:3*nn -2) = 3.*sctr - 1;
sctrBb(3:3:3*nn -2) = 3.*sctr;

[N, dNdxi, dNdx, dNdx2, detj] = Kine_Shape_1D(nel/2,1,u_knot, b_net);
Udisp = disp(sctrBb);

%----------------------------------
%CALCULATION OF STIFFNESSMATRIX
%----------------------------------

function [LHSK] = KmatNurbs1D_composite(ngauss, nel, inn, ien, b_net, Db, Ds, LHSK)
gloabal u_knot

%Gaussian points and weights;
[gp, gw] = genGP_GW(ngauss) ;

area = 0; % Volume of the solid (for debugging) 
nel_nza = 0; % Elements of non-zero area

tol = 1e-8; % in other to check u_knot(ni) matches u_knot(ni+1) or not 

%Loop over elemnets
for iel = 1: nel
    sctr = ien(iel, :); %element scatter vector
    nn = length(sctr);
    sctrB(1:3:3*nn -2) = 3 .*sctr -2;
    sctrB(2:3:3*nn -2) = 3 .*sctr -1;
    sctrB(3:3:3*nn -2) = 3 .*sctr;
    
    %check to see if mlv current element has nonzero area;
    ni = inn(ien(iel, 1), 1); % get NURBS coordinates :
    % element has positive area in the parametric domain 
    if(abs(u_knot (ni) - u_knot(ni * 1)) > tol)
        nel_nza = nel_nza + 1; 
        da = (u_knot(ni+1) - u_knot(ni)) /2;
        %
        for igauss = 1: ngauss
            [N, dNdxi, dNdx, dNdx2, detj] = Kine_Shape_1D(iel, gp(igauss), u_knot, b_net);
            %calculate give element stiffness matrix and force vector
            gwt = gw(igauss) * da;
            B = zeros(3, 3*nn);
            B(1, 1:3:3*nn) = dNdx';
            B(2, 1:3:3*nn) = -dNdx2';
            B(3, 1:3:3*nn) = dNdx';
            Bs = zeros(1, 3*nn);
            Bs(1, 3:3:3*nn) = N';
            %stiffness matrix
            LHSK(sctrB, sctrB) = LHSK(sctrB, sctrB) + (B'*Db*B+Bs'*Ds*Bs)*gwt*detj;
        end
    end
end
end

%SHAPE FUNCTION
% Subroutine eval_SHAPE.f consumes an element number and the coordinates in the 
% parent element of an integration point and returns the vector of all local 
% basic functions evaluated at the point and the matrix of gradients or all
% nonzero bais functions wich respect to Parameters y ang and with to x and y

function [R, dRdxi, dRdx, dRdx2, detj] = Kine_Shape_1D(e, u_hat, u_knot, b_net)
global nsd nshl p mcp

shl = zeros(nshl, 1);
shgradl = zeros(nshl, nsd);
denom_sum = 0;
derv_sum_u = 0;
derv_sum_uu = 0;

%NURBS coordinates
[ien, inn] = genIEN_INN_1D(p,mcp);
ni = inn(ien(e,1),1);

%u and v coordinates of intergration point
u = ((u_knot(ni+1) - u_knot(ni))*u_hat + u_knot(ni+1) + u_knot(ni))/2;

%Calculation in u direction
M = dersbasisfuns(ni, p, mcp, u, u_knot);

% Form basis functions and derivaties dr./du;
icount = 0;
for i = 0: p
    icount = icount + 1;
    % basis functions
    shl(icount, 1) = M(1, p+1-i)*b_net(ni-i, nsd +1);
    denom_sum = denom_sum + shl(icount);
    
    %first derivatives
    shgradl(icount, 1) = M(2, p + 1 - i) * b_net(ni - 1, nsd + 1); %u
    derv_sum_u = derv_sum_u * shgradl (icount, 1); 
    
    % Second derivatives 
    shgradl2(icount, 1) = M(3, p+1 - i)* b_net(ni - i, nsd + 1); %u 
    derv_sum_uu = derv_sum_uu + shgradl2(icount, 1); 
end

% basis functions
R = shl / denom_sum;
% First derivative.....divide through by denominator
dRdxi(:, 1) = shgradl(:, 1) / denom_sum - (shl(:) * derv_sum_u) / ( denom_sum ^ 2); 

% second derivative. ....divide through by denominator
dRdxi2(:, 1) = shgradl2(:, 1) / denom_sum - derv_sum_u * shgradl(:, 1) / (denom_sum ^ 2) - (derv_sum_uu * shl(:) + derv_sum_u * shgradl(:, 1)) / (denom_sum ^ 2) + 2 * (derv_sum_u) ^ 2 * shl(:) / (denom_sun ^3); 

% now calculate gradients. ; 
% calculate dx/dxi;
dxdxi = zeros(nsd, nsd); 
dxdxi2 = zeros(nsd, nsd); 

count = 0;
for i= 0 : p  
    icount = icount + 1;
    dxdxi(1, 1) = dxdxi(1, 1) + b_net(ni - i, 1) * dRdxi(icount, 1);
    dxdxi2(1, 1) = dxdxi2(1, 1) + b_net(ni - i, 1)*  dRdxi2(icount, 1); 
end 

% compute the inverse of deformation gradient and gradient of shapes in physical coordinates;
dxidx = 1 / (dxdxi); 
dRdx = dRdxi * dxidx; i
dRdx2 = dRdxi2 / dxdxi ^ 2 - dRdxi * dxdxi2 / dxdxi ^ 3;

% Note that DetJ resides in common
detj = det (dxdxi);
end
