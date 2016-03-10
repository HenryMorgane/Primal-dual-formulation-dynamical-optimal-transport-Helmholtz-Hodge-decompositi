% Primal-dual algorithm on two different staggered grids for the optimal 
% transport problem using Helmholtz-Hodge decomposition 

% Minimization of the functional
% min_(phi) J_h(K(phi)) + G(phi) where J_h(K(phi)) = J(K(phi)+grad(h))

addpath('toolbox/');
clear; close all; clc;

N = 32; P = 32; Q = 32; % size of the grid
d = [N,P,Q];
l = ones(3,1);

% Useful functions
mynorm = @(a)norm(a(:));
sum3 = @(a)sum(a(:));
mymin = @(x)min(min(min(x(:,:,:,1))));
dotp = @(a,b)sum(a(:).*b(:));
%--------------------------------------------------------------------------
% Linear operators
dx = @(V)diff(V,[],1);
dy = @(V)diff(V,[],2);
dz = @(V)diff(V,[],3);
% adjoint operator
dxS = @(V)cat(1,-V(1,:,:), -diff(V,[],1), V(end,:,:));
dyS = @(V)cat(2,-V(:,1,:), -diff(V,[],2), V(:,end,:));
dzS = @(V)cat(3,-V(:,:,1), -diff(V,[],3), V(:,:,end));

K  = @(X)pd_operator(dxS,dyS,dzS,X, +1);
KS = @(X)pd_operator(dxS,dyS,dzS,X, -1);

% Test for adjointness of the operators
U1 = rand(N,P,Q,3);
U2 = staggereds2(d); U2.M{1} = rand(N,P+1,Q+1); U2.M{2} = rand(N+1,P,Q+1); U2.M{3} = rand(N+1,P+1,Q);
if abs(dotp( U1, K(U2) ) - dotp_stag(KS(U1), U2))>1e-5
    warning('Adjointness problem');
    abs(dotp( U1, K(U2) ) - dotp_stag(KS(U1), U2))
end
%%--------------------------------------------------------------------------
% Boundary conditions
[Y,X] = meshgrid(linspace(0,1,P), linspace(0,1,N));

normalize = @(u)u/sum(u(:));
rho = 1e-12; % minimum density value
gaussian = @(a,b,sigmag)exp( -((X-a).^2+(Y-b).^2)/(2*sigmag^2) );
sigmag = 0.075;

rho0 = normalize( rho + gaussian(.35,.35,sigmag) );
rho1 = normalize( rho + gaussian(.65,.65,sigmag) );
epsilon = min(rho0(:));
 
% Display of rho_0 and rho_1
figure(1)
subplot(1,2,1); imageplot(rho0)
title('\rho_0')
subplot(1,2,2); imageplot(rho1)
title('\rho_1')
%--------------------------------------------------------------------------
% Poisson equation
u3Dbc = cat(3,-rho0, zeros(N,P,Q-2), rho1);
h=poisson3d_Neumann(u3Dbc*d(3)/l(3),l(1),l(2),l(3));
%--------------------------------------------------------------------------
% grad(h)
[gh1b,gh2b,gh3b] = Grad(h);

% Boundary conditions for grad(h)
gh = staggered(d);
gh.M{3} = cat(3,rho0,gh3b,rho1);
gh.M{2} = cat(2,zeros(N,1,Q),gh2b,zeros(N,1,Q));
gh.M{1} = cat(1,zeros(1,P,Q),gh1b,zeros(1,P,Q));
gradh = interp(gh);

%--------------------------------------------------------------------------
% Initial conditions
phi = staggereds2(d);
z = K(phi);
phit = phi;

% Operator norm
u = rand([d,3]);
u = u*(1/mynorm(u)); e = [];
for i=1:30
    v = K(KS(u));
    e(end+1) = dotp(u,v);
    u = v * (1/mynorm(v));
end

% Parameters
L = e(end);          % operator norm 
itmax = 10;        % nb of iterations
sigma = 100;         % first step of the algorithm
tau = 0.99/(sigma*L);% second step of the algorithm
theta = 1.;          % third step of the algorithm


% Memory
Jh = zeros(1,itmax);
Min = zeros(1,itmax);
divrhom = zeros(1,itmax);

% Primal-dual loop
tic
for i=1:itmax
    
    % Progress
    if rem(i,itmax/10)==0
        progressbar(i,itmax);
    end
    
    % (rho,m) on the centered grid
    rhom = K(phi) + gradh;
    
    % Functional, min(rho) and div(rho,m)
    Jh(i) = sum3(  sum(rhom(:,:,:,1:2).^2,4) ./ (max((rhom(:,:,:,3)),max(epsilon,1e-10)))  );
    Min(i) = mymin(rhom);
    % (rho,m) on the staggered grid
    rhoms = curls2(phi) + gh;
    divrhom(i) = Div(rhoms);
    
    % phi^{k-1}
    phim = phi;
    
    % z^{k+1}
    z = proxJ3D(z+sigma*(K(phit) + gradh));
    
    % phi^{k+1}
    phi = zeroboundary(phi - tau*KS(z),1);
    
    % tilde{phi^{k+1}}
    phit = phi + theta*(phi-phim);

end
toc

% Display of rho(.,t)
figure(2)
for i=1:Q
    imageplot(reshape(rhom(:,:,i,3),N,P));
    pause(0.2)
end

% Display
figure(3)
subplot(1,3,1); loglog(Jh)
title('Functional J_h(K(phi))')
xlabel('iterations i')
subplot(1,3,2); loglog(Min)
title('Minimal of \rho^i')
xlabel('iterations i')
subplot(1,3,3); loglog(divrhom)
title('Divergence of (\rho,m)')
xlabel('iterations i')
