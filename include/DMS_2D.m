function [ui,ei,crit] = DMS_2D(z,Lambda,choice)

%This algorithm minimizes
%
%           1/2*||z-u||^2 + beta*||(1-e).*Du||^2 + lambda*R(e)
%
%inputs:
%   - z in R^{n1xn2xm}: a noisy observation which is m components of
%   2D-grids n1xn2
%   - beta: hyperparameter associated to smoothness
%   - lambda: hyperparameter associated to contour detection 
%   - choice: structure to select the parameters of the problem
%       - edges: 'similar' (by default) for joint contour across components
% or 'distinct' for one contour per component
%   	- norm: 'L1' (by default) for l1-norm
%       - type: '2D' (by default) or '1D'
%   	- prior: 'gradient' (by default) or 'Laplacian'
%returns:
%   - u in R^{n1xn2xm}: the denoised image
%   - e in R^{n1xn2xM}: the contour where M is m for distinct edges and 1
%       otherwise
%
% Implementation Charles-Gérard LUCAS, ENS Lyon
% from 
% M. Foare, N. Pustelnik, and L. Condat: Semi-linearized
% proximal alternating minimization for a discrete Mumford–Shah model.
% IEEE Transactions on Image Processing. (2019)
% April 2021

% default parameters
if nargin == 2, choice = struct; end
if ~isfield(choice,'prior'), choice.prior = 'gradient'; end
if ~isfield(choice,'norm'), choice.norm = 'L1'; end
if ~isfield(choice,'edges'), choice.edges = 'similar'; end
if ~isfield(choice,'type'), choice.type = '2D'; end


beta = Lambda(1);lambda = Lambda(2);
% dim
[n1,n2,m] = size(z);

% version
if strcmp(choice.edges,'similar')
    l = 1;
elseif strcmp(choice.edges,'distinct')
    l = m;
end

% parameters
normD = sqrt(2);
ck = 1.01*beta*normD^2; 
dk = ck/1000;

% functionals
[L_prox,S_grad,R_prox,Psi] = functionals(beta,lambda,choice);

% initialization
u0 = z;
e0 = ones(n1,n2,2,l);
Psi_2 = Psi(u0,e0,z);
ui = u0;
ei = e0;
gap = 1;

k = 1;
while (gap > 1e-4)
    Psi_1 = Psi_2;
    
    % update u
    ui = L_prox(ui - beta*S_grad(ei,ui)/ck,z,1/ck);
      
    % update e
    if strcmp(choice.edges,'similar')
        Du2 = sum(D(ui).^2,4);
    elseif strcmp(choice.edges,'distinct')
        Du2 = D(ui).^2;
    end
    num = beta*Du2 + dk.*ei/2.;
    den = beta*Du2 + dk/2.;
    if strcmp(choice.edges,'similar')
        tau = lambda./(2*den);
    elseif strcmp(choice.edges,'distinct')
        tau = lambda./(2*den);
    end
    ei  = R_prox(num./den,tau);

    % compute error
    Psi_2 = Psi(ui,ei,z);    
    gap = abs(Psi_2 - Psi_1);
    crit.obj(k) = Psi_2;
    crit.gap(k) = abs(Psi_2 - Psi_1);
    k = k+1;
end

end

% functionals
function [L_prox,S_grad,R_prox,Psi] = functionals(beta,lambda,choice)
sum_all = @(M) sum(M(:));

L = @(u,z) sum_all((u-z).^2)/2;
L_prox = @(u,z,tau) (u+tau.*z)./(1+tau);

S = @(e,u) sum_all((repmat((1-e),[1 1 1 size(u,3)/size(e,4)]).*D(u)).^2);
S_grad = @(e,u) 2*Dstar(repmat((1-e).^2,[1 1 1 size(u,3)/size(e,4)]).*D(u));

if strcmp(choice.norm,'L1')
R = @(e) sum_all(abs(e));
R_prox = @(eta,tau) eta .* max(0,1 - tau./abs(eta));
elseif strcmp(choice.norm,'L12')
R = @(e) sum_all(sqrt(sum(e.^2,4)));
R_prox = @(eta,tau) eta .* max(0,1 - tau./sqrt(sum(eta.^2,4)));
end

Psi = @(u,e,z)  L(u,z) + beta*S(e,u) + lambda*R(e);
end

% finite difference operator
function y = D(u)
[n1,n2,m] = size(u);
y = zeros(n1,n2,2,m);
y(:,:,1,:) = [u(:,2:end,:)-u(:,1:end-1,:),zeros(n1,1,m)]./2.;
y(:,:,2,:) = [u(2:end,:,:)-u(1:end-1,:,:);zeros(1,n2,m)]./2.;
y = squeeze(y);
end 

function y = Dstar(x)
y1 = [x(:,1,1,:), x(:,2:end-1,1,:)-x(:,1:end-2,1,:), -x(:,end-1,1,:)]./2.;
y2 = [x(1,:,2,:); x(2:end-1,:,2,:)-x(1:end-2,:,2,:); -x(end-1,:,2,:)]./2.;
y = -(y1 + y2);
y = squeeze(y);
end 