function [sure,sugar,x] = sure_sugar_dms(data, choice, fdmc, lambda)
% Compute SURE and SUGAR from iteratively differentiated algorithm dslpalm
% for fixed regularization parameters lambda
% 
% Implementation C.G. LUCAS and B. PASCAL, ENS Lyon
% from 
% - Deledalle, C.-A. and Vaiter, S. and Fadili, J. and Peyré, G.: Stein
% Unbiased GrAdient estimator of the Risk (SUGAR) for multiple parameter 
% selection. SIAM J. Imaging Sci. (2014)
% and 
% -  C.-G. Lucas, B. Pascal, N. Pustelnik and P. Abry: Hyperparameter
% selection for the Discrete Mumford-Shah functional (2021)
% April 2021
%
%inputs:
%   - choice: structure to select the parameters of D-MS
%       - edges: 'similar' (by default) for joint contour across components
% or 'distinct' for one contour per component
%   	- norm: 'L1' (by default) for l1-norm
%       - type: '2D' (by default) or '1D'
%   	- prior: 'gradient' (by default) or 'Laplacian'
%   - param: structure to select the parameters of SUGAR D-MS
%       - edges: 'similar' (by default) for joint contour across components
% or 'distinct' for one contour per component
%   	- norm: 'L1' (by default) for l1-norm
%       - type: '2D' (by default) or '1D'
%   	- prior: 'gradient' (by default) or 'Laplacian'

av_sure = 0;
av_sugar = zeros(numel(lambda),1);

% Run differentiated discrete Mumford-Shah
if strcmp(choice.type,'1D')
    data = data';
elseif strcmp(choice.type,'2D')
    data = ipermute(data,[3,1,2]);
end

for i = 1:size(fdmc.delta,4)
    
    delta = fdmc.delta(:,:,:,i);
    [x,dx,crit] = dDMS(data,lambda(1),lambda(2),choice);
    [Ex,Edx,~] = dDMS(data + fdmc.eps*delta,lambda(1),lambda(2),choice);  
    param.lambda = lambda';
    
    % Compute SURE
    if numel(fdmc.sigma) == 1
        sure = sum((x(:)-data(:)).^2) + 2*fdmc.sigma.^2*sum((Ex(:)-x(:)).*delta(:))/fdmc.eps - fdmc.sigma.^2*numel(data);
    end  
    
    % Compute SUGAR
    sugar = zeros(numel(param.lambda),1);
    for nlbd = 1:numel(param.lambda)
        if numel(fdmc.sigma) == 1
            sugar(nlbd) = 2*sum(dx{nlbd}(:).*(x(:)-data(:))) + 2*fdmc.sigma^2*sum((Edx{nlbd}(:)-dx{nlbd}(:)).*delta(:))/fdmc.eps;
         end
    end
    
    % Averaged SURE and SUGAR
    av_sure = av_sure + sure;
    for nlbd = 1:numel(param.lambda)
        av_sugar(nlbd) = av_sugar(nlbd) + sugar(nlbd);
    end   
end    

    sure = av_sure/size(fdmc.delta,4);
    for nlbd = 1:numel(param.lambda)
        sugar(nlbd) = av_sugar(nlbd)/size(fdmc.delta,4);
    end
    
    crit_DMS = get_crits;
    if isfield(crit_DMS,'neval')
        crit_DMS.neval = crit_DMS.neval + 1;
    else
        crit_DMS.neval = 1;
    end
    
    crit_DMS.objective{crit_DMS.neval} = crit.obj;
    crit_DMS.gap{crit_DMS.neval} = crit.gap;
    crit_DMS.lambda{crit_DMS.neval} = lambda;
    crit_DMS.sure(crit_DMS.neval) = sure;
    crit_DMS.sugarnorm(crit_DMS.neval) = norm(sugar,'fro');
    set_crits(crit_DMS);
end