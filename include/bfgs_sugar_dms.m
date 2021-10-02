function  [lambda,crit] = bfgs_sugar_dms(data, param, choice)
    % Automated tuning of hyperparamters
    % Performs a quasi-Newton descent minimizing SURE
    %
    % Implementation C.G. LUCAS and B. PASCAL, ENS Lyon
    % April 2021
    
    % default parameters
    if nargin == 1, param = struct; end
    if nargin < 3, choice = struct; end
    if ~isfield(choice,'prior'), choice.prior = 'gradient'; end
    if ~isfield(choice,'norm'), choice.norm = 'L1'; end
    if ~isfield(choice,'edges'), choice.edges = 'similar'; end
    if ~isfield(choice,'type'), choice.type = '2D'; end
    if ~isfield(param,'R'), param.R = 1; end
    
    param.type = choice.type;
    param.prior = choice.prior;
    
    % Evaluate the standard deviation of the noise
    if strcmp(param.type,'1D')
        Nl = size(data,1);
        if isfield(param,'sigma')
            sigma = param.sigma;%.*ones(1,Nl);
        else
            sigma = zeros(1,Nl);
            for nl = 1:Nl
                [~,cH] = dwt(data(nl,:),'db1');
                C = abs(cH);
                sigma(nl) = median(C)/0.6745;
            end
        end
    elseif strcmp(param.type,'2D')
        data = ipermute(data,[2,3,1]);
        if isfield(param,'sigma')
            sigma = param.sigma;
        else
            Nl = size(data,1);
            sigma = zeros(1,Nl);
            for nl = 1:Nl
                [~,cH,cV,cD] = dwt2(squeeze(data(nl,:,:)),'db1');
                C = abs([cH(:) ; cV(:); cD(:)]);
                sigma = median(C)/0.6745;
            end
        end
    else
        error('Only (multivariate) signals of type 1D or images of type 2D allowed')
    end

    % Initialize the quasi-Newton algorithm
    if isfield(param,'fdmc')
        fdmc = param.fdmc;
    else
        fdmc.eps = 2*max(sigma(:))./numel(data(1,:,:))^.3;
        if strcmp(param.type,'1D')
            fdmc.delta = randn(size(data'));
        elseif strcmp(param.type,'2D')
            s = size(ipermute(data,[3,1,2]));
            if length(s)==2, s = [s 1]; end
            fdmc.delta = randn([s param.R]);
        end
        fdmc.sigma = sigma';
    end
    
    % Initial hyperparameter
    lambda_init = zeros(2,1);
    
    param.prior = choice.prior;
    D = @(z) opD(z,param);
    
    if isfield(param,'lambda')
        lambda_init(1) = param.lambda(1);
        lambda_init(2) = param.lambda(2);
    else
        if strcmp(param.type,'1D')
            Var = D(data');
        elseif strcmp(param.type,'2D')
            Var = D(ipermute(data,[3,1,2]));
        end
        lambda_init(1) = numel(data)*sigma(1)^2/(4*sum(Var(:).^2));
        lambda_init(2) = lambda_init(1)*sum(Var(:).^2)/(2*numel(data));
    end
    
    opts.x0 = lambda_init;
    % Initial Hessian 
    set_init(struct);
    set_crits(struct);
    [~,sugar] = sure_sugar_dms(data,choice,fdmc,lambda_init);
    if nargin==3 && isfield(param,'kappa')
        kappa = param.kappa;
    else
        kappa = 0.9;
    end
    opts.H0 = diag(abs(kappa*lambda_init)./abs(sugar));
    
    % Define SURE and SUGAR as a function
    sure_sugar_fun = @(lambda) sure_sugar_dms(data, choice, fdmc, lambda);
    
    % Run BFGS GRANSO algorithm
    opts.print_level = 1;
    opts.quadprog_info_msg = false;
    opts.prescaling_info_msg = false;
    opts.maxit = 20;
    opts.maxclocktime = 10*60;
    opts.halt_on_linesearch_bracket = true;
    
    soln = granso(2,sure_sugar_fun,@positivityConstraint,[],opts);
    lambda = soln.final.x;
    
    %init_DMS = get_init;
    %x = init_DMS.x;
    crit_DMS = get_crits;
    crit_DMS.lambda = cell2mat(crit_DMS.lambda);
    crit = crit_DMS;
    
    % Define positivity constraints
    function [ci,ci_grad] = positivityConstraint(lambda)
        % Impose that lambda >= 1e-2 * lambda_in (so that lambda > 0)
        ci      = 1e-2 * lambda_init - lambda;
        ci_grad = -eye(2);
    end
end