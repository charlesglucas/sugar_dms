function dprox = deprox(e,de,tau,type)
    
    if strcmp(type,'1D')
        dim = 2;
    elseif strcmp(type,'2D')
        dim = 4;
    end

    dot = @(e,de) sum(de.*e,dim);
    sigx = abs(e);
    sigx(sigx==0) = tau(sigx==0);
    dprox = (de - tau./sigx.*(de-dot(de,e)./sigx.^2.*e)).*(sigx>tau);
end