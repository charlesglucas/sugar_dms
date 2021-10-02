function y = opD(u,param)

if strcmp(param.type,'1D')
    if strcmp(param.prior,'gradient')
        y = [zeros(1,size(u,2));u(2:end,:)-u(1:end-1,:)]./2.;
    elseif strcmp(param.prior,'laplacian')
        y = [zeros(1,size(u,2));u(3:end,:)/4 - u(2:end-1,:)/2 + u(1:end-2,:)/4;zeros(1,size(u,2))];
    end
    y = squeeze(y);
elseif strcmp(param.type,'2D')
    [n1,n2,m] = size(u);
    y = zeros(n1,n2,2,m);
    y(:,:,1,:) = [u(:,2:end,:)-u(:,1:end-1,:),zeros(n1,1,m)]./2.;
    y(:,:,2,:) = [u(2:end,:,:)-u(1:end-1,:,:);zeros(1,n2,m)]./2.;
    y = squeeze(y);
end
end 