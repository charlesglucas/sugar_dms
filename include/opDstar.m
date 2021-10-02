function y = opDstar(x,param)

if strcmp(param.type,'1D')
if strcmp(param.prior,'gradient')
    y = [-x(2,:); x(2:end-1,:)-x(3:end,:); x(end,:)]./2.;
elseif strcmp(param.prior,'laplacian')
    y = ([x(2,:)/4;-x(2,:)/2 + x(3,:)/4; x(4:end-1,:)/4 - x(3:end-2,:)/2 ...
        + x(2:end-3,:)/4; x(end-2,:)/4 - x(end-1,:)/2;x(end-1,:)/4]);
end
y = squeeze(y);
elseif strcmp(param.type,'2D')
    y1 = [x(:,1,1,:), x(:,2:end-1,1,:)-x(:,1:end-2,1,:), -x(:,end-1,1,:)]./2.;
    y2 = [x(1,:,2,:); x(2:end-1,:,2,:)-x(1:end-2,:,2,:); -x(end-1,:,2,:)]./2.;
    y = -(y1 + y2);
    y = squeeze(y);
end
end 