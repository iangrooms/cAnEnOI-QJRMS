function OUT = RHS(~,Y,p)
% Computes RHS for my new L96 model
% Inputs:
%   p   struct with fields K, J, h, and F
%   Y   variables Y_{j,k}
% Output:
%   OUT dY/dt computed from new model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The following assumes p.K=41, and N=p.J*p.K

YHat = fft(Y(:));
XHat = [YHat(1:21);YHat(end-19:end)]*(1/p.J);
X = real(ifft(XHat));

nX = circshift(X,1).*(circshift(X,-1)-circshift(X,2));
NX = interpft(nX,p.J*41);

OUT = p.h*circshift(Y,-1).*(circshift(Y,1)-circshift(Y,-2)) - Y + p.F + NX;







% Y = reshape(IN,[p.J p.K]);
% X = repmat(mean(Y),p.J,1);
% y = Y-X;
% yv = sum(y.^2);
% 
% Xp  = circshift(X,[0 1]).*(circshift(X,[0 -1])-circshift(X,[0 2])) - X + p.F...
%     + (p.c*p.h/4)*repmat(2*yv+circshift(yv,[0 1])+circshift(yv,[0 -1]),p.J,1);
% 
% y = y(:);X = (2*X + circshift(X,[0 1]) + circshift(X,[0 -1]))/4;
% yp = circshift(y,-1).*(circshift(y,1)-circshift(y,-2)) - y - p.h*X(:).*y;
% yp = p.c*bsxfun(@plus,yp,-mean(yp));
% OUT = yp(:)+Xp(:);
