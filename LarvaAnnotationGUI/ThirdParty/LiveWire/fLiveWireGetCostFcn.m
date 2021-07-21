function dF = fLiveWireGetCostFcn(dImg, dWz, dWg, dWd)

if nargin < 2,
    dWz = 0.0;
    dWg = 1.0;
    dWd = 0.0;
end

% -------------------------------------------------------------------------
% Calculat the cost function

% The gradient strength cost Fg
dImg = double(dImg);
%[dY, dX] = gradient(dImg);
%dFg = sqrt(dX.^2 + dY.^2);
dFg = dImg;
dFg = 1 - dFg./max(dFg(:));

% The zero-crossing cost Fz
%lFz = ~edge(dImg, 'zerocross');

% The Sum:

%dF = dWz.*double(lFz)+ dWg.*dFg;
%test = 1;
dF = dFg;

%figure; imagesc(dF); colormap gray;
% 
% for i=1:20
%     for j=1:20
%         disp(num2str(dF(i, j)));
%     end
% end

% -------------------------------------------------------------------------