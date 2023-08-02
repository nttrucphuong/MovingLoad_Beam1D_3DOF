% CtrlPts = zeros(4, 2);
% 
% CtrlPts(1:3, 1) = [0; 0; 0];
% CtrlPts(1:3, 2) = [L; 0; 0];
% 
% CtrlPts(4, :) = 1;
% KntVect{1} = [0, 0, 1, 1];
% Dim = numel(KntVect);
% W = CtrlPts(4, :, :, :);
% CtrlPts3D = bsxfun(@rdivide, CtrlPts(1 : 3, :, :, :), W);
% NPts = size(CtrlPts);
% NCtrlPts = zeros(1, Dim);
% p = zeros(1, Dim);
% uqKntVect = cell(1, Dim);
% for i = 1 : Dim
%     NCtrlPts(i) = NPts(i + 1);
%     p(i) = numel(KntVect{i}) - NCtrlPts(i) - 1;
% end
q = x^2+3;
df = @(x) q;
for i = 1:3
    x = i;
    disp(i);
end
