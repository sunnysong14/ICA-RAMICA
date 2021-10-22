% o Liyan adds comments on Jan.6th 2016 Wednesday
% o codes from Haiping 
% 
% USAGE:
%   ???
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [posIdx] = ftrSel(vecYps, gndTX, sup)
% Take the projectins for training data and label to produce sorted basis!!!! 
% A general function, after projecting all,take traing idx

vecDim = size(vecYps, 1);

if sup == 0 
    % Sort by Variance
    % =======================================
    TVars=diag(vecYps*vecYps');
    [stTVars,posIdx]=sort(TVars,'descend');
else
    % Sort according to Fisher's discriminality
    % ==========================================
    classLabel = unique(gndTX);
    nClass = length(classLabel);    % Number of classes
    ClsIdxs = cell(nClass);
    Ns = zeros(nClass,1);
    for i = 1 : nClass
        ClsIdxs{i} = find(gndTX == classLabel(i));
        Ns(i) = length(ClsIdxs{i});
    end
    Ymean = mean(vecYps, 2);
    TSW = zeros(vecDim, 1);
    TSB = zeros(vecDim, 1);
    for i = 1 : nClass
        clsYp = vecYps(:, ClsIdxs{i});
        clsMean = mean(clsYp, 2);
        FtrDiff = clsYp - repmat(clsMean, 1, Ns(i));
        TSW = TSW + sum(FtrDiff .* FtrDiff, 2);
        meanDiff = clsMean - Ymean;
        TSB = TSB + Ns(i) * meanDiff .* meanDiff;
    end
    FisherRatio = TSB ./ TSW;
    [stRatio, posIdx] = sort(FisherRatio, 'descend');
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
