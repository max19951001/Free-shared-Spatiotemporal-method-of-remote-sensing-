function [out, outVec] = RMSE(ref,tar)
%--------------------------------------------------------------------------
% Root mean squared error (RMSE)
%
% USAGE
%   out = RMSE(ref,tar)
%
% INPUT
%   ref : reference HS data (rows,cols,bands)
%   tar : target HS data (rows,cols,bands)
%
% OUTPUT
%   out : RMSE (scalar)
%
%--------------------------------------------------------------------------
[rows,cols,bands] = size(ref);
outVec = zeros(1,bands);

for i = 1:bands
    outVec(1,i) = sum(sum((tar(:,:,i) - ref(:, :, i)).^2)/(rows*cols)).^0.5;
end

out = (sum(sum(sum((tar-ref).^2)))/(rows*cols*bands)).^0.5;