function [a,b,q] = guidedfilter_MS_low(I, p, r)
%   GUIDEDFILTER_MS   O(1) time implementation of guided filter using a MS image as the guidance.
%
%   - guidance image: I (should be a MS image)
%   - filtering input image: p (should be a gray-scale/single channel image)
%   - local window radius: r
%   - regularization parameter: eps

[hei, wid] = size(p);
dim=size(I,3);
N = boxfilter(ones(hei, wid), r); % the size of each local patch; N=(2r+1)^2 except for boundary pixels.

mean_p = boxfilter(p, r) ./ N;
for i=1:dim
    mean_I(:, :, i)= boxfilter(I(:, :, i), r) ./ N;
    mean_Ip(:, :, i)=boxfilter(I(:, :, i).*p, r) ./ N;
    cov_Ip(:, :, i)= mean_Ip(:, :, i)- mean_I(:, :, i).*mean_p;% covariance of (I, p) in each local patch.
end

% variance of I in each local patch: the matrix Sigma in Eqn (14).
% Note the variance in each local patch is a dim x dim symmetric matrix:
for i=1:dim
    for j=1:dim
        var_I(:,:,(i-1)*dim+j)=boxfilter(I(:, :, i).*I(:, :, j), r) ./ N - mean_I(:, :, i).*mean_I(:, :, j);  
    end
end

a = zeros(hei, wid, dim);
for y=1:hei
    for x=1:wid
        %for k=1:dim  Sigma(k,:)=var_I(y,x,k,:);  end
        Sigma=reshape(var_I(y,x,:),dim,dim);
        cov_Ip0=D3_D2(cov_Ip(y,x,:))';        
        a(y, x, :) = cov_Ip0 * inv(Sigma + eps * eye(dim)); % Eqn. (14) in the paper;
    end
end
b = mean_p - sum(a.* mean_I,3); % Eqn. (15) in the paper;

q=boxfilter(b, r);
for i=1:dim
    q=q+boxfilter(a(:,:,i), r).* I(:, :,i);
end
q=q./N;
end