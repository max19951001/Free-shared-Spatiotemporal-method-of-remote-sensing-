function nX = cml_fusion(LR_predictionDate, LR_earlierDate, HR_image, patchSize, useOverlap, shift, reg_param)

[lm,ln,b] = size(LR_predictionDate);
[m,n,~] = size(HR_image);
nX = zeros(m,n,b);

scale = fix(m/lm);
L = patchSize*scale;
N1 = fix(lm/patchSize);
N2 = fix(ln/patchSize);

for i=1:N1
    ir = (1:patchSize) + (i-1)*patchSize;
    if i==N1, ir = (i-1)*patchSize+1:lm; end
    hir = (1:L) + (i-1)*L;
    if i==N1, hir = (i-1)*L+1:m; end
    for j=1:N2
        jr = (1:patchSize) + (j-1)*patchSize;
        if j==N2, jr = (j-1)*patchSize+1:ln; end
        hjr = (1:L) + (j-1)*L;
        if j==N2, hjr = (j-1)*L+1:n; end
        
        x = LR_predictionDate(ir,jr,:);
        xc = LR_earlierDate(ir,jr,:);
        hrcolor = HR_image(hir,hjr,:);
        nX(hir,hjr,:) = cm_fusion(x, xc, hrcolor, reg_param);
    end
end
if useOverlap
    for xx = 1:shift:lm
        for yy = 1:shift:ln
            
            if(xx>=patchSize && yy>=patchSize && xx<lm-patchSize && yy<ln-patchSize)
                
                ir = xx-patchSize+1:xx+patchSize;
                jr = yy-patchSize+1:yy+patchSize;
                x = LR_predictionDate(ir,jr,:);
                xc = LR_earlierDate(ir,jr,:);
                hrcolor = HR_image(ir,jr,:);
                transform = cm_fusion(x,xc,hrcolor, reg_param);
                nX(xx:xx+1, yy:yy+1,1)=transform(patchSize:patchSize+1, patchSize:patchSize+1,1);
                
            end
        end
    end
end