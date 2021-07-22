function [result] =  Downscale_FractionChangeImages(CR_T1,CR_T2,FR_T1,bands,xH,yH,Factor)

xL = xH/Factor;
yL = yH/Factor;
%%%%%%%%%%%%%%%%%%%%%%%%%
FR_T1_ext = zeros(xH+4*Factor,yH+4*Factor,bands);      
FR_T1_ext(2*Factor+1:xH+2*Factor,2*Factor+1:yH+2*Factor,:)=FR_T1; 
CR_fractionDiff = (CR_T2 -  CR_T1);
FR_fractionDiff=zeros(xH,yH,bands);
for b=1:bands
      FR_fractionDiff(:,:,b) = SpatialInterpolation(CR_fractionDiff(:,:,b),Factor,'cubic');
end
DeltaF =  FR_fractionDiff;
DeltaF_ext = zeros(xH+4*Factor,yH+4*Factor,bands);      
DeltaF_ext(2*Factor+1:xH+2*Factor,2*Factor+1:yH+2*Factor,:)=DeltaF;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 20;
w = 1*Factor;
Matrix = zeros(Factor+1,Factor+1);
for i=1:Factor+1
    for j=1:Factor+1
         Matrix(i,j) =  sqrt(  (i-(Factor/2+1))^2 + (j-(Factor/2+1))^2  ) ;
    end
end

Dk = zeros(Factor+1,Factor+1);
for i=1:Factor+1
    for j=1:Factor+1
         Dk(i,j) =  1 + Matrix(i,j)/ (w/2);
    end
end
D=Dk;
ans = 1./D;
S = zeros(2*Factor+1,2*Factor+1);
T = zeros(2*Factor+1,2*Factor+1);
d = zeros(2*Factor+1,2*Factor+1);

predict = zeros(xH,yH,bands);
for i=2*Factor+1:2*Factor+xH
    for j=2*Factor+1:2*Factor+yH
        SD = FR_T1_ext( i-Factor/2:i+Factor/2,j-Factor/2:j+Factor/2,:  );
        for b=1:bands
            SD(:,:,b)  = SD(:,:,b) - FR_T1_ext(i,j,b) +0.0000000001;
        end
        SD=abs(SD);
        Sk = sum(SD,3);
        for b=1:bands        
            Sk_row = reshape(Sk,(Factor+1)*(Factor+1),1);
            tmp = sort(Sk_row);
            threhold_n = tmp(n+1);            
            SD_map = Sk<=threhold_n;
            C = D ;
            W = (  (1./C.* SD_map) / sum(sum(1./C.* SD_map))  );       
            predict(i-2*Factor,j-2*Factor,b) =   sum(sum(DeltaF_ext(i-Factor/2:i+Factor/2,j-Factor/2:j+Factor/2,b) .* W));            
        end
    end
end
result = predict;
