
function endmember_T1 = Estimate_endmember_T1(FR_T1_Probability,FR_T1,xH,yH,class,bands)

endmember_T1 = zeros(bands,class);% 
k=300;
location_fraction=zeros(xH,yH,class);

for c=1:class
    
    A=reshape(FR_T1_Probability(:,:,c),xH*yH,1);
    [fraction_sort, D]= sort(A,'descend');
    for n=1:k
        x_tmp=mod(D(n),xH);
        y_tmp=ceil(D(n)/xH);
        if x_tmp==0
            x_tmp=xH;
        end
        if fraction_sort(n)>0.01 
            location_fraction(x_tmp,y_tmp,c) = 1;
        end
    end
end
         
thresholdMatrix=zeros(xH,yH);
for c=1:class
    thresholdMatrix = thresholdMatrix+location_fraction(:,:,c);
end
thresholdMatrix=thresholdMatrix>0;

for b=1:bands
    
   x0 = ones(1,class)/class; 
   lb = ones(1,class) * min(min(FR_T1(:,:,b)));
   ub = ones(1,class) * max(max(FR_T1(:,:,b)));  
   options = optimoptions('fmincon','Display','off');  
   [x,fval] = fmincon(@(x) objfun_EndmT1(x,FR_T1_Probability,FR_T1(:,:,b),thresholdMatrix,xH,yH),x0,[],[],[],[],lb,ub,[],options);
   endmember_T1(b,:) = x;

end


