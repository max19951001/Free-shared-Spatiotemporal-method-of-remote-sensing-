function [fraction]=CLMM_Unmix(y,meanV,class,bands,xL,yL,n_iteration)
 fraction=zeros(xL,yL,class); 
 for i=1:xL
     for j=1:yL         
         fraction(i,j,:) = CLMM(class,meanV,squeeze(y(i,j,:)),n_iteration);
     end
 end