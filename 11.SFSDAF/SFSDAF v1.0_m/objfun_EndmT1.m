  function f = objfun_EndmT1(x,fraction,y,thresholdMatrix,xH,yH)
  total=0;    
 for i=1:xH
     for j=1:yH
         if thresholdMatrix(i,j)==1
            Temp = x*squeeze(fraction(i,j,:))-y(i,j);
            total = total + Temp.^2;
         end
     end
 end
 f=total;
 