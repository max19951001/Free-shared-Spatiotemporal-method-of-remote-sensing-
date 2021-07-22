  function f = objfun_deltaEndm(x,fraction,y_delta,thresholdMatrix,xL,yL)
  total=0;    
 for i=1:xL
     for j=1:yL
         if thresholdMatrix(i,j)==1
            Temp = x*squeeze(fraction(i,j,:))-y_delta(i,j);
            total = total + Temp.^2;
         end
     end
 end
 f=total  ;
 