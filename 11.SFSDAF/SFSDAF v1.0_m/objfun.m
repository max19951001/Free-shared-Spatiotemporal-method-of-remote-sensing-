function f = objfun(x,Endmember,PixelData)            
Temp = Endmember*x-PixelData;
 f = Temp'*Temp;
       