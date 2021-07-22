function ClassPoss = SpatialInterpolation(ClassFraction,Factor,Method)

[m,n] = size(ClassFraction);
Fraction = zeros(m+4,n+4);
Fraction(3:m+2,3:n+2)=ClassFraction;

[UsedLocationX,UsedLocationY] = meshgrid(Factor/2.0:Factor:Factor/2.0+(n+3)*Factor,Factor/2.0:Factor:Factor/2.0+(m+3)*Factor);        %%meshgrid(1:3,10:14)
[PredLocationX,PredLocationY]   = meshgrid(1.5+2*Factor:1:0.5+(n+2)*Factor,1.5+2*Factor:1:0.5+(m+2)*Factor);  

%% *************
if strcmp(Method,'bilinear')==1
ClassPoss = interp2(UsedLocationX,UsedLocationY,Fraction,PredLocationX,PredLocationY,'linear');
end
%% *************

%% *************
if strcmp(Method,'cubic')==1
ClassPoss = interp2(UsedLocationX,UsedLocationY,Fraction,PredLocationX,PredLocationY,'cubic');
end
%% *************

%% *************
if strcmp(Method,'spline')==1
ClassPoss = interp2(UsedLocationX,UsedLocationY,Fraction,PredLocationX,PredLocationY,'spline');
end
%% *************

%% *************
if strcmp(Method,'idw')==1
    ClassPoss = gIDW(UsedLocationX(:),UsedLocationY(:),Fraction(:),PredLocationX,PredLocationY,-2,'n',12);
end
%% *************

