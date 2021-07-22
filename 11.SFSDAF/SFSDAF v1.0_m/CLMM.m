function [Fraction] = CLMM(Class,Endmember,PixelData,n_iteration)
 
        x0 = ones(Class,1)*1.0/Class;    
        Aeq = ones(1,Class);
        beq = 1;
        lb = zeros(Class,1);
        ub = ones(Class,1);       
       
        options = optimoptions('fmincon','Display','off');
        [x,fval] = fmincon(@(x) objfun(x,Endmember,PixelData),x0,[],[],Aeq,beq,lb,ub,[],options);
    
        Fraction = x;

 


