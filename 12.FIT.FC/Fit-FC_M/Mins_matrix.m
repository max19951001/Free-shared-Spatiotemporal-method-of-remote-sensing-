%%%a function to find the N minmum values in the matrix
function [II,JJ]=Mins_matrix(A,n);
for i=1:n
    [ma, inda] = min(A(:));
    [r,c] = ind2sub(size(A), inda);
    II(i)=r;
    JJ(i)=c;
    A(r,c)=10^10;
end