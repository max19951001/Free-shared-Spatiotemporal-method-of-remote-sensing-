function Z_extend=Extend_cube(Z,w);
[a,b,c]=size(Z);
for i=1:c
    Z_extend1=[repmat(Z(:,1,i),[1,w]),Z(:,:,i),repmat(Z(:,end,i),[1,w])];%%%extend columns
    Z_extend(:,:,i)=[repmat(Z_extend1(1,:),[w,1]);Z_extend1;repmat(Z_extend1(end,:),[w,1])];%%%extend rows
end