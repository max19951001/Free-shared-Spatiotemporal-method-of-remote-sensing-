function Z_extend=Extend_plane(Z,w);
Z_extend1=[repmat(Z(:,1),[1,w]),Z,repmat(Z(:,end),[1,w])];%%%extend columns
Z_extend=[repmat(Z_extend1(1,:),[w,1]);Z_extend1;repmat(Z_extend1(end,:),[w,1])];%%%extend rows