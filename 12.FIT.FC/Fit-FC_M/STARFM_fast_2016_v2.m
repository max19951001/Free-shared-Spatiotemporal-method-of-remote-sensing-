%%%% Input: 
%%%% M0: available coarse MODIS image for predicted time, for only one band
%%%% M_cube and L_cube: available temporal neighboring image pairs of the same spatial size, for only one band
%%%% LR_cube and LNIR_cube: R and NIR images of Landsat
%%%% w: neignborhood window size (2*w+1)
%%%% N_S: number of similar pixels
%%%% A: a parameter defining the importance of spatial distance

%%%%Output: L0

function [L0,Lengthx]=STARFM_fast_2016_v2(M_cube,L_cube,M0,LR_cube,LNIR_cube,w,N_S,A);
[a,b,n]=size(L_cube);

M_cube_extend=Extend_cube(M_cube,w); clear M_cube;
L_cube_extend=Extend_cube(L_cube,w); clear L_cube;
M0_extend=Extend_plane(M0,w); clear M0;
%%%%red and NIR bands of Landsat (4 5 for Landsat 8 and 3 4 for TM/ETM+)
LR_cube_extend=Extend_cube(LR_cube,w); clear LR_cube;
LNIR_cube_extend=Extend_cube(LNIR_cube,w); clear LNIR_cube;

Simulated=zeros(a,b);
[a0,b0]=find(Simulated==0);La=length(a0);
aa0=a0+w;bb0=b0+w;

L0=zeros(a+2*w,b+2*w);
D=zeros(2*w+1,2*w+1);
for i=1:2*w+1
    for j=1:2*w+1
        D(i,j)=norm([i,j]-[w+1,w+1]);
    end
end
D=1+D/A;

for t=1:La
    i=aa0(t);
    j=bb0(t);
    W=zeros(2*w+1,2*w+1,n);%%%%weights for a (2*w+1)*(2*w+1) window of n images
    WC=zeros(2*w+1,2*w+1,n);
    
    Neighbors_L=L_cube_extend(i-w:i+w,j-w:j+w,:);
    Neighbors_M=M_cube_extend(i-w:i+w,j-w:j+w,:);
    Neighbors_M0=M0_extend(i-w:i+w,j-w:j+w);
    for k=1:n
        %%%%find similar pixels in the kth Landsat image
        Neighbors_LR=LR_cube_extend(i-w:i+w,j-w:j+w,k);
        Neighbors_LNIR=LNIR_cube_extend(i-w:i+w,j-w:j+w,k);
        Dif_R=abs(Neighbors_LR-LR_cube_extend(i,j,k)*ones(2*w+1,2*w+1));Dif_NIR=abs(Neighbors_LNIR-LNIR_cube_extend(i,j,k)*ones(2*w+1,2*w+1));
        Dif=Dif_R.^2+Dif_NIR.^2;
        [Location_II,Location_JJ]=Mins_matrix(Dif,N_S);

        S=abs(Neighbors_L(:,:,k)- Neighbors_M(:,:,k))+0.1^10;
        T=abs(Neighbors_M0- Neighbors_M(:,:,k))+0.1^10;
        
        for m=1:N_S
            ii=Location_II(m);
            jj=Location_JJ(m);
            W(ii,jj,k)=1/(D(ii,jj));
            %W(ii,jj,k)=1/(D(ii,jj)*S(ii,jj)*T(ii,jj));
            WC(ii,jj,k)=T(ii,jj);
        end
        Sum_WCk(k)=1/sum(sum(WC(:,:,k)));
    end
    Sum_W=sum(sum(D3_D2(W)));
    Sum_WC=sum(Sum_WCk);
    
    for k=1:n
        L0(i,j)=L0(i,j)+Sum_WCk(k)*L_cube_extend(i,j,k)/Sum_WC;
        for m=1:N_S
            ii=Location_II(m);
            jj=Location_JJ(m);
            L0(i,j)=L0(i,j)+W(ii,jj,k)*( Neighbors_M0(ii,jj)-Neighbors_M(ii,jj,k))/Sum_W;
        end
    end
end
L0=L0(w+1:end-w,w+1:end-w);