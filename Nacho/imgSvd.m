B=imread('/Users/isagues/Downloads/mandril_512x512.tiff');  
B=rgb2gray(B); 
doubleB=double(B); 
[U,S,V]=svd(doubleB);  
C=S; 
N=50; 
C(N+1:end,:)=0; 
C(:,N+1:end)=0; 
D=U*C*V'; 
imshow(uint8(D));
[[1:size(diag(S))(1)]' diag(S) cumsum(diag(S))/sum(diag(S))]