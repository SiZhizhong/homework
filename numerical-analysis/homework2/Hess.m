A=zeros(10,10)
for i=1:10
    for j=1:10
        if(i==j)
            A(i,j)=1.52*cos(i+1.2*j)
        else
            A(i,j)=sin(0.5*i+0.2*j)
        end
    end
end
[p,hessa]=hess(A);
[Q,R]=qr(hessa);
[V,D]=eig(A);
            