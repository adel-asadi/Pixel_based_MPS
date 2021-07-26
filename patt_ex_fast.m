function Tem=patt_ex_fast(TI, T, dim)

if dim==2
    
m=T(1,1); n=T(1,2);
%[m n]=size(T);
%p1=P(1,1); p2=P(1,2);
[k1 k2]=size(TI);
m1=(m-1)/2; n1=(n-1)/2;


Tem=[];
count=0;
for j=1:n
    for i=1:m
        ti=TI(i:i+k1-m,j:j+k2-n);
        Tem=[Tem reshape(ti,numel(ti),1)];
        count=count+1;
    end

end
else
if dim==3
m=T(1,1); n=T(1,2);o=T(1,3);
[k1 k2 k3]=size(TI);
m1=(m-1)/2; n1=(n-1)/2;o1=(o-1)/2;

Tem=[];
count=0;
for k=1:o
    for j=1:n
        for i=1:m
        ti=TI(i:i+k1-m,j:j+k2-n,k:k+k3-o);
        Tem=[Tem reshape(ti,numel(ti),1)];
        count=count+1;
        end
    end
end
end
end

