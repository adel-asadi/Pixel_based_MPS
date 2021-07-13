function [weight] = weight_gen(m,w_min,w_max)
s=m^2;
g=round(s/2);
gg=[w_min:((w_max-w_min)/(g-1)):w_max];
ggg=[w_max:-((w_max-w_min)/(g-2)):w_min];
gggg=[gg,ggg];
ss=sqrt(s);
ggggg=reshape(gggg,[ss ss]);
sss=zeros(ss,ss);
for i=1:ss
    for j=1:ss
        sss(i,j)=ggggg(j,i);
    end
end
k=zeros(ss,ss);
for i=1:ss
    for j=1:ss
        kk=sss(i,j);
        kkk=ggggg(i,j);
        k(i,j)=min(kk,kkk);
    end
end
weight=k(:);
end

