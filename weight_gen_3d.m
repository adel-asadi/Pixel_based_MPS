function [weight] = weight_gen_3d(m,w_min,w_max)
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
h=zeros(m,m,m);
hh=zeros(m,m,m);
hhh=zeros(m,m,m);
hhhh=zeros(m,m,m);
for i=1:m
    h(:,:,i)=k;
end
for i=1:m
    hh(:,i,:)=k;
end
for i=1:m
    hhh(i,:,:)=k;
end
for i=1:length(hhhh(:))
    a=h(:);
    aa=a(i);
    b=hh(:);
    bb=b(i);
    c=hhh(:);
    cc=c(i);
    hhhh(i)=min([aa bb cc]);
end
weight=hhhh(:);
end

