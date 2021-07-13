function [d1]=Eu_dist_v(x,y,p,WD,type,weight,out_dom)
s=length(x(1,:));
sss=find(x(1,:)==-1000);
ss=find(x(1,:)==-999);
if length(ss)==s
    d1=nan;
elseif out_dom==1 && (((length(ss))+(length(sss)))==s)
     d1=nan;
elseif ((length(ss)~=s) && (out_dom==0)) || (out_dom==1 && (((length(ss))+(length(sss)))~=s))
if WD==1 && type>1
    Q=[];
for i=1:type
    q=cat2bin(x,i);
    Q=[Q q];
end
Y=zeros(length(y(:,1)),(length(y(1,:))*3));
for j=1:length(y(:,1))
    vv=[];
    for i=1:type
        v=cat2bin(y(j,:),i);
        vv=[vv v];
    end
    Y(j,:)=vv;
end
V=Y;
sss=ss;
weightt=weight;
for i=1:(type-1)
    weightt=[weightt;weight];
    sss=[sss,(ss+((i*s)))];
end
Q(sss)=[];
V(:,sss)=[];
weightt(sss)=[];
d=(((weightt').*(abs(Q-V))).^p)';
d1=sum(d);
if p==2
    d1=sqrt(d1);
end
else
    Q=x;
    V=y;
    if out_dom==0
        Q(ss)=[];
        V(:,ss)=[];
        weight(ss)=[];
    elseif out_dom==1
        ssss=[ss,sss];
        Q(ssss)=[];
        V(:,ssss)=[];
        weight(ssss)=[];
    end
    d=(((weight').*(abs(Q-V))).^p)';
    if length(d(:,1))>1
        d1=sum(d);
    else
        d1=d;
    end
    if p==2
        d1=sqrt(d1);
    end
end
end
end

%s=length(x(1,:));
%ss=find(x(1,:)==-999);
%if sum(Q==-999)==s
%    d1=nan;
%else
%Q(ss)=[];
%V(:,ss)=[];
%weight(ss)=[];
%d=(((weight').*(abs(Q-V))).^p)';
%d=(((weight').*(Q~=V)).^p)';
%d1=sum(d)+sum(dd)+sum(ddd);
%d1=sum(d);
%d1 = sum(Q~=V);
%if p==2
%d1=sqrt(d1);
%end
%l=clusternum;t3=zeros(l,1);t2=ones(1,225);t3=Eu_dist_v(t2,centra,pwr,1,type,weight);