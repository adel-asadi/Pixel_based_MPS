function [XX]=Pixel_based_MPS_2D(TI,x2,grid,T,type,dim,pwr,paste,condi_uncondi,t1,t2,stationary,simu,mgmg,realization_mg,mg,Z,Tem,distfun,IDX,cntr,clusternum,km,W1,W2,weight,categorical,cdf_rnd,kc,out_dom,CC)
%%

p1=max(grid(:,1)); p2=max(grid(:,2));
lim=[min(grid(:,1)) max(grid(:,1)) min(grid(:,2)) max(grid(:,2))];
m=T(1,1); n=T(1,2);
m1=(m-1)/2; n1=(n-1)/2;

if distfun==1
    distf='euclidean';
elseif distfun==2
    distf='seuclidean';
elseif distfun==3
    distf='mahalanobis';
elseif distfun==4
    distf='cityblock';
elseif distfun==5
    distf='minkowski';
elseif distfun==6
    distf='chebychev';
elseif distfun==7
    distf='cosine';
elseif distfun==8
    distf='correlation';
elseif distfun==9
    distf='hamming';
elseif distfun==10
    distf='jaccard';
elseif distfun==11
    distf='spearman';
end

if mg==0 || mg==2
        X=[grid -999*ones(length(grid),1)];
end

if mg==2
    ij_ind1=[];
    ij_ind2=[];
    ij_ind3=[];
    for i=1:4:length(TI(:,1))
    for j=1:4:length(TI(1,:))
        ij_ind1=[ij_ind1;i,j];
    end
    end
    for i=3:4:length(TI(:,1))
    for j=3:4:length(TI(1,:))
        ij_ind2=[ij_ind2;i,j];
    end
    end
    for i=1:length(TI(:,1))
    for j=1:length(TI(1,:))
        ij_ind3=[ij_ind3;i,j];
    end
    end
    r1=rand(length(ij_ind1),1);
    r2=rand(length(ij_ind2),1);
    r3=rand(length(ij_ind3),1);
    [aa,b1]=sort(r1);
    [aa,b2]=sort(r2);
    [aa,b3]=sort(r3);
    bb={b1,b2,b3};
    ij_indd={ij_ind1,ij_ind2,ij_ind3};
end

r=rand(length(X),1);
[aa,b]=sort(r);

l=clusternum;
centra=zeros(clusternum,m*n);

for i=1:clusternum
centra(i,:)=cntr{i};
end

if (condi_uncondi==1 && mgmg==1) || (condi_uncondi==1 && (mg==0 || mg==2))
cpt=[];
    for i=1:length(x2)
        pi=find((x2(i,1)==X(:,1))& x2(i,2)==X(:,2));
        cpt=[cpt; pi];
    end
    X(cpt,3)=x2(:,3);
end

Y=(ones(lim(1,2)+(2*m1),lim(1,4)+(2*n1))).*-999;
for i=1:length(X(:,1))
Y((X(i,1)+m1),(X(i,2)+n1))=X(i,3);
end

%%

if stationary==true
    h=waitbar(0,'wait, calculation progressing ');
    if mg==0
    for k=1:length(b)
        if (X(b(k),3)==-999)
        c=X(b(k),1:2);
        YY=Y(c(1):(c(1)+(m-1)),c(2):(c(2)+(n-1)));
        t=(YY(:))';
        v=sum((t(:)==-999));
        if out_dom==1
            vv=sum((t(:)==-1000));
        end
        if v==m*n
                t4=randi(clusternum, [1 1]);
        elseif out_dom==1 && (v+vv)==m*n
            t4=randi(clusternum, [1 1]);
        else
            t2=t;
            t3=zeros(l,1);
            t3=Eu_dist_v(t2,centra,pwr,1,type,weight,out_dom);
            u1=t3;
            t44=find(u1==min(u1));
            t444=randi(length(t44));
            t4=t44(t444);
        end

        if categorical==2
            dstrb=cdf_rnd{t4};
            rn=randi(100);
            if rn<=dstrb(1,3)
                rnn=dstrb(1,1);
            else
                for i=1:length(dstrb(:,3))
                    if rn>dstrb(i,3) && rn<=dstrb(i+1,3)
                        rnn=dstrb(i+1,1);
                    end
                end
            end
        elseif (categorical==1) && ((kc(t4)>50)==1)
            dstrb=cdf_rnd{t4};
            rn=randi(100);
            if rn<=dstrb(1,3)
                rnn1=dstrb(1,1);
            %elseif rn>dstrb(1,3) && rn<=dstrb(end-1,3)
            else
                for i=1:(length(dstrb(:,3))-1)
                    if rn>dstrb(i,3) && rn<=dstrb(i+1,3)
                        rnn1=dstrb(i+1,1);
                    end
                end
            end
        else
        
        if type>=1
        t5=km{1,t4};
        else
        Cat=randi(type,[1 1]);
        ind=1;
        for i=1:type
            if Cat>=prob(i,t4)
            ind=ind+1;
            end
        end
        t5=km{1,t4}{1,ind};
        end

        t6=Tem(t5,:);
        ll=length(t6(:,1));
        
        if v==m*n
            t7=randi(length(t5), [1 1]);
        elseif out_dom==1 && (v+vv)==m*n
            t7=randi(length(t5), [1 1]);
        else
            t3=zeros(ll,1);
            t3=Eu_dist_v(t2,t6,pwr,1,type,weight,out_dom);
     
            u2=t3;
            t7=find(u2==min(u2));
        end
        
        r=randi(length(t7));
        t71=t6(t7(r),:);
        t8=reshape(t71,m,n);
        r=randi(length(t4));
        t711=CC{t4(r)};
        t811=reshape(t711,m,n);

        end
    
            if categorical==2
                X(b(k),3)=rnn;
                Y((c(1)+m1),(c(2)+n1))=rnn;
            elseif (categorical==1) && ((kc(t4)>50)==1)
                X(b(k),3)=rnn1;
                Y((c(1)+m1),(c(2)+n1))=rnn1;
            else
            X(b(k),3)=t8(m1+1,n1+1);
            Y((c(1)+m1),(c(2)+n1))=t8(m1+1,n1+1);
            end
        end
        waitbar(k/length(b),h) 
    end

    elseif mg==2
        for p=1:3
            bbb=bb{p};
            ij_ind=ij_indd{p};
        for k=1:length(bbb(:,1))
            b_ind=find((ij_ind(bbb(k),1)==X(:,1))& (ij_ind(bbb(k),2)==X(:,2)));
        if (X(b_ind,3)==-999)
        c=X(b_ind,1:2);
        %t=tmpldat(X,lim,c,D,0);
        YY=Y(c(1):(c(1)+(m-1)),c(2):(c(2)+(n-1)));
        t=(YY(:))';
        v=sum((t(:)==-999));
        if out_dom==1
            vv=sum((t(:)==-1000));
        end
        if v==m*n
            t4=randi(clusternum, [1 1]);
        elseif out_dom==1 && (v+vv)==m*n
            t4=randi(clusternum, [1 1]);
        else
            
            t2=t;
            t3=zeros(l,1);
            t3=Eu_dist_v(t2,centra,pwr,0,type,weight,out_dom);
            
            u1=t3;
            t44=find(u1==min(u1));
            t444=randi(length(t44));
            t4=t44(t444);
        end
        
        if categorical==2
            dstrb=cdf_rnd{t4};
            rn=randi(100);
            if rn<=dstrb(1,3)
                rnn=dstrb(1,1);
            else
                for i=1:length(dstrb(:,3))
                    if rn>dstrb(i,3) && rn<=dstrb(i+1,3)
                        rnn=dstrb(i+1,1);
                    end
                end
            end
        elseif (categorical==1) && ((kc(t4)>50)==1)
            dstrb=cdf_rnd{t4};
            rn=randi(100);
            if rn<=dstrb(1,3)
                rnn1=dstrb(1,1);
            %elseif rn>dstrb(1,3) && rn<=dstrb(end-1,3)
            else
                for i=1:(length(dstrb(:,3))-1)
                    if rn>dstrb(i,3) && rn<=dstrb(i+1,3)
                        rnn1=dstrb(i+1,1);
                    end
                end
            end
        else
        
        if type>=1
        t5=km{1,t4};
        else
        Cat=randi(type,[1 1]);
        ind=1;
        for i=1:type
            if Cat>=prob(i,t4)
            ind=ind+1;
            end
        end
        t5=km{1,t4}{1,ind};
        end

        t6=Tem(t5,:);
        ll=length(t6(:,1));
        
        if v==m*n
            t7=randi(length(t5), [1 1]);
        else
            t3=zeros(ll,1);
            t3=Eu_dist_v(t2,t6,pwr,1,type,weight,out_dom);
            u2=t3;
            t7=find(u2==min(u2));
        end
        
        r=randi(length(t7));
        t71=t6(t7(r),:);
        t8=reshape(t71,m,n);
        r=randi(length(t4));
        t711=CC{t4(r)};
        t811=reshape(t711,m,n);

        end
       
            if categorical==2
                X(b_ind,3)=rnn;
                Y((c(1)+m1),(c(2)+n1))=rnn;
            elseif (categorical==1) && ((kc(t4)>50)==1)
                X(b_ind,3)=rnn1;
                Y((c(1)+m1),(c(2)+n1))=rnn1;
            else
            X(b_ind,3)=t8(m1+1,n1+1);
            Y((c(1)+m1),(c(2)+n1))=t8(m1+1,n1+1);
            end
        waitbar(k/length(b),h) 
        end
        end
    end
    
close(h)

    end
    
    XX=X(:,3);
end