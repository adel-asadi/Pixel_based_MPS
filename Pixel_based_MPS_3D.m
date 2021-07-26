function [XX]=Pixel_based_MPS_3D(TI,x2,grid,T,type,dim,pwr,paste,condi_uncondi,t1,t2,t3,stationary,simu,mgmg,realization_mg,mg,Z,Tem,distfun,IDX,cntr,clusternum,km,W1,W2,weight,categorical,CC,cdf_rnd,topo,kc,ijk_ind_od,ijk_ind_odz)


p1=max(grid(:,1)); p2=max(grid(:,2)); p3=max(grid(:,3));
lim=[min(grid(:,1)) max(grid(:,1)) min(grid(:,2)) max(grid(:,2)) min(grid(:,3)) max(grid(:,3))];

m=T(1,1); n=T(1,2); o=T(1,3);
m1=(m-1)/2; n1=(n-1)/2; o1=(o-1)/2;

%%

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
    if topo==0
    X=[grid -999*ones(length(grid),1)];
    elseif topo==1
    X=[ijk_ind_od;ijk_ind_odz];
    end
end


if mg==2
    [length_i,length_j,length_k]=size(TI);
    ijk_ind1=[];
    ijk_ind2=[];
    ijk_ind3=[];
    for i=1:4:length_i
    for j=1:4:length_j
    for k=1:4:length_k
        ijk_ind1=[ijk_ind1;i,j,k];
    end
    end
    end
    for i=3:4:length_i
    for j=3:4:length_j
    for k=3:4:length_k
        ijk_ind2=[ijk_ind2;i,j,k];
    end
    end
    end
    for i=1:length_i
    for j=1:length_j
    for k=1:length_k
        ijk_ind3=[ijk_ind3;i,j,k];
    end
    end
    end
    r1=rand(length(ijk_ind1),1);
    r2=rand(length(ijk_ind2),1);
    r3=rand(length(ijk_ind3),1);
    [aa,b1]=sort(r1);
    [aa,b2]=sort(r2);
    [aa,b3]=sort(r3);
    bb={b1,b2,b3};
    ijk_indd={ijk_ind1,ijk_ind2,ijk_ind3};
else
r=rand(length(X),1);
[aa,b]=sort(r);
end


l=clusternum;

centra=zeros(clusternum,m*n*o);
for i=1:clusternum
centra(i,:)=cntr{i};
end


if condi_uncondi==1 && (mg==0 || mg==2)
cpt=[];
    for i=1:length(x2)
        pi=find((x2(i,1)==X(:,1))& x2(i,2)==X(:,2) & x2(i,3)==X(:,3));
        cpt=[cpt; pi];
    end
    X(cpt,4)=x2(:,4);
end


Y=(ones(lim(1,2)+(2*m1),lim(1,4)+(2*n1),lim(1,6)+(2*o1))).*-999;
for i=1:length(X(:,1))
Y((X(i,1)+m1),(X(i,2)+n1),(X(i,3)+o1))=X(i,4);
end
%%

if stationary==true
    if mg==0
    h=waitbar(0,'wait, calculation progressing ');
    for k=1:length(b)
        if (X(b(k),4)==-999)
        c=X(b(k),1:3);
        YY=Y(c(1):(c(1)+(m-1)),c(2):(c(2)+(n-1)),c(3):(c(3)+(o-1)));
        t=(YY(:))';
        v=sum((t(:)==-999));
        if topo==1
            vv=sum((t(:)==-1000));
        end
        if v==m*n*o
            t4=randi(clusternum, [1 1]);
        elseif topo==1 && (v+vv)==m*n*o
            t4=randi(clusternum, [1 1]);
        else
            t2=t;
            t3=zeros(l,1);
            t3=Eu_dist_v3(t2,centra,pwr,1, type,weight,topo);
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
            %elseif rn>dstrb(1,3) && rn<=dstrb(end-1,3)
            else
                for i=1:(length(dstrb(:,3))-1)
                    if rn>dstrb(i,3) && rn<=dstrb(i+1,3)
                        rnn=dstrb(i+1,1);
                    end
                end
            end
        elseif (categorical==1) && ((kc(t4)>100)==1)
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
        elseif categorical==3
            rnn3=centra(t4,:);
            t83=reshape(rnn3,m,n,o);
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
        
        if v==m*n*o
            t7=randi(length(t5), [1 1]);
        elseif topo==1 && (v+vv)==m*n*o
            t7=randi(length(t5), [1 1]);
        else
            t3=zeros(ll,1);
            t3=Eu_dist_v3(t2,t6,pwr,1,type,weight,topo);
            u2=t3;
            t7=find(u2==min(u2));
        end
        
        r=randi(length(t7));
        t71=t6(t7(r),:);
        t8=reshape(t71,m,n,o);
        r=randi(length(t4));
        t711=CC{t4(r)};
        t811=reshape(t711,m,n,o);
        end
        
            if categorical==11
                oo=t811(m1+1,n1+1,o1+1);
                ooo=rand(1);
            if ooo<=oo
                X(b(k),4)=1;
            else
                X(b(k),4)=0;
            end
            elseif categorical==2
                X(b(k),4)=rnn;
                Y((c(1)+m1),(c(2)+n1),(c(3)+o1))=rnn;
            elseif (categorical==1) && ((kc(t4)>100)==1)
                X(b(k),4)=rnn1;
                Y((c(1)+m1),(c(2)+n1),(c(3)+o1))=rnn1;
            elseif categorical==3
                X(b(k),4)=t83(m1+1,n1+1,o1+1);
                Y((c(1)+m1),(c(2)+n1),(c(3)+o1))=t83(m1+1,n1+1,o1+1);
            else
        X(b(k),4)=t8(m1+1,n1+1,o1+1);
        Y((c(1)+m1),(c(2)+n1),(c(3)+o1))=t8(m1+1,n1+1,o1+1);
            end
        end
        waitbar(k/length(b),h)
    end

elseif mg==2
      h=waitbar(0,'wait, calculation progressing ');
    for p=1:3
            bbb=bb{p};
            ijk_ind=ijk_indd{p};
        for k=1:length(bbb)
            b_ind=find((ijk_ind(bbb(k),1)==X(:,1))& (ijk_ind(bbb(k),2)==X(:,2))& (ijk_ind(bbb(k),3)==X(:,3)));
        if (X(b_ind,4)==-999)
        c=X(b_ind,1:3);
        YY=Y(c(1):(c(1)+(m-1)),c(2):(c(2)+(n-1)),c(3):(c(3)+(o-1)));
        t=(YY(:))';
        v=sum((t(:)==-999));
        if topo==1
            vv=sum((t(:)==3));
        end
        if v==m*n*o
            t4=randi(clusternum, [1 1]);
        elseif topo==1 && (v+vv)==m*n*o
            t4=randi(clusternum, [1 1]);
        else
            t2=t;
            t3=zeros(l,1);
            t3=Eu_dist_v3(t2,centra,pwr,1,type,weight,topo);
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
            %elseif rn>dstrb(1,3) && rn<=dstrb(end-1,3)
            else
                for i=1:(length(dstrb(:,3))-1)
                    if rn>dstrb(i,3) && rn<=dstrb(i+1,3)
                        rnn=dstrb(i+1,1);
                    end
                end
            end
        elseif (categorical==1) && ((kc(t4)>100)==1)
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
        
        if v==m*n*o
            t7=randi(length(t5), [1 1]);
        elseif topo==1 && (v+vv)==m*n*o
            t7=randi(length(t5), [1 1]);
        else
            t3=zeros(ll,1);
            t3=Eu_dist_v3(t2,t6,pwr,1, type,weight,topo);
            
            u2=t3;
            t7=find(u2==min(u2));
        end
        
        r=randi(length(t7));
        t71=t6(t7(r),:);
        t8=reshape(t71,m,n,o);
        r=randi(length(t4));
        t711=CC{t4(r)};
        t811=reshape(t711,m,n,o);
        end
        
            if categorical==11
                oo=t811(m1+1,n1+1,o1+1);
                ooo=rand(1);
            if ooo<=oo
                X(b(k),4)=1;
            else
                X(b(k),4)=0;
            end
            elseif categorical==2
                X(b_ind,4)=rnn;
                Y((c(1)+m1),(c(2)+n1),(c(3)+o1))=rnn;
            elseif (categorical==1) && ((kc(t4)>100)==1)
                X(b_ind,4)=rnn1;
                Y((c(1)+m1),(c(2)+n1),(c(3)+o1))=rnn1;
            else
        X(b_ind,4)=t8(m1+1,n1+1,o1+1);
        Y((c(1)+m1),(c(2)+n1),(c(3)+o1))=t8(m1+1,n1+1,o1+1);
            end
        end
        waitbar(k/length(bbb),h)
        end
    end
    end
end
    
%close(h)
XX= X(:,4);

end
