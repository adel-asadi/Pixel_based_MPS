function [km IDX cntr clusternum]=cluster3D_auto_dbscan(W2,Tem,T,type,stationary,epsilon,MinPts,distfun,topo,k_means,cluster_number)

if stationary==true
m=T(1,1);n=T(1,2);o=T(1,3);

if k_means==0

if topo==1
    D=pdist(W2);
else
    D=pdist(W2(:,4:end));
end
%IDX=dbscan(W2(:,4:end),epsilon,MinPts);
%IDX = dbscan(D,epsilon,MinPts,'Distance','precomputed');


if distfun==1
    %IDX = dbscan(W2(:,3:end),epsilon,MinPts); %Euclidean=default
    IDX = dbscan(D,epsilon,MinPts,'Distance','precomputed');
elseif distfun==2
    IDX = dbscan(W2(:,4:end),epsilon,MinPts,'Distance','seuclidean');
elseif distfun==3
    IDX = dbscan(W2(:,4:end),epsilon,MinPts,'Distance','mahalanobis');
elseif distfun==4
    IDX = dbscan(W2(:,4:end),epsilon,MinPts,'Distance','cityblock');
elseif distfun==5
    IDX = dbscan(W2(:,4:end),epsilon,MinPts,'Distance','minkowski','P',3);
elseif distfun==6
    IDX = dbscan(W2(:,4:end),epsilon,MinPts,'Distance','chebychev');
elseif distfun==7
    IDX = dbscan(W2(:,4:end),epsilon,MinPts,'Distance','cosine');
elseif distfun==8
    IDX = dbscan(W2(:,4:end),epsilon,MinPts,'Distance','correlation');
elseif distfun==9
    IDX = dbscan(W2(:,4:end),epsilon,MinPts,'Distance','hamming');
elseif distfun==10
    IDX = dbscan(W2(:,4:end),epsilon,MinPts,'Distance','jaccard');
elseif distfun==11
    IDX = dbscan(W2(:,4:end),epsilon,MinPts,'Distance','spearman');
end

clusternum=max(IDX);

if topo==1
    PlotClusterinResult(W2,IDX);
title(['DBSCAN Clustering (\epsilon = ' num2str(epsilon) ', MinPts = ' num2str(MinPts) ')']);
else
PlotClusterinResult(W2(:,4:end), IDX);
title(['DBSCAN Clustering (\epsilon = ' num2str(epsilon) ', MinPts = ' num2str(MinPts) ')']);
end

elseif k_means==1
    IDX = kmeans(Tem,cluster_number);
    clusternum=max(IDX);
end

%[IDX, C] = k_means(W2(:,3:end), clusternum);
v=ceil((m*n*o)/3);
%perform clustering algorithm to devide the training image data set
%for j=1:clusternum

if type>=1
jj=0;

if topo==1
    for j=1:clusternum
    t=find(IDX==j);
    if (isempty(t)~=1)
        jj=jj+1;
        cntr{jj}=mean(Tem(t,:),1);
        km{1,jj}=t;
    end
    end
else
for j=1:clusternum
    t=find(IDX==j);
    if (isempty(t)~=1)
        jj=jj+1;
        cntr{jj}=mean(Tem(t,:),1);
        km{1,jj}=t;
    end
end
end
%cntr

p1=[];
clusternum=jj;
elseif type==0
%if type>1
    jj=0;

for j=1:clusternum
    t=find(IDX==j);
    if (isempty(t)~=1)
        jj=jj+1;
        %cntr{jj}=mean(Tem(t,:),1);
        Q=[];        
        for i=1:type
            t1=find(Tem(t,v)==i);
            km{1,jj}{i}=t(t1,:);
            p(i,jj)=(length(t1)/length(t));
            q=cat2bin(Tem(t,:),i);
            Q=[Q q];
        end
        cntr{jj}=mean(Q,1);
        %end
    end
end
p1=cumsum(p);
clusternum=jj;
end
end
cntr;
end