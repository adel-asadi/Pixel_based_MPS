% % "MULTIPLE-POINT GEOSTATISTICS MATLAB CODE"

% Reference Paper:
% "Stochastic Embedding and Density-based Clustering of Image
% Patterns for Pixel-based Multiple-Point Geostatistical Simulation"
% Sumitted to the "Mathematical Geosciences" Journal.

% The codes are developed by Adel Asadi and Snehamoy Chatterjee at Michigan
% Technological University, as a novel pixel-based multiple-point geostatistical
% simulation method for stochastic modeling of Earth resources.
% Contact: "aasadi@mtu.edu"

%%

% Please follow this section's instructions to run the code for simulation
% of a single 2D variable. (Required Software: MATLAB 2019b or later)

tic

display('Step 1: Reading Parameters')

% Have your training image ready in the worksapce, and name it "TI".
%TI;
%load('TI_channel_2D_test.mat');

% Assign the number of simulations (realizations):
simu=1;

% Decide if you want to plot all of the realizations or only one of them.
% 0 : One realization (last) will be plotted.
% 1 : All of the realizations will be plotted.
plot_fig=0;

% Decide on whether your simulation is conditional or unconditional.
% 0 : unconditional simulation
% 1 : conditional simulation
condi_uncondi=0;

% Have your hard data (conditioning data) ready in the worksapce, and name it "x2".
% If no hard data are available, assign a random value to "x2", and do
% for unconditional simulation:
if condi_uncondi==0
    x2=1;
end

% Decide if you want to implement multiple-grid (multi-resolution)
% simulation approach.
% 0 : No multi-grid
% 2 : three-level coarse/medium/fine resolution multi-grid
mg=2; % Please assign "mg=0" if your simulation grid size is different that training image.

% Select the type of the TI variable via "type":
% 1 : continuous or bi-categorical TI
% Otherwise, assign number of categories to "type".
type=1;

% Template size determination for patterns extraction:
% 0 : manual selection by user
% 1 : automatic selection via entropy-based approach
auto_t=0;
% if "auto_t=0", assign the patterns search template size via "T".
T=[13 13];
% if "auto_t=0", assign the maximum value you want to consider.
% Recommended (can differ by case) : 40% of the size of the bigger axis of TI.
MAX_SIZE=41; % (odd value)

% Out-of-domain data availability:
% 0 : No Out-of-domain data
% 1 : yes
% If yes (out_dom=1), assign "-1000" to those nodes in the TI.
out_dom=0;

% End of user-defined parameters section.


%%

% "ADVANCED OPTIONS"

% PLEASE DO NOT CHANGE ANYTHING FORM THE FOLLOWINGS, UNLESS YOU ARE FAMILIAR
% WITH THE WHOLE CODE, SYNTAXES AND MPS SIMULATION PRECEDURE.

% Using Eppirical CDF (cumulative density distribution function - ccdf) for
% second step similarity detection (within the selected cluster), in order
% to increase the speed  by random sampling:
% 0 : No
% 1 : Yes
% 2 : Yes, but only for highly populated cluster (greater than 100 members)
categorical=1;
% This option only works for CATEGORICAL variables.

% Using "K-Means Clustering" on high-dimensional patterns database,
% INSTEAD OF applying "PCA and t-SNE" algorithms to map the data and subsequently
% cluster them by "DBSCAN" algorithm.
% 0 : "t-SNE && DBSCAN" - RECOMMENDED
% 1 : "K-Means Clustering" - Faster
k_means=0;
% "k-means=1"  provides faster simulation, but requires the number of
% clusters to be assigned by the user.
cluster_number=100; % Neglected if "k_means=0".
% If "k_means=0", assign the minimum number of patterns to form a cluster
% by DBSCAN. Recommnded value for 2D simulation is 2 (rule of thumb).
MinPts=2;

% Weight distribution for distance comparison of data event with patterns
% database:
w_min=0.9; % Minimum weight (far nodes within the template size)
w_max=1; % Maximum weight (nearest nodes within the template size)

% *** END OF USER-DEFINED PARAMETERS ***

%%

dim=2;

n1=length(TI(:,1));
n2=length(TI(1,:));
n3=1;

mgmg=0;
realization_mg=-999;
Z=TI;

grid=grille2d(1,1,n1,1,1,n2);

[t1,t2]=size(TI);
t3=1;

pwr=1;
paste=0;
stationary=true;

distfun=1;
% 4 : Manhattan (City-Block)
% 1 : Euclidean (L-2 Norm)
% 6 : 'chebychev'
% 12/5 : 'Minkovski-3'
% 3 : Mahalanobis
%1='euclidean'=default
%2='seuclidean'
%3='mahalanobis'
%4='cityblock'
%5='minkowski'
%6='chebychev'
%7='cosine'
%8='correlation'
%9='hamming'
%10='jaccard'
%11='spearman'
%12='Minkovski-3'
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

%%

if auto_t==1
ent = zeros(1,(MAX_SIZE-1)/2);
var1 = zeros(1,(MAX_SIZE-1)/2);
for i = 3:2:MAX_SIZE
    I2 = entropyfilt(TI,ones(i));
    ent((i-1)/2) = mean(mean(mean(I2)));
    var1((i-1)/2) = var(I2(:));
end
global l;
Tsize = 2*Select_Dimension_withLikelihood(ent) + 1;
T  = [Tsize Tsize];
end

%display('Template Size Determined')

%if ns==0
Tem=[];
Tem=patt_ex_fast(TI, T, dim);

if out_dom==1
    j=[];
    for i=1:length(Tem(:,1))
        f=find(Tem(i,:)==-1000);
        if isempty(f)==0
            j=[j;i];
        end
    end
     Tem(j,:)=[];
elseif out_dom==0
    j=[];
    for i=1:length(Tem(:,1))
        f=find(Tem(i,:)==-999);
        if isempty(f)==0
            j=[j;i];
        end
    end
     Tem(j,:)=[];
end

if out_dom==1
p=[0 0];
for i=1:n1
    for j=1:n2
            if TI(i,j)==-1000
                p=[p;i j];
            end
    end
end
p(1,:)=[];
p(:,3)=-1000;
x2=[x2;p];
%x2=p;
x2 = unique(x2(:,1:3), 'rows');
display('Out of Domain Removed From Data')
end

display('Patterns Extracted')
%%
%display('Stochastic Embedding Started for Coarse Grid')

if k_means==0
if length(Tem(1,:))>81
    %W1=tsne(Tem,'NumPCAComponents',81,'Standardize',true,'Distance',distf);
    W1t=tsne(Tem,'Standardize',true,'NumPCAComponents',81);
else
    W1t=tsne(Tem,'Distance',distf);
end
elseif k_means==1
    W1t=Tem;
end

display('Stochastic Embedding of Patterns Done')
%%

W1=normalize(W1t);
%%

%MinPts=dim;
if k_means==1
epsilon=1;
end

if k_means==0
if distfun==5
    [cIdx,cD] = knnsearch(W1,W1,'K',MinPts,'Distance',distf,'P',3);
else
    %[cIdx,cD] = knnsearch(W1,W1,'K',MinPts,'Distance',distf);
    [cIdx,cD] = knnsearch(W1,W1,'K',MinPts,'Distance','chebychev');
end
end
%%

if k_means==0
k_dist_mean=mean(cD(:,MinPts));
k_dist_std=std(cD(:,MinPts));
epsilon=k_dist_mean+(3*k_dist_std);
end
%%

node_coords=spatial_ext(T,t1,t2,1,dim);
node_coords=[node_coords(:,2) node_coords(:,1)];

if out_dom==0
W2=[node_coords,W1];
elseif out_dom==1
W2=W1;
end

[km,IDX,cntr,clusternum]=cluster_auto_dbscan(W2,Tem,T,type,stationary,epsilon,MinPts,distfun,k_means,out_dom,cluster_number);
CC=cntr;
%%
display('Clustering of Patterns Dataabase Done')

if ((categorical==2) || (categorical==1))
    for i=1:clusternum
        members=Tem(km{i},:);
    cent_ind=ceil(length(members(1,:))/2);
    cent_node=members(:,cent_ind);
    tbl=tabulate(cent_node);
    tbl(:,3)=cumsum(tbl(:,3));
    cdf_rnd{i}=tbl;
    end
[a b]=size(km);
kc=zeros(b,1);
for i=1:b
    kc(i)=length(km{i});
end
else
    cdf_rnd=11;
    kc=11;
end

if w_min~=w_max
if T(1)==T(2)
    weight=weight_gen(T(1),w_min,w_max);
else
    weight=ones(T(1)*T(2),1);
end
elseif w_min==w_max
    weight=ones(T(1)*T(2),1);
end

%%

display('Simulation Started')

if simu==1
for i=1:simu
    [XX]=Pixel_based_MPS_2D(TI,x2,grid,T,type,dim,pwr,paste,condi_uncondi,t1,t2,stationary,simu,mgmg,realization_mg,mg,Z,Tem,distfun,IDX,cntr,clusternum,km,W1,W2,weight,categorical,cdf_rnd,kc,out_dom,CC);
    realization_data{i}=reshape(XX(:,1),t1,t2);

end  
else
parfor i=1:simu
    [XX]=Pixel_based_MPS_2D(TI,x2,grid,T,type,dim,pwr,paste,condi_uncondi,t1,t2,stationary,simu,mgmg,realization_mg,mg,Z,Tem,distfun,IDX,cntr,clusternum,km,W1,W2,weight,categorical,cdf_rnd,kc,out_dom,CC);

    realization_data{i}=reshape(XX(:,1),t1,t2);
end
end

%%

if plot_fig==1
for i=1:simu
    figure
    imagesc(realization_data{i})
end
realization=realization_data{simu};
else
    realization=realization_data{simu};
    imagesc(realization)
end

display('Runtime:')

toc
