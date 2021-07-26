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
% of a single 3D variable. (Required Software: MATLAB 2019b or later)

tic

display('Step 1: Reading Parameters')

% Have your training image ready in the worksapce, and name it "TI".
%TI;
%load('TI_channel_3D_test.mat');

% Have your hard data (conditioning data) ready in the worksapce, and name it "x2".
% If no hard data are available, assign a random value to "x2", and do
% unconditional simulation.
%x2=1;

% Decide on whether your simulation is conditional or unconditional.
% 0 : unconditional simulation
% 1 : conditional simulation
condi_uncondi=1;
if condi_uncondi==0
    x2=1;
end

% Assign the number of simulations (realizations):
simu=1;

% Decide if you want to plot all of the realizations or only one of them.
% 0 : One realization (last) will be plotted.
% 1 : All of the realizations will be plotted.
plot_fig=1;

% Decide if you want to implement multiple-grid (multi-resolution)
% simulation approach.
% 0 : No multi-grid
% 2 : three-level coarse/medium/fine resolution multi-grid
mg=0; % Please assign "mg=0" if your simulation grid size is different that training image.

% Select the type of the TI variable via "type":
% 1 : continuous or bi-categorical TI
% Otherwise, assign number of categories to "type".
type=1;

% Template size determination for patterns extraction:
% 0 : manual selection by user
% 1 : automatic selection via entropy-based approach
auto_t=0;
% if "auto_t=0", assign the patterns search template size via "T".
T=[11 11 11];
% if "auto_t=0", assign the maximum value you want to consider.
% Recommended (can differ by case) : 40% of the size of the bigger axis of TI.
MAX_SIZE=25; % (odd value)

% Out-of-domain (above topography) data availability:
% 0 : No Out-of-domain data
% 1 : yes
% If yes (topo=1), assign "-1000" to those nodes in the TI.
topo=0;

% End of user-defined parameters section.

%%

% "ADVANCED OPTIONS"

% PLEASE DO NOT CHANGE ANYTHING FORM THE FOLLOWINGS, UNLESS YOU ARE FAMILIAR
% WITH THE WHOLE CODE, SYNTAXES AND MPS SIMULATION PRECEDURE.

% Weight distribution for distance comparison of data event with patterns
% database:
w_min=0.9; % Minimum weight (far nodes within the template size)
w_max=1; % Maximum weight (nearest nodes within the template size)

% Merging similar clusters (having exactly similar prototypes):
% 0 : No
% 1 : yes
simprot=0; % Only works for CATEGORICAL variables.

% Using Eppirical CDF (cumulative density distribution function - ccdf) for
% second step similarity detection (within the selected cluster), in order
% to increase the speed  by random sampling:
% 0 : No
% 1 : Yes
% 2 : Yes, but only for highly populated cluster (greater than 100 members)
categorical=0;
% This option only works for CATEGORICAL variables.

% Using "K-Means Clustering" on high-dimensional patterns database,
% INSTEAD OF applying "PCA and t-SNE" algorithms to map the data and subsequently
% cluster them by "DBSCAN" algorithm.
% 0 : "t-SNE && DBSCAN" - RECOMMENDED
% 1 : "K-Means Clustering"
k_means=0;
% "k-means=1"  provides faster simulation, but requires the number of
% clusters to be assigned by the user.
cluster_number=500; % Neglected if "k-means=0".
% If "k_means=0", assign the minimum number of patterns to form a cluster
% by DBSCAN. Recommnded values for 3D simulation are 2 (better accuracy)
% and 3 (faster - rule of thumb).
MinPts=3;

% *** END OF USER-DEFINED PARAMETERS ***

%%

dim=3;

[n1 n2 n3]=size(TI);

mgmg=0;
realization_mg=-999;
Z=TI;

grid=grille3d(1,n1,1,1,n2,1,1,n3,1);

[t1,t2,t3]=size(TI);

if auto_t==1
MAX_SIZE = 25;
ent = zeros(1,(MAX_SIZE-1)/2);
var1 = zeros(1,(MAX_SIZE-1)/2);
for i = 3:2:MAX_SIZE
    I2 = entropyfilt(TI,ones(i));
    ent((i-1)/2) = mean(mean(mean(I2)));
    var1((i-1)/2) = var(I2(:));
end
global l;
Tsize = 2*Select_Dimension_withLikelihood(ent) + 1;
T  = [Tsize Tsize Tsize];
end

display('Template Size Determined')
%%

Tem=[];
Tem=patt_ex_fast(TI, T, dim);

display('Patterns Extracted')

%%
if topo==1
TI_vec=TI(:);
TI_vec_od=find(TI_vec==-1000);
TI_vec(:)=0;
TI_vec(TI_vec_od)=1;
TI_mat_od=reshape(TI_vec,[n1,n2,n3]);
for i=1:n3
[row,col] = find(TI_mat_od(:,:,i));
[rowz,colz] = find(~TI_mat_od(:,:,i));
ijk_od{i}=[row,col,(ones(length(row),1).*i)];
ijk_odz{i}=[rowz,colz,(ones(length(rowz),1).*i)];
end
ijk_ind_od=[];
ijk_ind_odz=[];
for i=1:n3
ijk_ind_od=[ijk_ind_od;ijk_od{i}];
ijk_ind_odz=[ijk_ind_odz;ijk_odz{i}];
end
ijk_ind_od(:,4)=-1000;
ijk_ind_odz(:,4)=-999;
%x2=[x2;ijk_ind_od];
else
    ijk_ind_od=1;
    ijk_ind_odz=1;
end

display('Patterns with -1000 Removed')
display('Above Topography Removed From Data')

%%

pwr=1;
paste=0;
stationary=true;

distfun=1;

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

display('Stochastic Embedding Started for Coarse Grid')

if k_means==0
if length(Tem(1,:))>81
    optss = statset('MaxIter',1500);
    W1=tsne(Tem,'NumPCAComponents',81,'Perplexity',50,'Verbose',1,'options',optss);
    %W1=tsne(Tem,'NumPCAComponents',100,'Standardize',true);
    %W1=tsne(Tem,'NumPCAComponents',81);
else
    W1=tsne(Tem);
end
elseif k_means==1
    W1=Tem;
end

display('Stochastic Embedding Done for Coarse Grid')

W1=normalize(W1);

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

k_dist_mean=mean(cD(:,MinPts));

k_dist_std=std(cD(:,MinPts));

epsilon=k_dist_mean+(3*k_dist_std);
end

node_coords=spatial_ext(T,t1,t2,t3,dim);
node_coords=[node_coords(:,2) node_coords(:,1) node_coords(:,3)];

if topo==1
    W2=W1;
else
    W2=[node_coords,W1];
end
%%

out_dom=0;
[km IDX cntr clusternum]=cluster3D_auto_dbscan(W2,Tem,T,type,stationary,epsilon,MinPts,distfun,topo,k_means,cluster_number);

%%
if simprot==1
    for i=1:length(cntr)
        len{i}=find(cntr{i}~=0);
        leng(i)=length(len{i});
    end
    lengt=find(leng==0);
    clusternum=clusternum-length(lengt)+1;
    ccnn=cntr{lengt(1)};
    cntr(lengt)=[];
    cntr{length(cntr)+1}=ccnn;
    kkmm=km(lengt);
    km(lengt)=[];
    kkmmd=kkmm{1};
    for i=2:length(kkmm)
        kkmmd=[kkmmd;kkmm{i}];
    end
    km{length(km)+1}=kkmmd;
end

CC=cntr;

display('Clustering Done for Coarse Grid')

%%

if categorical==2 || categorical==1
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

%%

weight=weight_gen_3d(T(1),w_min,w_max);

%%

display('Simulation Started')

if simu==1
for i=1:simu
    [XX]=Pixel_based_MPS_3D(TI,x2,grid,T,type,dim,pwr,paste,condi_uncondi,t1,t2,t3,stationary,simu,mgmg,realization_mg,mg,Z,Tem,distfun,IDX,cntr,clusternum,km,W1,W2,weight,categorical,CC,cdf_rnd,topo,kc,ijk_ind_od,ijk_ind_odz);
    if topo==0
    realization_data{i}=reshape(XX(:,1),t1,t2,t3);
    end

end  
else
parfor i=1:simu
    [XX]=Pixel_based_MPS_3D(TI,x2,grid,T,type,dim,pwr,paste,condi_uncondi,t1,t2,t3,stationary,simu,mgmg,realization_mg,mg,Z,Tem,distfun,IDX,cntr,clusternum,km,W1,W2,weight,categorical,CC,cdf_rnd,topo,kc,ijk_ind_od,ijk_ind_odz);
    if topo==0
    realization_data{i}=reshape(XX(:,1),t1,t2,t3);
    end
end
end

%%
if plot_fig==1
for i=1:simu
    figure
    realization=realization_data{i};
    S = size(realization);
    [x,y,z] = ndgrid(1:S(1),1:S(2),1:S(3));
    scatter3(x(:),y(:),z(:),321,realization(:),'filled')
    colorbar
end
realization=realization_data{simu};
else
    realization=realization_data{simu};
    figure(i)
    S = size(realization);
    [x,y,z] = ndgrid(1:S(1),1:S(2),1:S(3));
    scatter3(x(:),y(:),z(:),321,realization(:),'filled')
end

%%
toc
