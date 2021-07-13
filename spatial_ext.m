function [ node_coords ] = spatial_ext(T,t1,t2,t3,dim)
%UNTITLED Summary of this function goes here
%   This function extracts the spatial coordinates of the nodes for each
%   pattern and stores them in a database whose indices correspond the the
%   patterns in the pattern database.


%3 dimensional case
if (dim==3);
    
    %calculating the space left on all sides of the grid that will not have
    %nodes for the patterns
    xspace=floor(T(1)/2);
    yspace=floor(T(2)/2);
    zspace=floor(T(3)/2);

    %gridding the entire space and then extracting only the coordinates which
    %correspond to the patterns
    [xcoords, ycoords, zcoords]=meshgrid(1:t2,1:t1,1:t3); %May return an error because of weird grid vector definitions, see line 40
    x_all=xcoords(xspace+1:t1-xspace,yspace+1:t2-yspace,zspace+1:t3-zspace);
    y_all=ycoords(xspace+1:t1-xspace,yspace+1:t2-yspace,zspace+1:t3-zspace);
    z_all=zcoords(xspace+1:t1-xspace,yspace+1:t2-yspace,zspace+1:t3-zspace);

    %reshaping the node coordinates into vectors
    xnodes=reshape(x_all,(t1-2*xspace)*(t2-2*yspace)*(t3-2*zspace),1);
    ynodes=reshape(y_all,(t1-2*xspace)*(t2-2*yspace)*(t3-2*zspace),1);
    znodes=reshape(z_all,(t1-2*xspace)*(t2-2*yspace)*(t3-2*zspace),1);
    
    %combining the coordinates into a single vector
    node_coords=[xnodes ynodes znodes];
    
    clear xspace yspace zspace xcoords ycoords zcoords x_all y_all z_all xnodes ynodes znodes
end

%2 dimensional case, exactly the same as above
if (dim==2)
    xspace=floor(T(1)/2);
    yspace=floor(T(2)/2);

    [ycoords, xcoords]=meshgrid(1:t2,1:t1);
    x_all=xcoords(xspace+1:t1-xspace,yspace+1:t2-yspace);
    y_all=ycoords(xspace+1:t1-xspace,yspace+1:t2-yspace);

    xnodes=reshape(x_all,(t1-2*xspace)*(t2-2*yspace),1);
    ynodes=reshape(y_all,(t1-2*xspace)*(t2-2*yspace),1);
     
    node_coords=[xnodes ynodes];
    
    clear xspace yspace xcoords ycoords x_all y_all xnodes ynodes

end

