function out=grille2d(xmin,pasx,nx,ymin,pasy,ny)
% grille2d : define a 2d grid
a=[kron(ones(ny,1),(1:nx)'*pasx) ,kron((1:ny)',ones(nx,1)*pasy)];
out(:,1)=xmin+a(:,1)-pasx;
out(:,2)=ymin+a(:,2)-pasy;
