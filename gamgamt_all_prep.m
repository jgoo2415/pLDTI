function gam0gam0t=gamgamt_all_prep(n,nx,ny,nz,nt,dnhat,xnhat,steps,h,y_array,B)

gam0gam0t=cell(nt,1); % save only a and b
for ind=1:nt
    gam0gam0t{ind}=zeros(6,6, steps+1);        
end
   
BBinv=inv(B'*B)*B';
factor=1/(((2*pi)^(1.5))*n*(h^3));

for q0=1:nt
    for s=1:(steps+1)
        u=xnhat{q0}(:,s);
        sig=zeros(48,48);
        
            for i=1:nx
                for j=1:ny
                    for k=1:nz
                        inner=(u(1)-i/nx)^2+(u(2)-j/ny)^2+(u(3)-k/nz)^2;
                        if inner<=16*h^2
                            kernel=exp(-0.5*(inner/h/h));
                   
                            yy=(y_array{q0}(nx*ny*(k-1)+nx*(j-1)+i, :)'-B*dnhat{q0}(:,s));
                            sig=sig+kernel*(yy*yy');
                        end
                    end
                end
            end
   
        gam0gam0t{q0}(:,:,s)=(BBinv)*sig*(BBinv)'*factor;
    end    
end


end 

