function [mat, dmdu, dm2du2]=dnhat_all(DD_array,n,nx,ny,nz,u,h,h1,h2)


mat=zeros(6,1);
dmdu=zeros(6,3);
dm2du2=zeros(6,1);

cut1=16*h^2;
cut2=1/(2*h^2);
cut1h2=16*h2^2;
cut1h1=16*h1^2;
cut2h1=1/(2*h1^2);
cut2h2=1/(2*h2^2);

for i=1:nx
    for j=1:ny
        for k=1:nz
            inner=(u(1)-i/nx)^2+(u(2)-j/ny)^2+(u(3)-k/nz)^2;
            if  inner<=cut1
                mat(:,1)=mat(:,1)+exp(-cut2*inner)*DD_array(:,nx*ny*(k-1)+nx*(j-1)+i);
            end
                
                    % first derivative 
            if inner<cut1h1
               factor1=exp(-inner*cut2h1)/h1;
               kerneldi=-(u(1)-i/nx)*factor1;
               kerneldj=-(u(2)-j/ny)*factor1;
               kerneldk=-(u(3)-k/nz)*factor1;
     
               dmdu(:,1)=dmdu(:,1)+kerneldi*DD_array(:,nx*ny*(k-1)+nx*(j-1)+i);
               dmdu(:,2)=dmdu(:,2)+kerneldj*DD_array(:,nx*ny*(k-1)+nx*(j-1)+i);
               dmdu(:,3)=dmdu(:,3)+kerneldk*DD_array(:,nx*ny*(k-1)+nx*(j-1)+i);
                  
            end
                
            if inner<cut1h2
               % second derivative
               dm2du2(:,1)=dm2du2(:,1)+(inner/h2/h2-3)*exp(-inner*cut2h2)*DD_array(:,nx*ny*(k-1)+nx*(j-1)+i);
            end
        end
    end
end


mat=mat/(((2*pi)^(1.5))*n*(h^3)); 
dmdu=dmdu/(((2*pi)^(1.5))*n*(h1^4));
dm2du2=dm2du2/(((2*pi)^(1.5))*n*(h2^5)); 
end