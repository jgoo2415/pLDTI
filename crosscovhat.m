function [crosshat]=crosscovhat(dnhat, dvdd, dddx, delta, steps)

crosshat=zeros(3,3, steps+1, steps+1);

factor2=1/(4*pi);

mat1=dnhat{1,1}(:,2);%time1
mat2=dnhat{2,1}(:,2);%time2
covpart=dvdd{1,1}(:,:,2)*(mat1*mat2')*dvdd{2,1}(:,:,2)';
crosshat(:,:,2,2)=delta*factor2*covpart;

    
for p=2:steps  
        mat1=dnhat{1,1}(:,p);
        mat2=dnhat{2,1}(:,p);
        covpart=dvdd{1,1}(:,:,p)*(mat1*mat2')*dvdd{2,1}(:,:,p)'; % 
        covpart1=dvdd{1,1}(:,:,p)*dddx{1,1}(:,:,p); % 
        covpart2=dvdd{2,1}(:,:,p)*dddx{2,1}(:,:,p);
        crosshat(:,:,p+1,p+1)=crosshat(:,:,p,p)+delta*covpart1*crosshat(:,:,p)+...
        delta*crosshat(:,:,p)*covpart2'+delta*factor2*covpart;  
end

for p=2:steps
    for k=1:p
        covpart1=dvdd{1,1}(:,:,p)*dddx{1,1}(:,:,p); % 
        crosshat(:,:,p+1,k)=crosshat(:,:,p,k)+delta*covpart1*crosshat(:,:,p,k);  
        crosshat(:,:,k,p+1)=crosshat(:,:,p+1,k);
    end
end



        
end