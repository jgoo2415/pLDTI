%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code Update: 07/07/2024
% Author: Juna Goo (junagoo@boisestate.edu)
% Reference: Goo, J., Sakhanenko, L. and Zhu, D. C. (2024) 
% Test for changes between paired integral curves with application 
% to Diffusion Tensor Imaging (not yet published)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation Example Code for time-invariant curves (H0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Theorem 1 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


nx=40; ny=40; nz=40; nt=2; n=nx*ny*nz; % sample size= # of spatial points
h=0.04; % bandwidth
beta=n*h^6; % beta
h1=(n/beta)^(-1/8); 
h2=(n/beta)^(-1/9);  

x0=[0.5*cos(0.3);0.5*sin(0.3);0.5]; % initial value

dir=48;% the number of gradient directions
[B, ~]=rotateb_new(dir); % create a 48*6 B-matrix
delta=0.01; steps=60; % Euler's method 
v0=[1;1;1]; % an initial eigenvector


a=ones(3,1);
chi_null2=zeros(100,1);

parpool('Processes',4);

tic
parfor (nsim=1:100,4)
rng(1234+nsim)

y_array=cell(nt,1); % each cell, a (nx*ny*nz)*48 matrix 
DD_array=cell(nt,1); % each cell, a 6*(nx*ny*nz) matrix 

for q0=1:nt
    y_array{q0,1}=init_template1b(nx,ny,nz,B); % generate the response vector
    DD_array{q0,1}=dtilda(y_array{q0,1},B);  % ols of the diffusion tensor D 
end

[xnhat, dnhat, dvdd, dddx, trH]=xnhat_all(n,nx,ny,nz,nt,DD_array,x0,delta,steps,h,h1,h2,v0); % estimate the integral curve, diffusion tensor estimator, etc
gam0gam0t=gamgamt_all_prep(n,nx,ny,nz,nt,dnhat,xnhat,steps,h,y_array(:,1),B); % estimate the noise tensor
[muhat,chat]=xnhatconf(dnhat,dvdd,dddx,gam0gam0t(:,1),trH,delta,steps,beta,nt); % estimated mean function and covariance function

ind=1; % at the 1st time point 
alpha=0.05; 
y=trackplot(xnhat{ind,1},muhat{ind,1},chat{ind,1},alpha,n,h,steps,'c'); % 95% confindence ellipsoids
hold on
ind=2; % at the 2nd time point 
y=trackplot(xnhat{ind,1},muhat{ind,1},chat{ind,1},alpha,n,h,steps,'m'); % 95% confindence ellipsoids

[crosshat]=crosscovhat(dnhat, dvdd, dddx, delta, steps);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Corollary 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wnhat=sqrt(n*h^2)*a'*(xnhat{2,1}-xnhat{1,1});
muvec=a'*(muhat{2,1}-muhat{1,1});

covhat=zeros((1+steps),(1+steps));
for j=1:(1+steps)
    for k=1:(1+steps)
        covhat(j,k)=a'*(chat{1,1}(:,:,j,k)+chat{2,1}(:,:,j,k)-2*crosshat(:,:,j,k))*a;
    end
end

diff=(wnhat(:,2:steps+1)-muvec(:,2:steps+1))';
P=covhat(2:steps+1,2:steps+1);
[U S V]=svds(P,2); % tsvd
chi_null2(nsim,1)=diff'*pinv(U*S*V')*diff;

end
toc

delete(gcp('nocreate'));

sum(chi_null2>0) % it should be 100 since the values of the chi-square test cannot be negative.
sum(chi_null2>chi2inv(0.95,2))/100 % type 1 error 
