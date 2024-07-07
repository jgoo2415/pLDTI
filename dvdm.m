function gv=dvdm(m,v,lambda)
% v is a 3 by 1 vector
% m=vech(D) is a 6 by 1 vector for a 3 by 3 symmetric diffusion tensor 
% D11=m(1,1), D12=D21=m(2,1), D13=D31=m(3,1), D22=m(4,1), D23=D32=m(5,1), D33=m(6,1)

gv=zeros(3,6); % dv/dvech(D)is a 3 by 6 matrix

z=diag(lambda)-[m(1,1) m(2,1) m(3,1); m(2,1) m(4,1) m(5,1); m(3,1) m(5,1) m(6,1)];
% pinvz is Moore-Penrose pseudoinverse of z 
p=pinv(z,0.1^8); % 0.1^8= tolerance

% gv=[dv(1,1)/dm(1,1) dv(1,1)/dm(1,2) ... dv(1,1)/dm(3,3);
%     dv(1,2)/dm(1,1) dv(1,2)/dm(1,2) ... dv(1,2)/dm(3,3);
%     dv(1,3)/dm(1,1) dv(1,3)/dm(1,2) ... dv(1,3)/dm(3,3)]

gv(1,1)=p(1,1)*v(1);
gv(1,2)=p(1,1)*v(2)+p(1,2)*v(1);
gv(1,3)=p(1,1)*v(3)+p(1,3)*v(1);
gv(1,4)=p(1,2)*v(2);
gv(1,5)=p(1,2)*v(3)+p(1,3)*v(2);
gv(1,6)=p(1,3)*v(3);

gv(2,1)=p(2,1)*v(1);
gv(2,2)=p(2,1)*v(2)+p(2,2)*v(1);
gv(2,3)=p(2,1)*v(3)+p(2,3)*v(1);
gv(2,4)=p(2,2)*v(2);
gv(2,5)=p(2,2)*v(3)+p(2,3)*v(2);
gv(2,6)=p(2,3)*v(3);

gv(3,1)=p(3,1)*v(1);
gv(3,2)=p(3,1)*v(2)+p(3,2)*v(1);
gv(3,3)=p(3,1)*v(3)+p(3,3)*v(1);
gv(3,4)=p(3,2)*v(2);
gv(3,5)=p(3,2)*v(3)+p(3,3)*v(2);
gv(3,6)=p(3,3)*v(3);
end


