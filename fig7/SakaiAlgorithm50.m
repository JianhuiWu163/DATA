close all;clc;clear all;
%*************************************************************************%
% Sakai algorithm: Flocking for Multirobots Without Distinguishing Robots and Obstacles
%*************************************************************************%
% ϵͳ��ʼ��
loop=3000;                      % ȷ��ѭ������
s=0.5;                          % sȡֵ��ΧΪ(0,1)
n=50;                           % ȷ�����������
m=50;                           % ȷ���쵼��ֱ��Ӱ����������Ŀ

filename = 'SA_50_16.xlsx';
for Dbeta=16:2:26
    for iters=1:100
        TestNum = [n,Dbeta,iters]
r=6;                           	% ȷ���������֪�뾶
ra=(1/s)*(sqrt(1+s*(r^2))-1);
dw=4;                        	% ȷ������Lattice����
dwa=(1/s)*(sqrt(1+s*(dw^2))-1);
h=0.2;                         	% �������h����ȡֵ��ΧΪ��0,1��

%%%%%%%%%%%%%%%%%%%%%%%%
beta1=5;   beta2=3;
c1=4;      c2=4;                % �쵼��Ӱ�캯��������c1,c2>0

step = 0.1;                     % ȷ������
dmin = inf*ones(1,loop);		% ���ݻ��棺�洢ÿ�ε�������С���
vmax = 5;

Xpbuff = zeros(n,loop);         % �洢λ������
Ypbuff = zeros(n,loop);
Xpt = zeros(1,loop);            % �洢λ������
Ypt = zeros(1,loop);
Xvbuff = zeros(n,loop);         % �洢�ٶ�����
Yvbuff = zeros(n,loop); 
% ��ʼ������λ�ú��ٶ����� 
p=zeros(2,n);
v=2*rand(2,n);
state=1;
i=1;
% InitRang=15;                    % �����ʼ��Χ

P0=zeros(2,n);
while state>=1
    if(i==1)
        P0(1,i)=15;             	% ��ʼ���������ʼλ������
        P0(2,i)=20;                 % ��ʼ���������ʼλ������
        i=i+1;
        if(i>n)
        	state=0;
    	end
    else
        P0(1,i)= 15 * rand + 5;     % ��ʼ���������ʼλ������
        P0(2,i)= 20 * rand + 10;     % ��ʼ���������ʼλ������
        RegenFlag=0;
        for j=1:i-1
            if(sqrt((P0(1,i)-P0(1,j))^2+(P0(2,i)-P0(2,j))^2)<=2)
                RegenFlag = 1;
            end
        end
        if(RegenFlag==0)
            i=i+1;
            if(i>n)
                state=0;
            end
        end
    end
end
p = P0;

Rip = zeros(3,n);
Sip = zeros(3,n);

% ��ʼ�쵼��λ�ú��ٶ����� 
pt(1,1)=80; pt(2,1)=80;   
vt(1,1)=0; vt(2,1)=0; 

% ��ʼ������Ŀ����λ�ú��ٶ����� 
pf=zeros(2,n);vf=zeros(2,n);

Xmax=100;
Ymax=Xmax;
%***************************************************************%
% ��ϰ�������
oRa = 8; 

% Dbeta = 26;
OX = 35;
Dd = 26;

oyk(1,1) = OX;                              oyk(2,1) = 21;
oyk(1,2) = OX;                              oyk(2,2) = oyk(2,1) + Dd;
oyk(1,3) = OX;                              oyk(2,3) = oyk(2,2) + Dd;

oyk(1,4) = OX + sqrt(Dbeta^2-(Dd/2)^2);     oyk(2,4) = oyk(2,2) + Dd/2;
oyk(1,5) = OX + sqrt(Dbeta^2-(Dd/2)^2);     oyk(2,5) = oyk(2,1) + Dd/2;

obsnum = 0;
dtheta = pi/300;
for i=1:5
    for oRad = 0:dtheta:2*pi-dtheta
        obsnum = obsnum + 1;
        op(1,obsnum) = oyk(1,i) + oRa * cos(oRad); 
        op(2,obsnum) = oyk(2,i) + oRa * sin(oRad); 
    end
end

obsnum = length(op(1,:));
ov=zeros(2,obsnum);
%***************************************************************%
% figure(1)
% set(gcf,'Position',[100 100 400 400]);%ͼƬ��С
% set(gca,'Position',[.10 .10 .85 .85]);%��������ռ����
% 
% plot(op(1,:),op(2,:),'k.');
% axis([0 Xmax 0 Ymax]);
% hold on;
% 
% plot(pt(1,1),pt(2,1),'g>','LineWidth',1,'MarkerSize',4);
% plot(p(1,:),p(2,:),'b*');
% pause(0.001);
% hold off;
%--------------------------------------------------------------%
FlyState = zeros(1,n);	%�����ʼʱ�������������֪��Χû���ϰ���
totalNum = obsnum + n;
oTotal0 = length(oyk(1,:)); 
oLen = obsnum/oTotal0;
oTotal = oTotal0 + n;             % �ϰ����ܸ�����  

for ld=1:loop
    dp = [op,p];	% �ϰ�������������
	dv = [ov,v];	% �ϰ�������������
    
    Xpt(1,ld) =  pt(1,1);
	Ypt(1,ld) =  pt(2,1);
	for i=1:n
		Xpbuff(i,ld) =  p(1,i);
		Ypbuff(i,ld) =  p(2,i);
		Xvbuff(i,ld) =  v(1,i);
		Yvbuff(i,ld) =  v(2,i);
    end 
    %----------------------------------------------------------------%
	% �����С���
    d_beta=zeros(n,totalNum);
    for i=1:n
        for k=1:totalNum
			if(i==k-obsnum)
				d_beta(i,k) = inf;	%������ų������屾��
			else
				d_beta(i,k) = sqrt((p(1,i)-dp(1,k))^2+(p(2,i)-dp(2,k))^2);
			end
			
            if(d_beta(i,k)<dmin(1,ld)) 
            	dmin(1,ld) = d_beta(i,k);
            end
        end
    end
    %----------------------------------------------------------------%
	% ���ϰ����Ӱ��
    % 1. ���干ʶ����A���ж�������֮����໥Ӱ��
    Ab0=zeros(n,totalNum);
    tData=zeros(n,totalNum);
    for i=1:n
        for k=1:totalNum
            if((d_beta(i,k)<=r)&&(i~=k-obsnum))
                Ab0(i,k)=1;
                
                if(dp(2,k)-p(2,i)>=0)
                    tData(i,k) = acos((dp(1,k)-p(1,i))/d_beta(i,k));
                else
                	tData(i,k) = 2*pi - acos((dp(1,k)-p(1,i))/d_beta(i,k));
                end
             	tData(i,k)=mod(tData(i,k),2*pi);    %������λ��0-2pi��Χ��
                    
                if(k>1)
                	for jj=1:k-1
                        if(abs(tData(i,k)-tData(i,jj))<=pi/360)
                            if(d_beta(i,k)<d_beta(i,jj))
                                Ab0(i,jj)=0;
                            else
                                Ab0(i,k)=0;
                                break;
                            end
                        end
                    end
                end 
            end
        end    
    end
    
    A_beta=zeros(n,oTotal);
    for i=1:n
        for k=1:oTotal
            if(k<=oTotal0)
                if(sum(Ab0(i,(k-1)*oLen+1:k*oLen))>0)
                    A_beta(i,k)=1;
                end
            else
                if(Ab0(i,oTotal0*oLen+(k-oTotal0))==1)
                    A_beta(i,k)=1;
                end
            end
        end
    end
%----------------------------------------------------------------%
	oI = eye(2);    % 2�׵�λ����
	pik = zeros(n,oTotal,2);
	vik = zeros(n,oTotal,2);
	oak = zeros(2,1);
	oPak = zeros(2,2);
    sik = zeros(2,1);
	ssT = zeros(2,2);
	for i=1:n
    	for k=1:oTotal
            if(A_beta(i,k)==1)
                if(k<=oTotal0)    % ���������ϰ����ͶӰλ�ú��ٶ�
                    ouik = oRa/sqrt((p(1,i)-oyk(1,k))^2+(p(2,i)-oyk(2,k))^2);
                    pik(i,k,1) = ouik * p(1,i) + (1-ouik) * oyk(1,k);	% �ϰ���ͶӰλ��
                    pik(i,k,2) = ouik * p(2,i) + (1-ouik) * oyk(2,k);         
                else
                    pik(i,k,1) = p(1,k-oTotal0);
                    pik(i,k,2) = p(2,k-oTotal0);
                end
                %---------------------------------------------------------------------------------------%
                if(sqrt((p(1,i)-pik(i,k,1))^2+(p(2,i)-pik(i,k,2))^2)~=0)
                    sik(1,1) = (p(1,i)-pik(i,k,1))/sqrt((p(1,i)-pik(i,k,1))^2+(p(2,i)-pik(i,k,2))^2);
                    sik(2,1) = (p(2,i)-pik(i,k,2))/sqrt((p(1,i)-pik(i,k,1))^2+(p(2,i)-pik(i,k,2))^2);
                else
                    sik(1,1) = 0;
                    sik(2,1) = 0;
                end
                ssT = sik*sik';
            	vik(i,k,1) = v(1,i) - (ssT(1,1)*v(1,i)+ssT(1,2)*v(2,i)) + (ssT(1,1)*vt(1,1)+ssT(1,2)*vt(2,1));
            	vik(i,k,2) = v(2,i) - (ssT(2,1)*v(1,i)+ssT(2,2)*v(2,i)) + (ssT(2,1)*vt(1,1)+ssT(2,2)*vt(2,1));
            end
        end
    end
    % ʵ��fya(z)
    % ����n(i,j)
    N_beta=zeros(n,oTotal,2);
    odik=zeros(n,oTotal);
    for i=1:n
        for k=1:oTotal
            if(A_beta(i,k)==1)
                odik(i,k) = sqrt((p(1,i)-pik(i,k,1))^2+(p(2,i)-pik(i,k,2))^2);       
                N_beta(i,k,1)=(pik(i,k,1)-p(1,i))/sqrt(1+s*odik(i,k)^2);  
                N_beta(i,k,2)=(pik(i,k,2)-p(2,i))/sqrt(1+s*odik(i,k)^2);
            end
        end
    end
    % ����da=||pj-pi||��
    da_beta=zeros(n,oTotal);
    for i=1:n
        for k=1:oTotal
            if(A_beta(i,k)==1)
                da_beta(i,k)=(1/s)*(sqrt(1+s*odik(i,k)^2)-1);
            end
        end
    end
    % ����fya(da)
    epsilon0 = 0.001;
    fya_beta=zeros(n,oTotal);
    ph_beta=zeros(n,oTotal);
    for i=1:n
        for k=1:oTotal
            if(A_beta(i,k)==1)
                z1_beta=da_beta(i,k)/dwa;
                if(z1_beta<h && z1_beta>=0)      
                    ph_beta(i,k)=1;
                elseif (z1_beta<=1 || z1_beta>=h)
                    ph_beta(i,k)=0.5*(1+cos(pi*((z1_beta-h)/(1-h))));
                else
                    ph_beta(i,k)=0;
                end
                fy(i,k)=-1/(da_beta(i,k)+epsilon0);
                fya_beta(i,k)=ph_beta(i,k)*fy(i,k)*A_beta(i,k);    % A(i,j)ʹi=jʱ��fya=0
            end
        end
    end
    %-----------------���������λ��Ӱ��------------------%
    % ��u11=fya*N��i,j��
    u11_beta=zeros(n,oTotal,2);
    for i=1:n
        for k=1:oTotal
            u11_beta(i,k,1) = fya_beta(i,k) * N_beta(i,k,1);
            u11_beta(i,k,2) = fya_beta(i,k) * N_beta(i,k,2);
        end
    end
    % ��λ�÷���u1
    u1_beta=zeros(2,n);
    for i=1:n
    	for k=1:oTotal
            u1_beta(1,i) = u1_beta(1,i) + u11_beta(i,k,1);
            u1_beta(2,i) = u1_beta(2,i) + u11_beta(i,k,2);
        end
    end
    %----------------����������ٶ�Ӱ��-------------------%
    %��u22=aij*(pj-pi)
    u22_beta=zeros(n,oTotal,2);
    ph_beta=zeros(n,oTotal);
    for i=1:n
        for k=1:oTotal
            if(A_beta(i,k)==1)
                z1_beta=da_beta(i,k)/ra;
                if(z1_beta<h && z1_beta>=0)      
                    ph_beta(i,k)=1;
                elseif (z1_beta<=1 || z1_beta>=h)
                    ph_beta(i,k)=0.5*(1+cos(pi*((z1_beta-h)/(1-h))));
                else
                    ph_beta(i,k)=0;
                end
            end
            
            u22_beta(i,k,1)=ph_beta(i,k) * (vik(i,k,1)-v(1,i)) * A_beta(i,k);
            u22_beta(i,k,2)=ph_beta(i,k) * (vik(i,k,2)-v(2,i)) * A_beta(i,k);
        end
    end
    u2_beta=zeros(2,n);
    for i=1:n
        for k=1:oTotal
            u2_beta(1,i) = u2_beta(1,i) + u22_beta(i,k,1);
            u2_beta(2,i) = u2_beta(2,i) + u22_beta(i,k,2);
        end
    end  
    
    uBeta=zeros(2,n);
    for i=1:2
		for j=1:n
			uBeta(i,j) = beta1*u1_beta(i,j)+beta2*u2_beta(i,j);
		end
	end
    %****************************************************************************%
    % ���쵼�ߵ�Ӱ��
    uGama=zeros(2,n);
    for i=1:n
        uGama(1,i) = -c1*(p(1,i)-pt(1,1))/sqrt(1+(p(1,i)-pt(1,1))^2+(p(2,i)-pt(2,1))^2)-c2*(v(1,i)-vt(1,1));
        uGama(2,i) = -c1*(p(2,i)-pt(2,1))/sqrt(1+(p(1,i)-pt(1,1))^2+(p(2,i)-pt(2,1))^2)-c2*(v(2,i)-vt(2,1));
    end
    % ����ٶ�u
    u=zeros(2,n);
    for i=1:2
    	for j=1:n
			u(i,j) = uBeta(i,j) + uGama(i,j);
		end
    end   
    % ����������м���
    for j=1:n       
        for i=1:2
            v(i,j) = v(i,j) + step * u(i,j);
        end
        %==============================================%
        for i=1:2
            p(i,j) = p(i,j) + step * v(i,j);
        end
    end
    %====================��̬��ͼ======================%
%     plot(op(1,:),op(2,:),'k.');
%     axis([0 Xmax 0 Ymax]);
%     hold on
%     
% %     plot(opfull(1,:),opfull(2,:),'.','color',[153/255 153/255 153/255]);
%     
%     plot(pt(1,1),pt(2,1),'g>','LineWidth',1,'MarkerSize',4);
%     plot(p(1,:),p(2,:),'b*');
% % 	for i=1:n
% %         [xCircle,yCircle] = circle(p(1,i),p(2,i),6,30);
% %         plot(p(1,i),p(2,i),'b*',xCircle,yCircle,'b:','LineWidth',0.5);
% %     end 
%     pause(0.01);
%     hold off;
end

%****************************************************************************************%
    for j=1:loop
        XlsDataBuff(iters,j) = dmin(1,j);
    end
end
filename(7)=fix(Dbeta/10)+48;
filename(8)=mod(Dbeta,10)+48;    
xlswrite(filename,XlsDataBuff);
end


