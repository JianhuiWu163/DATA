close all;clc;clear all;
%*************************************************************************%
% The proposed algorithm
%*************************************************************************%
% ϵͳ��ʼ��
loop=3000;                      % ȷ��ѭ������
s=0.5;                          % sȡֵ��ΧΪ(0,1)
n=50;                           % ȷ�����������
m=50;                           % ȷ���쵼��ֱ��Ӱ����������Ŀ

filename = 'Agent50_Data.xlsx';
for iters=1:100
    TestNum = [iters,n]
% ϵͳ��ʼ��
detectR =6;
r=6;                           	% ȷ���������֪�뾶
ra=(1/s)*(sqrt(1+s*(r^2))-1);
dw=4;                        	% ȷ������Lattice����
dwa=(1/s)*(sqrt(1+s*(dw^2))-1);
h=0.2;                         	% �������h����ȡֵ��ΧΪ��0,1��
%%%%%%%%%%%%%%%%%%%%%%%%
a_beta=10;  b_beta=100;         % 0<a_beta<=b_beta
beta1=1;    beta2=0.5;
c1=1;       c2=1;               % �쵼��Ӱ�캯��������c1,c2>0

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
        P0(1,i)=10;             	% ��ʼ���������ʼλ������
        P0(2,i)=20;                 % ��ʼ���������ʼλ������
        i=i+1;
        if(i>n)
        	state=0;
    	end
    else
        P0(1,i)= 10 * rand + 5;     % ��ʼ���������ʼλ������
        P0(2,i)= 30 * rand + 10;     % ��ʼ���������ʼλ������
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
pt(1,1)=60; pt(2,1)=20;   
vt(1,1)=0; vt(2,1)=0; 

% ��ʼ������Ŀ����λ�ú��ٶ����� 
pf=zeros(2,n);vf=zeros(2,n);

Xmax=100;
Ymax=Xmax;
%***************************************************************%
% ��ϰ�������
obsnum=0;
ImpRgb = imread('fig1.png');
Imp0 = rgb2gray(ImpRgb);
[Hlen,Wlen] = size(Imp0);
Imp = ones(Hlen,Wlen);
for i=1:Hlen
    Imp(Hlen+1-i,:)=Imp0(i,:);%���з�ת
end

for i=1:Hlen
    for j=1:Wlen
       if(Imp(i,j)==0) 
           obsnum = obsnum + 1;
           op(1,obsnum)=j*Xmax/Wlen;
           op(2,obsnum)=i*Ymax/Hlen;
       end
    end
end
ov=zeros(2,obsnum);
%--------------------------------------------------------------%
% figure(1)
% set(gcf,'Position',[100 100 400 400]);%ͼƬ��С
% set(gca,'Position',[.10 .10 .85 .85]);%��������ռ����
% 
% plot(op(1,:),op(2,:),'k.');hold on;
% axis([0 Xmax 0 Ymax]);
% plot(pt(1,1),pt(2,1),'g>','LineWidth',1,'MarkerSize',4);
% plot(p(1,:),p(2,:),'b*');
% pause(0.01);
% hold off;
%--------------------------------------------------------------%
FlyState = zeros(1,n);	%�����ʼʱ�������������֪��Χû���ϰ���
totalNum = obsnum + n;

SCount = zeros(1,n);
Dppt = zeros(3,n);
pn0 = p;    
vn0 = v;

test_num= zeros(1,n);
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
    %****************************************************************************%    
	% �ϰ�����������֮���������Ŀ���
	tData=zeros(n,totalNum);
	for i=1:n
		NlineFlag = 0;
		switch(FlyState(i))
			case 0	% ����Ŀ�������쵼��               
				if(pt(2,1)-p(2,i)>=0)
					alpha = acos((pt(1,1)-p(1,i))/sqrt((pt(1,1)-p(1,i))^2+(pt(2,1)-p(2,i))^2));
				else
					alpha = 2*pi - acos((pt(1,1)-p(1,i))/sqrt((pt(1,1)-p(1,i))^2+(pt(2,1)-p(2,i))^2));
                end
			
                Dppt(:,i) = zeros(3,1);
				num = 0;
				ExistObsFlag = 0;				
				for j=1:totalNum
					if(i==j-obsnum)
						rad_temp = inf;	%������ų������屾��
					else
						rad_temp = sqrt((dp(1,j)-p(1,i))^2 + (dp(2,j)-p(2,i))^2);
					end
					
					% ������λ
					if(rad_temp<=detectR)
						num = num + 1;
                        disData(i,num) = rad_temp;	% �洢����
						
						if(dp(2,j)-p(2,i)>=0)
							tData(i,num) = acos((dp(1,j)-p(1,i))/rad_temp) - alpha;
						else
							tData(i,num) = 2*pi - acos((dp(1,j)-p(1,i))/rad_temp) - alpha;
						end
						tData(i,num)=mod(tData(i,num),2*pi);    %������λ��0-2pi��Χ��
						
						pData(1,num) = dp(1,j);
						pData(2,num) = dp(2,j);
						
						% �޳�ɨ�費���ĵ�
						InvalidFlag = 0;
						if(num>1)
							for kk=1:num-1
								if(abs(tData(i,num)-tData(i,kk))<=pi/18000)
									if(disData(i,num)<disData(i,kk))
										tData(i,kk) = tData(i,num);
										pData(1,kk) = pData(1,num);
										pData(2,kk) = pData(2,num);
										InvalidFlag = 1;
										num = num - 1;
										break;
									end
								end
							end
                        end
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        if(InvalidFlag==0)
                            % �²�����򣺽��2����������Ŀ�����ת�˶�
                            if(sqrt((pt(1,1)-p(1,i))^2+(pt(2,1)-p(2,i))^2)<=sqrt((dp(1,j)-p(1,i))^2+(dp(2,j)-p(2,i))^2))
                                InvalidFlag = 1;
                            end
                            % �²�����򣬽��û�е���Ŀ��㣬���ﵽƽ��ʱ��ǿ��ת������״̬
                            if((sqrt((pt(1,1)-p(1,i))^2+(pt(2,1)-p(2,i))^2)>detectR)&&(sqrt((pn0(1,i)-p(1,i))^2+(pn0(2,i)-p(2,i))^2)<0.05))
                                SCount(i) = SCount(i) + 1;
                            else
                                SCount(i) = 0;
                            end
                            if(SCount(i)>= 5)
                                SCount(i) = 0;
                                ExistObsFlag = 1;
                            end
                        end
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
						if((InvalidFlag==0)&&((tData(i,num)>=359*pi/180)||(tData(i,num)<=pi/180)))	%���ϰ����赲
							ExistObsFlag = 1;
                        end
					end
				end
				
				if(ExistObsFlag==0)	%���ϰ����û���赲����Ŀ�������쵼��
					num = 0;				
				elseif((num>=1)&&(ExistObsFlag==1))	%���ϰ����赲����Ŀ�������쵼��
					Sip(1,i) = p(1,i);
					Sip(2,i) = p(2,i);
					Sip(3,i) = ld;
					FlyState(i) = 1;
				end
				
			case 1	% �������������쵼��
                if(pt(2,1)-p(2,i)>=0)
					alpha = acos((pt(1,1)-p(1,i))/sqrt((pt(1,1)-p(1,i))^2+(pt(2,1)-p(2,i))^2));
				else
					alpha = 2*pi - acos((pt(1,1)-p(1,i))/sqrt((pt(1,1)-p(1,i))^2+(pt(2,1)-p(2,i))^2));
                end
				%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % �²������:��ֹ����ʱ���ų��������ϰ���޷��ҵ�С�ڷ��ѵ��λ�ã�������ѭ��
                if(test_num(1,i)>=500)
                    test_num(1,i) = 0;
                    NlineFlag = 1;
                    FlyState(i) = 0;
                else
                    test_num(1,i) = test_num(1,i) + 1;
                end
%                 if(pt(2,1)-Sip(2,i)>=0)
% 					vtemp = mod(acos((pt(1,1)-Sip(1,i))/sqrt((pt(1,1)-Sip(1,i))^2+(pt(2,1)-Sip(2,i))^2))-alpha,2*pi);
% 				else
% 					vtemp = mod(2*pi-acos((pt(1,1)-Sip(1,i))/sqrt((pt(1,1)-Sip(1,i))^2+(pt(2,1)-Sip(2,i))^2))-alpha,2*pi);
%                 end
%                 
%                 if((vtemp>=330*pi/180)||(vtemp<=pi*30/180))
%                     Dppt(1,i) = Dppt(2,i); 
%                     Dppt(2,i) = Dppt(3,i); 
%                     Dppt(3,i) = sqrt((pt(1,1)-p(1,i))^2+(pt(2,1)-p(2,i))^2);
%                     if((Dppt(3,i)>Dppt(2,i))&&(Dppt(2,i)<Dppt(1,i)))
%                         NlineFlag = 1;
%                         FlyState(i) = 0;
%                     end
%                 end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				% �������Ƿ��С
                if(sqrt((pt(1,1)-p(1,i))^2 + (pt(2,1)-p(2,i))^2)<=sqrt((pt(1,1)-Sip(1,i))^2 + (pt(2,1)-Sip(2,i))^2))
					% �������ѵ�ʱ���������ϰ��ʹ����ϵx��������ָ��Ŀ�������쵼��
					NlineFlag = 1;
					FlyState(i) = 0;
                end
			
				num = 0;
				ExistObsFlag = 0;				
				for j=1:totalNum
					if(i==j-obsnum)
						rad_temp = inf;	%������ų������屾��
					else
						rad_temp = sqrt((dp(1,j)-p(1,i))^2 + (dp(2,j)-p(2,i))^2);
					end
					
					% ������λ
					if(rad_temp<=detectR)
						num = num + 1;
						disData(i,num) = rad_temp;	% �洢����
						
						if(dp(2,j)-p(2,i)>=0)
							tData(i,num) = acos((dp(1,j)-p(1,i))/rad_temp) - alpha;
						else
							tData(i,num) = 2*pi - acos((dp(1,j)-p(1,i))/rad_temp) - alpha;
						end
						tData(i,num)=mod(tData(i,num),2*pi);    %������λ��0-2pi��Χ��
						
						pData(1,num) = dp(1,j);
						pData(2,num) = dp(2,j);

						% �޳�ɨ�費���ĵ�
						InvalidFlag = 0;
						if(num>1)
							for kk=1:num-1
								if(abs(tData(i,num)-tData(i,kk))<=pi/18000)
									if(disData(i,num)<disData(i,kk))
										tData(i,kk) = tData(i,num);
										pData(1,kk) = pData(1,num);
										pData(2,kk) = pData(2,num);
										InvalidFlag = 1;
										num = num - 1;
										break;
									end
								end
							end
						end

						if(((tData(i,num)>=359*pi/180)||(tData(i,num)<=pi/180))&&(NlineFlag==1)&&(InvalidFlag==0))	% ���ϰ����赲
                            Sip(1,i) = p(1,i);
                            Sip(2,i) = p(2,i);
							Sip(3,i) = ld;
							ExistObsFlag = 1;
						end
					end
				end				
				if((ExistObsFlag==0)&&(NlineFlag==1))	%�ڵ����ߴ���û���ϰ����赲����Ŀ�������쵼��
					Rip(1,i) = p(1,i);
                    Rip(2,i) = p(2,i);
					Rip(3,i) = ld;
					num = 0;
				end
		end
		%------------------------------------------------------%
		NobjectFlag(i) = 0;
		if(num>=1)
			prFlag = 0;
			if(num==1)
				pf(1,i) = pData(1,num);
				pf(2,i) = pData(2,num);
			else
				[Res,index] =sort(tData(i,1:num),'ascend');		% �������У�Res��������index����ֵ
				cPNum = 20;	% �ϰ��������жϱ�׼
                if(num>cPNum)
                    tempRes = [Res(1,num-cPNum+2:num)-2*pi,Res];
                    tempNum = 0;
                    for k=2:num+cPNum-1
                        tempNum = tempNum + 1;
						tempAngle = tempRes(1,k)-tempRes(1,k-1);
                        if(tempAngle>=pi/6)
    						if(tempNum>=cPNum)		% ��Ѱ�������ϰ���
    							tempNum = 0; 
								pf(1,i) = pData(1,index(k-cPNum));
								pf(2,i) = pData(2,index(k-cPNum));
								prFlag = 1;
								break;
                            else
                                tempNum = 0; 
								prFlag = 0;
                            end
                        end
                    end
					if((prFlag==0)&&(tempNum~=0))
						pf(1,i) = pData(1,index(num));
						pf(2,i) = pData(2,index(num));
					elseif((prFlag==0)&&(tempNum==0))	% û����Ѱ�������ϰ���
						for k=2:num
							if((Res(1,k)-Res(1,k-1))>=pi/6)
								pf(1,i) = pData(1,index(k-1));
								pf(2,i) = pData(2,index(k-1));
								prFlag = 1;
								Sip(1,i) = p(1,i);
								Sip(2,i) = p(2,i);
								Sip(3,i) = ld;
								break;
							end
						end
						if((k==num)&&(prFlag==0))
							pf(1,i) = pData(1,index(k));
							pf(2,i) = pData(2,index(k));
						end
					end
                else
                    for k=2:num
                        if((Res(1,k)-Res(1,k-1))>=pi/6)
                            pf(1,i) = pData(1,index(k-1));
                            pf(2,i) = pData(2,index(k-1));
                            prFlag = 1;
                            Sip(1,i) = p(1,i);
                            Sip(2,i) = p(2,i);
                            Sip(3,i) = ld;
                            break;
                        end
                    end
					if((k==num)&&(prFlag==0))
						pf(1,i) = pData(1,index(k));
						pf(2,i) = pData(2,index(k));
					end
				end
            end     
%********************************************************************************************%
			rTemp = sqrt((pf(1,i)-p(1,i))^2+(pf(2,i)-p(2,i))^2);
			if(pf(2,i)-p(2,i)>=0)
				thetaTemp = acos((pf(1,i)-p(1,i))/rTemp);
			else
				thetaTemp = 2*pi - acos((pf(1,i)-p(1,i))/rTemp);
			end
			thetaTemp = mod(thetaTemp,2*pi);    %������λ��0-2pi��Χ��
			thetaTemp = thetaTemp + pi/4;
			
			pf(1,i) = p(1,i) + 1 * rTemp * cos(thetaTemp);
			pf(2,i) = p(2,i) + 1 * rTemp * sin(thetaTemp);
%********************************************************************************************%			
		else
			if((FlyState(i)==1)&&(NlineFlag==0))	%�����ϰ���ʱ�������屻�ų�����֪��Χû���ϰ�����������쵼�߱�����һʱ��λ��
				pf(1,i) = pf(1,i);
				pf(2,i) = pf(2,i);
			else
				pf(1,i) = pt(1,1);
				pf(2,i) = pt(2,1);
                NobjectFlag(i) = 1;
			end
        end
      
		% �����쵼���ٶȼ���
		vf(1,i) = (pf(1,i)-p(1,i))/step;
		vf(2,i) = (pf(2,i)-p(2,i))/step;
	end
    %***********************************************************************%
    %-----------------------------------------------------------------------%
	% �����С���
	Colli_beta=zeros(n,totalNum);
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
            % �ų���
            if(d_beta(i,k)<=2)
                Colli_beta(i,k) = 1;
            end
        end
    end
	%----------------------------------------------------------------%
	% ���ϰ����Ӱ��
    % 1. ���干ʶ����A���ж�������֮����໥Ӱ��
    A_beta=zeros(n,totalNum);
    for i=1:n
        for k=1:totalNum
            if((d_beta(i,k)<=r)&&(i~=k-obsnum))
                A_beta(i,k)=1;
            end
        end    
    end
    % ʵ��fya(z)
    % ����n(i,j)
    N_beta=zeros(n,totalNum,2);
    for i=1:n
        for k=1:totalNum
            N_beta(i,k,1)=(dp(1,k)-p(1,i))/sqrt(1+s*d_beta(i,k)^2);
            N_beta(i,k,2)=(dp(2,k)-p(2,i))/sqrt(1+s*d_beta(i,k)^2);
        end
    end
    % ����da=||pj-pi||��
    da_beta=zeros(n,totalNum);
    for i=1:n
        for k=1:totalNum
            da_beta(i,k)=(1/s)*(sqrt(1+s*d_beta(i,k)^2)-1);
        end
    end
    % ����fya(da)
    fya_beta=zeros(n,totalNum);
    ph_beta=zeros(n,totalNum);
    for i=1:n
        for k=1:totalNum
			z1_beta=da_beta(i,k)/ra;
			if(z1_beta<h && z1_beta>=0)      
				ph_beta(i,k)=1;
			elseif (z1_beta<=1 || z1_beta>=h)
				ph_beta(i,k)=0.5*(1+cos(pi*((z1_beta-h)/(1-h))));
			else
				ph_beta(i,k)=0;
			end
			z2_beta=da_beta(i,k)-dwa;	% dwa��������
			c_beta=abs(a_beta-b_beta)/sqrt(4*a_beta*b_beta);
			fy=((a_beta+b_beta)*((z2_beta+c_beta)/sqrt(1+(z2_beta+c_beta)^2))+(a_beta-b_beta))/2;
			fya_beta(i,k)=ph_beta(i,k)*fy*A_beta(i,k);    % A(i,j)ʹi=jʱ��fya=0
        end
    end
    %-----------------���������λ��Ӱ��------------------%
    % ��u11=fya*N��i,j��
    u11_beta=zeros(n,totalNum,2);
    for i=1:n
        for k=1:totalNum
			if((d_beta(i,k)>dw)&&(d_beta(i,k)<=r))
				tempC = 0.1;
			elseif(d_beta(i,k)<=2)
				tempC = 10000;
			else
				tempC = 1;
			end
            u11_beta(i,k,1) = tempC * fya_beta(i,k) * N_beta(i,k,1);
            u11_beta(i,k,2) = tempC * fya_beta(i,k) * N_beta(i,k,2);
        end
    end
    % ��λ�÷���u1
    u1_beta=zeros(2,n);
    for i=1:n
    	for k=1:totalNum
			if(i~=k-obsnum)
				u1_beta(1,i) = u1_beta(1,i) + u11_beta(i,k,1);
				u1_beta(2,i) = u1_beta(2,i) + u11_beta(i,k,2);
			end
        end
    end
    %----------------����������ٶ�Ӱ��-------------------%
    %��u22=aij*(pj-pi)
    u22_beta=zeros(n,totalNum,2);
    ovij=zeros(2,1);
    for i=1:n
        for k=1:totalNum
			if((i~=k-obsnum)&&(A_beta(i,k)==1))               
                if(FlyState(i)==0)  % �����ٶȵ�λ����һ�£�����ٶ�һ�²����������� 
                    %-------------------------------------------------------------%
                    SQvf = sqrt(vf(1,i)^2 + vf(2,i)^2);
                    if(SQvf==0)
                        SQ_vf1 = 0;
                        SQ_vf2 = 0;
                    else
                    	SQ_vf1 = vf(1,i)/SQvf;
                    	SQ_vf2 = vf(2,i)/SQvf;
                    end
                    %-------------------------------------------------------------%
                    SQv = sqrt(v(1,i)^2 + v(2,i)^2);
                    if(SQv==0)
                        SQ_v1 = 0;
                        SQ_v2 = 0;
                    else
                    	SQ_v1 = v(1,i)/SQv;
                        SQ_v2 = v(2,i)/SQv;
                    end
                    %-------------------------------------------------------------%
                    u22_beta(i,k,1) = ph_beta(i,k) * (SQ_vf1 - SQ_v1) * A_beta(i,k);
                	u22_beta(i,k,2) = ph_beta(i,k) * (SQ_vf2 - SQ_v2) * A_beta(i,k);            
                else                % �����ٶȶ��룬ʵ�ֿ��ٱ���
                    u22_beta(i,k,1) = ph_beta(i,k) * (vf(1,i) - v(1,i)) * A_beta(i,k);
                    u22_beta(i,k,2) = ph_beta(i,k) * (vf(2,i) - v(2,i)) * A_beta(i,k);
                end
			end
        end
    end
	
    u2_beta=zeros(2,n);
    for i=1:n
        for k=1:totalNum
            u2_beta(1,i) = u2_beta(1,i) + u22_beta(i,k,1);
            u2_beta(2,i) = u2_beta(2,i) + u22_beta(i,k,2);
        end
    end
    %****************************************************************************%
    uBeta=zeros(2,n);
    for i=1:2
		for j=1:n
            uBeta(i,j) = beta1 * u1_beta(i,j) + beta2 * u2_beta(i,j);
		end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ���쵼�ߵ�Ӱ��
    uGama=zeros(2,n);
    for i=1:n       
		uGama(1,i) = -c1 * (p(1,i)-pt(1,1)) - c2 * (v(1,i)-vt(1,1));
		uGama(2,i) = -c1 * (p(2,i)-pt(2,1)) - c2 * (v(2,i)-vt(2,1));
    end
    % ����ٶ�u
    u=zeros(2,n);
    for i=1:2
    	for j=1:n
            u(i,j) = uBeta(i,j) + uGama(i,j);
		end
    end   
    % ����������м���
    pn0 = p;    
    vn0 = v;
    
    for j=1:n       
        for i=1:2
            v(i,j) = v(i,j) + step * u(i,j);
        end
        vtemp = sqrt(v(1,j)^2+v(2,j)^2);
        if(vtemp>vmax) 
            v(1,j) = vmax * v(1,j)/vtemp;
            v(2,j) = vmax * v(2,j)/vtemp;    
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
%     plot(pt(1,1),pt(2,1),'g>','LineWidth',1,'MarkerSize',4);
%     plot(p(1,:),p(2,:),'b*'); 
% %     for i=1:n
% %         plot(p(1,i),p(2,i),'b*'); 
% %         text(p(1,i),p(2,i),num2str(i),'Color','red','Fontsize',7);
% %     end
%     pause(0.01);
%     hold off;
end

%****************************************************************************************%

for j=1:loop
	XlsDataBuff(iters,j) = dmin(1,j);
end
end	

xlswrite(filename,XlsDataBuff);
