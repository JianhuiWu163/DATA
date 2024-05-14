close all;clc;clear all;
%***************************************************************%
Xmax=100;
Ymax=Xmax;
r=6;

n=10;
loop=3000;                      % 确定循环周期
step = 0.1;
xt = 1:loop;
xt = xt*step;
%***************************************************************%
% 搭建障碍物区域
obsnum=0;
ImpRgb = imread('fig1.png');
Imp0 = rgb2gray(ImpRgb);
[Hlen,Wlen] = size(Imp0);
Imp = ones(Hlen,Wlen);
for i=1:Hlen
    Imp(Hlen+1-i,:)=Imp0(i,:);%排列反转
end

for i=1:Hlen
    for j=1:Wlen
       if(Imp(i,j)==0) 
           obsnum = obsnum + 1;
           op(1,obsnum)=j*Xmax/Wlen;
           op(2,obsnum)=i*Ymax/Hlen;
           op(3,obsnum)=2;
       end
    end
end
%***************************************************************%
n=10;
Data_10agents = xlsread('Agent10.xlsx');
for i=1:n
    px10(i,:) = Data_10agents(2*i-1,:);
    py10(i,:) = Data_10agents(2*i,:);
end
pz10 = 2*ones(1,length(px10(1,:)));

Xpt(1,:)=Data_10agents(2*n+1,:);
Ypt(1,:)=Data_10agents(2*n+2,:);
% dmin10(1,:)=Data_10agents(2*n+3,:);
%***************************************************************%
n=50;
Data_50agents = xlsread('Agent50.xlsx');
for i=1:n
    px50(i,:) = Data_50agents(2*i-1,:);
    py50(i,:) = Data_50agents(2*i,:);
end
pz50 = 2*ones(1,length(px50(1,:)));

% Xpt(1,:)=Data_50agents(2*n+1,:);
% Ypt(1,:)=Data_50agents(2*n+2,:);
% dmin50(1,:)=Data_50agents(2*n+3,:);
%**************************************************************************%
for j=1:loop
    n=10;
    Xcp10(1,j) = sum(px10(:,j))/n;
	Ycp10(1,j) = sum(py10(:,j))/n;
    
    n=50;
    Xcp50(1,j) = sum(px50(:,j))/n;
	Ycp50(1,j) = sum(py50(:,j))/n;
end
%**************************************************************************%
n=10;
figure(1)
% set(gcf,'Position',[100 100 330 200]);%图片大小
set(gcf,'Position',[100 100 165 100]);%图片大小
set(gca,'Position',[.13 .15 .8 .8]);%坐标轴所占比例

plot3(op(1,:),op(2,:),op(3,:),'.','color',[0/255,102/255,51/255]);hold on;
axis([0 100 0 100 0 5]);
set(gca,'xtick',[0:20:100],'ytick',[0:20:100]);
set(gca,'TickDir','in','FontName','Times New Roman','FontSize',9);
xlabel('x (m)','FontName','Times New Roman','FontSize',10,'position',[50,-15,0]);
ylabel('y (m)','FontName','Times New Roman','FontSize',10,'position',[-15,50,0]);

x1=0:0.1:100; y1=100*ones(1,length(x1));
y2=0:0.1:100; x2=100*ones(1,length(y2));
plot(x1,y1,'-','color',[153/255,153/255,153/255]);hold on;
plot(x2,y2,'-','color',[153/255,153/255,153/255]);hold on;
obsnum1 = 0;
for H=0:0.01:2-0.01
    for i=1:obsnum
    	op(3,i) = H; 
    end
    plot3(op(1,:),op(2,:),op(3,:),'.','color',[255/255,204/255,51/255]);hold on;
end

ax=gca;
ax.ZAxis.Visible='off';
view(-16,85);
%------------------------------------------------------------------------%
num=1;
count=0;
for i=1:n
    for j=1:n
        if((j>i)&&(sqrt((px10(i,num)-px10(j,num))^2 + (py10(i,num)-py10(j,num))^2)<=r))
            count=count+1;
            SdX1(1,count) = px10(i,num);
            SdY1(1,count) = py10(i,num);
            SdX1(2,count) = px10(j,num);
            SdY1(2,count) = py10(j,num);
        end
    end
end
[xNode1,yNode1] = AdjacencyLine(SdX1,SdY1);
count = length(xNode1(1,:));
for i=1:count
    line([xNode1(1,i),xNode1(2,i)],[yNode1(1,i),yNode1(2,i)],[pz10(1,num),pz10(1,num)],'color',[51/255 51/255 51/255]);
end
%----------------------------------------------------------------------%
num=loop;
count=0;
for i=1:n
    for j=1:n
        if((j>i)&&(sqrt((px10(i,num)-px10(j,num))^2 + (py10(i,num)-py10(j,num))^2)<=r))
            count=count+1;
            SdX2(1,count) = px10(i,num);
            SdY2(1,count) = py10(i,num);
            SdX2(2,count) = px10(j,num);
            SdY2(2,count) = py10(j,num);
        end
    end
end
[xNode2,yNode2] = AdjacencyLine(SdX2,SdY2);
count = length(xNode2(1,:));
for i=1:count
    line([xNode2(1,i),xNode2(2,i)],[yNode2(1,i),yNode2(2,i)],[pz10(1,num),pz10(1,num)],'color',[51/255 51/255 51/255]);
end
%----------------------------------------------------------------------%
for i=1:n
    plot3(px10(i,1),py10(i,1),pz10(1,1),'ro','MarkerSize',3);
    plot3(px10(i,2:5:800),py10(i,2:5:800),pz10(1,2:5:800),'r:');
    plot3(px10(i,loop),py10(i,loop),pz10(1,loop),'ro','MarkerSize',3);
end
%----------------------------------------------------------------------%
plot3(Xpt(1,1),Ypt(1,1),2.1,'g>','LineWidth',1,'MarkerSize',3);

hold off;

img =gcf;  %获取当前画图的句柄
% print(img, '-dpng', '-r600', './Agent10.png');
print(img, '-depsc', '-r600', './fig6a.eps');

%*************************************************************************************************%
n=50;
figure(2)
set(gcf,'Position',[100 100 165 100]);%图片大小
set(gca,'Position',[.13 .15 .8 .8]);%坐标轴所占比例

for i=1:obsnum
	op(3,i) = 2; 
end
plot3(op(1,:),op(2,:),op(3,:),'.','color',[0/255,102/255,51/255]);hold on;
axis([0 100 0 100 0 5]);
set(gca,'xtick',[0:20:100],'ytick',[0:20:100]);
set(gca,'TickDir','in','FontName','Times New Roman','FontSize',9);
xlabel('x (m)','FontName','Times New Roman','FontSize',10,'position',[50,-15,0]);
ylabel('y (m)','FontName','Times New Roman','FontSize',10,'position',[-15,50,0]);

x1=0:0.1:100; y1=100*ones(1,length(x1));
y2=0:0.1:100; x2=100*ones(1,length(y2));
plot(x1,y1,'-','color',[153/255,153/255,153/255]);hold on;
plot(x2,y2,'-','color',[153/255,153/255,153/255]);hold on;
obsnum1 = 0;
for H=0:0.01:2-0.01
    for i=1:obsnum
    	op(3,i) = H; 
    end
    plot3(op(1,:),op(2,:),op(3,:),'.','color',[255/255,204/255,51/255]);hold on;
end

ax=gca;
ax.ZAxis.Visible='off';
view(-16,85);
%------------------------------------------------------------------------%
num=1;
count=0;
for i=1:n
    for j=1:n
        if((j>i)&&(sqrt((px50(i,num)-px50(j,num))^2 + (py50(i,num)-py50(j,num))^2)<=r))
            count=count+1;
            SdX1(1,count) = px50(i,num);
            SdY1(1,count) = py50(i,num);
            SdX1(2,count) = px50(j,num);
            SdY1(2,count) = py50(j,num);
        end
    end
end
[xNode1,yNode1] = AdjacencyLine(SdX1,SdY1);
count = length(xNode1(1,:));
for i=1:count
    line([xNode1(1,i),xNode1(2,i)],[yNode1(1,i),yNode1(2,i)],[pz50(1,num),pz50(1,num)],'color',[51/255 51/255 51/255]);
end
%----------------------------------------------------------------------%
num=loop;
count=0;
for i=1:n
    for j=1:n
        if((j>i)&&(sqrt((px50(i,num)-px50(j,num))^2 + (py50(i,num)-py50(j,num))^2)<=r))
            count=count+1;
            SdX2(1,count) = px50(i,num);
            SdY2(1,count) = py50(i,num);
            SdX2(2,count) = px50(j,num);
            SdY2(2,count) = py50(j,num);
        end
    end
end
[xNode2,yNode2] = AdjacencyLine(SdX2,SdY2);
count = length(xNode2(1,:));
for i=1:count
    line([xNode2(1,i),xNode2(2,i)],[yNode2(1,i),yNode2(2,i)],[pz50(1,num),pz50(1,num)],'color',[51/255 51/255 51/255]);
end
%----------------------------------------------------------------------%
for i=1:n
    plot3(px50(i,1),py50(i,1),pz50(1,1),'ro','MarkerSize',3);
    plot3(px50(i,2:5:800),py50(i,2:5:800),pz50(1,2:5:800),'r:');
    plot3(px50(i,loop),py50(i,loop),pz50(1,loop),'ro','MarkerSize',3);
end
%----------------------------------------------------------------------%
plot3(Xpt(1,1),Ypt(1,1),2.1,'g>','LineWidth',1,'MarkerSize',3);

hold off;

img =gcf;  %获取当前画图的句柄
% print(img, '-dpng', '-r600', './Agent50.png')
print(img, '-depsc', '-r600', './fig6b.eps');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(3)
% set(gcf,'Position',[100 100 320 200]);%图片大小
% set(gca,'Position',[.14 .22 .82 .7]);%坐标轴所占比例
% 
% hp1 = plot(xt,Xcp10,'r',xt,Xcp50,'b',xt,Xpt,'g');hold on;
% axis([0 loop*step 0 120]);
% set(gca,'XTickMode','manual','Xtick',[0:loop*step/5:loop*step]);
% set(gca,'YTickMode','manual','Ytick',[0:20:120]);
% xlabel('t (s)','Fontname','times new Roman','Fontsize',10);
% ylabel('x (m)','Fontname','times new Roman','Fontsize',10);
% set(gca,'Fontname','times new Roman','Fontsize',9);%坐标刻度字体大小
% 
% hleg1=legend(hp1(1:3),'Center position of 10 \alpha-agents','Center position of 50 \alpha-agents','Position of a \gamma-agent');
% set(hleg1,'position', [0.62,.77,.1,.1],'Fontname','times new Roman','fontsize',9);
% legend boxoff;
% box off;
% hold off;
% 
% img =gcf;  %获取当前画图的句柄
% print(img, '-depsc', '-r600', './fig6c.eps');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(4)
% set(gcf,'Position',[100 100 320 200]);%图片大小
% set(gca,'Position',[.14 .22 .82 .7]);%坐标轴所占比例
% 
% hp1 = plot(xt,Ycp10,'r',xt,Ycp50,'b',xt,Ypt,'g');hold on;
% axis([0 loop*step 0 120]);
% set(gca,'XTickMode','manual','Xtick',[0:loop*step/5:loop*step]);
% set(gca,'YTickMode','manual','Ytick',[0:20:120]);
% xlabel('t (s)','Fontname','times new Roman','Fontsize',10);
% ylabel('y (m)','Fontname','times new Roman','Fontsize',10);
% set(gca,'Fontname','times new Roman','Fontsize',9);%坐标刻度字体大小
% box off;
% 
% hleg1=legend(hp1(1:3),'Center position of 10 \alpha-agents','Center position of 50 \alpha-agents','Position of a \gamma-agent');
% set(hleg1,'position', [0.62,.77,.1,.1],'Fontname','times new Roman','fontsize',9);
% legend boxoff;
% 
% img =gcf;  %获取当前画图的句柄
% print(img, '-depsc', '-r600', './fig6d.eps');

%***************************************************************%
% times = 100;
% xt = 1:times;
% dmin_agent10 = xlsread('Agent10_Data.xlsx');
% dmin_agent20 = xlsread('Agent20_Data.xlsx');
% dmin_agent30 = xlsread('Agent30_Data.xlsx');
% dmin_agent40 = xlsread('Agent40_Data.xlsx');
% dmin_agent50 = xlsread('Agent50_Data.xlsx');
% for i=1:times
%     dmin1(1,i) = min(dmin_agent10(i,:));
%     dmin2(1,i) = min(dmin_agent20(i,:));
%     dmin3(1,i) = min(dmin_agent30(i,:));
%     dmin4(1,i) = min(dmin_agent40(i,:));
%     dmin5(1,i) = min(dmin_agent50(i,:));
%     d5(1,i) = 1;
% end
% 
% figure(5)
% % set(gcf,'Position',[100 100 320 200]);%图片大小
% set(gcf,'Position',[100 100 280 160]);%图片大小
% set(gca,'Position',[.14 .22 .82 .7]);%坐标轴所占比例
% hp1 = plot(xt,dmin1,'r',xt,dmin2,'g',xt,dmin3,'b',xt,dmin4,'m',xt,dmin5,'c');
% axis([0 times 0.9 2.2]);
% 
% set(gca,'XTickMode','manual','Xtick',[0:times/5:times]);
% set(gca,'YTickMode','manual','Ytick',[1:0.3:2.2]);
% xlabel('Simulation times','Fontname','times new Roman','Fontsize',10);
% ylabel('Minimum distance (m)','Fontname','times new Roman','Fontsize',10);
% set(gca,'Fontname','times new Roman','Fontsize',9);%坐标刻度字体大小
% grid on;
% box off;
% 
% hleg1=legend(hp1(1:3),'10 \alpha-agents','20 \alpha-agents','30 \alpha-agents');
% hleg1.ItemTokenSize = [15,9];
% set(hleg1,'position', [0.47,.755,.1,.1],'Fontname','times new Roman','fontsize',9);
% legend boxoff;
% 
% ah=axes('position',get(gca,'position'),'visible','off');
% hleg2=legend(ah,hp1(4:5),'40 \alpha-agents','50 \alpha-agents');
% hleg2.ItemTokenSize = [15,9];
% set(hleg2,'position', [0.76,.798,.1,.1],'Fontname','times new Roman','fontsize',9);
% legend boxoff;
% 
% img =gcf;  %获取当前画图的句柄
% print(img, '-depsc', '-r600', './fig7.eps');
