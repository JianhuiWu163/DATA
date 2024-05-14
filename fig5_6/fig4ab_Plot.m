close all;clc;clear all;
%*************************************************************************%
% Sakai algorithm: Flocking for Multirobots Without Distinguishing Robots and Obstacles
%*************************************************************************%
% 搭建障碍物区域
oRa = 8; 

Dbeta = 16;
OX = 35;
Dd = 26;

oyk(1,1) = OX;                              oyk(2,1) = 21;
oyk(1,2) = OX;                              oyk(2,2) = oyk(2,1) + Dd;
oyk(1,3) = OX;                              oyk(2,3) = oyk(2,2) + Dd;

oyk(1,4) = OX + sqrt(Dbeta^2-(Dd/2)^2);     oyk(2,4) = oyk(2,1) + Dd/2;
oyk(1,5) = OX + sqrt(Dbeta^2-(Dd/2)^2);     oyk(2,5) = oyk(2,2) + Dd/2;

obsnum = 0;
dtheta = pi/300;
for i=1:5
    for oRad = 0:dtheta:2*pi-dtheta
        obsnum = obsnum + 1;
        op(1,obsnum) = oyk(1,i) + oRa * cos(oRad); 
        op(2,obsnum) = oyk(2,i) + oRa * sin(oRad); 
        op(3,obsnum) = 2; 
    end
end

Xmax=100;
Ymax=Xmax;
r=6;

n=10;
loop=3000;                      % 确定循环周期
step = 0.1;
xt = 1:loop;
xt = xt*step;

SakaiData = xlsread('SakaiAlgorithm_16.xlsx');
for i=1:n
    Spx10(i,:) = SakaiData(2*i-1,:);
    Spy10(i,:) = SakaiData(2*i,:);
end
Spz10 = 2*ones(1,length(Spx10(1,:)));

SXpt(1,:)=SakaiData(2*n+1,:);
SYpt(1,:)=SakaiData(2*n+2,:);
Sdmin10(1,:)=SakaiData(2*n+3,:);
%**************************************************************************%
figure(1)
set(gcf,'Position',[100 100 150 110]);%图片大小
set(gca,'Position',[.18 .15 .7 .8]);%坐标轴所占比例

plot3(op(1,:),op(2,:),op(3,:),'.','color',[0/255,102/255,51/255],'MarkerSize',3);hold on;
axis([0 100 0 100 0 5]);
set(gca,'xtick',[0:20:100],'ytick',[0:20:100]);
set(gca,'TickDir','in','FontName','Times New Roman','FontSize',7);
xlabel('x (m)','FontName','Times New Roman','FontSize',8,'position',[50,-18,0]);
ylabel('y (m)','FontName','Times New Roman','FontSize',8,'position',[-18,50,0]);

x1=0:0.1:100; y1=100*ones(1,length(x1));
y2=0:0.1:100; x2=100*ones(1,length(y2));
plot(x1,y1,'-','color',[153/255,153/255,153/255]);hold on;
plot(x2,y2,'-','color',[153/255,153/255,153/255]);hold on;
obsnum1 = 0;
for H=0:0.01:2-0.01
    for i=1:5
        for oRad = 0:dtheta:2*pi-dtheta
            obsnum1 = obsnum1 + 1;
            op1(1,obsnum1) = oyk(1,i) + (oRa-0.1) * cos(oRad); 
            op1(2,obsnum1) = oyk(2,i) + (oRa-0.1) * sin(oRad); 
            op1(3,obsnum1) = H; 
        end
    end
end

plot3(op1(1,:),op1(2,:),op1(3,:),'.','color',[255/255,204/255,51/255],'MarkerSize',1);hold on;
ax=gca;
ax.ZAxis.Visible='off';
view(-16,85);
%------------------------------------------------------------------------%
num=1;
count=0;
for i=1:n
    for j=1:n
        if((j>i)&&(sqrt((Spx10(i,num)-Spx10(j,num))^2 + (Spy10(i,num)-Spy10(j,num))^2)<=r))
            count=count+1;
            SdX1(1,count) = Spx10(i,num);
            SdY1(1,count) = Spy10(i,num);
            SdX1(2,count) = Spx10(j,num);
            SdY1(2,count) = Spy10(j,num);
        end
    end
end
[xNode1,yNode1] = AdjacencyLine(SdX1,SdY1);
count = length(xNode1(1,:));
for i=1:count
    line([xNode1(1,i),xNode1(2,i)],[yNode1(1,i),yNode1(2,i)],[Spz10(1,num),Spz10(1,num)],'color',[51/255 51/255 51/255]);
end
%----------------------------------------------------------------------%
num=loop;
count=0;
for i=1:n
    for j=1:n
        if((j>i)&&(sqrt((Spx10(i,num)-Spx10(j,num))^2 + (Spy10(i,num)-Spy10(j,num))^2)<=r))
            count=count+1;
            SdX2(1,count) = Spx10(i,num);
            SdY2(1,count) = Spy10(i,num);
            SdX2(2,count) = Spx10(j,num);
            SdY2(2,count) = Spy10(j,num);
        end
    end
end
[xNode2,yNode2] = AdjacencyLine(SdX2,SdY2);
count = length(xNode2(1,:));
for i=1:count
    line([xNode2(1,i),xNode2(2,i)],[yNode2(1,i),yNode2(2,i)],[Spz10(1,num),Spz10(1,num)],'color',[51/255 51/255 51/255]);
end
%----------------------------------------------------------------------%
for i=1:n
    plot3(Spx10(i,1),Spy10(i,1),Spz10(1,1),'ro','MarkerSize',1);
    plot3(Spx10(i,2:5:1800),Spy10(i,2:5:1800),Spz10(1,2:5:1800),'r:');
    plot3(Spx10(i,loop),Spy10(i,loop),Spz10(1,loop),'ro','MarkerSize',1);
end
%----------------------------------------------------------------------%
plot3(SXpt(1,1),SYpt(1,1),2.1,'g>','LineWidth',1,'MarkerSize',1);

% text(7,39,0,'\downarrow','Color','b','FontName','Times New Roman','FontSize',14);
% text(4,53,0,'\alpha-agent','Color','b','FontName','Times New Roman','FontSize',8);
% 
% text(20,83,0,'\rightarrow','Color','b','FontName','Times New Roman','FontSize',14);
% text(3,86.5,0,'Obstacle','Color','b','FontName','Times New Roman','FontSize',8);
% 
% text(76.5,73.5,2.1,'\uparrow','Color','b','FontName','Times New Roman','FontSize',14);
% text(73,69,0,'\gamma-agent','Color','b','FontName','Times New Roman','FontSize',8);

hold off;

img =gcf;  %获取当前画图的句柄
% print(img, '-dpng', '-r600', './SakaiAlgorithm_16.png')
print(img, '-depsc', '-r600', './fig4a.eps');
%*************************************************************************************************%

Data_10agents = xlsread('ProposedAlgorithm_16.xlsx');
for i=1:n
    px10(i,:) = Data_10agents(2*i-1,:);
    py10(i,:) = Data_10agents(2*i,:);
end
pz10 = 2*ones(1,length(px10(1,:)));

Xpt(1,:)=Data_10agents(2*n+1,:);
Ypt(1,:)=Data_10agents(2*n+2,:);
dmin10(1,:) = Data_10agents(2*n+3,:);
%**************************************************************************%
figure(2)
set(gcf,'Position',[100 100 150 110]);%图片大小
set(gca,'Position',[.18 .15 .7 .8]);%坐标轴所占比例

plot3(op(1,:),op(2,:),op(3,:),'.','color',[0/255,102/255,51/255],'MarkerSize',3);hold on;
axis([0 100 0 100 0 5]);
set(gca,'xtick',[0:20:100],'ytick',[0:20:100]);
set(gca,'TickDir','in','FontName','Times New Roman','FontSize',7);
xlabel('x (m)','FontName','Times New Roman','FontSize',8,'position',[50,-18,0]);
ylabel('y (m)','FontName','Times New Roman','FontSize',8,'position',[-18,50,0]);

x1=0:0.1:100; y1=100*ones(1,length(x1));
y2=0:0.1:100; x2=100*ones(1,length(y2));
plot(x1,y1,'-','color',[153/255,153/255,153/255]);hold on;
plot(x2,y2,'-','color',[153/255,153/255,153/255]);hold on;

plot3(op1(1,:),op1(2,:),op1(3,:),'.','color',[255/255,204/255,51/255],'MarkerSize',1);hold on;
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
            dX1(1,count) = px10(i,num);
            dY1(1,count) = py10(i,num);
            dX1(2,count) = px10(j,num);
            dY1(2,count) = py10(j,num);
        end
    end
end
[xNode1,yNode1] = AdjacencyLine(dX1,dY1);
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
            dX(1,count) = px10(i,num);
            dY(1,count) = py10(i,num);
            dX(2,count) = px10(j,num);
            dY(2,count) = py10(j,num);
        end
    end
end
[xNode,yNode] = AdjacencyLine(dX,dY);
count = length(xNode(1,:));
for i=1:count
    line([xNode(1,i),xNode(2,i)],[yNode(1,i),yNode(2,i)],[pz10(1,num),pz10(1,num)],'color',[51/255 51/255 51/255]);
end
%----------------------------------------------------------------------%
for i=1:n
    plot3(px10(i,1),py10(i,1),pz10(1,1),'ro','MarkerSize',1);
    plot3(px10(i,2:5:300),py10(i,2:5:300),pz10(1,2:5:300),'r:');
    plot3(px10(i,loop),py10(i,loop),pz10(1,loop),'ro','MarkerSize',1);
end
%----------------------------------------------------------------------%
plot3(Xpt(1,1),Ypt(1,1),2.1,'g>','LineWidth',1,'MarkerSize',1);

% text(7,39,0,'\downarrow','Color','b','FontName','Times New Roman','FontSize',14);
% text(4,53,0,'\alpha-agent','Color','b','FontName','Times New Roman','FontSize',8);
% 
% text(20,83,0,'\rightarrow','Color','b','FontName','Times New Roman','FontSize',14);
% text(3,86.5,0,'Obstacle','Color','b','FontName','Times New Roman','FontSize',8);
% 
% text(76.5,73.5,2.1,'\uparrow','Color','b','FontName','Times New Roman','FontSize',14);
% text(73,69,0,'\gamma-agent','Color','b','FontName','Times New Roman','FontSize',8);


hold off;

img =gcf;  %获取当前画图的句柄
% print(img, '-dpng', '-r600', './ProposedAlgorithm_16.png')
print(img, '-depsc', '-r600', './fig4b.eps');


