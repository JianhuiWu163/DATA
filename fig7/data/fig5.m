clc;close all;clear all;
% ----------------------------------------------------------------------%
filename = 'SA_10_16.xlsx';
for Dbeta=16:26
    filename(7)=fix(Dbeta/10)+48;
    filename(8)=mod(Dbeta,10)+48;  
    SA_10(Dbeta-15,:,:) = xlsread(filename);
end

for i=1:11
    for j=1:100
        dmin_sa10(i,j) = min(SA_10(i,j,:));
    end
end

for i=1:11
    Min_sa10(1,i) = min(dmin_sa10(i,:));
    Max_sa10(1,i) = max(dmin_sa10(i,:));
    Average_sa10(1,i) = mean(dmin_sa10(i,:));
end
% ----------------------------------------------------------------------%
filename = 'SA_50_16.xlsx';
for Dbeta=16:26
    filename(7)=fix(Dbeta/10)+48;
    filename(8)=mod(Dbeta,10)+48;  
    SA_50(Dbeta-15,:,:) = xlsread(filename);
end

for i=1:11
    for j=1:100
        dmin_sa50(i,j) = min(SA_50(i,j,:));
    end
end

for i=1:11
    Min_sa50(1,i) = min(dmin_sa50(i,:));
    Max_sa50(1,i) = max(dmin_sa50(i,:));
    Average_sa50(1,i) = mean(dmin_sa50(i,:));
end
% ----------------------------------------------------------------------%
filename = 'PA_10_16.xlsx';
for Dbeta=16:26
    filename(7)=fix(Dbeta/10)+48;
    filename(8)=mod(Dbeta,10)+48;  
    PA_10(Dbeta-15,:,:) = xlsread(filename);
end

for i=1:11
    for j=1:100
        dmin_pa10(i,j) = min(PA_10(i,j,:));
    end
end

for i=1:11
    Min_pa10(1,i) = min(dmin_pa10(i,:));
    Max_pa10(1,i) = max(dmin_pa10(i,:));
    Average_pa10(1,i) = mean(dmin_pa10(i,:));
end
% ----------------------------------------------------------------------%
filename = 'PA_50_16.xlsx';
for Dbeta=16:26
    filename(7)=fix(Dbeta/10)+48;
    filename(8)=mod(Dbeta,10)+48;  
    PA_50(Dbeta-15,:,:) = xlsread(filename);
end

for i=1:11
    for j=1:100
        dmin_pa50(i,j) = min(PA_50(i,j,:));
    end
end

for i=1:11
    Min_pa50(1,i) = min(dmin_pa50(i,:));
    Max_pa50(1,i) = max(dmin_pa50(i,:));
    Average_pa50(1,i) = mean(dmin_pa50(i,:));
end
% ----------------------------------------------------------------------%
% Average_pa10=[1.43308313222921,1.44480277501722,1.47305038541180,1.30751386061473,1.37506410476657,1.39692328377555,1.30621961934308,1.49611360901088,1.58241266716548,1.63151341549337,1.67619173986962];
% Average_pa50=[1.03320693810483,1.03189198286458,1.03171307017727,1.03099926626048,1.03316962219877,1.03598688958327,1.03111177161797,1.03343903601097,1.03845563672815,1.04169223235829,1.04890395982251];
% Average_sa10=[0.873985790663065,0.494474524094495,0.921638856721687,1.29908396401135,1.45525378633019,1.46928740050742,1.47997359055663,1.48823591554748,1.49933388594336,1.47451932845427,1.43525959603292];
% Average_sa50=[0.0197480623357118,0.441324904237252,0.179584998670045,0.114811956337403,0.0905733174908577,0.0744120040771516,0.0752254088948146,0.0732116441033098,0.0816107170127697,0.0765906967836158,0.0832616584906542];
% Min_pa10=[1.03135966606988,1.01462547670597,1.02420923562580,1.02339185448032,1.01603459154358,1.06171873726678,1.03950754539619,1.00663476507188,1.03745139241078,1.03745139241078,1.13932986580174];
% Min_pa50=[1.00497720191540,1.00439010416295,1.00226374475987,1.00270489513923,1.00094582265117,0.997056178273503,1.00146492703661,1.00325279538516,1.00169283139820,1.00322945494730,1.00331831642246];
% Min_sa10=[0.753547406201486,0.493394039428456,0.867547873158664,1.14217633885606,1.18543288451665,1.18500986014933,1.24439351468784,1.25248945028811,1.22758835166546,1.22758835166546,1.18580453228921];
% Min_sa50=[8.76669107098986e-05,0.0746901136115677,0.0172792851755248,0.0106027422588111,0.00888797505494625,0.0179345952184569,0.00465814803572305,0.00179291243601027,0.00594833646623793,0.00582688211490692,0.00659351098683391];
% Max_pa10=[1.87616858880738,1.88661735997278,1.88045562415200,1.69839285356007,1.84881738707535,1.86790905169831,1.84000322703564,1.90701536722314,2.14053390055852,2.07123574056795,2.08654605739567];
% Max_pa50=[1.11195041608066,1.09601248505032,1.09218846567149,1.07784066853879,1.10028100115970,1.12209250678579,1.09073860000802,1.08756477146562,1.13889616788579,1.12405546648878,1.16297481000885];
% Max_sa10=[1.05934853550699,0.496009909909634,0.986129215837018,1.48149054686472,1.68513978313083,1.69743097951952,1.74075487130760,1.74679678422615,1.75684703573056,1.67595167089020,1.67710914638153];
% Max_sa50=[0.0612057445775500,0.481827285308155,0.870056562628917,0.323900528872260,0.848757498565230,0.199417555109354,0.235391246675664,0.209446189171590,0.314082553898713,0.236152086231202,0.250369763175690];

% %-------------------------------------------------------------------------%
figure(1)
set(gcf,'Position',[100 100 150 110]);%图片大小
set(gca,'Position',[.17 .29 .74 .64]);%坐标轴所占比例

hp1 = plot([16:26],Average_sa10,'o:','LineWidth',1,'Markersize',1,'MarkerFaceColor',[255/255,102/255,0/255],'Color',[255/255,102/255,0/255]);hold on;
hp2 = plot([16:26],Average_sa50,'o:','LineWidth',1,'Markersize',1,'MarkerFaceColor',[0/255,255/255,102/255],'Color',[0/255,255/255,102/255]);
axis([16 26 0 3]);
set(gca,'xtick',[16:2:26],'ytick',[0:1:3],'FontName','Times New Roman','FontSize',7);
xlabel('d_\beta (m)','FontName','Times New Roman','FontSize',7);
ylabel('Minimum distance (m)','FontName','Times New Roman','FontSize',7);
grid on;
m1=fliplr(Max_sa10);%最大值反向
x1=[16:26 26:-1:16];
y1=[Min_sa10 m1];
h1=fill(x1,y1,[255/255,102/255,0/255]);
hold on;
set(h1,'edgealpha',0,'facealpha',0.2);

m2=fliplr(Max_sa50);%最大值反向
x2=[16:26 26:-1:16];
y2=[Min_sa50 m2];
h2=fill(x2,y2,[0/255,255/255,102/255]);
hold on;
set(h2,'edgealpha',0,'facealpha',0.2);

text(20.1,1.3,'\uparrow','Color','b','rotation',350);
text(17,0.75,'Average value','Color','b','FontName','Times New Roman','FontSize',7);
text(23,1.05,'\uparrow','Color','b','rotation',10);
text(22.3,0.5,'Min value','Color','b','FontName','Times New Roman','FontSize',7);
text(19.5,1.96,'\downarrow','Color','b','rotation',-10);
text(18,2.4,'Max value','Color','b','FontName','Times New Roman','FontSize',7);

hleg1=legend([hp1;hp2],'N=10','N=50');
hleg1.ItemTokenSize = [15,9];
set(hleg1,'position', [0.71,.78,.1,.1],'Fontname','times new Roman','fontsize',7);
legend boxoff;
box off;

img =gcf;  %获取当前画图的句柄
print(img, '-depsc', '-r600', './fig5a.eps');
% %-------------------------------------------------------------------------%
figure(2)
set(gcf,'Position',[100 100 150 110]);%图片大小
set(gca,'Position',[.17 .29 .74 .64]);%坐标轴所占比例

hp1 = plot([16:26],Average_pa10,'o:','LineWidth',1,'Markersize',1,'MarkerFaceColor',[255/255,102/255,0/255],'Color',[255/255,102/255,0/255]);hold on;
hp2 = plot([16:26],Average_pa50,'o:','LineWidth',1,'Markersize',1,'MarkerFaceColor',[0/255,255/255,102/255],'Color',[0/255,255/255,102/255]);
axis([16 26 0 3]);
set(gca,'xtick',[16:2:26],'ytick',[0:1:3],'FontName','Times New Roman','FontSize',7);
xlabel('d_\beta (m)','FontName','Times New Roman','FontSize',7);
ylabel('Minimum distance (m)','FontName','Times New Roman','FontSize',7);
grid on;

m1=fliplr(Max_pa10);%最大值反向
x1=[16:26 26:-1:16];
y1=[Min_pa10 m1];
h1=fill(x1,y1,[255/255,102/255,0/255]);
hold on;
set(h1,'edgealpha',0,'facealpha',0.2);

m2=fliplr(Max_pa50);%最大值反向
x2=[16:26 26:-1:16];
y2=[Min_pa50 m2];
h2=fill(x2,y2,[0/255,255/255,102/255]);
hold on;
set(h2,'edgealpha',0,'facealpha',0.2);

% text(22,1.05,'\uparrow','Color','b','rotation',30);
% text(22,0.75,'Min','Color','b','FontName','Times New Roman','FontSize',8);
% text(20.1,1.8,'\downarrow','Color','b','rotation',-30);
% text(19.7,2.1,'Average','Color','b','FontName','Times New Roman','FontSize',8);
% text(21.8,2,'\downarrow','Color','b','rotation',-30);
% text(21.7,2.25,'Max','Color','b','FontName','Times New Roman','FontSize',8);

hleg1=legend([hp1;hp2],'N=10','N=50');
hleg1.ItemTokenSize = [15,9];
set(hleg1,'position', [0.71,.78,.1,.1],'Fontname','times new Roman','fontsize',7);
legend boxoff;
box off;

img =gcf;  %获取当前画图的句柄
print(img, '-depsc', '-r600', './fig5b.eps');

