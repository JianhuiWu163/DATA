%===========================================================================%
% LineX(1,:):线条起始节点的X轴坐标   LineY(1,:):线条起始节点的Y轴坐标
% LineX(2,:):线条终止节点的X轴坐标   LineY(2,:):线条终止节点的Y轴坐标
% 该函数暂时没有考虑三维坐标，三维坐标时下面的polyxpoly需要替换掉
%===========================================================================%
function [xNode,yNode] = AdjacencyLine(LineX,LineY)

count = length(LineX(1,:));

num=0;
for i=1:count
    ConnectFlag = 0;
    for j=1:count
        if(j~=i)
            [xi,yi]=polyxpoly(LineX(:,i),LineY(:,i),LineX(:,j),LineY(:,j));%线i与线j是否相交
            if(~isempty(xi)) % 非空集:表示两线相交
                if(((xi==LineX(1,i))&&(yi==LineY(1,i)))||((xi==LineX(2,i))&&(yi==LineY(2,i)))||((xi==LineX(1,j))&&(yi==LineY(1,j)))||((xi==LineX(2,j))&&(yi==LineY(2,j))))
                    ConnectFlag = 0;
                else
                    if(sqrt((LineX(1,i)-LineX(2,i))^2+(LineY(1,i)-LineY(2,i))^2) >= sqrt((LineX(1,j)-LineX(2,j))^2+(LineY(1,j)-LineY(2,j))^2))
                        ConnectFlag = 1;
                        break;
                    end
                end
            end
        end
    end
    
    if(ConnectFlag==0)
        num = num + 1;
        xNode(1,num)=LineX(1,i);
        xNode(2,num)=LineX(2,i);
        yNode(1,num)=LineY(1,i);
        yNode(2,num)=LineY(2,i);
    end
end
