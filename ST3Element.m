classdef ST3Element<handle
    % 常应变单元的三角形几何条件
    properties (SetAccess=private,GetAccess=public)
        Tri=[];
        TriEdgeBoundary=[];
        %%
        Points=[];
        Node_Boundary=[];
        NodeAttachments=[];
        %%
        Edges=[];%
        EdgeAttachments=[];
        EdgeBoundary=[]; % 1 is boundary edge; 0 is inner edge.
        %%
        % 边中点的编号ID,对应边的序号，按6节点单元的排序方式放置。
    end
    properties(Dependent)
        Triangulation;
        Tri_Area;
        BoundaryEdge;
    end
    methods
        function obj=ST3Element(T,P)
            %T3_LINEARELEMENTGEOMETRY 构造此类的实例
            %   此处显示详细说明
            if nargin==0
                obj.Tri=[];
            else
                obj.Points=P;
                obj.Tri=T;
                obj.ConstructGeometry(T,P);
            end
        end
        function ConstructGeometry(obj,T,P)
            %% 返回三角形单元的必要拓扑关系
            %
            %%
            trinum=size(T,1);
            TriID=1:trinum;
            TriID=repmat(TriID,3,1); % copy 3 rows
            TriID=TriID(:);
            %% Get Edges
            Edge=[T(:,1),T(:,2),T(:,2),T(:,3),T(:,3),T(:,1)]';
            Edge=Edge(:);
            Edge=reshape(Edge,2,[]);
            Edge=Edge';
            Temp1=Edge;
            Edge=sort(Edge,2); % sort node order for each edge
            [edges,ia,ic]=unique(Edge,'rows'); % the repeat edges
            obj.Edges=Temp1(ia,:); % 保证边界边是逆时针排列
            %% Get Edges Attached triangle elements
            Temp2=accumarray(ic,TriID,[],@(x) {x});
            obj.EdgeAttachments=Temp2;
            %% Get boundary information of edges
            A=cellfun(@(x) 2-length(x),Temp2); % A(i)=2-length(Temp{i});Note: length(Temp{i})=1 or length(Temp{i})=2;% A(i)=1: length(x)=1 boudary edge.;% A(i)=0: length(x)=2 inner edge.
            obj.EdgeBoundary=A; % 1 boundary 0 inner
            %% Get boundary information of elements' edges
            Temp3=A(ic);             % index of edge in triangle elements
            Temp4=reshape(Temp3,3,[]);
            obj.TriEdgeBoundary=Temp4';  % 1 boundary; 0 inner edge
            %% Get Edge and Edge Attached triangle elements
            BEdge=edges(A==1,:);
            obj.Node_Boundary=unique(BEdge(:));
            %% GetNode Attached triangle elements
            Temp5=T';
            NodeID=Temp5(:);
            [~,~,ic]=unique(NodeID);
            obj.NodeAttachments=accumarray(ic,TriID,[],@(x) {x});  %  NodeAttachments{i}=TriID(ic==i);;
        end
        function Triangulation=get.Triangulation(obj)
            T=obj.Tri;
            P=obj.Points;
            Triangulation=triangulation(T,P);
        end
        function Tri_Area=get.Tri_Area(obj)
            T=obj.Tri;
            P=obj.Points;
            x=P(:,1);
            y=P(:,2);
            Tri_x=x(T);
            Tri_y=y(T);
            Tri_Area=polyarea(Tri_x,Tri_y,2);% 可GPU
        end
        function BoundaryEdge=get.BoundaryEdge(obj)
            TR=obj.Triangulation;
            BoundaryEdge=freeBoundary(TR);
        end
        function PlotElements(obj)
            T=obj.Tri;
            P=obj.Points;
            TR=triangulation(T,P);
            triplot(TR,'k');
            axis equal;
            hold off
        end
        function ChangePoints(obj,P)
            obj.Points=P;
        end
    end
end

