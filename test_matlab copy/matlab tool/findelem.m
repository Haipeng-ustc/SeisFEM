function findelem(node,elem,range,varargin)
%% FINDELEM highlights some elements
%
%    FINDELEM(node,elem,range) finds all elements whose indices are in the
%    range array by displaying these elements in yellow.
%
%    FINDELEM(node,elem) finds all elements.
%   
%    FINDELEMnode,elem,range,'noindex') skip the display of indices.
%
%    FINDELEM(node,elem,range,'param','value','param','value'...) allows
%    additional patch param/value pairs to highlight the elements.
%    
% Example:
%     node = [1,0; 1,1; 0,1; -1,1; -1,0; -1,-1; 0,-1; 0,0];
%     elem = [1,2,8; 3,8,2; 8,3,5; 4,5,3; 7,8,6; 5,6,8];
%     subplot(1,2,1);
%     showmesh(node,elem);
%     findelem(node,elem,1,'index','FaceColor','r','MarkerSize',24);
%     subplot(1,2,2);
%     showmesh(node,elem);
%     findelem(node,elem);
%
%   See also findelem3, findnode3, findedge.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

hold on

if (nargin==2) || isempty(range) || (ischar(range) && strcmp(range,'all'))
    range = (1:size(elem,1))'; 
end
if islogical(range)
    range = find(range); 
end
if size(range,2)>size(range,1)
    range = range'; 
end
center=(node(elem(range,1),:)+node(elem(range,2),:)+node(elem(range,3),:))/3;
if length(range) < size(elem,1)
    x = reshape(node(elem(range,:),1),size(range,1), size(elem,2))';
    y = reshape(node(elem(range,:),2),size(range,1), size(elem,2))';
    h = patch(x,y,'y');
    if nargin > 3
        if strcmp(varargin{1},'noindex') || strcmp(varargin{1},'index')
           if size(varargin,2)>=2
                set(h,varargin{2:end});
           end
        else
            set(h,varargin{:});
        end
    end
end
if (nargin <=3) || ~(strcmp(varargin{1},'noindex'))
    if size(node,2) == 2
        plot(center(:,1),center(:,2),'o','LineWidth',1,'MarkerEdgeColor','k',...
         'MarkerFaceColor','y','MarkerSize',18);    
        text(center(:,1)-0.02,center(:,2),int2str(range),'FontSize',12,...
        'FontWeight','bold','Color','k');
    elseif size(node,2) == 3 % surface mesh
        plot3(center(:,1),center(:,2),center(:,3),'o','LineWidth',1,'MarkerEdgeColor','k',...
         'MarkerFaceColor','y','MarkerSize',18);    
        text(center(:,1)-0.02,center(:,2),center(:,3),int2str(range),'FontSize',12,...
        'FontWeight','bold','Color','k');        
    end
end
hold off