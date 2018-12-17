%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This functions takes in an array of events and display them. Color
%%% coding depends on user's option. The function doesn't erase or reset
%%% the graph.
%%%=============================================================================
%%% events      |   The trajectory wished to be displayed
%%% varagin     |   Optional variables. The first one is color option.
%%%             |   'energy' or a 1x3 array of hsv values.
%%%             |   2-nd/3-rd ones are acid/excitation positions
%%%             |   4-th one is the max energy
%%%             |   5-th colormap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function traceHandle = displayTrajectory(events,varargin)
    tic
    %%% Taking care of the custom variables
    nArgCustom  =   nargin - 1;
    
    %%% Color option
    if nArgCustom >= 1
        colorOption     =   varargin{1};
    else
        colorOption     =   'energy';   
    end
    
    if nArgCustom >= 3
        acid_xyz        =   varargin{2};
        acid_act_e_xyz  =   varargin{3};
    else
        acid_xyz        =   [];
        axid_act_e_xyz  =   [];
    end
    
    if nArgCustom >=4
        eMax = varargin{4};
    else 
        eMax = events{1}.Ein;
    end
    
    if nArgCustom >=5
        colorMap = varargin{5};
    else 
        colorMap = jet(256);
    end
          
    
    %% Plotting the quivers       
    tic
    for eventIterator = 1:length(events)
        this = events{eventIterator};
        qX = this.xyz_init(1);
        qY = this.xyz_init(2);
        qZ = this.xyz_init(3);
        qU = this.xyz(1);
        qV = this.xyz(2);
        qW = this.xyz(3);
        
        traceHandle = quiver3(qX,qY,qZ,qU-qX,qV-qY,qW-qZ,0,'-','LineWidth',2,...
        'Color',colorOptionWrapper(colorOption,this,eMax,colorMap));
        hold on;
    end  
    
    if length(acid_xyz) >=1
        quiver3(acid_act_e_xyz(:,1),acid_act_e_xyz(:,2),acid_act_e_xyz(:,3),...
            acid_xyz(:,1) - acid_act_e_xyz(:,1),...
            acid_xyz(:,2) - acid_act_e_xyz(:,2),...
            acid_xyz(:,3) - acid_act_e_xyz(:,3),...
            0,'--');
        scatter3(acid_xyz(:,1),acid_xyz(:,2),acid_xyz(:,3),'o')            
    end
    toc
    
end

function color = colorOptionWrapper(colorOption,varargin)
    if nargin >= 3
        eMax = varargin{2};
    else
        eMax = 80;
    end
    
    if length(colorOption) ==   3
        color = hsv2rgb(colorOption);
    elseif strcmp(colorOption,'energy') && nargin >= 2
        color = varargin{3}(ceil(min(varargin{1}.Ein/eMax,1)*size(varargin{3},1)),:);
    elseif strcmp(colorOption,'eLoss') && nargin >= 2
        %color =hsv2rgb([1-varargin{1}.Eloss/eMax,1,1]);
        color = varargin{3}(ceil(min(varargin{1}.Eloss/eMax,1)*size(varargin{3},1)),:);
    else
        color = colorOption;
    end
end