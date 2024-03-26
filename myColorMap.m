function Cmap = myColorMap(Colour_List,Num_Steps)
%% myColorMap creats smooth costmised unlimited colours colormaps
% Colour_List is the colour matrix
% Num_Steps is a +ve scalar
% Copyright 2016 Dr Ahmed Abass a.abass@liverpool.ac.uk 
% -------------------------------------------------------------------------
% Example
% Colour_01 = [0 0 0];% black
% Colour_02 = [0 0 1];% blue
% Colour_03 = [0 1 0];% green
% Colour_04 = [1 1 0];% yellow
% Colour_05 = [1 0 0];% red
% Colour_06 = [1 1 1];% white
% Colour_List = cat(1,Colour_01,Colour_02,Colour_03,Colour_04,Colour_05,Colour_06); % List of colours in one matrix
% Number_of_Steps = 50;
% my_Cmap = myColorMap(Colour_List,Number_of_Steps);
% figure
% surf(peaks(100),'edgecolor','none')
% shading interp
% colormap(my_Cmap);
% colorbar
% grid off
% set(gca,'Color','w')
% set(gcf,'Color','w')
% -------------------------------------------------------------------------
if size(Colour_List,2) ~= 3
    error(ErrorMessage)
elseif (max(Colour_List(:)) >1) || (min(Colour_List(:)) <0)
    error('Colors must have values in [0,1].')
end
if nargin == 1
    Num_Steps = 100;
elseif Num_Steps < size(Colour_List,1)
    Num_Steps = size(Colour_List,1);
else
    Num_Steps = round(Num_Steps);
end
if nargin<= 2
    X0 = transpose(linspace(0,1,size(Colour_List,1)));
    Y1 = Colour_List(:,1);
    Y2 = Colour_List(:,2);
    Y3 = Colour_List(:,3);
    X1 = transpose(linspace(0,1,Num_Steps));
    Y1 = interp1(X0,Y1,X1,'PCHIP');
    Y2 = interp1(X0,Y2,X1,'PCHIP');
    Y3 = interp1(X0,Y3,X1,'PCHIP');
    Cmap = [Y1 Y2 Y3];
else
    error('correct input format. Type ''help myColorMap'' for help.')
end
% Copyright 2016 Dr Ahmed Abass a.abass@liverpool.ac.uk 