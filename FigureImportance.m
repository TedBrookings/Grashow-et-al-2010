function FigureImportance(importance,...
                           prop_names,...
                           param_names,...
                           r_max,...
                           importance_se)

%FigureImportance   Makes a figure plotting values as the size of
%                    circles in a grid.
%   --- Mangled version of Adam Taylor's original code ---
%       Only change:  doesn't make a new figure
%   figure_importance(importance,prop_names,param_names) makes a
%   figure plotting each element of importance as a circle, with the
%   area proportional to the value of that element.  Properties are
%   shown as rows, and parameters as columns.  importance should be 
%   n_properties x n_parameters.  prop_names should be a cell array 
%   of length n_properties, and param_names a cell array of length
%   n_parameters.
%
%   figure_importance(importance,prop_names,param_names,r_max) 
%   makes the circles such that an element with importance equal to
%   one will have a circle of radius r_max
%   
%   figure_importance(importance,prop_names,param_names,r_max,...
%                     importance_se) 
%   plots the standard error of the importance estimtes as lighter-
%   colored rings around the circles
%
%   Copyright 2009 Adam L. Taylor
                         
% args
if nargin<4 || isempty(r_max)
  r_max=0.5;  
end
if nargin<5 
  importance_se=[];
end
                      
% get dims
[n_props,n_params]=size(importance);
                      
% declare some colors
param_clrs_red=[ [1 0 0] ; ...
                 [150 0 200]/255 ; ...
                 [0 0 1] ; ...
                 [0 200 255]/255 ; ...
                 [0 200 0]/255 ; ...
                 [200 200 0]/255 ; ...
                 [255 100 0]/255 ; ...
                 [0 0 0] ];
param_clrs_red=repmat(param_clrs_red,[3 1]);              

% make the figures
%figure;
% set_figure_size([8 4]);
%This was uncommented in Adam's original code, but it causes
%   problems with docked figures:
%set_figure_size([10 5]);
axes;
set_axes_size([7.5 3.5]);
for i=1:n_props
  y=i;
  prop_name=prop_names{i};
  % draw the background bars
  if mod(i,2)==1
    patch([n_params+0.5 n_params+0.5 0.5 0.5],...
          [y+0.5 y-0.5 y-0.5 y+0.5],...
          0.9*[1 1 1],...
          'edgecolor','none');
  end
  for j=1:n_params
    x=j;
    f_this=sqrt(importance(i,j));
    if f_this>0
       if ~isempty(importance_se)
         f_this_se=sqrt(importance(i,j)+2*importance_se(i,j));
         square(j,y,...
                r_max*f_this_se,...
                0.5*[1 1 1]+0.5*param_clrs_red(j,:));
       end
      circle(x,y,...
             r_max*f_this,...
             param_clrs_red(j,:));
    end
  end
end
box on;
h_axis = gca;
set(h_axis,'layer','top');
xlim([0.5 n_params+0.5]);
ylim([0.5 n_props+0.5]);
set(h_axis,'xtick',1:n_params);
set(h_axis,'xticklabel',param_names);
%set(gca,'xticklabel',{});
set(h_axis,'ytick',1:n_props);
set(h_axis,'yticklabel',prop_names);
%set(gca,'yticklabel',{});
set(h_axis,'ydir','reverse');
set(h_axis,'ticklength',[0 0]);
set(h_axis,'dataaspectratio',[1 1 1]);
set(h_axis,'FontName', 'Arial')
set(h_axis,'FontSize', 12)
