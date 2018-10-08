function cbar = colorbar_my(ax)
% COLORBAR MY - add colorbar, keep axes in same place
% h = colobar_my
%
% label it too: set(get(h,'Ylabel'),'string','Log_1_0(K_z)')



if nargin < 1, ax = gca; end
pos = get(ax,'position');
cbar = colorbar('peer',ax);
set(ax,'position',pos)
axes(ax)

% move colorbar to just outside of right axis
cb_pos = get(cbar,'position');
set(cbar,'position',[pos(1)+pos(3) cb_pos(2:end)])

end