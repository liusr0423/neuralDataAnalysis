function figscroll(src,event)
%function enables us to scroll left and right using the arrow keys
ax = get(src, 'CurrentAxes');
x_orig = get(ax, 'XLim');
x_step = (x_orig(1)-x_orig(2))/2;
 
switch event.Key 
    case 'leftarrow'
        x_new = x_orig + x_step;
    case 'rightarrow'
        x_new = x_orig - x_step;
end
 
set(ax,'XLim',x_new);