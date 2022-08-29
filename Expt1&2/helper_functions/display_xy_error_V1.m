function display_xy_error_V1(x,y, xe, ye, bar_col)


for i=1:length(x)
   
    plot( [x(i),x(i)], [y(i)-ye, y(i)+ye], 'color', bar_col);
    plot( [x(i)-xe,x(i)+xe], [y(i), y(i)], 'color', bar_col);
     
end


end