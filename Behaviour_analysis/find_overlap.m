function overlap = find_overlap( int_series1, int_series2)


    nx = size(int_series1,1);
    ny = size(int_series2,1);
    
    if nx < ny
        int_series1 = cat(1, int_series1, nan(ny-nx,2));
    else
        int_series2 = cat(1, int_series2, nan(nx-ny,2));
    end
    
    xy = cat(3, int_series1, int_series2 );
    
    int_start = 1; int_end = 2;
    i1 = 1; i2 = 1;
    indx(1) = 1; indx(2) = 1;
    if xy( indx(1), int_start, 1) < xy( indx(2), int_start, 2)
        early_int = 1;
        late_int = 2;
    else
        early_int = 2;
        late_int = 1;
    end
    
    overlap = [];
    while indx(1) <= nx && indx(2) <= ny
       if xy( indx(early_int), int_start, early_int ) > xy( indx(late_int), int_start, late_int) 
          early_int = late_int;
          late_int = 1 + 1*(early_int == 1 );
       end
       if xy( indx( late_int), int_start, late_int ) <=   xy( indx( early_int), int_end, early_int )
           % overlap exists because x_start <= y_start <= x_end
           if xy( indx(late_int), int_end, late_int) <= xy( indx(early_int), int_end, early_int)
              overlap = cat(1, overlap, reshape( xy( indx(late_int), :, late_int ), [1,2]) );
              indx(late_int) = indx(late_int) +1;
           else
              overlap = cat(1, overlap, [xy( indx(late_int), int_start, late_int ), xy( indx(early_int), int_end, early_int )] );
              indx(early_int) = indx(early_int) + 1;
           end
       else
          indx(early_int) = indx(early_int) + 1;
       end
    end
end