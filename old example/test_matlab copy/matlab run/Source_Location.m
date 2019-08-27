function  [Source_node] = Source_Location( x, z, Source_x, Source_z, Node_num, ave_size)
    flag = 0;
    for i = 1 : Node_num
        if (abs(x(i) - Source_x) <= ave_size  && abs (z(i) - Source_z) <= ave_size )
            Source_node = i;
            flag = flag + 1;
        end 
    end 
    if (flag == 1 )
       % fprintf('Find the Seismic Source Location at Node: %d \n', Source_node );
    else if ( flag == 0 )
        fprintf('Warning: The source location cannot be found and Check the Source_x and Source_z\n');
        pause
    else
        %fprintf( 'Warning: Find More Than One Sources! \n');
    end
    end
    %fprintf( 'Find the Seismic Source Location at (x, z): (%f, %f) \n', x(Source_node), z(Source_node) );

end 

 