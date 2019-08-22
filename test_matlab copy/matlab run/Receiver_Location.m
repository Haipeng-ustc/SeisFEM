function  [Rec_node] = Receiver_Location( x, z, Rec_x, Rec_z, Node_num, ave_size)
    flag = 0;
    for i = 1 : Node_num
        if (abs(x(i) - Rec_x) <= ave_size  && abs (z(i) - Rec_z) <= ave_size )
            Rec_node = i;
            flag = flag + 1;
        end 
    end 
    if (flag == 1 )
       % fprintf('Find the receiver Location at Node: %d \n', Rec_node );
    else if ( flag == 0 )
        fprintf('Warning: The receiver location cannot be found and Check the Rec_x and Rec_z!\n');
        pause
    else
      %  fprintf( 'Warning: Find more than one receivers! \n');
    end
    end
    %fprintf( 'Find the receiver location at (x, z): (%f, %f) \n', x(Rec_node), z(Rec_node) );

end 

 