function [ data_after_dump, dump_position ] = selectdata(data,  percent_element)
% data: our object
% percent_element: the percent for the data we want to dump
num_elements = sum(sum( data ~= 0));
dump_elements = floor( num_elements * percent_element);

pos_of_elements = find( data ~= 0);
dump = datasample( pos_of_elements, dump_elements, 'Replace', false );
temp = data;
temp(dump) = 0;
data_after_dump = temp;
dump_position = dump;
end
