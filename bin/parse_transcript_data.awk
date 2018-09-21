BEGIN{start=0;};$1 ~ /series_matrix_table_end/{start=0;}; { if(start==1){print;} }; $1 ~ /series_matrix_table_begin/{start=1;};
