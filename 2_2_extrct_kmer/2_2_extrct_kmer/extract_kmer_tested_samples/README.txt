For paf 29
valid_start_idx <- (start_pos >= 20 & start_pos <= 2666)
valid_end_idx <- (end_pos >= 20 & end_pos <= 2666)

For other sample 
valid_start_idx <- (start_pos >= 20 & start_pos <= 2666 & start_pos != 234)
valid_end_idx <- (end_pos >= 20 & end_pos <= 2666 & end_pos != 237)