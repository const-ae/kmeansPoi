## code to prepare `log_lookup_table` dataset goes here


LOG_LOOKUP_VECTOR <- numeric(length = get_LOG_LOOKUP_TABLE_EXPONENTIAL_LENGTH() * get_LOG_LOOKUP_TABLE_SIGNIFCAND_LENGTH())
for(exp_idx in seq(1, get_LOG_LOOKUP_TABLE_EXPONENTIAL_LENGTH())-1L){
  for(sig_idx in seq(1, get_LOG_LOOKUP_TABLE_SIGNIFCAND_LENGTH())-1L){
    eq_double <- fuse_ints_for_double(exp_idx, sig_idx)
    LOG_LOOKUP_VECTOR[(exp_idx * get_LOG_LOOKUP_TABLE_SIGNIFCAND_LENGTH() + sig_idx) + 1] <- log(eq_double)
    stopifnot(get_exponent_as_integer_from_double(eq_double) == exp_idx,
              get_significand_as_integer_from_double(eq_double) == sig_idx)
  }
}


LOG_LOOKUP_VECTOR2 <- numeric(length = get_LOG_LOOKUP_TABLE_EXPONENTIAL_LENGTH() * get_LOG_LOOKUP_TABLE_SIGNIFCAND_LENGTH())
offset <- get_LOG_LOOKUP_TABLE_EXPONENTIAL_LENGTH()/2
for(exp_idx in seq(-get_LOG_LOOKUP_TABLE_EXPONENTIAL_LENGTH()/2, get_LOG_LOOKUP_TABLE_EXPONENTIAL_LENGTH()/2-1)){
  for(sig_idx in seq(1, get_LOG_LOOKUP_TABLE_SIGNIFCAND_LENGTH())-1L){
    eq_double <- fuse_ints_for_double(exp_idx, sig_idx)
    LOG_LOOKUP_VECTOR2[((exp_idx + offset) * get_LOG_LOOKUP_TABLE_SIGNIFCAND_LENGTH() + sig_idx) + 1] <- log(eq_double)
    stopifnot(get_exponent_as_integer_from_double(eq_double) == exp_idx,
              get_significand_as_integer_from_double(eq_double) == sig_idx)
  }
}

usethis::use_data(LOG_LOOKUP_VECTOR, LOG_LOOKUP_VECTOR2, internal = TRUE, overwrite = TRUE)

