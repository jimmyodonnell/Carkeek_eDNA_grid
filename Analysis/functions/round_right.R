# the function 'round' in base R does things that are not expected by a typical person.

round_right <- function(x) trunc(x+0.5)

# example:
some_numbers <- 0:10 + 0.5
round(some_numbers)
round_right(some_numbers)
# which looks more intuitive to you?