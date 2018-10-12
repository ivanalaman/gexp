latin <- function(n, nrand = 20) {
             x = matrix(LETTERS[1:n], n, n)
             x = t(x)
             for (i in 2:n) x[i, ] = x[i, c(i:n, 1:(i - 1))]
             if (nrand > 0) {
               for (i in 1:nrand) {
                 x = x[sample(n), ]
                 x = x[, sample(n)]
               }
             }
             x
           }    
