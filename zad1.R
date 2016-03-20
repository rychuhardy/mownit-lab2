gauss_jordan <- function(A, B) { # Regular gauss-jordan elimination
  
  if(nrow(A) == ncol(A) && nrow(A)==nrow(B) && ncol(B)==1) { # Validate input arguments
    A = cbind(A,B) # Join A and B for easier manipulation
    for(k in 1:(nrow(A))) { # Perform elimination
      M = A[k,k]
      for(i in 1:nrow(A)) {
        if( i!=k) {
          N = A[i,k]
          for(j in 1:ncol(A)) {
            A[i,j] = A[i,j] - N / M * A[k,j]
          }
        }
      }
    }  
  }
}

ppivot_gj <- function(A, B) { # Gauss-Jordan elimination with partial pivoting
  if(nrow(A) == ncol(A) && nrow(A)==nrow(B) && ncol(B)==1) { # Validate input arguments
    A = cbind(A,B) # Join A and B for easier manipulation
    for(k in 1:(nrow(A))) { # Perform elimination
      m = which.max(abs(A[k:nrow(A), k])) # Choose the row with largest element in the k-th column
      
      A[c(m,k),] <- A[c(k,m),] # Swap rows so that the chosen maximum is on diagonal
      
      M = A[k,k]
      for(i in 1:nrow(A)) {
        if( i!=k) {
          N = A[i,k]
          for(j in 1:ncol(A)) {
            A[i,j] = A[i,j] - N / M * A[k,j]
          }
        }
      }
      
    }
  }
}

fpivot_gj <- function(A,B) 
{
  C <- c() # Used for remembering swaps between columns. 
  #They are pairs (1,x1), (2,x2) ...  (k, xk) so there is no need to remeber first value
  
  if(nrow(A) == ncol(A) && nrow(A)==nrow(B) && ncol(B)==1) { # Validate input arguments
    A = cbind(A,B) # Join A and B for easier manipulation
    for(k in 1:(nrow(A))) { # Perform elimination
      m = which.max(abs(A[k:nrow(A), k:(ncol(A)-1)])) # Find max elem index (its value is like creating matrix with byrow=false)
      
      m <- m+(k-1)*(nrow(A)+1) # convert m to index of whole matrix A and not its subset from which.max argument
      col <- if(m %/% (ncol(A)-1) == m / (ncol(A)-1)) m / (ncol(A)-1) else m %/% (ncol(A)-1) + 1
      row <- if( m %% (ncol(A)-1) == 0) (ncol(A)-1) else m %% (ncol(A)-1)
      
      #Perform row and col swapping
      A[c(row, k),] <- A[c(k,row),]
      A[,c(col, k)] <- A[,c(k,col)]
      
      C <- append(C, col) # Remember which column is swapped
      
      M = A[k,k]
      for(i in 1:nrow(A)) {
        if( i!=k) {
          N = A[i,k]
          for(j in 1:ncol(A)) {
            A[i,j] = A[i,j] - N / M * A[k,j]
          }
        }
      }
      
      
    }
    # Put columns to their original places (in reversed order)
    for (k in nrow(A):1) {
      A[,c(k, C[k])] <- A[,c(C[k],k)]
    }
    A
  }
}
