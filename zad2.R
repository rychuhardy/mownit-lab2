lu__ <- function(A) # Returns list with L and U
{ 
  if(ncol(A)==nrow(A)) { # Validate argument
    N <- nrow(A)
    
    L <- matrix(rep(0,N*N), nrow=N, ncol=N)
    U <- matrix(rep(0,N*N), nrow=N, ncol=N)
    
    #Factorisation
    for(i in 1:N) {
      for(j in i:N) {
        if(i > 0) {
        U[i,j] = A[i,j] - sum(L[i,1:(i-1)]*U[1:(i-1),j])
        }
        else {
          U[i,j] = A[i,j] - sum(L[i,1:1]*U[1:1,j])
        }
      }
      if(i<N) { # Because for i==N for loop would generate decreasing sequence
        for(j in (i+1):N) {
          if(i > 0) {
          L[j,i] = (A[j,i]-sum(L[j,1:(i-1)]*U[1:(i-1),i]))/U[i,i]
          }
          else {
            L[j,i] = (A[j,i]-sum(L[j,1:1]*U[1:1,i]))/U[i,i]
          }
        }
      }
    }
  }
  return (list(L,U))
}