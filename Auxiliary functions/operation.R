operation <- function(op, A, nodes) {
  x <- nodes[1]
  y <- nodes[2]

  if(op == 1) {
    A[x,y] = 1
    return(A)
  }

  if(op == 2) {
    A[x,y] = 0
    return(A)
  }

  if(op == 3) {
    A[x,y] = 0
    A[y,x] = 1
    return(A)
  }
}
