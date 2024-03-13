


Uhex <- matrix(0,6,2)
for(i in 1:6){
  Uhex[i,1] <- cos(i*60*pi/180) 
  Uhex[i,2] <- sin(i*60*pi/180)
}
Uhex <- zapsmall(Uhex)


Uras <- matrix(0,4,2)
for(i in 1:4){
  Uras[i,1] <- sqrt(2)*cos(i*90*pi/180) 
  Uras[i,2] <- sqrt(2)sin(i*90*pi/180)
}
Uras <- zapsmall(Uras)


Uq <- matrix(0,8,2)
for(i in 1:8){
  Uq[i,1] <- cos(i*45*pi/180) 
  Uq[i,2] <- sin(i*45*pi/180)
}
Uq <- zapsmall(Uq)


m <- rnorm(2, 0, 3)
# m <- c(1,0)
# U <- Uras
U <- Uq
Um <- U %*% m
# Um <- pmax(0,Um)
# U <- U/1.5

xxx <- matrix(0, 1, 2)
for(i in 1:nrow(U)){
  # print(U[i,] %*% t(U[i,]))
  xxx <- xxx + Um[i] * t(U[i,])
}
zapsmall(xxx/4) - m


xxx <- matrix(0, 2, 2)
for(i in c(1:8)){
  print( outer(U[i,], U[i,]))
  xxx <- xxx + outer(U[i,] , U[i,])
}
xxx




   