fmmnum <- function(s){
  return(3*s-1)
}

## 2s (m/sd), s-1 (transition)

fmmsinnum <- function(s){
  return(5*s - 3)
}

## 2s (m/sd), 3s-3 (transition; ie: )

hmmnum <- function(s){
  return(s^2 + 2*s -1)
}

hmmsinnum <- function(s){
  return(3*s^2 - 1)
}

hmmhourlynum <- function(s){
  return(s-1 + 24*s*(s-1) + 2*s)
}

hmmhourlynum(3:4)
