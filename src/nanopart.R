library(mvtnorm)

# Read one hourly file
get.raw <- function(file.name) {
  data <- readBin(file.name, "raw", 10e6)
  df <- as.data.frame(t(matrix(as.numeric(data), nrow=10)))
  two.compl <- function(x, y) { (x + 2^8*y) -> z; ifelse(z < 2^15, z, z - 2^16) }
  out <- data.frame(
    second = two.compl(df$V1, df$V2),
    u      = two.compl(df$V3, df$V4) / 100.0,
    v      = two.compl(df$V5, df$V6) / 100.0,
    w      = two.compl(df$V7, df$V8) / 100.0,
    t      = two.compl(df$V9, df$V10) / 100.0
  )
  out <- out[out$second < 5000,]
  out$vel <- sqrt(out$u^2 + out$v^2)
  out$dir <- 180/pi * atan2(-out$u, -out$v)
  out$dir[out$dir < 0] <- out$dir[out$dir < 0] + 360.0
  n <- length(out$second)
  out$second <- 3600*(0:(n-1))/n
  return(out)
}

# Reattribute time stamp, id number of data is fair
check.file <- function(d, freq=10, threshold.percent=0.1) {
  n <- length(d$second)
  ok <- (abs(n-freq*3600)/(freq*36) < threshold.percent)
  return(ok)
}


met.proc <- function(file.name, time.step=1, time.width=1) {

  d <- get.raw(file.name)
  if(!check.file(d)) {
    print("File is not 0.1% within expected 36000 records length")
    return(NULL)
  }

  time.instants <- seq(from=0, to=3600, by=time.step)
  n <- length(time.instants)
  u <- numeric(n)
  v <- numeric(n)
  uu <- numeric(n)
  vv <- numeric(n)
  uv <- numeric(n)
  for(i in 1:n) {
    t.i <- time.instants[i]
    print(t.i)
    time.idx <- which((t.i-time.width/2. <= d$second) & (d$second <= t.i+time.width/2.))
    k <- length(time.idx)
    if(k > 1) {
      u[i]  <- mean(d$u[time.idx])
      v[i]  <- mean(d$v[time.idx])
      uu[i] <- var(d$u[time.idx])
      vv[i] <- var(d$v[time.idx])
      uv[i] <- cov(d$u[time.idx], d$v[time.idx])
    }
    else {
      u[i]  <- NA
      v[i]  <- NA
      uu[i] <- NA
      vv[i] <- NA
      uv[i] <- NA
    }
  }

  met.set <- data.frame(
    tm = time.instants,
    u.m  = u,
    v.m  = v,
    u.v  = uu,
    v.v  = vv,
    uv.c = uv
  )

  return(met.set)

}


simulate.conc <- function(met.set, max.time, x.rec, y.rec, rec.width=2, max.seconds=7200, parts.per.step=10, debug=FALSE) {

  # Reserve particle space
  n <- min(c(length(met.set$tm), max.time))
  conc <- numeric(n)
  print(n)
  max.points <- max.seconds * parts.per.step
  x <- numeric(max.points)
  y <- numeric(max.points)

  # Iterate over meteo steps
  part.next <- 1
  part.count <- 0
  for(i in 1:n) {

    # Release new particles
    x[part.next:(part.next+parts.per.step)] <- 0
    y[part.next:(part.next+parts.per.step)] <- 0
    part.next <- part.next + parts.per.step
    if(part.next > max.points) {
      part.next <- 1
    }
    part.count <- part.count + parts.per.step
    if(part.count > max.points) {
      part.count <- max.points
    }

    # Get meteo data, and use them to compose the random velocities
    vel   <- c(met.set$u.m[i], met.set$v.m[i])
    sigma <- cbind( c(met.set$u.v[i], met.set$uv.c[i]), c(met.set$uv.c[i], met.set$v.v[i]))

    # Apply shift (note: 1s period!)
    if(debug) {print(i)}
    shifts <- rmvnorm(part.count, mean=vel, sigma=sigma)
    part.range <- 1:part.count
    x[part.range] <- x[part.range] + shifts[,1]
    y[part.range] <- y[part.range] + shifts[,2]

    # Count particles in range
    conc[i] <- sum(sqrt((x-x.rec)^2+(y-y.rec)^2) < rec.width)

  }

  out <- list(
    x = x,
    y = y,
    conc = conc
  )
  return(out)

}


summarize.pos <- function(model.out) {
  x.m <- mean(model.out$x)
  y.m <- mean(model.out$y)
  dev <- mean(sqrt((model.out$x - x.m)^2 + (model.out$y - y.m)^2))
  ang <- atan2(x.m, y.m)
  out <- list(
    x.m = x.m,
    y.m = y.m,
    dev = dev,
    ang = ang
  )
  return(out)
}


summarize.met <- function(met.set, rec.dist=100) {
  u   <- mean(met.set$u.m)
  v   <- mean(met.set$v.m)
  rad <- atan2(u, v)
  deg <- 180*rad/pi
  rec.x <- rec.dist*sin(rad)
  rec.y <- rec.dist*cos(rad)
  out <- list(rad=rad, deg=deg, x=rec.x, y=rec.y)
  return(out)
}


spatial.plot <- function(pos, file.name, edge.len=500, edge.center.x=250, edge.center.y=0) {
  edge.int.x <- c(edge.center.x-edge.len/2, edge.center.x+edge.len/2)
  edge.int.y <- c(edge.center.y-edge.len/2, edge.center.y+edge.len/2)
  png(file=file.name, width=8, height=8, units="in", res=96*4)
  plot(pos$x, pos$y, cex=0.1, xlim=edge.int.x, ylim=edge.int.y, xlab="", ylab="")
  dev.off()
}


time.plot <- function(cnt, file.name) {
  png(file=file.name, width=8, height=6, units="in", res=96*4)
  plot(cnt$conc, type="l", xlab="Time (s)", ylab="Concentration (nominal)")
  dev.off()
}


