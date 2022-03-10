set.seed(123456)
gs=list(s1=list(g1=sample(1000,abs(floor(100*rnorm(1)))),
                  g2=sample(1000,abs(floor(100*rnorm(1))))),
        s2=list(g1=sample(1000,abs(floor(100*rnorm(1)))),
                  g2=sample(1000,abs(floor(100*rnorm(1))))),
        s3=list(g1=sample(1000,abs(floor(100*rnorm(1)))),
                  g2=sample(1000,abs(floor(100*rnorm(1))))))

df <- data.frame()

for (i in seq_along(gs)) {
  l <- gs[[i]]
  for (j in seq_along(l)) {
    ll <- l[[j]]
    len <- length(ll)
    dfl <-
      cbind(
        sample = rep(names(gs)[i], len),
        gene = rep(names(l)[j], len),
        data = ll
      )
    df <- rbind(df, dfl)
  }
}

# OR:
deg <- gs
deg_list <- lapply(names(deg), function(y) {
  tmp <- deg[[y]]
  data.frame(group = paste(y, unlist(lapply(names(tmp), function(x) {
    rep(x, length(tmp[[x]]))
  })), sep = '_') ,
  gene = unlist(tmp))
})
group_g<- do.call(rbind, deg_list)
