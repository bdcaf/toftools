library(dplyr)
library(purrr)

canMerge <- function(sp1, sp2, max_dif = 10) { (sp2[[1,'starts']] - sp1[[1,'ends']]) < max_dif }
mergeSp <- function(sp1, sp2){
  sp1$ends <- sp2$ends
  sp1
}

mergeFun <- function(acc, cur, tib, max_dif){
  if (canMerge(acc, cur, max_dif)){
        tibout <- tib
        acc2 <- mergeSp(acc, cur)
  } else {
        tibout <- rbind(tib, acc)
        acc2 <- cur
  }
  list(acc = acc2, tib = tibout)
}

TibAccRed <- function(f, tibin, max_dif=10){
  init <- data_frame()
  acc <- tibin[1,]
  tibin <- tibin[-1,]

  while (nrow(tibin) > 0){
        cur <- tibin[1,]
        tibin <- tibin[-1,]

        val <- f(acc, cur, init, max_dif)
        acc <- val$acc
        init <- val$tib
  }

  rbind(init, acc)
}

#dd <- TibAccRed(mergeFun, only_sig)


