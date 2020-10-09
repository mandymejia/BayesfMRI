# Internal function that fins connected components on a mesh. It  takes a vector of indices for
# nodes in a mesh and returns a list of the connected components in that vector
find.connected <- function(ind,mesh){
  if(class(mesh) != "inla.mesh"){
    stop("mesh should be of class inla.mesh")
  }

  sets <- list()
  ind.add <- ind
  k=1
  while(length(ind.add)>0){
    set.k <- c(ind.add[1])
    ind.add <- ind.add[-1]
    #check which other nodes that are connected to the first
    tri.k <- which(rowSums(mesh$graph$tv==set.k[1])>0)
    rem.con <- which(ind.add %in% c(mesh$graph$tv[tri.k,]))
    while(length(rem.con)>0){
      #add connected nodes
      set.k <- c(set.k,ind.add[rem.con])
      ind.add <- ind.add[-rem.con]
      #check if further nodes are connected to the added
      tc <- 0*mesh$graph$tv
      for(i in 1:length(set.k)){
        tc <- tc+(mesh$graph$tv==set.k[i])
      }
      tri.k <- which(rowSums(tc)>0)
      rem.con <- which(ind.add %in% c(mesh$graph$tv[tri.k,]))
    }
    sets[[k]] <- set.k
    k=k+1
  }
  return(sets)
}

# Internal function that find the smalles excursion set. If res is is of class excurobj, the
# continuous function is not used to calculate the area of the connected components. Instead the
# function uses the approximation that the contribution from each element
# \phi_i is \int \phi_i(s)ds = C_ii.
# The argument area.el contains these integrals. The argument is thus only used if res is of class
# excurobj.  area.limit is a limit for the area of the sets. NULL is returned if all sets have areas
# larger than area.limit. If factor>0 the set is expanded by an amount given by factor*area.limit.

#' Find smallest activation
#' 
#' @param res Undocumented
#' @param mesh Undocumented
#' @param area.limit Undocumented
#' @param factor Undocumented
#' @param area.el Undocumented
#' 
#' @import sp
#' 
find.smallest.activation <- function(res,mesh,area.limit,factor,area.el){

  # Check to see that the INLA package is installed
  if (!requireNamespace("grDevices", quietly = TRUE)) {
    stop("`find.smallest.activation` requires the `grDevices` package. Please install it.", call. = FALSE)
  }

  if (class(res) == "excurobj") {

    ind.E <- which(res$E == 1)
    ind.rem <- NULL
    areas <- NULL
    if(length(ind.E)>0){
      #extract areas of connected components
      s <- find.connected(ind.E,mesh)
      areas <- rep(0,length(s))
      for (i in 1:length(s)) {
        areas[i] <- sum(area.el[s[[i]]])
      }

      #find components to remove

      if (sum(areas < area.limit) > 0) {
        ind.rem <- s[[which.min(areas)]]
      }
    }

  } else {
    #extract areas of connected components
    areas <- unlist(lapply(res$M["1"]@polygons[[1]]@Polygons,function(x) x@area))

    #find components to remove
    ind.rem <- NULL
    if (sum(areas < area.limit) > 0) {
      ind <- which.min(areas)
      set.rem <- res$M["1"]
      set.rem@polygons[[1]]@Polygons <- set.rem@polygons[[1]]@Polygons[ind]
      mesh.locs <- data.frame(x = mesh$loc[,1],y = mesh$loc[,2])
      sp::coordinates(mesh.locs) <- ~ x + y
      ov <- sp::over(mesh.locs,set.rem)

      center.points <- mesh.locs[!is.na(ov)]
      p.list <- NULL
      for (i in 1:length(center.points)) {
        p1 <- center.points[i,]
        p.list <- rbind(p.list,cbind(c(p1@coords[1,1] - sqrt(area.limit)*factor,
                                       p1@coords[1,1],
                                       p1@coords[1,1] + sqrt(area.limit)*factor,
                                       p1@coords[1,1]),
                                     c(p1@coords[1,2],
                                       p1@coords[1,2] - sqrt(area.limit)*factor,
                                       p1@coords[1,2],
                                       p1@coords[1,2] + sqrt(area.limit)*factor)))

      }
      ch <- grDevices::chull(p.list)
      p.coords <- p.list[c(ch, ch[1]), ]
      p <- sp::Polygon(p.coords)
      sps = sp::SpatialPolygons(list(sp::Polygons(list(p),1)))
      ov <- sp::over(mesh.locs,sps)
      ind.rem <- unique(which(!is.na(ov)))
    }
  }
  if(is.null(areas)){
    min.area <- NULL
  } else {
    min.area <- min(areas)
  }
  return(list(ind.rem = ind.rem, smallest.area = min.area))
}

#' Calculation of excursion sets with the removal of small areas
#'
#' @description Function similar to the excursions function, but which also takes a mesh as input
#' and computes modified excursion sets where all regions with areas smaller than area.limit are
#' removed in the internal calculations.
#' @param alpha Error probability for the excursion set.
#' @param u Excursion or contour level.
#' @param mu Expectation vector.
#' @param Q Precision matrix.
#' @param type Type of region:
#'  \itemize{
#'     \item{'>' }{positive excursion region}
#'     \item{'<' }{negative excursion region}
#'     \item{'!=' }{contour avoiding region}
#'     \item{'=' }{contour credibility region}}
#' @param n.iter Number or iterations in the MC sampler that is used for approximating probabilities. The default value is 10000.
#' @param Q.chol The Cholesky factor of the precision matrix (optional).
#' @param F.limit The limit value for the computation of the F function. F is set to NA for all nodes where F<1-F.limit. Default is F.limit = \code{alpha}.
#' @param vars Precomputed marginal variances (optional).
#' @param rho Marginal excursion probabilities (optional). For contour regions, provide \eqn{P(X>u)}.
#' @param method Method for handeling the latent Gaussian structure:
#'  \itemize{
#'       \item{'EB' }{Empirical Bayes (default)}
#'       \item{'QC' }{Quantile correction, rho must be provided if QC is used.}}
#' @param ind Indices of the nodes that should be analysed (optional).
#' @param max.size Maximum number of nodes to include in the set of interest (optional).
#' @inheritParams verbose_Param_direct_FALSE
#' @inheritParams max.threads_Param
#' @inheritParams seed_Param
#' @param area.limit Positive number. All connected excursion sets with an area smaller than this
#' numbere are removed.
#' @param use.continuous Logical parameter indicating whether the areas of the excursion sets
#' should be calculated using the \code{continuous} function in \code{excursions}. If FALSE, the
#' function uses the approximation that the area for each node is the integral of the FEM basis function \eqn{\phi_i}.
#' @param factor Non-negative number. If larger than zero, each set is expanded by a small amount
#' that is proportional to this argument.
#' @param plot.progress Logical parameter that indicates whether the results should be plotted.
#' @param mesh The mesh on which the model is defined.
#'
#' @return If \code{use.continuous = FALSE}, an item of class \code{excurobj}. Otherwise a list with
#' the same elements as the output of \code{continuous}.
#' @export
#' @import excursions
#' @importFrom INLA inla.mesh.projector
#'
excursions.no.spurious <- function(alpha,
                                   u,
                                   mu,
                                   Q,
                                   type,
                                   n.iter=10000,
                                   Q.chol,
                                   F.limit,
                                   vars,
                                   rho=NULL,
                                   method='EB',
                                   ind,
                                   max.size,
                                   verbose=FALSE,
                                   max.threads=0,
                                   seed,
                                   area.limit,
                                   use.continuous = FALSE,
                                   factor = 0, #1/2 for removing area.limit around each point
                                   plot.progress = FALSE,
                                   mesh)

{
  if (method == 'QC') {
    qc = TRUE
  } else if (method == 'EB') {
    qc = FALSE
  } else {
    stop('only EB and QC methods are supported.')
  }
  if (missing(alpha))
    stop('Must specify error probability')

  if (missing(mesh))
    stop('Must specify mesh')

  if (missing(area.limit))
    stop('Must specify area limit')

  if (missing(u))
    stop('Must specify level')

  if (missing(mu)) {
    stop('Must specify mean value')
  } else {
    mu <- excursions:::private.as.vector(mu)
  }
  if (missing(Q) && missing(Q.chol))
    stop('Must specify a precision matrix or its Cholesky factor')

  if (missing(type))
    stop('Must specify type of excursion set')

  if (qc && missing(rho))
    stop('rho must be provided if QC is used.')

  if (missing(F.limit)) {
    F.limit = alpha
  } else {
    F.limit = max(alpha,F.limit)
  }

  if (!missing(Q.chol) && !is.null(Q.chol)) {
    ## make the representation unique (i,j,v) and upper triangular
    Q = excursions:::private.as.dgTMatrix(excursions:::private.as.dtCMatrixU(Q.chol))
    is.chol = TRUE
  } else {
    ## make the representation unique (i,j,v)
    Q = excursions:::private.as.dgTMatrix(Q)
    is.chol = FALSE
  }

  if (missing(vars)) {
    if (is.chol) {
      vars <- excursions.variances(L = Q)
    } else {
      vars <- excursions.variances(Q = Q)
    }
  } else {
    vars <- excursions:::private.as.vector(vars)
  }

  if (!missing(rho))
    rho <- excursions:::private.as.vector(rho)

  if (!missing(ind))
    ind <- excursions:::private.as.vector(ind)

  if (verbose)
    cat("Calculate marginals\n")
  if(is.null(rho)){
    marg <- excursions:::excursions.marginals(type = type, vars = vars,
                                              mu = mu, u = u, QC = qc)

  } else {
    marg <- excursions:::excursions.marginals(type = type, rho = rho,vars = vars,
                                              mu = mu, u = u, QC = qc)

  }


  if (missing(max.size)) {
    m.size = length(mu)
  } else {
    m.size = max.size
  }
  if (!missing(ind)) {
    if (is.logical(ind)) {
      indices = ind
      if (missing(max.size)) {
        m.size = sum(ind)
      } else {
        m.size = min(sum(ind),m.size)
      }
    } else {
      indices = rep(FALSE,length(mu))
      indices[ind] = TRUE
      if (missing(max.size)) {
        m.size = length(ind)
      } else {
        m.size = min(length(ind),m.size)
      }
    }
  } else {
    indices = rep(TRUE,length(mu))
  }

  n = length(mu[ind])
  ii = which(rho[1:n] > 0)
  if (length(ii) == 0) i = n + 1 else i = min(ii)

  F = Fe  = E = G = rep(0,n)
  F = marg$rho

  E[F > 1 - alpha] = 1

  if (type == '=') {
    F = 1 - F
  }

  if (type == "<") {
    G[mu > u] = 1
  } else {
    G[mu >= u] = 1
  }

  M = rep(-1,n)
  if (type == "<") {
    M[E == 1] = 0
  } else if (type == ">") {
    M[E == 1] = 1
  } else if (type == "!=" || type == "=") {
    M[E == 1 & mu > u] = 1
    M[E == 1 & mu < u] = 0
  }

  if (missing(ind) || is.null(ind)) {
    ind <- seq_len(n)
  } else if (is.logical(ind)) {
    ind <- which(ind)
  }

  res.marg <- list(F = F,
                 G = G,
                 M = M,
                 E = E,
                 mean = mu,
                 vars = vars,
                 rho = marg$rho,
                 meta = (list(calculation = "excursions",
                              type = type,
                              level = u,
                              F.limit = 1,
                              alpha = alpha,
                              n.iter = n.iter,
                              method = method,
                              ind = ind,
                              Fe = Fe,
                              call = match.call())))
  class(res.marg) <- "excurobj"
  res.marg <- crop.lists(res.marg)
  if (verbose)
    cat("Calculate sizes\n")
  if(use.continuous){
    sets.marg <- excursions::continuous(ex = res.marg, geometry = mesh, alpha = alpha)
    areas <- find.smallest.activation(sets.marg,mesh,area.limit = area.limit,factor = factor)
  } else {
    #area.el <- diag(inla.fmesher.smorg(mesh$loc,mesh$graph$tv, fem = 0, output = list("c0"))$c0)
    area.el <- inla.fmesher.smorg(mesh$loc,mesh$graph$tv, fem = 0, output = list("c0"))$c0@x
    areas <- find.smallest.activation(res.marg,mesh,area.limit = area.limit,factor = factor,
                                      area.el= area.el)
  }

  ind.rem <- areas$ind.rem

  rho <- marg$rho
  rho.crop = rho[ind]
  if(use.continuous == FALSE && plot.progress){
    proj <- inla.mesh.projector(mesh, dims = c(100,100))
  }
  ###########
  if(is.null(ind.rem)){
    cat("Calculated excursion set\n")
    res.exc <-  excursions(mu = mu,
                           Q = Q,
                           alpha = alpha,
                           type = type,
                           u = u,
                           rho = rho,
                           ind = ind,
                           n.iter = n.iter,
                           Q.chol = Q.chol,
                           F.limit = F.limit,
                           method = method,
                           vars = vars,
                           max.size,
                           verbose = FALSE,
                           max.threads = max.threads,
                           seed = seed)
    res.exc <- crop.lists(res.exc)
    if (verbose)
      cat("Calculate sizes\n")
    if(use.continuous){
      sets <- excursions::continuous(ex = res.exc, geometry = mesh, alpha = alpha)
      areas <- find.smallest.activation(sets, mesh, area.limit = area.limit, factor = factor)
    } else {
      areas <- find.smallest.activation(res.exc, mesh, area.limit = area.limit, factor = factor,
                                        area.el= area.el)
    }
    ind.rem <- areas$ind.rem
    cat(sprintf("smallest area = %f (limit = %f)\n", areas$smallest.area, area.limit))
    if (plot.progress) {
      # if(use.continuous){
      #   plot(mesh, vertex.color = rgb(0.8, 0.8, 0.8), draw.segments = FALSE,
      #        edge.color = rgb(0.8, 0.8, 0.8), main = " ")
      #   plot(sets$M["1"], col = "red", xlim = range(mesh$loc[,1]),
      #        ylim = range(mesh$loc[,2]), add = TRUE)
      # } else {
      #   E1 <- res.exc$E
      #   E1[ind.rem] <- 0
      #   image(proj$x, proj$y, inla.mesh.project(proj, field = E1),asp = 1)
      # }

    }
  }

  #############
  while (!is.null(ind.rem)) {
    cat(sprintf("smallest area = %f (limit = %f)\n", areas$smallest.area,area.limit))
    rho.crop[ind.rem] <- 0
    rho[ind] <- rho.crop
    if (verbose)
      cat("Update excursion set\n")
    res.exc <-  excursions(mu = mu,
                           Q = Q,
                           alpha = alpha,
                           type = type,
                           u = u,
                           rho = rho,
                           ind = ind,
                           n.iter = n.iter,
                           Q.chol = Q.chol,
                           F.limit = F.limit,
                           method = method,
                           vars = vars,
                           max.size,
                           verbose = FALSE,
                           max.threads = max.threads,
                           seed = seed)
    res.exc <- crop.lists(res.exc)
    if (verbose)
      cat("Calculate sizes\n")
    if(use.continuous){
      sets <- excursions::continuous(ex = res.exc, geometry = mesh, alpha = alpha)
      areas <- find.smallest.activation(sets, mesh, area.limit = area.limit, factor = factor)
    } else {
      areas <- find.smallest.activation(res.exc, mesh, area.limit = area.limit, factor = factor,
                                        area.el= area.el)
    }
    ind.rem <- areas$ind.rem
    if (plot.progress) {
      # if(use.continuous){
      #   plot(mesh, vertex.color = rgb(0.8, 0.8, 0.8), draw.segments = FALSE,
      #        edge.color = rgb(0.8, 0.8, 0.8), main = " ")
      #   plot(sets$M["1"], col = "red", xlim = range(mesh$loc[,1]),
      #        ylim = range(mesh$loc[,2]), add = TRUE)
      # } else {
      #   E1 <- res.exc$E
      #   E1[ind.rem] <- 0
      #   image(proj$x, proj$y, inla.mesh.project(proj, field = E1),asp = 1)
      # }

    }
  }
  if(use.continuous){
    return(sets)
  } else {
    return(res.exc)
  }

}

#' Calculation of excursion sets with the removal of small areas
#'
#' @description Function similar to the excursions.inla function, but which also takes a mesh as input
#' and computes modified excursion sets where all regions with areas smaller than area.limit are
#' removed in the internal calculations.
#' @param result.inla Result object from INLA call.
#' @param stack The stack object used in the INLA call.
#' @param name The name of the component for which to do the calculation. This argument should
#' only be used if a stack object is not provided, use the tag argument otherwise.
#' @param tag The tag of the component in the stack for which to do the calculation. This argument
#' should only be used if a stack object is provided, use the name argument otherwise.
#' @param u Excursion or contour level.
#' @param alpha Error probability for the excursion set.
#' @param type Type of region:
#'  \itemize{
#'     \item{'>' }{positive excursion region}
#'     \item{'<' }{negative excursion region}
#'     \item{'!=' }{contour avoiding region}
#'     \item{'=' }{contour credibility region}}
#' @param n.iter Number or iterations in the MC sampler that is used for approximating probabilities. The default value is 10000.
#' @param F.limit The limit value for the computation of the F function. F is set to NA for all nodes where F<1-F.limit. Default is F.limit = \code{alpha}.
#' @param u.link If u.link is TRUE, u is assumed to be in the scale of the data and is then transformed to the scale of the linear predictor (default FALSE).
#' @param method Method for handeling the latent Gaussian structure:
#'  \itemize{
#'       \item{'EB' }{Empirical Bayes (default)}
#'       \item{'QC' }{Quantile correction, rho must be provided if QC is used.}}
#' @param ind Indices of the nodes that should be analysed (optional).
#' @param max.size Maximum number of nodes to include in the set of interest (optional).
#' @inheritParams verbose_Param_direct_FALSE
#' @inheritParams max.threads_Param
#' @inheritParams seed_Param
#' @param mesh The mesh on which the model is defined.
#' @param area.limit Positive number. All connected excursion sets with an area smaller than this
#' numbere are removed.
#' @param use.continuous Logical parameter indicating whether the areas of the excursion sets
#' should be calculated using the \code{continuous} function in \code{excursions}. If FALSE, the
#' function uses the approximation that the area for each node is the integral of the FEM basis function \eqn{\phi_i}.
#' @param plot.progress Logical parameter that indicates whether the results should be plotted.
#'
#' @return If \code{use.continuous = FALSE}, an item of class \code{excurobj}. Otherwise a list with
#' the same elements as the output of \code{continuous}.
#' @export
excursions.inla.no.spurious <- function(result.inla,
                            stack,
                            name=NULL,
                            tag=NULL,
                            u,
                            alpha=1,
                            type,
                            n.iter=10000,
                            F.limit,
                            u.link = FALSE,
                            method,
                            ind=NULL,
                            max.size,
                            verbose=FALSE,
                            max.threads=0,
                            seed=NULL,
                            mesh,
                            area.limit,
                            use.continuous = FALSE,
                            plot.progress = FALSE)


{
  if (!requireNamespace("INLA", quietly = TRUE))
    stop('This function requires the `INLA` package (see www.r-inla.org/download)')
  if (missing(result.inla))
    stop('Must supply INLA result object')
  if (missing(method)) {
    cat('No method selected, using EB\n')
    method = 'EB'
  }
  if (missing(mesh))
    stop('Must specify mesh')
  if (missing(area.limit))
    stop('Must specify area.limit')
  if (missing(u))
    stop('Must specify level u')
  if (missing(type))
    stop('Must specify type of excursion')

  if (result.inla$.args$control.compute$config == FALSE)
    stop('INLA result must be calculated using control.compute$config=TRUE')

  if (missing(F.limit)) {
    F.limit = alpha
  } else {
    F.limit = max(alpha,F.limit)
  }

  n = length(result.inla$misc$configs$config[[1]]$mean)

  #Get indices for the component of interest in the configs
  ind.stack <- excursions:::inla.output.indices(result.inla,
                                                name = name,
                                                stack = stack,
                                                tag = tag)
  n.out = length(ind.stack)
  #Index vector for the nodes in the component of interest
  ind.int <- seq_len(n.out)

  #ind is assumed to contain indices within the component of interest
  if (!missing(ind) && !is.null(ind)) {
    ind <- excursions:::private.as.vector(ind)
    ind.int <- ind.int[ind]
    ind.stack <- ind.stack[ind]
  }
  ind = ind.stack


  # If u.link is TRUE, the limit is given in linear scale
  # then transform to the scale of the linear predictor
  u.t = rho = rep(0,n)
  if (u.link == TRUE) {
    links = result.inla$misc$linkfunctions$names[
      result.inla$misc$linkfunctions$link]
    u.tmp = sapply(ind, function(i) excursions:::private.link.function(u,links[i]))
    u.t[ind] = u.tmp
  } else {
    u.t = u
  }


  #Calculate marginal probabilities
  #If stack and tag are provided, we are interested in the predictor
  #Otherwise, we are interested in some of the effects

  if (verbose)
    cat('Calculating marginal probabilities\n')

  #Are we interested in a random effect?
  random.effect = FALSE
  if (!missing(name) && !is.null(name) && name != "APredictor" && name != "Predictor") {
    random.effect = TRUE
    if (is.null(result.inla$marginals.random))
      stop('INLA result must be calculated using return.marginals.random=TRUE if excursion sets to be calculated for a random effect of the model')
  }

  if (!random.effect && is.null(result.inla$marginals.linear.predictor))
    stop('INLA result must be calculated using return.marginals.linear.predictor=TRUE if excursion sets are to be calculated for the linear predictor.')

  if (random.effect) {
    rho.ind <- sapply(1:length(ind), function(j) excursions:::inla.get.marginal(ind.int[j],
                                                                                u = u,
                                                                                result = result.inla,
                                                                                effect.name = name,
                                                                                u.link = u.link,
                                                                                type = type))
  } else {
    rho.ind <- sapply(1:length(ind), function(j) excursions:::inla.get.marginal(ind[j],
                                                                                u = u,
                                                                                result = result.inla,
                                                                                u.link = u.link,
                                                                                type = type))
  }
  rho[ind] = rho.ind

  n.theta = result.inla$misc$configs$nconfig
  for (i in 1:n.theta) {
    config = excursions:::private.get.config(result.inla,i)
    if (config$lp == 0)
      break
  }

  if (verbose)
    cat('Calculating excursion function using the ', method, ' method\n')

  if (method == 'EB' || method == 'QC') {
    res <- excursions.no.spurious(alpha = alpha,
                                  u = 0,
                                  mu = config$mu - u.t,
                                  Q = config$Q,
                                  mesh = mesh,
                                  area.limit = area.limit,
                                  type = type,
                                  method = method,
                                  use.continuous = use.continuous,
                                  F.limit = F.limit,
                                  vars = config$vars,
                                  rho = rho,
                                  ind = ind,
                                  n.iter = n.iter,
                                  max.threads = max.threads,
                                  seed = seed,
                                  verbose = verbose,
                                  plot.progress = plot.progress)

  } else {
   stop('Only implemented for EB and QC methods so far')
  }
  return(res)
}

# Internal function that crops output of the excursions function.
crop.lists <- function(in.list)
  {
  out <- in.list
  out$F <- out$F[in.list$meta$ind]
  out$G <- out$G[in.list$meta$ind]
  out$M <- out$M[in.list$meta$ind]
  out$E <- out$E[in.list$meta$ind]
  out$mean <- out$mean[in.list$meta$ind]
  out$vars <- out$vars[in.list$meta$ind]
  out$rho <- out$rho[in.list$meta$ind]
  out$meta$Fe <- out$meta$Fe[in.list$meta$ind]
  out$meta$ind <-NULL
  return(out)
}
