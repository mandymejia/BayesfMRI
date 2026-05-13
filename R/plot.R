#' S3 method: use \code{\link[ciftiTools]{view_xifti}} to plot a \code{"BGLM"} object
#'
#' @param x An object of class "BGLM"
#' @param Bayes \code{TRUE} for plotting Bayesian results, \code{FALSE} for plotting
#' classical GLM results. Default: \code{NULL}, which will use the Bayesian results
#' if available and the classical results if not.
#' @param idx Which field should be plotted? Give the numeric indices or the
#'  names. \code{NULL} (default) will show all fields. This argument overrides
#'  the \code{idx} argument to \code{\link[ciftiTools]{view_xifti}}.
#' @param title If NULL, the field names associated with idx will be used.
#' @param session Which session should be plotted? \code{NULL} (default) will
#'  use the first.
#' @param zlim Overrides the \code{zlim} argument for
#'  \code{\link[ciftiTools]{view_xifti}}. Default: \code{c(-1, 1)}.
#' @param ... Additional arguments to \code{\link[ciftiTools]{view_xifti}}
#'
#' @method plot BGLM
#'
#' @importFrom ciftiTools view_xifti
#' @export
#'
#' @return Result of the call to \code{ciftiTools::view_cifti}.
#'
plot.BGLM <- function(x, Bayes=NULL, idx=NULL, title=NULL, session=NULL, zlim=c(-1, 1), ...){

  # Method
  if (is.null(Bayes)) {
    method <- if (!is.null(x$estimate_xii$Bayes[[1]])) { "Bayes" } else { "classical" }
  } else if (isTRUE(Bayes) || isFALSE(Bayes)) {
    method <- if (isTRUE(Bayes)) { "Bayes" } else { "classical" }
  } else {
    stop('`Bayes` must be `TRUE`, `FALSE` or `NULL`.')
  }
  # if (is.null(x$estimate_xii[[method]])) {
  #   stop(paste("Method", gsub("betas_", "", method, fixed=TRUE), "does not exist."))
  # }

  # Session
  if (is.null(session)) {
    if (length(x$estimate_xii[[method]]) > 1) { message("Plotting the first session.") }
    session <- 1
  } else if (is.numeric(session)) {
    stopifnot(length(session)==1)
    stopifnot(session %in% seq(length(x$session_names)))
  }
  the_xii <- x$estimate_xii[[method]][[session]]
  if (is.null(the_xii)) { stop(paste("Session", session, "does not exist.")) }

  # Column index
  if (is.null(idx)) {
    idx <- seq_len(ncol(do.call(rbind, the_xii$data)))
  } else if (is.character(idx)) {
    idx <- match(idx, the_xii$meta$cifti$names)
  }

  # Names
  idx_names <- the_xii$meta$cifti$names[idx]

  # Title(s)
  if(is.null(title)){
    title <- idx_names
  }

  # Plot
  ciftiTools::view_xifti(the_xii, idx=idx, title=title, zlim=zlim, fname_suffix = idx_names, ...)
}

#' S3 method: use \code{\link[ciftiTools]{view_xifti}} to plot a \code{"act_BGLM"} object
#'
#' @param x An object of class "act_BGLM"
#' @param idx Which field should be plotted? Give the numeric indices or the
#'  names. \code{NULL} (default) will show all fields. This argument overrides
#'  the \code{idx} argument to \code{\link[ciftiTools]{view_xifti}}.
#' @param title If NULL, the field names associated with idx will be used.
#' @param session Which session should be plotted? \code{NULL} (default) will
#'  use the first.
#' @param ... Additional arguments to \code{\link[ciftiTools]{view_xifti}}
#'
#' @method plot act_BGLM
#'
#' @importFrom ciftiTools view_xifti
#' @export
#'
#' @return Result of the call to \code{ciftiTools::view_cifti_surface}.
#'
plot.act_BGLM <- function(x, idx=NULL, title=NULL, session=NULL, ...){

  # Session
  if (is.null(session)) {
    if (length(x$activations_xii) > 1) { message("Plotting the first session.") }
    session <- 1
  } else if (is.numeric(session)) {
    stopifnot(length(session)==1)
    stopifnot(session %in% seq(length(x$activations_xii)))
  }
  the_xii <- x$activations_xii[[session]]
  if (is.null(the_xii)) { stop(paste("Session", session, "does not exist.")) }

  # Column index
  if (is.null(idx)) {
    idx <- seq_len(ncol(do.call(rbind, the_xii$data)))
  } else if (is.character(idx)) {
    idx <- match(idx, the_xii$meta$cifti$names)
  }

  # Names
  idx_names <- the_xii$meta$cifti$names[idx]

  # Title(s)
  if(is.null(title)){
    title <- idx_names
  }

  # Values and colors
  vals <- unique(c(as.matrix(the_xii)))
  vals <- setdiff(vals, 0)
  #for a single level (e.g., only a single gamma value was given to activations()), use red to color activations
  if(length(vals) == 1){
    the_xii <- convert_to_dlabel(the_xii, colors = "red",
                                 labels = rownames(the_xii$meta$cifti$labels[[1]]))
  }

  # Plot
  ciftiTools::view_xifti(the_xii, idx=idx, title=title, fname_suffix = idx_names, ...)
}

#' S3 method: use \code{\link[ciftiTools]{view_xifti}} to plot a \code{"BGLM2"} object
#'
#' @param x An object of class \code{"BGLM2"}.
#' @param idx Which contrast should be plotted? Give the numeric indices or the
#'  names. \code{NULL} (default) will show all contrasts for \code{stat = "contrasts"},
#'  \code{stat = "activations"}, and \code{stat = "nested_activations"}.
#'  For \code{stat = "activation_levels"}, \code{idx} must specify exactly one contrast.
#' @param stat Estimates of the \code{"contrasts"} (default), the single-threshold
#'  \code{"activations"}, the \code{"nested_activations"}, or the per-threshold
#'  \code{"activation_levels"}.
#' @param level For \code{stat = "activation_levels"} only, the integer index (or
#'  indices) of the activation level(s) to plot within the selected contrast.
#'  Default: \code{NULL}, which shows all levels for that contrast.
#' @param alpha For \code{stat = "activation_levels"} only, the alpha threshold
#'  value(s) to plot within the selected contrast. Specify either \code{level} or
#'  \code{alpha}, but not both.
#' @param zlim Overrides the \code{zlim} argument for
#'  \code{\link[ciftiTools]{view_xifti}}. Default: \code{c(-1, 1)}.
#' @param ... Additional arguments to \code{\link[ciftiTools]{view_xifti}}.
#'
#' @method plot BGLM2
#'
#' @importFrom ciftiTools view_xifti
#' @export
#'
#' @return Result of the call to \code{ciftiTools::view_xifti()}.
#'
plot.BGLM2 <- function(
    x,
    idx = NULL,
    stat = c(
      "contrasts",
      "activations",
      "nested_activations",
      "activation_levels",
      "activated_contrast",
      "activated_level_contrast"
    ),
    level = NULL,
    alpha = NULL,
    zlim = c(-1, 1),
    ...
) {
  stat <- match.arg(stat)

  if (!is.null(level) && !is.null(alpha)) {
    stop("Specify only one of `level` or `alpha`.")
  }

  ## Shortcut:
  ## If user asks for nested_activations + level/alpha,
  ## interpret it as plotting exact activation levels for one contrast.
  if (stat == "nested_activations" && (!is.null(level) || !is.null(alpha))) {
    stat <- "activation_levels"
  }

  ## `level` and `alpha` are only meaningful for activation-level plots
  ## or activated-level contrast plots.
  if (
    stat == "activated_contrast" &&
    (!is.null(level) || !is.null(alpha))
  ) {
    stop(
      "`level` and `alpha` are not used when `stat = 'activated_contrast'`. ",
      "Use `stat = 'activated_level_contrast'` if you want to mask beta maps ",
      "by a specific nested activation level."
    )
  }

  ## Helper: resolve contrast index by numeric index or contrast name
  resolve_contrast_idx <- function(idx, available_names, what = "contrast") {
    if (is.null(idx)) return(NULL)

    if (is.character(idx)) {
      idx_match <- match(idx, available_names)
      if (anyNA(idx_match)) {
        bad <- idx[is.na(idx_match)]
        stop(
          "Unknown ", what, " name(s): ",
          paste(bad, collapse = ", "),
          ". Available names are: ",
          paste(available_names, collapse = ", ")
        )
      }
      return(idx_match)
    }

    idx <- as.integer(idx)
    if (anyNA(idx) || any(idx < 1L) || any(idx > length(available_names))) {
      stop(
        "`idx` is out of range. Valid indices are 1:",
        length(available_names), "."
      )
    }
    idx
  }

  ## Exact level maps
  if (stat == "activation_levels") {
    xii_list <- x$activation_levels_xii

    if (is.null(xii_list)) {
      stop(
        "No activation-level maps in `BGLM2` object. ",
        "Re-run `BayesGLM2()` with `alpha_nested`."
      )
    }

    contrast_names <- names(xii_list)
    if (is.null(contrast_names)) {
      contrast_names <- x$contrast_estimate_xii$meta$cifti$names
      names(xii_list) <- contrast_names
    }

    if (is.null(idx)) {
      stop(
        "For exact activation levels, please specify exactly one contrast in `idx`."
      )
    }

    if (length(idx) != 1L) {
      stop(
        "For exact activation levels, `idx` must specify exactly one contrast."
      )
    }

    contrast_idx <- resolve_contrast_idx(idx, contrast_names, what = "contrast")
    the_xii <- xii_list[[contrast_idx]]

    if (is.null(the_xii)) {
      stop("No activation-level maps found for the selected contrast.")
    }

    chosen_contrast <- contrast_names[contrast_idx]
    aa <- x$BayesGLM2_results$alpha_grid[[contrast_idx]]

    if (!is.null(alpha)) {
      alpha <- as.numeric(alpha)
      if (anyNA(alpha)) {
        stop("`alpha` must be numeric.")
      }

      match_alpha <- function(a, grid) {
        tol <- sqrt(.Machine$double.eps)
        hits <- which(abs(grid - a) < tol)
        if (length(hits) == 0L) return(NA_integer_)
        hits[1]
      }

      level_idx <- vapply(alpha, match_alpha, integer(1), grid = aa)

      if (anyNA(level_idx)) {
        bad <- alpha[is.na(level_idx)]
        stop(
          "Requested alpha value(s) not found for contrast `", chosen_contrast, "`: ",
          paste(bad, collapse = ", "),
          ". Available alpha values are: ",
          paste(aa, collapse = ", ")
        )
      }
    } else if (!is.null(level)) {
      level_idx <- as.integer(level)
      if (anyNA(level_idx) || any(level_idx < 1L) || any(level_idx > length(aa))) {
        stop(
          "`level` must contain integer indices between 1 and ", length(aa),
          " for contrast `", chosen_contrast, "`."
        )
      }
    } else {
      level_idx <- seq_along(aa)
    }

    ciftiTools::view_xifti(the_xii, idx = level_idx, zlim = zlim, ...)
    return(invisible(NULL))
  }

  ## Helper: mask beta maps by an activation-like xifti object
  mask_beta_by_activation <- function(beta_xii, mask_xii, beta_idx, mask_idx) {
    out_xii <- beta_xii

    for (part in names(out_xii$data)) {
      beta_part <- out_xii$data[[part]]
      mask_part <- mask_xii$data[[part]]

      if (is.null(beta_part) || is.null(mask_part)) next

      active_mask <- !is.na(mask_part[, mask_idx]) & mask_part[, mask_idx] != 0

      beta_part[!active_mask, beta_idx] <- NA_real_
      out_xii$data[[part]] <- beta_part
    }

    out_xii
  }

  ## Activated contrast:
  ## Use single-threshold activation maps to mask contrast beta maps.
  if (stat == "activated_contrast") {
    beta_xii <- x$contrast_estimate_xii
    act_xii  <- x$activations_xii

    if (is.null(beta_xii)) {
      stop("No contrast beta maps found in `x$contrast_estimate_xii`.")
    }

    if (is.null(act_xii)) {
      stop(
        "No single-threshold activation maps found in `x$activations_xii`. ",
        "Re-run `BayesGLM2()` with single-threshold activation settings, ",
        "or use `stat = 'contrasts'` to plot unmasked beta maps."
      )
    }

    beta_names <- beta_xii$meta$cifti$names
    act_names  <- act_xii$meta$cifti$names

    if (is.null(beta_names) || is.null(act_names)) {
      stop("Both beta maps and activation maps must have contrast names.")
    }

    available_names <- beta_names

    if (is.null(idx)) {
      idx <- seq_along(available_names)
    } else {
      idx <- resolve_contrast_idx(idx, available_names, what = "contrast")
    }

    masked_beta_xii <- beta_xii

    for (ii in idx) {
      contrast_name <- beta_names[ii]
      jj <- match(contrast_name, act_names)

      if (is.na(jj)) {
        stop(
          "No matching activation map found for contrast `",
          contrast_name,
          "`."
        )
      }

      masked_beta_xii <- mask_beta_by_activation(
        beta_xii = masked_beta_xii,
        mask_xii = act_xii,
        beta_idx = ii,
        mask_idx = jj
      )

      masked_beta_xii$meta$cifti$names[ii] <- paste0(
        contrast_name,
        " beta masked by activation"
      )
    }

    ciftiTools::view_xifti(masked_beta_xii, idx = idx, zlim = zlim, ...)
    return(invisible(NULL))
  }

  ## Activated level contrast:
  ## Use one nested activation level to mask the corresponding contrast beta map.
  if (stat == "activated_level_contrast") {
    beta_xii <- x$contrast_estimate_xii
    xii_list <- x$activation_levels_xii

    if (is.null(beta_xii)) {
      stop("No contrast beta maps found in `x$contrast_estimate_xii`.")
    }

    if (is.null(xii_list)) {
      stop(
        "No activation-level maps found in `x$activation_levels_xii`. ",
        "Re-run `BayesGLM2()` with `alpha_nested`."
      )
    }

    contrast_names <- names(xii_list)

    if (is.null(contrast_names)) {
      contrast_names <- beta_xii$meta$cifti$names
      names(xii_list) <- contrast_names
    }

    beta_names <- beta_xii$meta$cifti$names

    if (is.null(beta_names)) {
      stop("The contrast beta maps must have contrast names.")
    }

    if (is.null(idx)) {
      stop(
        "For `stat = 'activated_level_contrast'`, ",
        "please specify exactly one contrast in `idx`."
      )
    }

    if (length(idx) != 1L) {
      stop(
        "For `stat = 'activated_level_contrast'`, ",
        "`idx` must specify exactly one contrast."
      )
    }

    contrast_idx <- resolve_contrast_idx(idx, beta_names, what = "contrast")
    chosen_contrast <- beta_names[contrast_idx]

    level_contrast_idx <- match(chosen_contrast, contrast_names)

    if (is.na(level_contrast_idx)) {
      stop(
        "No activation-level maps found for contrast `",
        chosen_contrast,
        "`."
      )
    }

    level_xii <- xii_list[[level_contrast_idx]]

    if (is.null(level_xii)) {
      stop(
        "No activation-level map object found for contrast `",
        chosen_contrast,
        "`."
      )
    }

    aa <- x$BayesGLM2_results$alpha_grid[[level_contrast_idx]]

    if (!is.null(level) && !is.null(alpha)) {
      stop("Specify only one of `level` or `alpha`.")
    }

    if (is.null(level) && is.null(alpha)) {
      stop(
        "For `stat = 'activated_level_contrast'`, ",
        "please specify exactly one nested level using `level` or `alpha`."
      )
    }

    if (!is.null(alpha)) {
      alpha <- as.numeric(alpha)

      if (length(alpha) != 1L || anyNA(alpha)) {
        stop("`alpha` must be exactly one numeric value.")
      }

      match_alpha <- function(a, grid) {
        tol <- sqrt(.Machine$double.eps)
        hits <- which(abs(grid - a) < tol)
        if (length(hits) == 0L) return(NA_integer_)
        hits[1]
      }

      level_idx <- match_alpha(alpha, aa)

      if (is.na(level_idx)) {
        stop(
          "Requested alpha value not found for contrast `",
          chosen_contrast,
          "`: ",
          alpha,
          ". Available alpha values are: ",
          paste(aa, collapse = ", ")
        )
      }

      level_label <- paste0("alpha=", alpha)
    } else {
      level_idx <- as.integer(level)

      if (
        length(level_idx) != 1L ||
        is.na(level_idx) ||
        level_idx < 1L ||
        level_idx > length(aa)
      ) {
        stop(
          "`level` must be exactly one integer between 1 and ",
          length(aa),
          " for contrast `",
          chosen_contrast,
          "`."
        )
      }

      level_label <- paste0("level=", level_idx, ", alpha=", aa[level_idx])
    }

    masked_beta_xii <- mask_beta_by_activation(
      beta_xii = beta_xii,
      mask_xii = level_xii,
      beta_idx = contrast_idx,
      mask_idx = level_idx
    )

    masked_beta_xii$meta$cifti$names[contrast_idx] <- paste0(
      chosen_contrast,
      " beta masked by nested activation ",
      level_label
    )

    ciftiTools::view_xifti(
      masked_beta_xii,
      idx = contrast_idx,
      zlim = zlim,
      ...
    )

    return(invisible(NULL))
  }

  ## Other stats
  stat_name <- switch(
    stat,
    contrasts = "contrast_estimate_xii",
    activations = "activations_xii",
    nested_activations = "nested_activations_xii"
  )

  the_xii <- x[[stat_name]]

  if (is.null(the_xii)) {
    if (stat_name == "activations_xii") {
      stop(
        "No single-threshold activations in `BGLM2` object. ",
        "Use `stat = 'nested_activations'` for nested mode, or re-run `BayesGLM2()` with single-threshold activation settings."
      )
    }
    if (stat_name == "nested_activations_xii") {
      stop(
        "No nested activations in `BGLM2` object. ",
        "Re-run `BayesGLM2()` with `alpha_nested`."
      )
    }
    stop("Requested statistic is not available in this `BGLM2` object.")
  }

  ## Single-threshold activation labels
  if (stat_name == "activations_xii") {
    excur_type <- x$BayesGLM2_results$excursion_type
    contrast_names <- the_xii$meta$cifti$names

    if (!is.null(excur_type) && !identical(excur_type, "none")) {
      if (length(excur_type) == 1L) {
        excur_type <- rep(excur_type, length(contrast_names))
      }
      if (length(excur_type) == length(contrast_names)) {
        names(the_xii$meta$cifti$labels) <- paste0(
          contrast_names, ", '", excur_type, "'"
        )
      }
    }
  }

  ## Nested activation labels
  if (stat_name == "nested_activations_xii") {
    alpha_grid <- x$BayesGLM2_results$alpha_grid
    contrast_names <- the_xii$meta$cifti$names

    if (!is.null(alpha_grid) && length(alpha_grid) == length(contrast_names)) {
      names(the_xii$meta$cifti$labels) <- vapply(
        seq_along(contrast_names),
        function(i) {
          aa <- alpha_grid[[i]]
          paste0(
            contrast_names[i],
            " (nested: ",
            paste(paste0("alpha=", aa), collapse = " > "),
            ")"
          )
        },
        character(1)
      )
    }
  }

  available_names <- the_xii$meta$cifti$names

  if (is.null(idx)) {
    idx <- seq_along(available_names)
  } else {
    idx <- resolve_contrast_idx(idx, available_names, what = "contrast")
  }

  ciftiTools::view_xifti(the_xii, idx = idx, zlim = zlim, ...)
}

#' S3 method: use \code{\link[ciftiTools]{view_xifti}} to plot a \code{"prev_BGLM"} object
#'
#' @param x An object of class "prev_BGLM"
#' @param idx Which task should be plotted? Give the numeric indices or the
#'  names. \code{NULL} (default) will show all tasks. This argument overrides
#'  the \code{idx} argument to \code{\link[ciftiTools]{view_xifti}}.
#' @param session Which session should be plotted? \code{NULL} (default) will
#'  use the first.
#' @param drop_zeros Color locations without any activation across all results
#'  (zero prevalence) the same color as the medial wall? Default: \code{NULL} to
#'  drop the zeros if only one \code{idx} is being plotted.
#' @param colors,zlim See \code{\link[ciftiTools]{view_xifti}}.
# Here, the defaults are overrided to use the Viridis \code{"plasma"} color scale between
#  \code{1/nA} and 1, where \code{nA} is the number of results in \code{x}.
#' @param ... Additional arguments to \code{\link[ciftiTools]{view_xifti}}
#'
#' @method plot prev_BGLM
#'
#' @importFrom ciftiTools view_xifti
#' @importFrom fMRItools is_1
#' @export
#'
#' @return Result of the call to \code{ciftiTools::view_cifti_surface}.
#'
plot.prev_BGLM <- function(
  x, idx=NULL, session=NULL,
  drop_zeros=NULL, colors="plasma",
  #zlim=c(round(1/x$n_results-.005, 2), 1), ...){
  zlim=c(0, 1), ...){

  # Session
  if (is.null(session)) {
    if (length(x$prev_xii) > 1) { message("Plotting the first session.") }
    session <- 1
  } else if (is.numeric(session)) {
    stopifnot(length(session)==1)
    stopifnot(session %in% seq(length(x$prev_xii)))
  }
  the_xii <- x$prev_xii[[session]]
  if (is.null(the_xii)) { stop(paste("Session", session, "does not exist.")) }

  # Column index
  if (is.null(idx)) {
    idx <- seq_len(ncol(do.call(rbind, the_xii$data)))
  } else if (is.character(idx)) {
    idx <- match(idx, the_xii$meta$cifti$names)
  }

  if (is.null(drop_zeros)) { drop_zeros <- length(idx) == 1 }
  stopifnot(is_1(drop_zeros, "logical"))
  if (drop_zeros) {
    the_xii <- move_to_mwall(the_xii, 0)
  }

  # Plot
  ciftiTools::view_xifti(the_xii, idx=idx, colors=colors, zlim=zlim, ...)
}
