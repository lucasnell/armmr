

new_armmMod <- function(.stan, .call, .hmc, .x_means_sds, .y_means_sds, .stan_data) {

    stopifnot(inherits(.call, "call"))
    stopifnot(inherits(.hmc, "logical"))
    stopifnot(is.null(.x_means_sds) || inherits(.x_means_sds, "data.frame"))
    stopifnot(is.null(.y_means_sds) || inherits(.y_means_sds, "data.frame"))

    # So it doesn't show the whole function if using do.call:
    if (.call[1] != as.call(quote(armm()))) {
        .call[1] <- as.call(quote(armm()))
    }

    armm_obj <- structure(list(stan = .stan, call = .call,
                              hmc = .hmc,
                              x_means_sds = .x_means_sds,
                              y_means_sds = .y_means_sds,
                              stan_data = .stan_data),
                         class = "armmMod")

    return(armm_obj)

}


#' Print an `armmMod` object.
#'
#' @export
#'
#' @noRd
#'
print.armmMod <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    summary(x, digits = digits, ...)
}

#' Summary for an `armmMod` object.
#'
#' Note that `Lower` and `Upper` fields are from 95% HPDIs.
#'
#' @export
#'
#' @noRd
#'
summary.armmMod <- function(object,
                            digits = max(3, getOption("digits") - 3),
                            ...) {

    cat("Autoregressive mixed model\n")
    cat(sprintf("  method:         %s\n",
                ifelse(object$hmc, "Hamiltonian Monte Carlo", "Direct optimization")))
    cat(sprintf("  family:         %s\n", family(object)))
    form <- if (inherits(object$call$formula, "formula")) object$call$formula else {
        eval(object$call$formula, parent.frame(1L))
    }
    cat("  formula:       ", paste(trimws(deparse(form)), collapse = " "), "\n")
    cat("  data:          ", paste(trimws(deparse(object$call$data)), collapse = " "))
    scale_strs <- c()
    if (!is.null(object$x_means_sds)) scale_strs <- c(scale_strs, "scaled X")
    if (!is.null(object$y_means_sds)) scale_strs <- c(scale_strs, "scaled Y")
    if (length(scale_strs) > 0) {
        cat(sprintf(" << %s >>\n", paste(scale_strs, collapse = " & ")))
    } else cat("\n")
    cat("  observations:  ", nobs(object), "\n")
    if (!is.null(object$call$rstan_control)) {
        cat("  rstan options: ", gsub("list\\(|\\)", "",
                                      deparse(object$call$rstan_control)), "\n")
    }
    cat(sprintf("  obs. error:     %s\n", "sig_obs" %in% names(object$stan)))
    cat("------\n")
    cat("Autoregressive parameters:\n")
    print(autoreg(object), digits = digits)

    # cat("------\n")
    # cat("Random effects:\n")
    # print(ranef(object), digits = digits)

    cat("------\n")
    cat("Fixed effects:\n")
    print(fixef(object), digits = digits)


    # z <- rstan::extract(object$stan, "log_lik_sum")[[1]]
    # armmr:::hpdi(z, 0.68)

}



# LEFT OFF HERE ----
# ranef.armmMod <- function(object, ...)
# object$stan





#' Extract autoregressive parameters
#'
#' @param object A fitted model with class `armmMod`
#' @param \dots Additional arguments, ignored for method compatibility.
#'
#' @return AR parameters.
#'
#' @export
autoreg <- function(object, ...) {
    UseMethod("autoreg")
}
#' @export
#' @describeIn autoreg Autoregressive parameters for an `armmMod` object
autoreg.armmMod <- function(object, ...) {
    if (is.null(eval(object$call$ar_form))) return(NULL)
    # Add `+0` to get all levels of groups:
    new_form <- as.formula(paste0(deparse(eval(object$call$ar_form)), "+0"))
    ar_names <- colnames(model.matrix(new_form, eval(object$call$data)))
    phis <- rstan::extract(object$stan, "phi")[[1]]
    cis <- lapply(1:ncol(phis), function(i) armmr:::hpdi(phis[,i], 0.95))
    autoreg_df <- data.frame(Median = apply(phis, 2, median),
                           Lower = sapply(cis, function(x) x[["lower"]]),
                           Upper = sapply(cis, function(x) x[["upper"]]))
    rownames(autoreg_df) <- ar_names
    return(autoreg_df)
}


#' Extract fixed-effects estimates from an `armmMod` object
#'
#' @param object A fitted model with class `armmMod`.
#' @param ... Ignored.
#'
#' @return A dataframe of fixed-effects estimates.
#'
#' @export
fixef.armmMod <- function(object, ...) {

    A <- rstan::extract(object$stan, "alpha")[[1]]

    cis <- lapply(1:ncol(A), function(i) armmr:::hpdi(A[,i], 0.95))

    fixef_df <- data.frame(Median = apply(A, 2, median),
                           Lower = sapply(cis, function(x) x[["lower"]]),
                           Upper = sapply(cis, function(x) x[["upper"]]))

    rownames(fixef_df) <- colnames(object$stan_data$x)

    return(fixef_df)
}


#' Residuals of `armmMod` objects
#'
#' Getting different types of residuals for `armmMod` objects.
#'
#' @param object A fitted model with class `armmMod`.
#' @param type Type of residuals, currently only `"response"` is programmed.
#' @param \dots Additional arguments, ignored for method compatibility.
#' @return A vector of residuals.
#' @export
residuals.communityPGLMM <- function(object,
                                     type = "response",
                                     ...) {
    if (family(object) == "normal"){
        y <- object$stan_data$y
        mu <- fitted(object)
        res <- y - mu
        return(res)
    }

    stop("\nresiduals only for normal so far.")
}


#' Fitted values for an `armmMod` object
#'
#' @param object A fitted model with class `armmMod`
#' @param \dots Additional arguments, ignored for method compatibility.
#' @return Fitted values.
#' @export
fitted.armmMod <- function(object, ...) {
    if (!object$hmc) {
        stop("\n`fitted.armmMod` method not yet implemented for direct ",
             "optimization", call. = FALSE)
    }
    apply(rstan::extract(object$stan, "y_pred")[[1]], 2, median)
}


#' Extracting the model frame from an `armmMod` object
#'
#' @inheritParams stats::model.frame
#' @method model.frame armmMod
#'
#' @export
model.frame.armmMod <- function(formula, ...) {
    model.frame(formula$formula, formula$data)
}


#' Number of observations in an `armmMod` object
#'
#' @inheritParams stats::nobs
#' @export
nobs.armmMod <- function(object, ...) {
    return(object$stan_data$n_obs)
}


#' Family for an `armmMod` object
#'
#' @inheritParams stats::family
#'
#' @export
family.armmMod <- function(object, ...) {
    fam <- ifelse(is.null(object$call$distr), formals(armm)[["distr"]],
                 object$call$distr)
    return(fam)
}

