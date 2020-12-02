

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
#' @method print armmMod
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
#' @method summary armmMod
#'
#' @noRd
#'
summary.armmMod <- function(object,
                            digits = max(3, getOption("digits") - 3),
                            ...) {

    cat("Autoregressive mixed model\n")
    cat(sprintf("  method:         %s\n",
                ifelse(object$hmc, "Hamiltonian Monte Carlo",
                       "Direct optimization")))
    cat(sprintf("  family:         %s\n", family(object)))
    form <- if (inherits(object$call$formula, "formula")) {
        object$call$formula
    } else {
        eval(object$call$formula, parent.frame(1L))
    }
    cat("  formula:       ", paste(trimws(deparse(form)), collapse = " "), "\n")
    cat("  data:          ", paste(trimws(deparse(object$call$data)),
                                   collapse = " "))
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
    LL <- median(rstan::extract(object$stan, "log_lik_sum")[[1]])
    cat("Median posterior logLik:", LL, "\n")
    cat("------\n")

    AR <- autoreg(object)
    if (!is.null(AR)) {
        cat("Autoregressive parameters:\n")
        rownames(AR) <- paste0(" ", rownames(AR))
        print(AR, digits = digits)
    } else cat("No autoregressive parameters\n")

    if (any(grepl("^sig_beta", names(object$stan)))) {
        cat("------\n")
        cat("Random effects:\n")
        print_sigma_betas(object, digits = digits)
    } else cat("------\nNo random effects\n")

    cat("------\n")
    cat("Fixed effects:\n")
    FE <- fixef(object)
    rownames(FE) <- paste0(" ", rownames(FE))
    print(FE, digits = digits)

    invisible(NULL)

}




}


# #' @noRd
# #'
# #' @export
# ranef <- function(object, ...) {
#     UseMethod("ranef")
# }
# LEFT OFF HERE ----
# ranef.armmMod <- function(object, ...)
# object$stan

# object$stan_data$b_groups
#
# names(object$stan)[grepl("beta", names(object$stan))]




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
#' @param method A single string, for either using quantiles (`"quantile"`)
#'     or HPDI (`"hpdi"`) to compute the standard errors.
#'     Defaults to `"quantile"`.
#' @method autoreg armmMod
#' @describeIn autoreg Autoregressive parameters for an `armmMod` object
autoreg.armmMod <- function(object,
                            method = c("quantile", "hpdi"),
                            ...) {

    if (is.null(eval(object$call$ar_form))) return(NULL)

    method <- match.arg(tolower(method), c("quantile", "hpdi"))
    # Add `+0` to get all levels of groups:
    new_form <- as.formula(paste0(deparse(eval(object$call$ar_form)), "+0"))
    ar_names <- colnames(model.matrix(new_form, eval(object$call$data)))
    phis <- rstan::extract(object$stan, "phi")[[1]]

    if (method == "quantile") {
        upper <- unname(sapply(1:ncol(phis), function(i) quantile(phis[,i], 0.84)))
        lower <- unname(sapply(1:ncol(phis), function(i) quantile(phis[,i], 0.16)))
    } else {
        ints <- lapply(1:ncol(phis), function(i) hpdi(phis[,i], 0.68))
        upper <- sapply(ints, function(x) x[["upper"]])
        lower <- sapply(ints, function(x) x[["lower"]])
    }
    SEs <- 0.5 * (upper - lower)

    autoreg_df <- data.frame(Median = apply(phis, 2, median),
                             `Std.Error` = SEs)
    rownames(autoreg_df) <- ar_names
    return(autoreg_df)
}





#' @name fixef
#' @title Extract fixed-effects estimates from an `armmMod` object
#' @aliases fixef fixed.effects fixef.armmMod
#' @docType methods
#' @param object A fitted model with class `armmMod`.
#' @param method A single string, for either using quantiles (`"quantile"`)
#'     or HPDI (`"hpdi"`) to compute the standard errors.
#'     Defaults to `"quantile"`.
#' @param ... Ignored.
#' @return A dataframe of fixed-effects estimates. Standard errors are
#'     based on either quantiles or HPDIs, and are half the width of
#'     the 68% uncertainty interval.
#' @importFrom lme4 fixef
#' @export fixef
#' @method fixef armmMod
#' @export
fixef.armmMod <- function(object, method = c("quantile", "hpdi"), ...) {

    method <- match.arg(tolower(method), c("quantile", "hpdi"))

    A <- rstan::extract(object$stan, "alpha")[[1]]

    if (method == "quantile") {
        upper <- unname(sapply(1:ncol(A), function(i) quantile(A[,i], 0.84)))
        lower <- unname(sapply(1:ncol(A), function(i) quantile(A[,i], 0.16)))
    } else {
        ints <- lapply(1:ncol(A), function(i) hpdi(A[,i], 0.68))
        upper <- sapply(ints, function(x) x[["upper"]])
        lower <- sapply(ints, function(x) x[["lower"]])
    }
    SEs <- 0.5 * (upper - lower)

    fixef_df <- data.frame(Median = apply(A, 2, median),
                           `Std.Error` = SEs)

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
#'
#' @method residuals armmMod
#'
#' @export
residuals.armmMod <- function(object,
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
#'
#' @method fitted armmMod
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
#'
#' @export
model.frame.armmMod <- function(formula, ...) {
    model.frame(formula$formula, formula$data)
}


#' Number of observations in an `armmMod` object
#'
#' @inheritParams stats::nobs
#' @method nobs armmMod
#' @export
nobs.armmMod <- function(object, ...) {
    return(object$stan_data$n_obs)
}


#' Family for an `armmMod` object
#'
#' @inheritParams stats::family
#' @method family armmMod
#'
#' @export
family.armmMod <- function(object, ...) {
    fam <- ifelse(is.null(object$call$family), formals(armm)[["family"]],
                 object$call$family)
    return(fam)
}

