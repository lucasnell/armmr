context("test-rnd_names.R")


test_that("proper rnd_names output", {


    expand <- function(...) expand.grid(..., KEEP.OUT.ATTRS = FALSE,
                                        stringsAsFactors = FALSE)

    form <- y ~ x1 * x2 + x3 + (x1 * x2 | g1 + g2) + (x3 | g1) + (1 | g2) + x4
    time_form <- ~ t | tg + g1
    ar_form <- ~ g1
    y_scale <- ~ g1

    set.seed(1)
    x1_coef <- 1.5
    x3_coefs <- runif(10, 1, 5)
    data <- data.frame(g1 = factor(rep(1:5, each = 20)),
                       g2 = factor(rep(2:1, each = 50)),
                       y = rnorm(100),
                       x1 = runif(100),
                       x2 = factor(rep(1:4, 25)),
                       x3 = rnorm(100),
                       x4 = rnorm(100),
                       t = rep(1:10, 10),
                       tg = factor(rep(1:10, each = 10)))
    data$y <- data$y + data$x1 * x1_coef
    data$y <- data$y + data$x3 *
        sapply(as.integer(interaction(data$g1, data$g2)),
               function(i) x3_coefs[i])

    # Suppressing warnings bc I'm only interested in the output, not a
    # properly fitting model
    mod <- suppressWarnings(armm(form, time_form, ar_form, y_scale, data,
                                 rstan_control = list(chains = 1, iter = 40,
                                                      verbose = FALSE,
                                                      show_messages = FALSE)))

    rnd_names <- rbind(
        expand(c("g2"), "(Intercept)"),
        expand(c("g1", "g2"), c("x1", "x22", "x23", "x24")),
        expand(c("g1"), "x3"),
        expand(c("g1", "g2"), c("x1:x22", "x1:x23", "x1:x24")))
    rnd_names <- rnd_names[,2:1]
    colnames(rnd_names) <- c("Name", "Groups")


    rnd_lvl_names <- armmr:::f_apply(split(rnd_names, 1:nrow(rnd_names)),
                                     function(x) {
                                         lvls <- levels(eval(parse(text = x$Groups),
                                                             data))
                                         expand(Name = x$Name,
                                                Groups = x$Groups,
                                                Level = lvls)
                                     }, rbind)
    rownames(rnd_lvl_names) <- NULL



    expect_equal(rnd_names, mod$rnd_names)
    expect_equal(rnd_lvl_names, mod$rnd_lvl_names)
})
