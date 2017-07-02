## calculate similarity of daily discharge patterns

#' Get Dry Weather Flow Patterns
#'
#' @param x The discharge time series (xts)
#' @param h Height of Cluster (defaults to 0.1)
#' @param separate_weekend logical. Should the dry weather flow separately analyzed?
#' @param make_equidistant logical. Should the time series be made equidistant?
#' @param plot logical. Should the sequences be plotted?
#' @param fig.opts USe additional fig.ops (currently legend.position only)
#' @param verbose logical. Need more information?
#' @return A list with dwf pattern (median) and sequences.
#' @export get_dwf_pattern
#' @importFrom magrittr "%>%"
get_dwf_pattern <- function(x,
                            h = 0.1,
                            separate_weekend = FALSE,
                            make_equidistant = FALSE,
                            plot = TRUE,
                            fig.opts = NULL,
                            verbose = FALSE) {


  # take the original discharge time series and create a list for weekdays and
  # weekends
  list_of_daily_xts <- .prepare_daily_series(x = x,
                                             separate_weekend = separate_weekend,
                                             make_equidistant = make_equidistant)

  # perform algorithm for weekdays and weekends separately
  dwf_seq <- lapply(list_of_daily_xts, function(weektype) {

    # take list of xts objects
    t <- lapply(weektype, zoo::coredata)

    # remove days which have no variance
    t <- t[which(unlist(lapply(t, stats::sd)) != 0)]

    ## create a matrix and zscale it
    t_scaled <- do.call(cbind, t) %>%
      scale

    # perform SBD distance calculation and hierarchical clusters
    hcl <- proxy::dist(t_scaled,
                       method = "SBD",
                       by_rows = FALSE) %>%
      base::as.matrix(.) %>%
      stats::as.dist(.) %>%
      stats::hclust(d = ., method = "complete")

    # cut the hclut tree based on height
    ct <- stats::cutree(hcl, h = h)

    # get elements within the most densed cluster aka most similar sequences
    dwf_seq <- which(ct == Mode(ct))

    return(dwf_seq)

  })

  # How many sequences are in the most densed cluster?
  if (verbose) message(paste("No of sequences in the most densed cluster:",
                             paste(sapply(dwf_seq, length), collapse = ",")))

  # subset list based on dwf_seq
  list_of_xts_pattern <- mapply("[", list_of_daily_xts, dwf_seq, SIMPLIFY = FALSE)

  # compute mean dry weather flow
  list_of_median_xts_pattern <- lapply(list_of_xts_pattern, function(x) {
    ret <- cbind(matrixStats::colMedians(t(do.call(cbind, lapply(x, zoo::coredata)))))
    return(ret)
  })

  # return a list with dry weather flow days and median
  res <- list(dwf_sequences = list_of_xts_pattern,
              dwf_median = list_of_median_xts_pattern)


  if (plot) .plot_pattern(x = res,
                          fig.opts = fig.opts)

  return(res)

}

#' @keywords internal
#' @importFrom magrittr "%>%"
.plot_pattern <- function(x, fig.opts = NULL) {

  if (is.null(fig.opts)) {
    fig.opts <- list(title = "NN",
                     legend.position = "right")
  }

  panel_list <- lapply(x$dwf_sequences, function(x) {

      # extract the index for plotting
      # we take the index of the first element and assume it to be the same for
      # the others
      index <- strftime(zoo::index(x[[1]]),
                        format = "%H:%M:%S") %>%
        # add a dummy date
        paste("2016-12-12", .) %>%
        # make posixct again
        as.POSIXct(., tz = xts::tzone(x[[1]]))

      x <- lapply(x, zoo::coredata)

      # create data.frame to be ggplotable
      m <- data.frame(index, do.call(cbind, x))

      colmed <- matrixStats::colMedians(x = t(m[,-1]))

      colnames(m) <- c("index", seq_len(length(x)))
      dfcm <- data.frame(index = index, colmed = colmed)

      gg_obj <- m %>%
        # wide to long
        tidyr::gather(id, discharge, -index) %>%
        # change id to factor to be sortable
        dplyr::mutate(id = factor(id,
                                  levels = unique(as.numeric(id)))) %>%
        ggplot2::ggplot(., ggplot2::aes(x = index,
                                        y = discharge,
                                        color = id)) +
        ggplot2::geom_line() +
        ggplot2::geom_line(data = dfcm,
                           mapping = ggplot2::aes(x = index, y = colmed),
                           size = 2,
                           linetype = 1,
                           color = "gray45") +
        # scale x labels
        ggplot2::scale_x_datetime(labels = scales::date_format("%H:%M")) +
        ggplot2::labs(x = "",
                      y = "discharge (l/s)",
                      subtitle = paste("no of sequences:", length(x))) +
        #.udtools_plot_theme() +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = fig.opts$legend.position)

      return(gg_obj)

    })

  # change labs manually
  for (i in 1:length(panel_list)) {
    tmp <- panel_list[[i]]$labels$subtitle
    panel_list[[i]]$labels$subtitle <- paste0(names(panel_list)[i], " (", tmp, ")")
  }

  # make grobs
  grob_list <- lapply(panel_list, ggplot2::ggplotGrob)

  # set same widths
  widths <- do.call(grid::unit.pmax, lapply(grob_list, "[[", "widths"))
  grob_list <- lapply(grob_list, function(g) {g$widths <- widths; g})

  # new page
  grid::grid.newpage()

  if (!is.null(fig.opts$top)) {
    t <- paste()
  }

  # assemble grobs
  p <- gridExtra::arrangeGrob(grobs = grob_list, top = "Analysis of Dry Weather Flow Pattern",
                              ncol = 1,
                              heights = rep(4,length(grob_list)))

  # plot
  grid::grid.draw(p)

}

#' @param x The xts object to prepare
#' @param separate_weekend If weekends should be separated
#' @param make_equidistant Equidistance required?
#' @keywords internal
.prepare_daily_series <- function(x, separate_weekend = FALSE, make_equidistant = TRUE) {

  # we need the following packages:
  # dtwclust, magrittr

  # set zero to NA
  x[x < 0] <- NA

  # remove NA
  x <- stats::na.omit(x)

  # error handling
  if (nrow(x) == 0) {
    warning(paste(names(x), "only NA"))
    return(NULL)
  }

  if (separate_weekend) {

    # create a list of xts objects distinguished by weekday and weekend
    list_of_xts <- separate_xts(xts = x,
                                interval = "wday",
                                index = list(1:5, c(0,6)),
                                names = c("weekdays", "weekend"))

  } else {

    # create a list of xts objects with no further segmentation
    list_of_xts <- separate_xts(xts = x,
                                interval = "wday",
                                index = list(0:6),
                                names = c("weekdays and weekend"))

  }

  # create now daily xts objects
  list_of_daily_series <- lapply(list_of_xts,
                                 FUN = function(x) separate_xts(xts = x,
                                                                interval = "day",
                                                                index = unique(xts::.indexday(x))))
  # need to make equidistant?
  if (make_equidistant) {
    list_of_daily_series <- list_of_daily_series %>%
      purrr::map(., ~ purrr::map(., ~ udtools::make_equidistant(., mode = 3, by = "15 mins", maxgap = 4, all = c(FALSE, TRUE))))
  }

  # remove all days with too little data (e.g. less than 24, 48, 96 values (depends on the frequency...))
  # first estimate median periodicity and find nu
  m_period <- 60 / as.numeric(mode_periodicity(x, units = "mins")) * 24

  # now use R's Higher-Order Function to filter the list
  list_of_daily_series_cleaned <- lapply(list_of_daily_series,
                                         function(x) Filter(function(y) length(y) == m_period, x))


  # if not sequences are generated, create an equistant one but give warning
  if (max(unlist(lapply(list_of_daily_series_cleaned, length))) == 0) {
    warning("sequences needed to be harmonized.")

    list_of_daily_series <- list_of_daily_series %>%
      purrr::map(., ~ purrr::map(., ~ udtools::make_equidistant(., mode = 3, by = "15 mins", maxgap = 4, all = c(FALSE, TRUE))))

    # remove all days with too little data (e.g. less than 24, 48, 96 values (depends on the frequency...))
    # first estimate median periodicity and find nu
    m_period <- 60 / as.numeric(mode_periodicity(x, units = "mins")) * 24

    # now use R's Higher-Order Function to filter the list
    list_of_daily_series_cleaned <- lapply(list_of_daily_series,
                                           function(x) Filter(function(y) length(y) == m_period, x))

  }


  return(list_of_daily_series_cleaned)

}
