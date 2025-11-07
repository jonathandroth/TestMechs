#' @title This function processes the inputted dataframe to remove missing values
#' @param df A data.frame
#' @param d Name of the treatment variable in the df
#' @param m Name of the mediator variable
#' @param y Name of the outcome variable
#' @param w (Optional) Name of weighting variable. If null, equal weights are used
#' @param reg_formula (Optional) A fixest regression formula. If provided, we remove observations that have missing covariates for the regression
#' @export
remove_missing_from_df <- function(df, d, m, y, w = NULL, reg_formula = NULL){

  missing_d <- is.na(df[[d]])
  missing_m <- !complete.cases(df[,m]) #allow for vector-valued m
  missing_y <- is.na(df[[y]])

  df <- df[!(missing_d | missing_m | missing_y), , drop = FALSE]

  if (!is.null(w)) {
    missing_w <- is.na(df[[w]])
    df <- df[!missing_w, , drop = FALSE]
  }

  if(!is.null(reg_formula)){
    mod <- fixest::feols(reg_formula, data = df, notes = FALSE)

    if(!is.null(mod$obs_selection$obsRemoved)){
      df <- df[mod$obs_selection$obsRemoved , , drop = FALSE]
    }
  }

  return(df)
}
