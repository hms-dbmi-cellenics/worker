#' generateErrorMessage
#'
#' @param code error code.
#' @param user_message message to display to a user.
#'
#' @return
#' @export
#'
#' @examples
generateErrorMessage <- function(code, user_message) {
  return(paste0(code, ":|:", user_message))
}

#' extractErrorList
#'
#' @param error_message generated with generateErrorMessage.
#'
#' @return sepparates an error message into code and user message.
#' @export
#'
#' @examples
extractErrorList <- function(error_message) {
  error_string <- unlist(strsplit(error_message, ":|:", fixed = TRUE))

  is.expected <- !is.na(error_string[2])
  error_code <- ifelse(is.expected, error_string[1], error_codes$UNHANDLED_ERROR)
  user_message <- ifelse(is.expected, error_string[2], error_message)

  return(
    list(
      error_code = error_code,
      user_message = user_message
    )
  )
}

#' format response
#'
#' @param data
#' @param error
#'
#' @return formats response for the UI.
#' @export
#'
#' @examples
formatResponse <- function(data, error) {
  return(
    list(
      data = data,
      error = error
    )
  )
}
