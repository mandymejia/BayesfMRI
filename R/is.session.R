#' Title
#'
#' @param sess 
#' @return
#' @export
#'
#' @examples
is.session <- function(sess){
                                        # a valid "session" object is a list with precisely the following fields
                                        # - BOLD: T \times V matrix of BOLD responses, rows are time points, columns are voxels
                                        # - design: T \times K matrix, where colnames holds the stimuli names
                                        # - nuisance: T \times J matrix, where colnames holds the nuisance regressor names

                                        # these common dimensions must match across the objects
                                        # - T: The number of time points
                                        # - V: the number of voxels
                                        # - K: the number of stimuli (regressors)
                                        # - J: the number of nuisance regressors
                                        #
                                        # the checks that we need to run are
                                        # 1. are there precisely three fields?
                                        # 2. do the field names match
                                        # 3. are each fields of the correct type
                                        # 4. do the dimensions match

    ## check number of fields
    if(length(sess) != 3){stop('I expected the session to have 3 fields, but it does not.')}

    ## check identities of fields
    fields = c('BOLD','design','nuisance')
    if(! all.equal(names(sess),fields)){stop(
                                            paste0('You are missing the following fields',setdiff(names(sess),fields)))}

    ## check each field's type
    if(! (is.numeric(sess$BOLD))){stop('I expected BOLD to be numeric, but it is not')}
    if(! (is.matrix(sess$BOLD))){stop('I expected BOLD to be a matrix, but it is not')}

    if(! (is.matrix(sess$design))){stop('I expected design to be a matrix, but it is not')}
    if(! (is.matrix(sess$design))){stop('I expected design to be a matrix, but it is not')}

    if(! (is.matrix(sess$nuisance))){stop('I expected nuisance to be a matrix, but it is not')}
    if(! (is.matrix(sess$nuisance))){stop('I expected nuisance to be a matrix, but it is not')}

    ## check the dimensions of each field: T
    if(nrows(BOLD) != nrows(design)){stop("BOLD and design don't have the same number of rows (time points)")}
    if(nrows(BOLD) != nrows(nuisance)){stop("BOLD and nuisance don't have the same number of rows (time points)")}
    if(nrows(design) != nrows(nuisance)){stop("design and nuisance don't have the same number of rows (time points)")}

    
    return(TRUE)    
    

}
