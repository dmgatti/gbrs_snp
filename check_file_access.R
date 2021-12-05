################################################################################
# Given a file path and a permission level, check that the file exists and
# has the correct permissoins (i.e. read/write).
#
# Daniel Gatti
# dan.gatti@jax.org
# 2021-12-04
################################################################################

# Arguments:
# filename: string containing the full path to the file to test.
# mode: integer indicating existance, read, write. From file.access().
#       0: test for existence.
#       1: test for execute permission.
#       2: test for write permission.
#       4: test for read permission. 
# Returns: nothing. Stops execution if is not found and has correct permission.
check_file_access = function(filename, mode) {

  if(!file.exists(filename)) {
    stop(paste('ERROR: file not found:', filename))
  }
 
  # 0 indicates success from file.access()
  if(file.access(filename, mode) != 0) {
    stop(paste('ERROR: file does not have correct permissions: ', filename))
  }

} # check_file_access()
