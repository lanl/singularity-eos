#------------------------------------------------------------------------------
# placeholder
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# override the default 'message' command
# useful for QoL improvements like suppressing
# output and coloring.
#------------------------------------------------------------------------------
function(message)
  if (NOT MESSAGE_QUIET)
    _message(${ARGN})
  endif()
endfunction()
