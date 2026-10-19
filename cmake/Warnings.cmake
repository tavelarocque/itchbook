function(itchbook_set_warnings target)
  if(MSVC)
    target_compile_options(${target} PRIVATE /W4)
    return()
  endif()
  target_compile_options(${target} PRIVATE
    -Wall -Wextra -Wpedantic
    -Wshadow -Wconversion -Wsign-conversion
    -Wold-style-cast -Wnon-virtual-dtor -Woverloaded-virtual
    -Wdouble-promotion -Wimplicit-fallthrough
  )
  if(ITCHBOOK_WERROR)
    target_compile_options(${target} PRIVATE -Werror)
  endif()
endfunction()
