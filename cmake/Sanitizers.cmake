# ITCHBOOK_SANITIZE=address,undefined adds the sanitizers to every target in
# the build. UBSan is made fatal so CI fails on the first report instead of
# printing a warning and carrying on.
if(ITCHBOOK_SANITIZE)
  string(REPLACE "," ";" _san_list "${ITCHBOOK_SANITIZE}")
  set(_san_flags "-fsanitize=${ITCHBOOK_SANITIZE}" -fno-omit-frame-pointer -g)
  if("undefined" IN_LIST _san_list)
    list(APPEND _san_flags -fno-sanitize-recover=all)
  endif()
  add_compile_options(${_san_flags})
  add_link_options("-fsanitize=${ITCHBOOK_SANITIZE}")
  message(STATUS "itchbook: sanitizers enabled: ${ITCHBOOK_SANITIZE}")
endif()
