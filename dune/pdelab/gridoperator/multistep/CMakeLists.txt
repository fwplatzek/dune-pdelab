set(gridoperatormultistepdir ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/gridoperator/multistep)

set(gridoperatormultistep_HEADERS
  localassembler.hh
  )

install(FILES ${gridoperatormultistep_HEADERS} DESTINATION ${gridoperatormultistepdir})
