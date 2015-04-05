set(gridoperatorimexonestepdir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/gridoperator/imexonestep)

set(gridoperatorimexonestep_HEADERS
        enginebase.hh
        jacobianengine.hh
        localassembler.hh
        patternengine.hh
        prestageengine.hh
        residualengine.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${gridoperatorimexonestep_HEADERS} DESTINATION ${gridoperatorimexonestepdir})
