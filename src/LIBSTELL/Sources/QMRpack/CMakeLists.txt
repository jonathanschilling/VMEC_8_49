
list (APPEND vmec_sources
	${CMAKE_CURRENT_SOURCE_DIR}/daxpby.f
	${CMAKE_CURRENT_SOURCE_DIR}/drandn.f
	${CMAKE_CURRENT_SOURCE_DIR}/dutfx.f
)
set (vmec_sources "${vmec_sources}" PARENT_SCOPE)

