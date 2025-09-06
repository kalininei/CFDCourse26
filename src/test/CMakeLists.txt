set (SRC
	cfd26_test.cpp
	grid_test.cpp
	mat_test.cpp
	geom_test.cpp
		
	utils/filesystem.cpp
	utils/vecmat.cpp
)

set (HEADERS
	cfd26_test.hpp

	utils/filesystem.hpp
	utils/vecmat.hpp
)

add_executable(cfd26_test ${SRC} ${HEADERS})
target_link_libraries(cfd26_test PRIVATE ${CFD_TARGET})
target_link_libraries(cfd26_test PRIVATE Catch2::Catch2)
target_include_directories(cfd26_test PRIVATE ${CMAKE_SOURCE_DIR}/src/test)
target_compile_definitions(cfd26_test PUBLIC TEST_DIRECTORY="${CMAKE_SOURCE_DIR}/test_data/")
