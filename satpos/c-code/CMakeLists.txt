cmake_minimum_required(VERSION 3.0)
enable_language(C)
project(satpos)

set(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_BUILD_TYPE Release)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu99 -lm -Wall -Wextra -lcheck")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Og -ggdb -Werror -Wno-error=unused-result -Wno-error=unused-parameter -Wno-error=unused-variable -Wno-error=unused-function")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O2")

add_executable(
	main.out
	src/io/file.c
	src/io/csv_reader.c
	src/coordinate.c
	src/coordinate_factory.c
	src/util/array.c
	src/util/str.c
	
	test/io/file_test.c
	test/io/csv_reader_test.c
	test/coordinate_test.c
	test/coordinate_factory_test.c
	test/util/array_test.c
	test/util/str_test.c
	test/test.c
)
