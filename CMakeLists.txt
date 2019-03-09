cmake_minimum_required(VERSION 3.0)

include(${CMAKE_CURRENT_SOURCE_DIR}/../../max-api/script/max-pretarget.cmake)

include_directories( 
	"${C74_INCLUDES}"
)

include_directories(
        "${CMAKE_CURRENT_SOURCE_DIR}/mi"
)

set(MARBLES_SRC
        mi/marbles/resources.cc
        mi/marbles/random/random_generator.h
        mi/marbles/random/random_stream.h
        mi/marbles/random/t_generator.cc
        mi/marbles/random/x_y_generator.cc
        mi/marbles/note_filter.h
        mi/marbles/ramp/ramp_extractor.cc
        mi/stmlib/dsp/units.cc
        mi/marbles/random/output_channel.cc
        mi/marbles/random/lag_processor.cc
        mi/marbles/random/quantizer.cc
)

set(MIMARBLESLIB_SRC 
       ${MARBLES_SRC}
)

add_library(MIMarblesLib ${MIMARBLESLIB_SRC} )

add_library( 
	${PROJECT_NAME} 
	MODULE
	${PROJECT_NAME}.cpp
)

target_link_libraries(${PROJECT_NAME} MIMarblesLib)

include(${CMAKE_CURRENT_SOURCE_DIR}/../../max-api/script/max-posttarget.cmake)