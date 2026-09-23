message(STATUS "Testing whether a double is evaluated as a double")
try_compile(HAVE_IEEE_DOUBLE_EVALUATION
    ${PROJECT_BINARY_DIR}/config
    ${PROJECT_SOURCE_DIR}/config/excess-precision.cpp)
if(HAVE_IEEE_DOUBLE_EVALUATION)
    message(STATUS "Testing whether a double is evaluated as a double -- Success")
else()
    message(STATUS "Testing whether a double is evaluated as a double -- Failed (excess precision, x87?)")
endif()
