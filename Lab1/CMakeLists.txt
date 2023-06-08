cmake_minimum_required(VERSION 3.6)
project(untitled)

set(CMAKE_CXX_STANDARD 11)

#Set your directories.  The dot representes the root application folder.
#Thus my the path to my domain folder:
set(APP_DOMAIN ./domain)
set(APP_CALC ./calc)

#Include the directories and now your cpp files will recognize your headers
include_directories(${APP_DOMAIN})
include_directories(${APP_CALC})

set(
        SOURCE_FILES
        main.cpp
        ${APP_DOMAIN}/person.h
        ${APP_DOMAIN}/person.cpp
        ${APP_CALC}/add.h
        ${APP_CALC}/add.cpp
        mm2.c mm3.c mm5_4.c mm5_8.c mm7.c mm8.c mm9.c)
add_executable(untitled ${SOURCE_FILES})

#CHEARS!