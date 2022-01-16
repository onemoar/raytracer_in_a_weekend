@echo off

mkdir build
pushd build
rem set linker_flags=/subsystem:console /entry:WinMainCRTStartup

cl -O2 -Oi -GL -fp:fast -fp:except-  /D NDEBUG /D _WINDOWS -MT -GS- -Gy -Gd -wd4068 -wd4530 -wd4577 -nologo ..\main.cpp user32.lib gdi32.lib /link %linker_flags%

popd