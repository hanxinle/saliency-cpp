@echo off
rem Build saliency.exe with VS 2022 MSVC (x64 Release).
call "C:\Program Files\Microsoft Visual Studio\2022\Enterprise\VC\Auxiliary\Build\vcvarsall.bat" x64 >nul
if errorlevel 1 (
    echo vcvarsall.bat failed - check VS 2022 install path.
    exit /b 1
)
cl /O2 /std:c++17 /EHsc /W3 /utf-8 /nologo /Fe:saliency.exe main.cpp salmap.cpp wta.cpp segment.cpp
if errorlevel 1 exit /b 1
echo Build OK: saliency.exe
